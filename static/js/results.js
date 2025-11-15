// version1.0
// ========================================
// RESULTS.JS - Results Page JavaScript
// ========================================

// GLOBAL selected epitopes storage
var selectedEpitopes = {
    bcell: [],
    mhc1: [],
    mhc2: []
};

$(document).ready(function() {
    console.log('Results page loaded, initializing...');
    
    // Initialize DataTables with minimal config first
    var bcellTable = $('#bcellTable').DataTable({
        pageLength: 10,
        lengthMenu: [10, 25, 50, 100],
        order: [[1, 'desc']],
        scrollX: true,
        autoWidth: false,
        drawCallback: function() {
            console.log('B-cell table redrawn, restoring states');
            attachCheckboxHandlers();
        }
    });
    
    var mhc1Table = $('#mhc1Table').DataTable({
        pageLength: 10,
        lengthMenu: [10, 25, 50, 100],
        order: [[1, 'desc']],
        scrollX: true,
        autoWidth: false,
        drawCallback: function() {
            console.log('MHC-I table redrawn, restoring states');
            attachCheckboxHandlers();
        }
    });
    
    var mhc2Table = $('#mhc2Table').DataTable({
        pageLength: 10,
        lengthMenu: [10, 25, 50, 100],
        order: [[1, 'desc']],
        scrollX: true,
        autoWidth: false,
        drawCallback: function() {
            console.log('MHC-II table redrawn, restoring states');
            attachCheckboxHandlers();
        }
    });
    
    // Initial attachment
    attachCheckboxHandlers();
    
    // Function to attach handlers and restore states
    function attachCheckboxHandlers() {
        // Remove any existing handlers first
        $('.bcell-checkbox, .mhc1-checkbox, .mhc2-checkbox').off('change');
        
        // B-cell checkboxes
        $('.bcell-checkbox').each(function() {
            var checkbox = $(this);
            var epitope = checkbox.attr('data-epitope');
            
            // Restore state if selected
            if (selectedEpitopes.bcell.indexOf(epitope) !== -1) {
                checkbox.prop('checked', true);
            }
            
            // Attach handler
            checkbox.on('change', function() {
                handleCheckboxChange('bcell', epitope, $(this).is(':checked'));
            });
        });
        
        // MHC-I checkboxes
        $('.mhc1-checkbox').each(function() {
            var checkbox = $(this);
            var epitope = checkbox.attr('data-epitope');
            
            // Restore state if selected
            if (selectedEpitopes.mhc1.indexOf(epitope) !== -1) {
                checkbox.prop('checked', true);
            }
            
            // Attach handler
            checkbox.on('change', function() {
                handleCheckboxChange('mhc1', epitope, $(this).is(':checked'));
            });
        });
        
        // MHC-II checkboxes
        $('.mhc2-checkbox').each(function() {
            var checkbox = $(this);
            var epitope = checkbox.attr('data-epitope');
            
            // Restore state if selected
            if (selectedEpitopes.mhc2.indexOf(epitope) !== -1) {
                checkbox.prop('checked', true);
            }
            
            // Attach handler
            checkbox.on('change', function() {
                handleCheckboxChange('mhc2', epitope, $(this).is(':checked'));
            });
        });
    }
    
    // Handle checkbox change
    function handleCheckboxChange(type, epitope, isChecked) {
        console.log('Checkbox changed:', type, epitope, isChecked);
        
        if (isChecked) {
            // Add to selected if not already there
            if (selectedEpitopes[type].indexOf(epitope) === -1) {
                selectedEpitopes[type].push(epitope);
                console.log('Added:', epitope, 'to', type);
            }
        } else {
            // Remove from selected
            var index = selectedEpitopes[type].indexOf(epitope);
            if (index !== -1) {
                selectedEpitopes[type].splice(index, 1);
                console.log('Removed:', epitope, 'from', type);
            }
        }
        
        updateSelectionDisplay();
    }
    
    // Update selection display
    function updateSelectionDisplay() {
        $('#bcellCount').text(selectedEpitopes.bcell.length);
        $('#mhc1Count').text(selectedEpitopes.mhc1.length);
        $('#mhc2Count').text(selectedEpitopes.mhc2.length);
        
        // Update proceed button state
        var totalSelected = selectedEpitopes.bcell.length + 
                          selectedEpitopes.mhc1.length + 
                          selectedEpitopes.mhc2.length;
        
        if (totalSelected > 0) {
            $('#proceedBtn').prop('disabled', false).addClass('pulse-animation');
        } else {
            $('#proceedBtn').prop('disabled', true).removeClass('pulse-animation');
        }
        
        console.log('Selection updated:', selectedEpitopes);
    }
    
    // Select all handlers
    $('#selectAllBcell').on('change', function() {
        var isChecked = $(this).is(':checked');
        bcellTable.$('.bcell-checkbox', {"page": "current"}).each(function() {
            var epitope = $(this).attr('data-epitope');
            $(this).prop('checked', isChecked);
            handleCheckboxChange('bcell', epitope, isChecked);
        });
    });
    
    $('#selectAllMhc1').on('change', function() {
        var isChecked = $(this).is(':checked');
        mhc1Table.$('.mhc1-checkbox', {"page": "current"}).each(function() {
            var epitope = $(this).attr('data-epitope');
            $(this).prop('checked', isChecked);
            handleCheckboxChange('mhc1', epitope, isChecked);
        });
    });
    
    $('#selectAllMhc2').on('change', function() {
        var isChecked = $(this).is(':checked');
        mhc2Table.$('.mhc2-checkbox', {"page": "current"}).each(function() {
            var epitope = $(this).attr('data-epitope');
            $(this).prop('checked', isChecked);
            handleCheckboxChange('mhc2', epitope, isChecked);
        });
    });
    
    // Proceed to construct page
    window.proceedToConstruct = function() {
        var totalSelected = selectedEpitopes.bcell.length + 
                          selectedEpitopes.mhc1.length + 
                          selectedEpitopes.mhc2.length;
        
        if (totalSelected === 0) {
            alert('Please select at least one epitope before proceeding');
            return;
        }
        
        console.log('Proceeding with selections:', selectedEpitopes);
        
        // Save selections to session before screening
        $.ajax({
            url: '/save_selected_epitopes',
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify(selectedEpitopes),
            success: function(response) {
                if (response.success) {
                    // Redirect to screening page
                    window.location.href = '/screening';
                } else {
                    alert('Error saving selections: ' + response.error);
                }
            },
            error: function(xhr, status, error) {
                alert('Error: ' + error);
            }
        });
    };
    
    // Download button animation
    $('.btn-download').on('click', function() {
        $(this).find('i').addClass('fa-spin');
        setTimeout(function() {
            $(this).find('i').removeClass('fa-spin');
        }, 1000);
    });
    
    // Initial display update
    updateSelectionDisplay();
});

// Add CSS for animations
var style = document.createElement('style');
style.textContent = `
    .pulse-animation {
        animation: pulse-glow 1.5s ease-in-out infinite;
    }
    
    @keyframes pulse-glow {
        0%, 100% {
            box-shadow: 0 0 10px rgba(0, 255, 255, 0.5);
        }
        50% {
            box-shadow: 0 0 20px rgba(0, 255, 255, 1);
        }
    }
`;
document.head.appendChild(style);
