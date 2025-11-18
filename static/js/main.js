// version1.0
// ========================================
// MAIN.JS - Index Page JavaScript
// ========================================

// Workflow selection functions
function showPredictionWorkflow() {
    $('.workflow-selection').fadeOut(300, function() {
        $('#predictionFormCard').fadeIn(300);
    });
}

function showEpitopeWorkflow() {
    $('.workflow-selection').fadeOut(300, function() {
        $('#epitopeFormCard').fadeIn(300);
    });
}

function resetWorkflow() {
    $('#predictionFormCard, #epitopeFormCard').fadeOut(300, function() {
        $('.workflow-selection').fadeIn(300);
    });
}

$(document).ready(function() {
    // Prediction form submission handler
    $('#predictionForm').on('submit', function(e) {
        e.preventDefault();
        
        // Validate form
        if (!validateForm()) {
            return;
        }
        
        // Get form data
        const formData = new FormData(this);
        
        // Show loading overlay
        showLoading();
        
        // Start polling for status
        startStatusPolling();
        
        // Submit form via AJAX
        $.ajax({
            url: '/predict',
            type: 'POST',
            data: formData,
            processData: false,
            contentType: false,
            success: function(response) {
                if (response.success) {
                    // Success - redirect to results page
                    updateLoadingText('All predictions complete! Redirecting...');
                    $('#progressBar').css('width', '100%');
                    setTimeout(function() {
                        window.location.href = '/results';
                    }, 1000);
                } else {
                    // Error
                    hideLoading();
                    stopStatusPolling();
                    showError(response.error || 'An error occurred during prediction');
                }
            },
            error: function(xhr, status, error) {
                hideLoading();
                stopStatusPolling();
                showError('Server error: ' + error);
            }
        });
    });

    // Epitope form submission handler
    $('#epitopeForm').on('submit', function(e) {
        e.preventDefault();
        
        // Get epitope data
        const bcellEpitopes = $('#bcell_epitopes').val().trim().split('\n').filter(e => e.trim());
        const mhc1Epitopes = $('#mhc1_epitopes').val().trim().split('\n').filter(e => e.trim());
        const mhc2Epitopes = $('#mhc2_epitopes').val().trim().split('\n').filter(e => e.trim());
        
        // Validate at least one epitope
        if (bcellEpitopes.length === 0 && mhc1Epitopes.length === 0 && mhc2Epitopes.length === 0) {
            showError('Please enter at least one epitope');
            return;
        }
        
        // Save epitopes to session and redirect to screening
        $.ajax({
            url: '/save_selected_epitopes',
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({
                bcell: bcellEpitopes,
                mhc1: mhc1Epitopes,
                mhc2: mhc2Epitopes
            }),
            success: function(response) {
                if (response.success) {
                    // Redirect to screening page
                    window.location.href = '/screening';
                } else {
                    showError('Error saving epitopes: ' + response.error);
                }
            },
            error: function(xhr, status, error) {
                showError('Server error: ' + error);
            }
        });
    });
    
    // Form validation
    function validateForm() {
        const fasta = $('#fasta_sequence').val().trim();
        const mhc1 = $('#mhc1_alleles').val().trim();
        const mhc2 = $('#mhc2_alleles').val().trim();
        
        if (!fasta) {
            showError('Please enter a protein sequence in FASTA format');
            return false;
        }
        
        if (!fasta.startsWith('>')) {
            showError('FASTA sequence must start with ">" followed by a header');
            return false;
        }
        
        if (!mhc1) {
            showError('Please enter MHC-I alleles');
            return false;
        }
        
        if (!mhc2) {
            showError('Please enter MHC-II alleles');
            return false;
        }
        
        return true;
    }
    
    // Show loading overlay
    function showLoading() {
        $('#loadingOverlay').addClass('show');
        $('#progressBar').css('width', '0%');
    }
    
    // Hide loading overlay
    function hideLoading() {
        $('#loadingOverlay').removeClass('show');
    }
    
    // Update loading text
    function updateLoadingText(text) {
        $('#loadingText').text(text);
    }
    
    // Status polling variables
    let statusPollingInterval = null;
    
    // Start polling for prediction status
    function startStatusPolling() {
        updateLoadingText('Initializing predictions...');
        $('#progressBar').css('width', '0%');
        
        statusPollingInterval = setInterval(function() {
            $.ajax({
                url: '/check_prediction_status',
                type: 'GET',
                success: function(status) {
                    updateProgressFromStatus(status);
                },
                error: function() {
                    // Continue polling even on error
                }
            });
        }, 3000); // Poll every 3 seconds
    }
    
    // Stop polling
    function stopStatusPolling() {
        if (statusPollingInterval) {
            clearInterval(statusPollingInterval);
            statusPollingInterval = null;
        }
    }
    
    // Update progress based on actual file status
    function updateProgressFromStatus(status) {
        console.log('Status update:', status);
        let progress = 0;
        let message = 'Initializing predictions...';
        
        if (!status.bcell_done) {
            progress = 10;
            message = 'Running B-cell epitope prediction (NetBCE)...';
        } else if (!status.mhc1_done) {
            progress = 40;
            message = 'Running MHC-I epitope prediction (NetMHCPan)...';
        } else if (!status.mhc2_done) {
            progress = 70;
            message = 'Running MHC-II epitope prediction (NetMHCIIPan)...';
        } else {
            progress = 95;
            message = 'Finalizing results...';
            stopStatusPolling();
        }
        
        console.log('Progress:', progress, 'Message:', message);
        updateLoadingText(message);
        $('#progressBar').css('width', progress + '%');
    }
    
    // Simulate progress bar
    function simulateProgress() {
        let progress = 0;
        const steps = [
            { progress: 20, text: 'Running NetBCE prediction...' },
            { progress: 50, text: 'Running MHC-I prediction...' },
            { progress: 80, text: 'Running MHC-II prediction...' },
            { progress: 100, text: 'Finalizing results...' }
        ];
        
        let currentStep = 0;
        
        const interval = setInterval(function() {
            if (currentStep < steps.length) {
                progress = steps[currentStep].progress;
                updateLoadingText(steps[currentStep].text);
                $('#progressBar').css('width', progress + '%');
                currentStep++;
            } else {
                clearInterval(interval);
            }
        }, 3000);
    }
    
    // Show error message
    function showError(message) {
        alert('Error: ' + message);
    }
    
    // Add animation to cards on scroll
    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.style.opacity = '1';
                entry.target.style.transform = 'translateY(0)';
            }
        });
    }, {
        threshold: 0.1
    });
    
    // Observe all animated elements
    document.querySelectorAll('.fade-in, .slide-up').forEach(el => {
        observer.observe(el);
    });
    
    // Button hover effect
    $('.btn-construct').on('mouseenter', function() {
        $(this).find('i').addClass('fa-spin');
    }).on('mouseleave', function() {
        $(this).find('i').removeClass('fa-spin');
    });
    
    // Auto-resize textareas
    $('textarea').on('input', function() {
        this.style.height = 'auto';
        this.style.height = (this.scrollHeight) + 'px';
    });
});
