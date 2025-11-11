// ========================================
// CONSTRUCT.JS - Construct Page JavaScript
// ========================================

$(document).ready(function() {
    // Retrieve selected epitopes from sessionStorage
    const selectedEpitopes = JSON.parse(sessionStorage.getItem('selectedEpitopes')) || {
        bcell: [],
        mhc1: [],
        mhc2: []
    };
    
    // Display selection summary
    displaySelectionSummary();
    
    // Display selected epitopes
    displaySelectedEpitopes();
    
    // Adjuvant card selection handler
    // Adjuvant card selection handler
    $('.adjuvant-card').on('click', function() {
        $('.adjuvant-card').removeClass('selected');
        $(this).addClass('selected');
        $(this).find('input[type="radio"]').prop('checked', true);
    });
    
    // Design card selection handler
    $('.design-card').on('click', function() {
        $('.design-card').removeClass('selected');
        $(this).addClass('selected');
        $(this).find('input[type="radio"]').prop('checked', true);
    });
    
    // Build construct button handler
    $('#buildConstructBtn').on('click', function() {
        buildConstruct();
    });
    
    // Copy button handlers
    $('#copyBtn, #copySequenceBtn').on('click', function() {
        copyToClipboard();
    });
    
    // ProtParam button handler
    $('#protparamBtn').on('click', function() {
        analyzeProtParam();
    });
    
    // 3D Structure prediction button handler
    $('#predict3DBtn').on('click', function() {
        predict3DStructure();
    });
    
    // ProDy NMA button handler
    $('#runProDyNMABtn').on('click', function() {
        window.open('/prody_nma', '_blank');
    });
    
    // Display selection summary
    function displaySelectionSummary() {
        $('#bcellCountDisplay').text(selectedEpitopes.bcell.length);
        $('#mhc1CountDisplay').text(selectedEpitopes.mhc1.length);
        $('#mhc2CountDisplay').text(selectedEpitopes.mhc2.length);
        
        // Animate counts
        animateCounts();
    }
    
    // Display selected epitopes
    function displaySelectedEpitopes() {
        // B-cell epitopes
        const bcellContainer = $('#bcellEpitopesList');
        bcellContainer.empty();
        if (selectedEpitopes.bcell.length === 0) {
            bcellContainer.html('<span class="text-muted">No B-cell epitopes selected</span>');
        } else {
            selectedEpitopes.bcell.forEach(epitope => {
                bcellContainer.append(`<span class="epitope-badge">${epitope}</span>`);
            });
        }
        
        // MHC-I epitopes
        const mhc1Container = $('#mhc1EpitopesList');
        mhc1Container.empty();
        if (selectedEpitopes.mhc1.length === 0) {
            mhc1Container.html('<span class="text-muted">No MHC-I epitopes selected</span>');
        } else {
            selectedEpitopes.mhc1.forEach(epitope => {
                mhc1Container.append(`<span class="epitope-badge">${epitope}</span>`);
            });
        }
        
        // MHC-II epitopes
        const mhc2Container = $('#mhc2EpitopesList');
        mhc2Container.empty();
        if (selectedEpitopes.mhc2.length === 0) {
            mhc2Container.html('<span class="text-muted">No MHC-II epitopes selected</span>');
        } else {
            selectedEpitopes.mhc2.forEach(epitope => {
                mhc2Container.append(`<span class="epitope-badge">${epitope}</span>`);
            });
        }
    }
    
    // Animate counts
    function animateCounts() {
        $('.epitope-summary-card h4').each(function() {
            const finalValue = parseInt($(this).text());
            const element = $(this);
            let currentValue = 0;
            
            const interval = setInterval(function() {
                if (currentValue < finalValue) {
                    currentValue++;
                    element.text(currentValue);
                } else {
                    clearInterval(interval);
                }
            }, 50);
        });
    }
    
    // Build construct
    function buildConstruct() {
        // Get selected adjuvant
        const selectedAdjuvant = $('input[name="adjuvant"]:checked').val();
        
        if (!selectedAdjuvant) {
            alert('Please select an adjuvant before building the construct');
            return;
        }
        
        // Get selected construct design
        const selectedDesign = $('input[name="construct_design"]:checked').val();
        
        if (!selectedDesign) {
            alert('Please select a construct design');
            return;
        }
        
        // Check if at least one epitope is selected
        const totalEpitopes = selectedEpitopes.bcell.length + 
                            selectedEpitopes.mhc1.length + 
                            selectedEpitopes.mhc2.length;
        
        if (totalEpitopes === 0) {
            alert('No epitopes selected. Please go back and select epitopes.');
            return;
        }
        
        // Show loading animation
        const button = $('#buildConstructBtn');
        button.prop('disabled', true);
        button.html('<i class="fas fa-spinner fa-spin"></i> Building Construct...');
        
        // Prepare data
        const data = {
            bcell_epitopes: selectedEpitopes.bcell,
            mhc1_epitopes: selectedEpitopes.mhc1,
            mhc2_epitopes: selectedEpitopes.mhc2,
            adjuvant: selectedAdjuvant,
            design: selectedDesign
        };
        
        // Send request to backend
        $.ajax({
            url: '/construct',
            type: 'POST',
            contentType: 'application/json',
            data: JSON.stringify(data),
            success: function(response) {
                if (response.success) {
                    displayConstruct(response);
                } else {
                    alert('Error: ' + response.error);
                }
                
                // Reset button
                button.prop('disabled', false);
                button.html('<i class="fas fa-hammer"></i> Build Vaccine Construct');
            },
            error: function(xhr, status, error) {
                alert('Server error: ' + error);
                
                // Reset button
                button.prop('disabled', false);
                button.html('<i class="fas fa-hammer"></i> Build Vaccine Construct');
            }
        });
    }
    
    // Display construct
    function displayConstruct(response) {
        // Show result section
        $('#constructResult').slideDown(500);
        
        // Update length
        $('#constructLength').text(response.length);
        
        // Display clean sequence
        $('#sequenceContent').text(response.clean_sequence);
        
        // Display schematic
        displaySchematic(response.construct);
        
        // Scroll to result
        $('html, body').animate({
            scrollTop: $('#constructResult').offset().top - 100
        }, 500);
    }
    
    // Display schematic with colored parts
    function displaySchematic(constructString) {
        const schematicContainer = $('#schematicDisplay');
        schematicContainer.empty();
        
        // Split by hyphens and color-code
        const parts = constructString.split('-');
        let html = '';
        
        parts.forEach((part, index) => {
            let className = '';
            
            // Determine part type
            if (index === 0) {
                className = 'adjuvant';
            } else if (part === 'EAAAK' || part === 'AAY' || part === 'GPGPG' || part === 'KK') {
                className = 'linker';
            } else if (selectedEpitopes.bcell.includes(part)) {
                className = 'bcell';
            } else if (selectedEpitopes.mhc1.includes(part)) {
                className = 'mhc1';
            } else if (selectedEpitopes.mhc2.includes(part)) {
                className = 'mhc2';
            }
            
            html += `<span class="construct-part ${className}">${part}</span>`;
            
            if (index < parts.length - 1) {
                html += '<span class="construct-part linker">-</span>';
            }
        });
        
        schematicContainer.html(html);
        
        // Add legend
        const legend = `
            <div class="mt-3" style="font-size: 0.9rem;">
                <strong>Legend:</strong>
                <span class="construct-part adjuvant">Adjuvant</span>
                <span class="construct-part linker">Linker</span>
                <span class="construct-part bcell">B-Cell</span>
                <span class="construct-part mhc1">MHC-I</span>
                <span class="construct-part mhc2">MHC-II</span>
            </div>
        `;
        schematicContainer.append(legend);
    }
    
    // Copy to clipboard
    function copyToClipboard() {
        const sequence = $('#sequenceContent').text();
        
        if (!sequence) {
            alert('No sequence to copy');
            return;
        }
        
        // Copy to clipboard
        navigator.clipboard.writeText(sequence).then(function() {
            // Show success toast
            showToast();
        }).catch(function(err) {
            // Fallback method
            const textarea = document.createElement('textarea');
            textarea.value = sequence;
            document.body.appendChild(textarea);
            textarea.select();
            document.execCommand('copy');
            document.body.removeChild(textarea);
            
            showToast();
        });
    }
    
    // Show toast notification
    function showToast() {
        const toast = $('#copyToast');
        toast.addClass('show');
        
        setTimeout(function() {
            toast.removeClass('show');
        }, 3000);
    }
    
    // Analyze with ProtParam
    function analyzeProtParam() {
        const sequence = $('#sequenceContent').text().replace(/-/g, '');
        
        if (!sequence) {
            alert('No sequence found to analyze');
            return;
        }
        
        // Show loading
        $('#protparamBtn').html('<i class="fas fa-spinner fa-spin"></i> Analyzing...');
        $('#protparamBtn').prop('disabled', true);
        
        $.ajax({
            url: '/protparam',
            type: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({ sequence: sequence }),
            success: function(response) {
                if (response.success) {
                    displayProtParamResults(response);
                    $('#protparamResult').slideDown();
                    
                    // Scroll to results
                    $('html, body').animate({
                        scrollTop: $('#protparamResult').offset().top - 100
                    }, 1000);
                } else {
                    alert('Error: ' + response.error);
                }
                
                // Reset button
                $('#protparamBtn').html('<i class="fas fa-chart-line"></i> Analyze with ProtParam');
                $('#protparamBtn').prop('disabled', false);
            },
            error: function(xhr, status, error) {
                alert('Analysis failed: ' + error);
                $('#protparamBtn').html('<i class="fas fa-chart-line"></i> Analyze with ProtParam');
                $('#protparamBtn').prop('disabled', false);
            }
        });
    }
    
    // Display ProtParam results
    function displayProtParamResults(data) {
        // Basic properties
        $('#param_mw').text(data.molecular_weight.toFixed(2) + ' Da');
        $('#param_length').text(data.sequence_length + ' aa');
        $('#param_pi').text(data.isoelectric_point.toFixed(2));
        
        // Physicochemical properties
        $('#param_arom').text(data.aromaticity.toFixed(4));
        $('#param_gravy').text(data.gravy.toFixed(4));
        $('#param_instab').text(data.instability_index.toFixed(2));
        $('#param_stability').text(data.stability);
        
        // Color code stability
        if (data.stability === 'Stable') {
            $('#param_stability').css('color', '#00ff88');
        } else {
            $('#param_stability').css('color', '#ff3366');
        }
        
        // Secondary structure
        $('#param_helix').text((data.helix_fraction * 100).toFixed(2) + '%');
        $('#param_turn').text((data.turn_fraction * 100).toFixed(2) + '%');
        $('#param_sheet').text((data.sheet_fraction * 100).toFixed(2) + '%');
        
        // Extinction coefficients
        $('#param_ext_red').text(data.extinction_coefficient_reduced + ' M⁻¹ cm⁻¹');
        $('#param_ext_ox').text(data.extinction_coefficient_oxidized + ' M⁻¹ cm⁻¹');
        
        // Amino acid composition
        const aaBody = $('#aaCompositionBody');
        aaBody.empty();
        
        // Sort amino acids by count (descending)
        const sortedAA = Object.entries(data.amino_acid_composition)
            .sort((a, b) => b[1] - a[1]);
        
        sortedAA.forEach(([aa, percent]) => {
            const count = data.amino_acid_counts[aa] || 0;
            if (count > 0) {
                aaBody.append(`
                    <tr>
                        <td><strong>${aa}</strong></td>
                        <td>${count}</td>
                        <td>${(percent * 100).toFixed(2)}%</td>
                    </tr>
                `);
            }
        });
    }
    
    // Add hover effects to adjuvant cards
    $('.adjuvant-card').hover(
        function() {
            $(this).find('i').addClass('fa-bounce');
        },
        function() {
            $(this).find('i').removeClass('fa-bounce');
        }
    );
    
    // Add animation to epitope badges
    $('.epitope-badge').each(function(index) {
        $(this).css({
            'animation': 'fade-in 0.5s ease',
            'animation-delay': (index * 0.05) + 's'
        });
    });
    
    // 3D Structure Prediction Function
    function predict3DStructure() {
        const sequence = $('#sequenceContent').text().replace(/-/g, '');
        
        if (!sequence) {
            alert('No construct sequence available! Please build the construct first.');
            return;
        }
        
        // Show loading state
        $('#predict3DBtn').prop('disabled', true).html('<i class="fas fa-spinner fa-spin"></i> Predicting...');
        
        // Hide structure result initially
        $('#structure3DResult').hide();
        
        $.ajax({
            url: '/predict_structure',
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({
                sequence: sequence
            }),
            success: function(response) {
                if (response.success) {
                    // Show structure viewer first
                    $('#structure3DResult').show();
                    
                    // Scroll to structure immediately
                    $('html, body').animate({
                        scrollTop: $('#structure3DResult').offset().top - 100
                    }, 300);
                    
                    // Render 3D structure using 3Dmol.js (after a short delay to ensure DOM is ready)
                    setTimeout(function() {
                        render3DStructure(response.pdb_content);
                    }, 100);
                } else {
                    alert('Error: ' + response.error);
                }
            },
            error: function(xhr, status, error) {
                alert('Structure prediction failed: ' + error);
            },
            complete: function() {
                $('#predict3DBtn').prop('disabled', false).html('<i class="fas fa-cube"></i> Predict 3D Structure');
            }
        });
    }
    
    // Render 3D Structure using 3Dmol.js
    function render3DStructure(pdbContent) {
        // Clear previous viewer
        $('#viewer3d').empty();
        
        // Create 3Dmol viewer with optimized settings
        let element = $('#viewer3d');
        let config = { 
            backgroundColor: '#1a0028',
            antialias: false,  // Disable antialiasing for better performance
            disableFog: true
        };
        let viewer = $3Dmol.createViewer(element, config);
        
        // Add model
        viewer.addModel(pdbContent, 'pdb');
        
        // Set simple cartoon style only (no surface for performance)
        viewer.setStyle({}, {
            cartoon: {
                color: 'spectrum',
                thickness: 0.6
            }
        });
        
        // Center and zoom
        viewer.zoomTo();
        viewer.zoom(0.8);
        
        // Render once
        viewer.render();
        
        // NO auto-rotation by default (user can rotate manually)
        // Uncomment next line if you want slow auto-rotation:
        // viewer.spin('y', 0.3);
        
        // Add controls info
        setTimeout(function() {
            if (!$('.viewer-controls').length) {
                element.append(`
                    <div class="viewer-controls" style="position: absolute; top: 10px; right: 10px; background: rgba(0,0,0,0.85); padding: 12px; border-radius: 8px; color: white; font-size: 11px; border: 1px solid #00ffff; backdrop-filter: blur(5px);">
                        <div style="margin-bottom: 8px; color: #00ffff; font-weight: bold; font-size: 12px;">
                            <i class="fas fa-hand-pointer"></i> Controls
                        </div>
                        <div style="margin: 4px 0;"><i class="fas fa-mouse" style="color: #00ffff; width: 16px;"></i> Left Click: Rotate</div>
                        <div style="margin: 4px 0;"><i class="fas fa-mouse" style="color: #00ffff; width: 16px;"></i> Scroll: Zoom</div>
                        <div style="margin: 4px 0;"><i class="fas fa-mouse" style="color: #00ffff; width: 16px;"></i> Right Click: Pan</div>
                    </div>
                `);
            }
        }, 100);
    }
    
    // Redirect back if no epitopes selected
    if (selectedEpitopes.bcell.length === 0 && 
        selectedEpitopes.mhc1.length === 0 && 
        selectedEpitopes.mhc2.length === 0) {
        alert('No epitopes selected. Redirecting to results page...');
        window.location.href = '/results';
    }
});

// Add custom CSS for selected adjuvant card
const style = document.createElement('style');
style.textContent = `
    .adjuvant-card.selected {
        border-color: var(--accent-cyan) !important;
        background: rgba(0, 255, 255, 0.15) !important;
        box-shadow: 0 0 30px rgba(0, 255, 255, 0.5);
    }
    
    .adjuvant-card.selected::after {
        content: '✓';
        position: absolute;
        top: 10px;
        right: 10px;
        background: var(--accent-cyan);
        color: var(--primary-dark);
        width: 30px;
        height: 30px;
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: bold;
        font-size: 1.2rem;
    }
`;
document.head.appendChild(style);
