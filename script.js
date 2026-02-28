// ==========================================
// UPDATE QUICK SUMMARY (non-floating)
// ==========================================
function updateQuickSummary(results) {
    const quickSummary = document.getElementById('quickSummary');
    const quickSummaryContent = document.getElementById('quickSummaryContent');
    
    const systemType = document.getElementById('systemType').value;
    const isBritish = systemType.includes('ethanol');
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
    
    quickSummaryContent.innerHTML = `
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>Stages:</span>
                <strong>${results.stages}</strong>
            </div>
        </div>
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>Feed Stage:</span>
                <strong>${results.feed_stage}</strong>
            </div>
        </div>
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>R min:</span>
                <strong>${results.RMin}</strong>
            </div>
        </div>
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>Qc:</span>
                <strong>${results.QcKW} kW</strong>
            </div>
        </div>
    `;
    
    quickSummary.style.display = 'block';
}

// ==========================================
// UPDATE ENTHALPY TABLE
// ==========================================
function updateEnthalpyTable(results) {
    const systemType = document.getElementById('systemType').value;
    const isBritish = systemType.includes('ethanol');
    
    let html = '';
    if (results.stage_compositions && results.stage_compositions.length > 0) {
        results.stage_compositions.forEach((stage, i) => {
            const idxLiq = Math.round(stage.x * 199);
            const idxVap = Math.round(stage.y * 199);
            const HL = results.HL_curve && results.HL_curve[idxLiq] ? 
                     results.HL_curve[idxLiq].toFixed(0) : 'N/A';
            const HV = results.HV_curve && results.HV_curve[idxVap] ? 
                     results.HV_curve[idxVap].toFixed(0) : 'N/A';
            
            html += `<tr><td>Stage ${i+1}</td><td>${HL}</td><td>${HV}</td></tr>`;
        });
    } else {
        html = '<tr><td colspan="3" class="text-center text-muted">No enthalpy data available</td></tr>';
    }
    document.getElementById('enthalpyBody').innerHTML = html;
}

// ==========================================
// MODIFY DISPLAY RESULTS
// ==========================================
function displayResults(results) {
    createPlot(results);
    updateQuickSummary(results);
    updateEnthalpyTable(results);
    
    const systemType = document.getElementById('systemType').value;
    const isBritish = systemType.includes('ethanol');
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
    const powerUnit = 'kW';
    
    const summaryHtml = `
        <tr><td>Distillate Flow Rate (D)</td><td>${results.D} kmol/hr</td></tr>
        <tr><td>Bottoms Flow Rate (W)</td><td>${results.W} kmol/hr</td></tr>
        <tr><td>Δ_R</td><td>(${results.xDeltaR.toFixed(3)}, ${results.HDeltaR.toFixed(2)} ${enthalpyUnit})</td></tr>
        <tr><td>Δ_S</td><td>(${results.xDeltaS.toFixed(3)}, ${results.HDeltaS.toFixed(2)} ${enthalpyUnit})</td></tr>
        <tr><td>Condenser Duty (Qc)</td><td>${results.QcKW} ${powerUnit}</td></tr>
        <tr><td>Reboiler Duty (Qr)</td><td>${results.QrKW} ${powerUnit}</td></tr>
        <tr><td>Δ_R min</td><td>(${results.xD.toFixed(3)}, ${results.QPrimeMin.toFixed(2)} ${enthalpyUnit})</td></tr>
        <tr><td>Δ_S min</td><td>(${results.xB.toFixed(3)}, ${results.QDoublePrimeMin.toFixed(2)} ${enthalpyUnit})</td></tr>
        <tr><td>Minimum Reflux Ratio</td><td>${results.RMin}</td></tr>
        <tr><td>Number of Stages</td><td>${results.stages}</td></tr>
        <tr><td>Feed Stage</td><td>${results.feed_stage}</td></tr>
    `;
    document.getElementById('summaryBody').innerHTML = summaryHtml;
    
    let stagesRows = '';
    if (results.stage_compositions && results.stage_compositions.length > 0) {
        results.stage_compositions.forEach((stage, i) => {
            stagesRows += `<tr><td>Stage ${i+1}</td><td>${stage.x.toFixed(4)}</td><td>${stage.y.toFixed(4)}</td></tr>`;
        });
    } else {
        stagesRows = '<tr><td colspan="3" class="text-center text-muted">No stage data available</td></tr>';
    }
    document.getElementById('stagesBody').innerHTML = stagesRows;
    
    document.getElementById('exportBtn').disabled = false;
    currentResults = results;
    
    showToast('✅ Calculation completed successfully!', 'success');
}

// ==========================================
// MODIFY LOAD ETHANOL-WATER DATASET
// ==========================================
function loadEthanolWaterDataset(pressure) {
    const dataset = ETHANOL_WATER_DATA[pressure];
    if (!dataset) return;
    
    console.log(`Loading ${dataset.name} dataset...`);
    
    // Update input fields with British units data
    document.getElementById('xData').value = JSON.stringify(dataset.x);
    document.getElementById('yData').value = JSON.stringify(dataset.y);
    document.getElementById('Hl').value = JSON.stringify(dataset.Hl);
    document.getElementById('Hv').value = JSON.stringify(dataset.Hv);
    
    // Store temperature data for reference
    document.getElementById('TData').value = JSON.stringify(dataset.T);
    
    // Update dataset info panel
    document.getElementById('datasetName').textContent = dataset.name;
    
    let azeoInfo = '';
    if (dataset.azeotrope) {
        azeoInfo = ` | Azeotrope: x = ${dataset.azeotrope.x.toFixed(3)} at ${dataset.azeotrope.T.toFixed(1)}°F`;
    }
    
    document.getElementById('datasetDetails').innerHTML = `
        <i class="fas fa-database me-1"></i> ${dataset.x.length} points
        <i class="fas fa-fire ms-2 me-1"></i> Hl: ${dataset.Hl[0]}-${dataset.Hl[dataset.Hl.length-1]} BTU/lbmole
        <i class="fas fa-temperature-high ms-2 me-1"></i> T: ${dataset.T[0].toFixed(1)}-${dataset.T[dataset.T.length-1].toFixed(1)}°F
        ${azeoInfo}
    `;
    
    // Update system info display
    const systemInfo = document.getElementById('systemInfo');
    const systemInfoText = document.getElementById('systemInfoText');
    systemInfo.style.display = 'block';
    
    systemInfoText.innerHTML = `
        <b>${dataset.name}</b><br>
        📊 ${dataset.x.length} data points<br>
        🔥 Enthalpy: <b>${dataset.units.enthalpy}</b><br>
        🌡️ Temperature: <b>${dataset.units.temperature}</b>
    `;
    
    // Update unit display
    updateUnitDisplay(`ethanol-${pressure}`);
    
    // Show dataset info panel
    document.getElementById('datasetInfoPanel').style.display = 'block';
    
    // Update preview
    updatePreview();
    
    // Show success message
    showToast(`✅ Loaded ${dataset.name} (British Units)`, 'success');
}

// ==========================================
// MODIFY SYSTEM TYPE CHANGE HANDLER
// ==========================================
document.getElementById('systemType').addEventListener('change', function() {
    const systemType = this.value;
    currentSystem = systemType;
    
    console.log(`System changed to: ${systemType}`);
    
    // Hide quick summary when switching systems
    document.getElementById('quickSummary').style.display = 'none';
    
    // Reset system info display
    const systemInfo = document.getElementById('systemInfo');
    const datasetInfoPanel = document.getElementById('datasetInfoPanel');
    
    if (systemType === 'user-defined') {
        systemInfo.style.display = 'none';
        datasetInfoPanel.style.display = 'none';
        
        // Restore example data
        document.getElementById('xData').value = '[0, 0.08, 0.18, 0.25, 0.49, 0.65, 0.79, 0.91, 1.0]';
        document.getElementById('yData').value = '[0, 0.28, 0.43, 0.51, 0.73, 0.83, 0.90, 0.96, 1.0]';
        document.getElementById('Hl').value = '[24.3, 24.1, 23.2, 22.8, 22.05, 21.75, 21.7, 21.6, 21.4]';
        document.getElementById('Hv').value = '[61.2, 59.6, 58.5, 58.1, 56.5, 55.2, 54.4, 53.8, 53.3]';
        
        updateUnitDisplay('user-defined');
        showToast('📝 Switched to User Defined mode (any units)', 'info');
    } 
    else if (systemType === 'ethanol-water-1atm') {
        loadEthanolWaterDataset('1atm');
    }
    else if (systemType === 'ethanol-water-76mmHg') {
        loadEthanolWaterDataset('76mmHg');
    }
    
    updatePreview();
});

// ==========================================
// ADD AUTO-SCROLL TO RESULTS
// ==========================================
function scrollToResults() {
    const resultsSection = document.querySelector('.results-section');
    if (resultsSection) {
        resultsSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
}

// Modify the displayResults function to include auto-scroll
const originalDisplayResults = displayResults;
displayResults = function(results) {
    originalDisplayResults(results);
    setTimeout(scrollToResults, 100); // Small delay to ensure plot is rendered
};
