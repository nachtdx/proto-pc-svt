// script.js
// ==========================================
// STATE & VARIABLES GLOBAL
// ==========================================
let pyodide = null;
let pyodideReady = false;
let currentResults = null;

// ==========================================
// PARSE ARRAY STRING
// ==========================================
function parseArrayString(str) {
    try {
        let cleanStr = str.replace(/[\[\]']/g, '').trim();
        if (cleanStr === '') return [];
        return cleanStr.split(',').map(x => {
            let num = parseFloat(x.trim());
            return isNaN(num) ? null : num;
        }).filter(x => x !== null);
    } catch (e) {
        return [];
    }
}

// ==========================================
// VALIDATE ARRAY
// ==========================================
function validateArray(arr, name) {
    if (!Array.isArray(arr) || arr.length === 0) {
        return {valid: false, error: `${name} tidak boleh kosong`};
    }
    if (arr.some(isNaN)) {
        return {valid: false, error: `${name} harus berupa angka`};
    }
    if (name.includes('x') || name.includes('y')) {
        if (arr.some(v => v < 0 || v > 1)) {
            return {valid: false, error: `${name} harus antara 0 dan 1`};
        }
    }
    return {valid: true};
}

// ==========================================
// UPDATE PREVIEW TABLE
// ==========================================
function updatePreview() {
    const xData = parseArrayString(document.getElementById('xData').value);
    const yData = parseArrayString(document.getElementById('yData').value);
    const Hl = parseArrayString(document.getElementById('Hl').value);
    const Hv = parseArrayString(document.getElementById('Hv').value);
    
    const lengths = [xData.length, yData.length, Hl.length, Hv.length];
    const maxLength = Math.max(...lengths);
    
    let html = '';
    for (let i = 0; i < maxLength; i++) {
        html += '<tr>';
        html += `<td>${i + 1}</td>`;
        html += `<td>${i < xData.length ? xData[i].toFixed(3) : '-'}</td>`;
        html += `<td>${i < yData.length ? yData[i].toFixed(3) : '-'}</td>`;
        html += `<td>${i < Hl.length ? Hl[i].toFixed(2) : '-'}</td>`;
        html += `<td>${i < Hv.length ? Hv[i].toFixed(2) : '-'}</td>`;
        html += '</tr>';
    }
    
    document.getElementById('previewBody').innerHTML = html;
    
    const warning = document.getElementById('previewWarning');
    if (new Set(lengths).size > 1) {
        if (!warning) {
            const div = document.createElement('div');
            div.id = 'previewWarning';
            div.className = 'alert alert-warning mt-2';
            div.innerHTML = '<i class="fas fa-exclamation-triangle"></i> Panjang array tidak sama!';
            document.querySelector('#previewTable').after(div);
        }
    } else if (warning) {
        warning.remove();
    }
}

// ==========================================
// LOAD CALCULATOR CODE
// ==========================================
async function loadCalculatorCode() {
    try {
        // Kode Python yang akan dijalankan di Pyodide
        const pythonCode = `
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import root_scalar
import json

# ==========================================
# INTERPOLATION FUNCTIONS
# ==========================================
def linear_interpolate(x, x_val, y_val):
    if not x_val or not y_val or len(x_val) != len(y_val):
        return float('nan')
    if x <= x_val[0]:
        return y_val[0]
    if x >= x_val[-1]:
        return y_val[-1]
    f = interp1d(x_val, y_val, kind='linear', fill_value='extrapolate')
    return float(f(x))

def cubic_interpolate(x, x_val, y_val):
    if not x_val or not y_val or len(x_val) != len(y_val):
        return float('nan')
    if x < x_val[0] or x > x_val[-1]:
        return float('nan')
    cs = CubicSpline(x_val, y_val, extrapolate=False)
    return float(cs(x))

# ==========================================
# CALCULATE STAGES
# ==========================================
def calculate_stages(xD, xB, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data):
    y = xD
    stages = 0
    stage_points = [{'x': xD, 'y': cubic_interpolate(xD, data['yData'], data['Hv'])}]
    tie_lines = []
    construction_lines = []
    stage_compositions = []
    in_rectifying = True
    feed_stage = 0
    error = ""

    while stages < 100:
        # Find x_n from y_n using equilibrium data
        def find_x(x):
            return cubic_interpolate(x, data['xData'], data['yData']) - y
        sol = root_scalar(find_x, bracket=[0, 1], method='bisect')
        if not sol.converged:
            error = f"Stage {stages + 1}: Invalid liquid composition."
            break
        x_n = sol.root

        if x_n <= xB:
            HLxB = linear_interpolate(xB, data['xData'], data['Hl'])
            if np.isfinite(HLxB):
                stage_points.append({'x': xB, 'y': HLxB})
                tie_lines.append({'x': [xB, y], 'y': [HLxB, cubic_interpolate(y, data['yData'], data['Hv'])]})
                stage_compositions.append({'x': xB, 'y': y})
                stages += 1
            break

        HLx_n = linear_interpolate(x_n, data['xData'], data['Hl'])
        if not np.isfinite(HLx_n):
            error = f"Stage {stages + 1}: Failed to interpolate liquid enthalpy."
            break
        stage_points.append({'x': x_n, 'y': HLx_n})

        HVy_n = cubic_interpolate(y, data['yData'], data['Hv'])
        if not np.isfinite(HVy_n):
            error = f"Stage {stages + 1}: Failed to interpolate vapor enthalpy."
            break
        tie_lines.append({'x': [x_n, y], 'y': [HLx_n, HVy_n]})
        stage_compositions.append({'x': x_n, 'y': y})
        stages += 1

        if in_rectifying and x_n <= zF:
            in_rectifying = False
            feed_stage = stages

        xDelta = xDeltaR if in_rectifying else xDeltaS
        HDelta = HDeltaR if in_rectifying else HDeltaS

        if abs(x_n - xDelta) < 1e-6:
            error = f"Stage {stages + 1}: Composition too close to difference point."
            break
        slope = (HDelta - HLx_n) / (xDelta - x_n)
        if not np.isfinite(slope):
            error = f"Stage {stages + 1}: Invalid slope calculation."
            break

        def find_y(y):
            return cubic_interpolate(y, data['yData'], data['Hv']) - (HDelta + slope * (y - xDelta))
        sol = root_scalar(find_y, bracket=[xB, xD], method='bisect')
        if not sol.converged:
            yMin = max(0, y - 0.2)
            yMax = min(1, y + 0.2)
            sol = root_scalar(find_y, bracket=[yMin, yMax], method='bisect')
            if not sol.converged:
                error = f"Stage {stages + 1}: Failed to find valid y_{stages + 2}."
                break
        yNext = sol.root
        HVyNext = cubic_interpolate(yNext, data['yData'], data['Hv'])
        if not np.isfinite(HVyNext):
            error = f"Stage {stages + 1}: Failed to interpolate vapor enthalpy."
            break
        stage_points.append({'x': yNext, 'y': HVyNext})

        if in_rectifying:
            def find_x_end(x):
                return linear_interpolate(x, data['xData'], data['Hl']) - (HDelta + slope * (x - xDelta))
            sol = root_scalar(find_x_end, bracket=[0, 1], method='bisect')
            if not sol.converged:
                error = f"Stage {stages + 1}: Failed to find liquid line intersection."
                break
            xEnd = sol.root
            HEnd = linear_interpolate(xEnd, data['xData'], data['Hl'])
        else:
            xEnd = yNext
            HEnd = HVyNext
        if not np.isfinite(HEnd):
            error = f"Stage {stages + 1}: Invalid construction line endpoint."
            break
        construction_lines.append({'x': [xDelta, xEnd], 'y': [HDelta, HEnd]})
        y = yNext

    if stages == 0 or len(stage_points) < 2:
        error = "Failed to calculate stages. Check input data."
    
    return {
        'stages': stages,
        'feed_stage': feed_stage,
        'tie_lines': tie_lines,
        'construction_lines': construction_lines,
        'stage_compositions': stage_compositions,
        'error': error
    }

# ==========================================
# MAIN CALCULATION
# ==========================================
def calculate(data, params):
    # Extract params
    zF = params['zF']
    F = params['F']
    xD = params['xD']
    xB = params['xB']
    q = params['q']
    R = params['R']
    
    # Calculate enthalpies
    yF = cubic_interpolate(zF, data['xData'], data['yData'])
    HLzF = linear_interpolate(zF, data['xData'], data['Hl'])
    HVzF = cubic_interpolate(yF, data['yData'], data['Hv'])
    HF = q * HLzF + (1 - q) * HVzF
    HD = linear_interpolate(xD, data['xData'], data['Hl'])
    HW = linear_interpolate(xB, data['xData'], data['Hl'])
    
    if not all(np.isfinite([yF, HLzF, HVzF, HF, HD, HW])):
        return {'error': "Interpolation failed. Ensure compositions are within valid range."}
    
    # Material balance
    D = F * (zF - xB) / (xD - xB)
    W = F - D
    
    # Condenser duty
    HVxD = cubic_interpolate(xD, data['yData'], data['Hv'])
    Qc = D * (HVxD - HD) * (R + 1)
    QcKW = Qc * 0.27778
    xDeltaR = xD
    HDeltaR = HD + Qc / D
    
    # Stripping difference point
    xDeltaS = xB
    slope = (HDeltaR - HF) / (xDeltaR - zF)
    HDeltaS = HF + slope * (xDeltaS - zF)
    
    # Reboiler duty
    Qr = W * (HW - HDeltaS)
    QrKW = Qr * 0.27778
    
    # Minimum reflux
    yFMin = cubic_interpolate(zF, data['xData'], data['yData'])
    HVyF = cubic_interpolate(yFMin, data['yData'], data['Hv'])
    slopeMin = (HVyF - HF) / (yFMin - zF)
    QPrimeMin = HF + slopeMin * (xD - zF)
    QDoublePrimeMin = HF + slopeMin * (xB - zF)
    RMin = (QPrimeMin - HVxD) / (HVxD - HD)
    
    # Calculate stages
    stage_results = calculate_stages(
        xD, xB, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data
    )
    
    if stage_results['error']:
        return {'error': stage_results['error']}
    
    # Prepare curves for plotting
    x_range = np.linspace(0, 1, 200).tolist()
    HL_curve = [linear_interpolate(xi, data['xData'], data['Hl']) for xi in x_range]
    HV_curve = [cubic_interpolate(xi, data['yData'], data['Hv']) for xi in x_range]
    y_equilibrium = [cubic_interpolate(xi, data['xData'], data['yData']) for xi in x_range]
    
    # Axis limits
    y_values = data['Hl'] + data['Hv'] + [HF, HD, HW, HDeltaR, HDeltaS, QPrimeMin, QDoublePrimeMin, HVyF]
    y_values = [y for y in y_values if np.isfinite(y)]
    yMin, yMax = min(y_values) - 10, max(y_values) + 10
    
    # Results
    results = {
        'D': round(D, 2),
        'W': round(W, 2),
        'xDeltaR': round(xDeltaR, 3),
        'HDeltaR': round(HDeltaR, 2),
        'xDeltaS': round(xDeltaS, 3),
        'HDeltaS': round(HDeltaS, 2),
        'QcKW': round(QcKW, 2),
        'QrKW': round(QrKW, 2),
        'QPrimeMin': round(QPrimeMin, 2),
        'QDoublePrimeMin': round(QDoublePrimeMin, 2),
        'RMin': round(RMin, 2),
        'stages': stage_results['stages'],
        'feed_stage': stage_results['feed_stage'],
        'stage_compositions': stage_results['stage_compositions'],
        'x_range': x_range,
        'HL_curve': HL_curve,
        'HV_curve': HV_curve,
        'y_equilibrium': y_equilibrium,
        'yMin': yMin,
        'yMax': yMax,
        'tie_lines': stage_results['tie_lines'],
        'construction_lines': stage_results['construction_lines'],
        'HF': HF,
        'zF': zF,
        'xD': xD,
        'xB': xB,
        'yFMin': yFMin,
        'HVyF': HVyF
    }
    
    return results

# Wrapper untuk JS
def calculate_from_js(xData, yData, Hl, Hv, zF, F, xD, xB, q, R):
    try:
        data = {
            'xData': [float(x) for x in xData],
            'yData': [float(y) for y in yData],
            'Hl': [float(h) for h in Hl],
            'Hv': [float(h) for h in Hv]
        }
        
        params = {
            'zF': float(zF),
            'F': float(F),
            'xD': float(xD),
            'xB': float(xB),
            'q': float(q),
            'R': float(R)
        }
        
        result = calculate(data, params)
        
        # Convert numpy arrays to lists
        if isinstance(result, dict):
            for key, value in result.items():
                if isinstance(value, np.ndarray):
                    result[key] = value.tolist()
                elif isinstance(value, list):
                    for i, item in enumerate(value):
                        if isinstance(item, dict):
                            for k, v in item.items():
                                if isinstance(v, np.ndarray):
                                    item[k] = v.tolist()
        
        return json.dumps(result)
        
    except Exception as e:
        return json.dumps({'error': str(e)})
`;
        
        // Jalankan kode Python di Pyodide
        pyodide.runPython(pythonCode);
        console.log('✅ Calculator code loaded!');
        
    } catch (error) {
        console.error('Error loading calculator code:', error);
        throw error;
    }
}

// ==========================================
// RUN CALCULATION
// ==========================================
async function runCalculation(inputData) {
    if (!pyodideReady) {
        throw new Error('Pyodide belum siap. Tunggu sebentar...');
    }
    
    try {
        // Panggil Python function
        const resultJson = pyodide.runPython(`
            calculate_from_js(
                ${JSON.stringify(inputData.xData)},
                ${JSON.stringify(inputData.yData)},
                ${JSON.stringify(inputData.Hl)},
                ${JSON.stringify(inputData.Hv)},
                ${inputData.zF},
                ${inputData.F},
                ${inputData.xD},
                ${inputData.xB},
                ${inputData.q},
                ${inputData.R}
            )
        `);
        
        return JSON.parse(resultJson);
        
    } catch (error) {
        console.error('Error running calculation:', error);
        throw error;
    }
}

// ==========================================
// CREATE PLOT
// ==========================================
function createPlot(results) {
    const stageColors = ['#FF6B6B', '#4ECDC4', '#FF9F1C', '#6A4C93', '#2E86AB', 
                         '#A23B72', '#F18F01', '#2D6A4F', '#9E2A2B', '#540D6E'];
    
    // Buat figure dengan subplots
    const traces = [];
    
    // ===== DIAGRAM H-x-y =====
    // Saturated liquid curve
    traces.push({
        x: results.x_range,
        y: results.HL_curve,
        mode: 'lines',
        name: 'Saturated Liquid',
        line: {color: '#2E86AB', width: 4},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // Saturated vapor curve
    traces.push({
        x: results.x_range,
        y: results.HV_curve,
        mode: 'lines',
        name: 'Saturated Vapor',
        line: {color: '#A23B72', width: 4},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // xD line
    traces.push({
        x: [results.xD, results.xD],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'x<sub>D</sub>',
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // xB line
    traces.push({
        x: [results.xB, results.xB],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'x<sub>B</sub>',
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // zF line
    traces.push({
        x: [results.zF, results.zF],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'z<sub>F</sub>',
        line: {color: '#2D6A4F', width: 3, dash: 'dash'},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // Difference points
    traces.push({
        x: [results.xDeltaR, results.xDeltaS],
        y: [results.HDeltaR, results.HDeltaS],
        mode: 'markers+text',
        name: 'Difference Points',
        marker: {color: '#F97316', size: 14, symbol: 'star', line: {color: 'white', width: 1}},
        text: ['Δ<sub>R</sub>', 'Δ<sub>S</sub>'],
        textposition: ['top center', 'bottom center'],
        textfont: {size: 14, color: '#F97316'},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // Operating line
    traces.push({
        x: [results.xDeltaR, results.zF, results.xDeltaS],
        y: [results.HDeltaR, results.HF, results.HDeltaS],
        mode: 'lines+markers',
        name: 'Operating Line',
        line: {color: '#0A9396', width: 3},
        marker: {size: 8, color: '#0A9396', line: {color: 'white', width: 1}},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // Minimum reflux line
    traces.push({
        x: [results.xB, results.zF, results.yFMin, results.xD],
        y: [results.QDoublePrimeMin, results.HF, results.HVyF, results.QPrimeMin],
        mode: 'lines+markers',
        name: 'Minimum Reflux Line',
        line: {color: '#E9C46A', width: 2.5, dash: 'dash'},
        marker: {size: 8, color: '#E9C46A', line: {color: 'white', width: 1}},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // Stage tie lines
    results.tie_lines.forEach((tie, i) => {
        const color = stageColors[i % stageColors.length];
        traces.push({
            x: tie.x,
            y: tie.y,
            mode: 'lines',
            name: `Stage ${i+1}`,
            line: {color: color, width: 3},
            xaxis: 'x',
            yaxis: 'y'
        });
    });
    
    // ===== DIAGRAM VLE =====
    // Equilibrium curve
    traces.push({
        x: results.x_range,
        y: results.y_equilibrium,
        mode: 'lines',
        name: 'Equilibrium Curve',
        line: {color: '#1E1E1E', width: 3.5},
        xaxis: 'x2',
        yaxis: 'y2'
    });
    
    // y = x line
    traces.push({
        x: [0, 1],
        y: [0, 1],
        mode: 'lines',
        name: 'y = x',
        line: {color: '#6C757D', width: 2, dash: 'dash'},
        xaxis: 'x2',
        yaxis: 'y2'
    });
    
    // Stage points in VLE
    results.stage_compositions.forEach((stage, i) => {
        const color = stageColors[i % stageColors.length];
        traces.push({
            x: [stage.x, stage.y],
            y: [stage.y, stage.y],
            mode: 'markers+lines',
            name: `Stage ${i+1} VLE`,
            line: {color: color, width: 2},
            marker: {color: color, size: 8, symbol: ['circle', 'diamond']},
            xaxis: 'x2',
            yaxis: 'y2',
            showlegend: false
        });
    });
    
    // Layout
    const layout = {
        title: {
            text: 'Ponchon–Savarit Diagram',
            font: {size: 20, family: 'Arial Black', color: '#1E1E1E'},
            x: 0.5
        },
        grid: {
            rows: 2,
            columns: 1,
            pattern: 'independent',
            roworder: 'top to bottom'
        },
        xaxis: {
            title: 'Mole Fraction (x or y)',
            range: [0, 1],
            tickformat: '.2f',
            gridcolor: '#E0E0E0'
        },
        yaxis: {
            title: 'Enthalpy (MJ/kmol)',
            range: [results.yMin, results.yMax],
            gridcolor: '#E0E0E0'
        },
        xaxis2: {
            title: 'Mole Fraction (x or y)',
            range: [0, 1],
            tickformat: '.2f',
            gridcolor: '#E0E0E0'
        },
        yaxis2: {
            title: 'y (Vapor Fraction)',
            range: [0, 1],
            tickformat: '.2f',
            gridcolor: '#E0E0E0'
        },
        height: 800,
        showlegend: true,
        legend: {
            x: 1.02,
            y: 1,
            font: {size: 10},
            bgcolor: 'rgba(255,255,255,0.9)',
            bordercolor: '#1E1E1E',
            borderwidth: 1
        },
        hovermode: 'closest'
    };
    
    Plotly.newPlot('plotDiv', traces, layout, {responsive: true});
}

// ==========================================
// DISPLAY RESULTS
// ==========================================
function displayResults(results) {
    // Create plot
    createPlot(results);
    
    // Summary table
    const summaryHtml = `
        <tr><td>Distillate Flow Rate (D)</td><td>${results.D} kmol/hr</td></tr>
        <tr><td>Bottoms Flow Rate (W)</td><td>${results.W} kmol/hr</td></tr>
        <tr><td>Δ_R</td><td>(${results.xDeltaR}, ${results.HDeltaR} MJ/kmol)</td></tr>
        <tr><td>Δ_S</td><td>(${results.xDeltaS}, ${results.HDeltaS} MJ/kmol)</td></tr>
        <tr><td>Condenser Duty (Qc)</td><td>${results.QcKW} kW</td></tr>
        <tr><td>Reboiler Duty (Qr)</td><td>${results.QrKW} kW</td></tr>
        <tr><td>Minimum Reflux Ratio</td><td>${results.RMin}</td></tr>
        <tr><td>Number of Stages</td><td>${results.stages}</td></tr>
        <tr><td>Feed Stage</td><td>${results.feed_stage}</td></tr>
    `;
    document.getElementById('summaryBody').innerHTML = summaryHtml;
    
    // Stages table
    let stagesHtml = '';
    results.stage_compositions.forEach((stage, i) => {
        stagesHtml += `<tr>
            <td>Stage ${i + 1}</td>
            <td>${stage.x.toFixed(3)}</td>
            <td>${stage.y.toFixed(3)}</td>
        </tr>`;
    });
    document.getElementById('stagesBody').innerHTML = stagesHtml;
    
    // Enable export button
    document.getElementById('exportBtn').disabled = false;
    currentResults = results;
}

// ==========================================
// INITIALIZE PYODIDE (VERSI 0.23.4)
// ==========================================
async function initPyodide() {
    const loadingDiv = document.getElementById('pyodide-loading');
    const progressBar = document.getElementById('loading-progress');
    const statusText = document.getElementById('loading-status');
    const statusDiv = document.getElementById('pyodide-status');
    const calculateBtn = document.getElementById('calculateBtn');
    
    loadingDiv.style.display = 'flex';
    
    try {
        // Step 1: Load Pyodide core
        statusText.innerText = 'Loading Pyodide core...';
        progressBar.style.width = '20%';
        progressBar.innerText = '20%';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading Python...';
        
        pyodide = await loadPyodide({
            indexURL: "https://cdn.jsdelivr.net/pyodide/v0.23.4/full/"
        });
        
        // Step 2: Load packages
        statusText.innerText = 'Loading numpy...';
        progressBar.style.width = '40%';
        progressBar.innerText = '40%';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading numpy...';
        await pyodide.loadPackage('numpy');
        
        statusText.innerText = 'Loading scipy...';
        progressBar.style.width = '60%';
        progressBar.innerText = '60%';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading scipy...';
        await pyodide.loadPackage('scipy');
        
        statusText.innerText = 'Loading pandas...';
        progressBar.style.width = '80%';
        progressBar.innerText = '80%';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading pandas...';
        await pyodide.loadPackage('pandas');
        
        // Step 3: Load calculator
        statusText.innerText = 'Initializing calculator...';
        progressBar.style.width = '90%';
        progressBar.innerText = '90%';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading calculator...';
        await loadCalculatorCode();
        
        // Selesai
        progressBar.style.width = '100%';
        progressBar.innerText = '100%';
        statusText.innerText = 'Ready!';
        statusDiv.innerHTML = '<i class="fas fa-check-circle text-success"></i> Python environment ready!';
        pyodideReady = true;
        calculateBtn.disabled = false;
        
        setTimeout(() => {
            loadingDiv.style.display = 'none';
        }, 500);
        
        console.log('✅ Pyodide siap digunakan!');
        
    } catch (error) {
        console.error('Error loading Pyodide:', error);
        statusText.innerText = 'Error loading Python environment';
        statusDiv.innerHTML = `<i class="fas fa-exclamation-circle text-danger"></i> Error: ${error.message}`;
        loadingDiv.style.display = 'none';
    }
}

// ==========================================
// EVENT LISTENERS
// ==========================================

// Preview update
['xData', 'yData', 'Hl', 'Hv'].forEach(id => {
    const element = document.getElementById(id);
    if (element) {
        element.addEventListener('input', updatePreview);
    }
});

// Custom q handling
document.getElementById('q').addEventListener('change', function() {
    document.getElementById('customQDiv').style.display = 
        this.value === 'custom' ? 'block' : 'none';
});

// Calculate button
document.getElementById('calculateBtn').addEventListener('click', async function() {
    const loading = document.getElementById('loading');
    const btn = this;
    
    // Parse arrays
    const xData = parseArrayString(document.getElementById('xData').value);
    const yData = parseArrayString(document.getElementById('yData').value);
    const Hl = parseArrayString(document.getElementById('Hl').value);
    const Hv = parseArrayString(document.getElementById('Hv').value);
    
    // Validate
    const validations = [
        validateArray(xData, 'xData'),
        validateArray(yData, 'yData'),
        validateArray(Hl, 'Hl'),
        validateArray(Hv, 'Hv')
    ];
    
    for (let v of validations) {
        if (!v.valid) {
            alert('Error: ' + v.error);
            return;
        }
    }
    
    // Check lengths
    if (xData.length !== yData.length || 
        xData.length !== Hl.length || 
        xData.length !== Hv.length) {
        alert('Error: Semua array harus memiliki panjang yang sama!');
        return;
    }
    
    // Get q value
    let q = document.getElementById('q').value;
    if (q === 'custom') {
        q = parseFloat(document.getElementById('customQ').value);
    } else {
        q = parseFloat(q);
    }
    
    // Prepare input data
    const inputData = {
        xData: xData,
        yData: yData,
        Hl: Hl,
        Hv: Hv,
        zF: parseFloat(document.getElementById('zF').value),
        F: parseFloat(document.getElementById('F').value),
        xD: parseFloat(document.getElementById('xD').value),
        xB: parseFloat(document.getElementById('xB').value),
        q: q,
        R: parseFloat(document.getElementById('R').value)
    };
    
    // Disable button and show loading
    btn.disabled = true;
    loading.style.display = 'block';
    
    try {
        // Run calculation
        const results = await runCalculation(inputData);
        
        if (results.error) {
            alert('Error: ' + results.error);
        } else {
            displayResults(results);
        }
        
    } catch (error) {
        alert('Error: ' + error.message);
    } finally {
        btn.disabled = false;
        loading.style.display = 'none';
    }
});

// Export button
document.getElementById('exportBtn').addEventListener('click', function() {
    if (!currentResults) return;
    
    // Create CSV content
    let csv = 'Stage,x (Liquid),y (Vapor)\n';
    currentResults.stage_compositions.forEach((stage, i) => {
        csv += `${i+1},${stage.x.toFixed(4)},${stage.y.toFixed(4)}\n`;
    });
    
    // Download CSV
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'ponchon_savarit_results.csv';
    a.click();
    window.URL.revokeObjectURL(url);
});

// Initial preview
updatePreview();

// Initialize Pyodide
initPyodide();
