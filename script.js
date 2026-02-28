// ==========================================
// STATE & VARIABLES GLOBAL
// ==========================================
let pyodide = null;
let pyodideReady = false;
let currentResults = null;
let currentSystem = 'user-defined';

// ==========================================
// ETHANOL-WATER DATASET (BRITISH UNITS)
// ==========================================
const ETHANOL_WATER_DATA = {
    '1atm': {
        name: 'Ethanol-Water at 1 atm',
        x: [0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 1.0],
        y: [0, 0.175, 0.350, 0.505, 0.655, 0.745, 0.805, 0.845, 0.880, 0.910, 0.935, 0.960, 0.975, 0.990, 1.0],
        Hl: [1392, 1380, 1362, 1329, 1288, 1255, 1226, 1202, 1185, 1168, 1157, 1146, 1140, 1134, 1128],
        Hv: [3506, 3465, 3425, 3351, 3265, 3197, 3140, 3088, 3042, 3002, 2968, 2939, 2922, 2910, 2905],
        T: [212.0, 203.9, 197.2, 189.5, 184.0, 180.1, 178.7, 177.3, 176.0, 174.9, 173.8, 173.1, 172.8, 172.7, 172.6],
        units: { enthalpy: 'BTU/lbmole', temperature: '°F', duty: 'BTU/hr' },
        azeotrope: { x: 0.895, T: 172.7 }
    },
    '76mmHg': {
        name: 'Ethanol-Water at 76 mmHg',
        x: [0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 1.0],
        y: [0, 0.192, 0.377, 0.527, 0.713, 0.746, 0.771, 0.794, 0.822, 0.912, 0.942, 0.959, 0.978, 0.990, 1.0],
        Hl: [1037, 1019, 1002, 974, 928, 893, 865, 842, 819, 802, 791, 779, 773, 768, 762],
        Hv: [2778, 2745, 2704, 2635, 2544, 2478, 2426, 2380, 2338, 2302, 2274, 2248, 2234, 2225, 2221],
        T: [212.0, 203.4, 197.2, 189.2, 184.5, 181.7, 179.6, 177.8, 176.2, 174.3, 173.0, 172.8, 172.7, 172.8, 173.0],
        units: { enthalpy: 'BTU/lbmole', temperature: '°F', duty: 'BTU/hr' },
        azeotrope: { x: 0.86, T: 172.8 }
    }
};

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
        console.error('Error parsing array:', e);
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
// SHOW TOAST MESSAGE
// ==========================================
function showToast(message, type = 'info') {
    const toastContainer = document.getElementById('toastContainer');
    if (!toastContainer) return;
    
    const toastId = 'toast_' + Date.now();
    const toast = document.createElement('div');
    toast.id = toastId;
    toast.className = `toast align-items-center text-white bg-${type} border-0`;
    toast.setAttribute('role', 'alert');
    toast.setAttribute('aria-live', 'assertive');
    toast.setAttribute('aria-atomic', 'true');
    
    toast.innerHTML = `
        <div class="d-flex">
            <div class="toast-body">
                ${message}
            </div>
            <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
        </div>
    `;
    
    toastContainer.appendChild(toast);
    const bsToast = new bootstrap.Toast(toast, { delay: 3000 });
    bsToast.show();
    
    toast.addEventListener('hidden.bs.toast', () => {
        toast.remove();
    });
}

// ==========================================
// UPDATE UNIT DISPLAY
// ==========================================
function updateUnitDisplay(systemType) {
    const isBritish = systemType.includes('ethanol');
    const hlUnit = document.getElementById('hlUnit');
    const hvUnit = document.getElementById('hvUnit');
    const previewUnitBadge = document.getElementById('previewUnitBadge');
    const unitNote = document.getElementById('unitNote');
    const flowUnit = document.getElementById('flowUnit');
    const flowNote = document.getElementById('flowNote');
    const dutyUnit = document.getElementById('dutyUnit');
    
    if (isBritish) {
        if (hlUnit) {
            hlUnit.textContent = 'BTU/lbmole';
            hlUnit.className = 'unit-badge bg-primary text-white';
        }
        if (hvUnit) {
            hvUnit.textContent = 'BTU/lbmole';
            hvUnit.className = 'unit-badge bg-primary text-white';
        }
        if (previewUnitBadge) {
            previewUnitBadge.textContent = 'BTU/lbmole';
            previewUnitBadge.className = 'unit-badge bg-primary text-white';
        }
        if (unitNote) unitNote.innerHTML = '<b>British Units:</b> Enthalpy in BTU/lbmole, Flow in lbmol/hr, Duty in BTU/hr';
        if (flowUnit) {
            flowUnit.textContent = 'lbmol/hr';
            flowUnit.className = 'unit-badge bg-success text-white';
        }
        if (flowNote) flowNote.innerHTML = 'Flow rate in lbmol/hr (British units)';
        if (dutyUnit) dutyUnit.textContent = 'BTU/hr';
        
        document.getElementById('xData').readOnly = true;
        document.getElementById('yData').readOnly = true;
        document.getElementById('Hl').readOnly = true;
        document.getElementById('Hv').readOnly = true;
    } else {
        if (hlUnit) {
            hlUnit.textContent = 'any units';
            hlUnit.className = 'unit-badge';
        }
        if (hvUnit) {
            hvUnit.textContent = 'any units';
            hvUnit.className = 'unit-badge';
        }
        if (previewUnitBadge) {
            previewUnitBadge.textContent = 'any units';
            previewUnitBadge.className = 'unit-badge';
        }
        if (unitNote) unitNote.innerHTML = 'User Defined: bebas menggunakan satuan apapun';
        if (flowUnit) {
            flowUnit.textContent = 'kmol/hr';
            flowUnit.className = 'unit-badge bg-success text-white';
        }
        if (flowNote) flowNote.innerHTML = 'Flow rate in kmol/hr (SI units)';
        if (dutyUnit) dutyUnit.textContent = 'kW';
        
        document.getElementById('xData').readOnly = false;
        document.getElementById('yData').readOnly = false;
        document.getElementById('Hl').readOnly = false;
        document.getElementById('Hv').readOnly = false;
    }
}

// ==========================================
// LOAD ETHANOL-WATER DATASET
// ==========================================
function loadEthanolWaterDataset(pressure) {
    const dataset = ETHANOL_WATER_DATA[pressure];
    if (!dataset) return;
    
    console.log(`Loading ${dataset.name} dataset...`);
    
    document.getElementById('xData').value = JSON.stringify(dataset.x);
    document.getElementById('yData').value = JSON.stringify(dataset.y);
    document.getElementById('Hl').value = JSON.stringify(dataset.Hl);
    document.getElementById('Hv').value = JSON.stringify(dataset.Hv);
    document.getElementById('TData').value = JSON.stringify(dataset.T);
    
    // Update dataset info panel
    const datasetInfoPanel = document.getElementById('datasetInfoPanel');
    const datasetName = document.getElementById('datasetName');
    const datasetDetails = document.getElementById('datasetDetails');
    
    if (datasetInfoPanel) datasetInfoPanel.style.display = 'block';
    if (datasetName) datasetName.textContent = dataset.name;
    
    let azeoInfo = '';
    if (dataset.azeotrope) {
        azeoInfo = ` | Azeotrope: x = ${dataset.azeotrope.x.toFixed(3)} at ${dataset.azeotrope.T.toFixed(1)}°F`;
    }
    
    if (datasetDetails) {
        datasetDetails.innerHTML = `
            <i class="fas fa-database me-1"></i> ${dataset.x.length} points
            <i class="fas fa-fire ms-2 me-1"></i> Hl: ${dataset.Hl[0]}-${dataset.Hl[dataset.Hl.length-1]} BTU/lbmole
            <i class="fas fa-temperature-high ms-2 me-1"></i> T: ${dataset.T[0].toFixed(1)}-${dataset.T[dataset.T.length-1].toFixed(1)}°F
            ${azeoInfo}
        `;
    }
    
    const systemInfo = document.getElementById('systemInfo');
    const systemInfoText = document.getElementById('systemInfoText');
    if (systemInfo) systemInfo.style.display = 'block';
    if (systemInfoText) {
        systemInfoText.innerHTML = `
            <b>${dataset.name}</b><br>
            📊 ${dataset.x.length} data points<br>
            🔥 Enthalpy: <b>${dataset.units.enthalpy}</b><br>
            🌡️ Temperature: <b>${dataset.units.temperature}</b>
        `;
    }
    
    updateUnitDisplay(`ethanol-${pressure}`);
    updatePreview();
    showToast(`✅ Loaded ${dataset.name} (British Units)`, 'success');
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
    
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish = systemType.includes('ethanol');
    
    let html = '';
    for (let i = 0; i < maxLength; i++) {
        html += '<tr>';
        html += `<td>${i + 1}</td>`;
        html += `<td>${i < xData.length ? xData[i].toFixed(3) : '-'}</td>`;
        html += `<td>${i < yData.length ? yData[i].toFixed(3) : '-'}</td>`;
        
        if (i < Hl.length) {
            html += `<td>${isBritish ? Hl[i].toFixed(0) : Hl[i].toFixed(2)}</td>`;
        } else {
            html += '<td>-</td>';
        }
        
        if (i < Hv.length) {
            html += `<td>${isBritish ? Hv[i].toFixed(0) : Hv[i].toFixed(2)}</td>`;
        } else {
            html += '<td>-</td>';
        }
        
        html += '</tr>';
    }
    
    const previewBody = document.getElementById('previewBody');
    if (previewBody) previewBody.innerHTML = html;
    
    const warning = document.getElementById('previewWarning');
    if (new Set(lengths).size > 1) {
        if (!warning) {
            const div = document.createElement('div');
            div.id = 'previewWarning';
            div.className = 'alert alert-warning mt-2';
            div.innerHTML = '<i class="fas fa-exclamation-triangle"></i> Panjang array tidak sama!';
            document.querySelector('#previewTable')?.after(div);
        }
    } else if (warning) {
        warning.remove();
    }
}

// ==========================================
// LOAD CALCULATOR CODE (PYTHON)
// ==========================================
async function loadCalculatorCode() {
    try {
        const pythonCode = `
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import root_scalar
import json

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

def calculate_stages(xD, xB, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data):
    y = xD
    stages = 0
    tie_lines = []
    construction_lines = []
    stage_compositions = []
    in_rectifying = True
    feed_stage = 0
    error = ""

    # Stage 1: dari xD
    stage_compositions.append({'x': xD, 'y': xD})
    
    while stages < 20:
        # Cari x dari kurva kesetimbangan
        def find_x(x):
            return cubic_interpolate(x, data['xData'], data['yData']) - y
        
        try:
            sol = root_scalar(find_x, bracket=[0, 1], method='bisect')
            if not sol.converged:
                break
            x_n = sol.root
        except:
            break

        # Cek apakah sudah mencapai bottom
        if x_n <= xB:
            stage_compositions.append({'x': xB, 'y': y})
            stages += 1
            break

        # Hitung entalpi untuk tie line
        HLx_n = linear_interpolate(x_n, data['xData'], data['Hl'])
        HVy_n = cubic_interpolate(y, data['yData'], data['Hv'])
        
        tie_lines.append({'x': [x_n, y], 'y': [HLx_n, HVy_n]})
        stage_compositions.append({'x': x_n, 'y': y})
        stages += 1

        if in_rectifying and x_n <= zF:
            in_rectifying = False
            feed_stage = stages

        xDelta = xDeltaR if in_rectifying else xDeltaS
        HDelta = HDeltaR if in_rectifying else HDeltaS

        if abs(x_n - xDelta) < 1e-6:
            break
            
        slope = (HDelta - HLx_n) / (xDelta - x_n)
        
        # Cari y berikutnya
        def find_y(y_val):
            return cubic_interpolate(y_val, data['yData'], data['Hv']) - (HDelta + slope * (y_val - xDelta))
        
        try:
            sol = root_scalar(find_y, bracket=[xB, xD], method='bisect')
            y = sol.root
        except:
            break
        
        # Cari x_end untuk construction line
        if in_rectifying:
            def find_x_end(x):
                return linear_interpolate(x, data['xData'], data['Hl']) - (HDelta + slope * (x - xDelta))
            try:
                sol = root_scalar(find_x_end, bracket=[0, 1], method='bisect')
                xEnd = sol.root
                HEnd = linear_interpolate(xEnd, data['xData'], data['Hl'])
                construction_lines.append({'x': [xDelta, xEnd], 'y': [HDelta, HEnd]})
            except:
                pass

    return {
        'stages': stages,
        'feed_stage': feed_stage,
        'tie_lines': tie_lines,
        'construction_lines': construction_lines,
        'stage_compositions': stage_compositions,
        'error': error
    }

def calculate(data, params):
    zF = params['zF']
    F = params['F']
    xD = params['xD']
    xB = params['xB']
    q = params['q']
    R = params['R']
    
    yF = cubic_interpolate(zF, data['xData'], data['yData'])
    HLzF = linear_interpolate(zF, data['xData'], data['Hl'])
    HVzF = cubic_interpolate(yF, data['yData'], data['Hv'])
    HF = q * HLzF + (1 - q) * HVzF
    HD = linear_interpolate(xD, data['xData'], data['Hl'])
    HW = linear_interpolate(xB, data['xData'], data['Hl'])
    
    if not all(np.isfinite([yF, HLzF, HVzF, HF, HD, HW])):
        return {'error': "Interpolation failed. Check if compositions are within data range."}
    
    D = F * (zF - xB) / (xD - xB)
    W = F - D
    
    HVxD = cubic_interpolate(xD, data['yData'], data['Hv'])
    Qc = D * (HVxD - HD) * (R + 1)
    QcDuty = round(Qc, 2)
    
    xDeltaR = xD
    HDeltaR = HD + Qc / D
    
    xDeltaS = xB
    slope = (HDeltaR - HF) / (xDeltaR - zF)
    HDeltaS = HF + slope * (xDeltaS - zF)
    
    Qr = W * (HW - HDeltaS)
    QrDuty = round(Qr, 2)
    
    yFMin = cubic_interpolate(zF, data['xData'], data['yData'])
    HVyF = cubic_interpolate(yFMin, data['yData'], data['Hv'])
    slopeMin = (HVyF - HF) / (yFMin - zF)
    QPrimeMin = HF + slopeMin * (xD - zF)
    QDoublePrimeMin = HF + slopeMin * (xB - zF)
    RMin = (QPrimeMin - HVxD) / (HVxD - HD)
    
    stage_results = calculate_stages(
        xD, xB, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data
    )
    
    if stage_results['error']:
        return {'error': stage_results['error']}
    
    x_range = np.linspace(0, 1, 200).tolist()
    HL_curve = [linear_interpolate(xi, data['xData'], data['Hl']) for xi in x_range]
    HV_curve = [cubic_interpolate(xi, data['yData'], data['Hv']) for xi in x_range]
    y_equilibrium = [cubic_interpolate(xi, data['xData'], data['yData']) for xi in x_range]
    
    y_values = data['Hl'] + data['Hv'] + [HF, HD, HW, HDeltaR, HDeltaS, QPrimeMin, QDoublePrimeMin, HVyF]
    y_values = [y for y in y_values if np.isfinite(y)]
    yMin, yMax = min(y_values) - 10, max(y_values) + 10
    
    return {
        'D': round(D, 2), 'W': round(W, 2),
        'xDeltaR': round(xDeltaR, 3), 'HDeltaR': round(HDeltaR, 2),
        'xDeltaS': round(xDeltaS, 3), 'HDeltaS': round(HDeltaS, 2),
        'QcDuty': QcDuty, 'QrDuty': QrDuty,
        'QPrimeMin': round(QPrimeMin, 2), 'QDoublePrimeMin': round(QDoublePrimeMin, 2),
        'RMin': round(RMin, 2), 
        'stages': stage_results['stages'],
        'feed_stage': stage_results['feed_stage'],
        'stage_compositions': stage_results['stage_compositions'],
        'tie_lines': stage_results['tie_lines'],
        'construction_lines': stage_results['construction_lines'],
        'x_range': x_range, 'HL_curve': HL_curve, 'HV_curve': HV_curve,
        'y_equilibrium': y_equilibrium, 'yMin': yMin, 'yMax': yMax,
        'HF': HF, 'zF': zF, 'xD': xD, 'xB': xB,
        'yFMin': yFMin, 'HVyF': HVyF
    }

def calculate_from_js(xData, yData, Hl, Hv, zF, F, xD, xB, q, R):
    try:
        data = {
            'xData': [float(x) for x in xData],
            'yData': [float(y) for y in yData],
            'Hl': [float(h) for h in Hl],
            'Hv': [float(h) for h in Hv]
        }
        params = {'zF': float(zF), 'F': float(F), 'xD': float(xD), 
                  'xB': float(xB), 'q': float(q), 'R': float(R)}
        return json.dumps(calculate(data, params))
    except Exception as e:
        return json.dumps({'error': str(e)})
`;
        
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
    if (!pyodideReady) throw new Error('Pyodide belum siap.');
    try {
        return JSON.parse(pyodide.runPython(`
            calculate_from_js(
                ${JSON.stringify(inputData.xData)},
                ${JSON.stringify(inputData.yData)},
                ${JSON.stringify(inputData.Hl)},
                ${JSON.stringify(inputData.Hv)},
                ${inputData.zF}, ${inputData.F}, ${inputData.xD}, 
                ${inputData.xB}, ${inputData.q}, ${inputData.R}
            )
        `));
    } catch (error) {
        console.error('Error running calculation:', error);
        throw error;
    }
}

// ==========================================
// CREATE PLOT - PONCHON-SAVARIT YANG BENAR
// ==========================================
function createPlot(results) {
    console.log('Creating plot with stage compositions:', results.stage_compositions);
    console.log('Tie lines:', results.tie_lines);
    
    const stageColors = [
        '#FF6B6B', '#4ECDC4', '#FF9F1C', '#6A4C93', '#2E86AB', 
        '#A23B72', '#F18F01', '#2D6A4F', '#9E2A2B', '#540D6E'
    ];
    
    const traces = [];
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish = systemType.includes('ethanol');
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
    
    // ========== KURVA DASAR H-x-y ==========
    traces.push({
        x: results.x_range,
        y: results.HL_curve,
        mode: 'lines',
        name: 'Saturated Liquid',
        line: {color: '#2E86AB', width: 4},
        xaxis: 'x', 
        yaxis: 'y'
    });
    
    traces.push({
        x: results.x_range,
        y: results.HV_curve,
        mode: 'lines',
        name: 'Saturated Vapor',
        line: {color: '#A23B72', width: 4},
        xaxis: 'x', 
        yaxis: 'y'
    });
    
    // ========== GARIS VERTIKAL ==========
    traces.push({
        x: [results.xD, results.xD],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'x<sub>D</sub>',
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x', 
        yaxis: 'y'
    });
    traces.push({
        x: [results.xD, results.xD],
        y: [0, 1],
        mode: 'lines',
        showlegend: false,
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x2', 
        yaxis: 'y2'
    });
    
    traces.push({
        x: [results.xB, results.xB],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'x<sub>B</sub>',
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x', 
        yaxis: 'y'
    });
    traces.push({
        x: [results.xB, results.xB],
        y: [0, 1],
        mode: 'lines',
        showlegend: false,
        line: {color: '#6C757D', width: 2.5, dash: 'dash'},
        xaxis: 'x2', 
        yaxis: 'y2'
    });
    
    traces.push({
        x: [results.zF, results.zF],
        y: [results.yMin, results.yMax],
        mode: 'lines',
        name: 'z<sub>F</sub>',
        line: {color: '#2D6A4F', width: 3, dash: 'dash'},
        xaxis: 'x', 
        yaxis: 'y'
    });
    traces.push({
        x: [results.zF, results.zF],
        y: [0, 1],
        mode: 'lines',
        showlegend: false,
        line: {color: '#2D6A4F', width: 3, dash: 'dash'},
        xaxis: 'x2', 
        yaxis: 'y2'
    });
    
    // ========== DIFFERENCE POINTS ==========
    traces.push({
        x: [results.xDeltaR, results.xDeltaS],
        y: [results.HDeltaR, results.HDeltaS],
        mode: 'markers+text',
        name: 'Difference Points',
        marker: {color: '#F97316', size: 14, symbol: 'star'},
        text: ['Δ<sub>R</sub>', 'Δ<sub>S</sub>'],
        textposition: ['top center', 'bottom center'],
        xaxis: 'x', 
        yaxis: 'y'
    });
    
    // ========== OPERATING LINE ==========
    traces.push({
        x: [results.xDeltaR, results.zF, results.xDeltaS],
        y: [results.HDeltaR, results.HF, results.HDeltaS],
        mode: 'lines+markers',
        name: 'Operating Line',
        line: {color: '#0A9396', width: 3},
        marker: {size: 8, color: '#0A9396'},
        xaxis: 'x', 
        yaxis: 'y'
    });
    
    // ========== MINIMUM REFLUX LINE ==========
    if (results.yFMin && results.HVyF) {
        traces.push({
            x: [results.xB, results.zF, results.yFMin, results.xD],
            y: [results.QDoublePrimeMin, results.HF, results.HVyF, results.QPrimeMin],
            mode: 'lines+markers',
            name: 'Minimum Reflux Line',
            line: {color: '#E9C46A', width: 2.5, dash: 'dash'},
            marker: {size: 8, color: '#E9C46A'},
            xaxis: 'x', 
            yaxis: 'y'
        });
    }
    
    // ========== VLE CURVE ==========
    traces.push({
        x: results.x_range,
        y: results.y_equilibrium,
        mode: 'lines',
        name: 'Equilibrium Curve',
        line: {color: '#1E1E1E', width: 3.5},
        xaxis: 'x2', 
        yaxis: 'y2'
    });
    
    // ========== y = x LINE ==========
    traces.push({
        x: [0, 1],
        y: [0, 1],
        mode: 'lines',
        name: 'y = x',
        line: {color: '#6C757D', width: 2, dash: 'dash'},
        xaxis: 'x2', 
        yaxis: 'y2'
    });
    
    // ========== TIE LINES (DIAGRAM ATAS) ==========
    if (results.tie_lines?.length) {
        results.tie_lines.forEach((tie, i) => {
            const color = stageColors[i % stageColors.length];
            
            // Tie line
            traces.push({
                x: tie.x,
                y: tie.y,
                mode: 'lines',
                name: `Stage ${i+1}`,
                line: {color: color, width: 3},
                legendgroup: `stage_${i+1}`,
                xaxis: 'x', 
                yaxis: 'y'
            });
            
            // Markers
            traces.push({
                x: [tie.x[0], tie.x[1]],
                y: [tie.y[0], tie.y[1]],
                mode: 'markers',
                marker: {
                    color: [color, color],
                    size: [10, 12],
                    symbol: ['circle', 'diamond'],
                    line: {color: 'white', width: 1}
                },
                showlegend: false,
                legendgroup: `stage_${i+1}`,
                xaxis: 'x', 
                yaxis: 'y'
            });
        });
    }
    
    // ========== STAGE TRACING (DIAGRAM BAWAH) ==========
    if (results.stage_compositions?.length > 1) {
        for (let i = 0; i < results.stage_compositions.length - 1; i++) {
            const currentStage = results.stage_compositions[i];
            const nextStage = results.stage_compositions[i + 1];
            const color = stageColors[i % stageColors.length];
            
            // Garis vertikal: dari (x, y) ke (x, x)
            traces.push({
                x: [currentStage.x, currentStage.x],
                y: [currentStage.y, currentStage.x],
                mode: 'lines',
                line: {color: color, width: 2.5},
                showlegend: false,
                legendgroup: `stage_${i+1}`,
                xaxis: 'x2', 
                yaxis: 'y2'
            });
            
            // Garis horizontal: dari (x, x) ke (nextX, nextX)
            traces.push({
                x: [currentStage.x, nextStage.x],
                y: [currentStage.x, nextStage.x],
                mode: 'lines',
                line: {color: color, width: 2.5},
                showlegend: false,
                legendgroup: `stage_${i+1}`,
                xaxis: 'x2', 
                yaxis: 'y2'
            });
            
            // Titik untuk current stage
            traces.push({
                x: [currentStage.x],
                y: [currentStage.y],
                mode: 'markers',
                marker: {color: color, size: 10, symbol: 'circle'},
                showlegend: false,
                legendgroup: `stage_${i+1}`,
                xaxis: 'x2', 
                yaxis: 'y2'
            });
            
            // Titik untuk next stage (hanya sekali)
            if (i === results.stage_compositions.length - 2) {
                traces.push({
                    x: [nextStage.x],
                    y: [nextStage.y],
                    mode: 'markers',
                    marker: {color: color, size: 10, symbol: 'circle'},
                    showlegend: false,
                    legendgroup: `stage_${i+2}`,
                    xaxis: 'x2', 
                    yaxis: 'y2'
                });
            }
        }
    }
    
    // ========== GARIS PROYEKSI VERTIKAL ==========
    if (results.stage_compositions?.length) {
        results.stage_compositions.forEach((stage, i) => {
            if (i === 0) return; // Skip stage 0 (xD) karena sudah diwakili
            
            const color = stageColors[(i-1) % stageColors.length];
            const idxLiq = Math.round(stage.x * 199);
            const idxVap = Math.round(stage.y * 199);
            const H_liq = results.HL_curve?.[idxLiq];
            const H_vap = results.HV_curve?.[idxVap];
            
            if (H_liq && H_vap) {
                // Proyeksi liquid
                traces.push({
                    x: [stage.x, stage.x],
                    y: [H_liq, results.yMin],
                    mode: 'lines',
                    showlegend: false,
                    line: {color: color, width: 1.5, dash: 'dot'},
                    legendgroup: `stage_${i}`,
                    xaxis: 'x', 
                    yaxis: 'y'
                });
                
                // Proyeksi vapor
                traces.push({
                    x: [stage.y, stage.y],
                    y: [H_vap, results.yMin],
                    mode: 'lines',
                    showlegend: false,
                    line: {color: color, width: 1.5, dash: 'dot'},
                    legendgroup: `stage_${i}`,
                    xaxis: 'x', 
                    yaxis: 'y'
                });
            }
        });
    }
    
    // ========== LAYOUT ==========
    const layout = {
        title: {
            text: '<b>Ponchon–Savarit Diagram: Ethanol-Water System</b>',
            font: {size: 18},
            x: 0.5
        },
        grid: {
            rows: 2,
            columns: 1,
            pattern: 'independent',
            roworder: 'top to bottom'
        },
        margin: { l: 70, r: 130, t: 50, b: 70 },
        xaxis: {
            domain: [0.1, 0.9],
            title: '',
            range: [0, 1],
            tickformat: '.2f',
            showline: true,
            linecolor: '#1E1E1E',
            mirror: true
        },
        yaxis: {
            domain: [0.5, 0.95],
            title: `<b>Enthalpy (${enthalpyUnit})</b>`,
            range: [results.yMin, results.yMax],
            showline: true,
            linecolor: '#1E1E1E',
            mirror: true
        },
        xaxis2: {
            domain: [0.1, 0.9],
            title: '<b>Mole Fraction (x or y)</b>',
            range: [0, 1],
            tickformat: '.2f',
            showline: true,
            linecolor: '#1E1E1E',
            mirror: true
        },
        yaxis2: {
            domain: [0.05, 0.45],
            title: '<b>Vapor Fraction (y)</b>',
            range: [0, 1],
            showline: true,
            linecolor: '#1E1E1E',
            mirror: true
        },
        legend: {
            x: 1.02,
            y: 1,
            xanchor: 'left',
            yanchor: 'top',
            font: {size: 9},
            bgcolor: 'rgba(255,255,255,0.9)',
            bordercolor: '#1E1E1E',
            borderwidth: 1
        },
        height: 700,
        width: document.querySelector('.main-content')?.clientWidth - 40 || 1000
    };

    Plotly.newPlot('plotDiv', traces, layout, {responsive: true});
    
    window.addEventListener('resize', () => {
        const containerWidth = document.querySelector('.main-content')?.clientWidth || 900;
        Plotly.relayout('plotDiv', {width: containerWidth - 40});
    });
}
// ==========================================
// UPDATE QUICK SUMMARY
// ==========================================
function updateQuickSummary(results) {
    const quickSummary = document.getElementById('quickSummary');
    const quickSummaryContent = document.getElementById('quickSummaryContent');
    
    if (!quickSummary || !quickSummaryContent) return;
    
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish = systemType.includes('ethanol');
    const flowUnit = isBritish ? 'lbmol/hr' : 'kmol/hr';
    const dutyUnit = isBritish ? 'BTU/hr' : 'kW';
    
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
                <strong>${results.QcDuty?.toLocaleString() || 0} ${dutyUnit}</strong>
            </div>
        </div>
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>D:</span>
                <strong>${results.D} ${flowUnit}</strong>
            </div>
        </div>
        <div class="col-md-3 col-sm-6">
            <div class="summary-item">
                <span>W:</span>
                <strong>${results.W} ${flowUnit}</strong>
            </div>
        </div>
    `;
    
    quickSummary.style.display = 'block';
}

// ==========================================
// UPDATE ENTHALPY TABLE
// ==========================================
function updateEnthalpyTable(results) {
    const enthalpyBody = document.getElementById('enthalpyBody');
    if (!enthalpyBody) return;
    
    let html = '';
    if (results.stage_compositions?.length) {
        results.stage_compositions.forEach((stage, i) => {
            const idxLiq = Math.round(stage.x * 199);
            const idxVap = Math.round(stage.y * 199);
            const HL = results.HL_curve?.[idxLiq]?.toFixed(0) || 'N/A';
            const HV = results.HV_curve?.[idxVap]?.toFixed(0) || 'N/A';
            
            html += `<tr><td>Stage ${i+1}</td><td>${HL}</td><td>${HV}</td></tr>`;
        });
    } else {
        html = '<tr><td colspan="3" class="text-center text-muted">No enthalpy data available</td></tr>';
    }
    enthalpyBody.innerHTML = html;
}

// ==========================================
// DISPLAY RESULTS
// ==========================================
function displayResults(results) {
    console.log('Displaying results:', results);
    
    createPlot(results);
    updateQuickSummary(results);
    updateEnthalpyTable(results);
    
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish = systemType.includes('ethanol');
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
    const flowUnit = isBritish ? 'lbmol/hr' : 'kmol/hr';
    const dutyUnit = isBritish ? 'BTU/hr' : 'kW';
    
    // Summary table
    const summaryBody = document.getElementById('summaryBody');
    if (summaryBody) {
        const summaryHtml = `
            <tr><td>Distillate Flow (D)</td><td>${results.D} ${flowUnit}</td></tr>
            <tr><td>Bottoms Flow (W)</td><td>${results.W} ${flowUnit}</td></tr>
            <tr><td>Δ_R</td><td>(${results.xDeltaR}, ${results.HDeltaR} ${enthalpyUnit})</td></tr>
            <tr><td>Δ_S</td><td>(${results.xDeltaS}, ${results.HDeltaS} ${enthalpyUnit})</td></tr>
            <tr><td>Condenser Duty (Qc)</td><td>${results.QcDuty?.toLocaleString() || 0} ${dutyUnit}</td></tr>
            <tr><td>Reboiler Duty (Qr)</td><td>${results.QrDuty?.toLocaleString() || 0} ${dutyUnit}</td></tr>
            <tr><td>Δ_R min</td><td>(${results.xD}, ${results.QPrimeMin} ${enthalpyUnit})</td></tr>
            <tr><td>Δ_S min</td><td>(${results.xB}, ${results.QDoublePrimeMin} ${enthalpyUnit})</td></tr>
            <tr><td>Minimum Reflux Ratio</td><td>${results.RMin}</td></tr>
            <tr><td>Number of Stages</td><td>${results.stages}</td></tr>
            <tr><td>Feed Stage</td><td>${results.feed_stage}</td></tr>
        `;
        summaryBody.innerHTML = summaryHtml;
    }
    
    // Stages table
    const stagesBody = document.getElementById('stagesBody');
    if (stagesBody) {
        let stagesRows = '';
        if (results.stage_compositions?.length) {
            results.stage_compositions.forEach((stage, i) => {
                stagesRows += `<tr><td>Stage ${i+1}</td><td>${stage.x.toFixed(4)}</td><td>${stage.y.toFixed(4)}</td></tr>`;
            });
        } else {
            stagesRows = '<tr><td colspan="3" class="text-center text-muted">No stage data available</td></tr>';
        }
        stagesBody.innerHTML = stagesRows;
    }
    
    const exportBtn = document.getElementById('exportBtn');
    if (exportBtn) exportBtn.disabled = false;
    
    currentResults = results;
    showToast('✅ Calculation completed!', 'success');
}

// ==========================================
// INITIALIZE PYODIDE
// ==========================================
async function initPyodide() {
    const loadingDiv = document.getElementById('pyodide-loading');
    const progressBar = document.getElementById('loading-progress');
    const statusText = document.getElementById('loading-status');
    const statusDiv = document.getElementById('pyodide-status');
    const calculateBtn = document.getElementById('calculateBtn');
    
    if (!loadingDiv || !progressBar || !statusText || !statusDiv || !calculateBtn) return;
    
    loadingDiv.style.display = 'flex';
    
    try {
        statusText.innerText = 'Loading Pyodide core...';
        progressBar.style.width = '20%'; 
        progressBar.innerText = '20%';
        
        pyodide = await loadPyodide({
            indexURL: "https://cdn.jsdelivr.net/pyodide/v0.23.4/full/"
        });
        
        statusText.innerText = 'Loading numpy...';
        progressBar.style.width = '40%'; 
        progressBar.innerText = '40%';
        await pyodide.loadPackage('numpy');
        
        statusText.innerText = 'Loading scipy...';
        progressBar.style.width = '60%'; 
        progressBar.innerText = '60%';
        await pyodide.loadPackage('scipy');
        
        statusText.innerText = 'Loading pandas...';
        progressBar.style.width = '80%'; 
        progressBar.innerText = '80%';
        await pyodide.loadPackage('pandas');
        
        statusText.innerText = 'Initializing calculator...';
        progressBar.style.width = '90%'; 
        progressBar.innerText = '90%';
        await loadCalculatorCode();
        
        progressBar.style.width = '100%'; 
        progressBar.innerText = '100%';
        statusText.innerText = 'Ready!';
        statusDiv.innerHTML = '<i class="fas fa-check-circle text-success"></i> Python ready!';
        pyodideReady = true;
        calculateBtn.disabled = false;
        
        showToast('✅ Pyodide loaded!', 'success');
        setTimeout(() => loadingDiv.style.display = 'none', 500);
        
    } catch (error) {
        console.error('Pyodide error:', error);
        statusDiv.innerHTML = `<i class="fas fa-exclamation-circle text-danger"></i> Error: ${error.message}`;
        loadingDiv.style.display = 'none';
    }
}

// ==========================================
// EVENT LISTENERS
// ==========================================
document.addEventListener('DOMContentLoaded', function() {
    // Preview updates
    ['xData', 'yData', 'Hl', 'Hv'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('input', updatePreview);
    });
    
    // System type change
    const systemType = document.getElementById('systemType');
    if (systemType) {
        systemType.addEventListener('change', function() {
            const type = this.value;
            currentSystem = type;
            
            const quickSummary = document.getElementById('quickSummary');
            if (quickSummary) quickSummary.style.display = 'none';
            
            const systemInfo = document.getElementById('systemInfo');
            const datasetInfoPanel = document.getElementById('datasetInfoPanel');
            
            if (type === 'user-defined') {
                if (systemInfo) systemInfo.style.display = 'none';
                if (datasetInfoPanel) datasetInfoPanel.style.display = 'none';
                
                document.getElementById('xData').value = '[0, 0.08, 0.18, 0.25, 0.49, 0.65, 0.79, 0.91, 1.0]';
                document.getElementById('yData').value = '[0, 0.28, 0.43, 0.51, 0.73, 0.83, 0.90, 0.96, 1.0]';
                document.getElementById('Hl').value = '[24.3, 24.1, 23.2, 22.8, 22.05, 21.75, 21.7, 21.6, 21.4]';
                document.getElementById('Hv').value = '[61.2, 59.6, 58.5, 58.1, 56.5, 55.2, 54.4, 53.8, 53.3]';
                
                updateUnitDisplay('user-defined');
                showToast('📝 User Defined mode', 'info');
            } else if (type === 'ethanol-water-1atm') {
                loadEthanolWaterDataset('1atm');
            } else if (type === 'ethanol-water-76mmHg') {
                loadEthanolWaterDataset('76mmHg');
            }
            
            updatePreview();
        });
    }
    
    // q selector
    const qSelect = document.getElementById('q');
    if (qSelect) {
        qSelect.addEventListener('change', function() {
            const customQDiv = document.getElementById('customQDiv');
            if (customQDiv) {
                customQDiv.style.display = this.value === 'custom' ? 'block' : 'none';
            }
        });
    }
    
    // Calculate button
    const calculateBtn = document.getElementById('calculateBtn');
    if (calculateBtn) {
        calculateBtn.addEventListener('click', async function() {
            const loading = document.getElementById('loading');
            const btn = this;
            
            if (!pyodideReady) {
                alert('Pyodide belum siap. Tunggu...');
                return;
            }
            
            const xData = parseArrayString(document.getElementById('xData').value);
            const yData = parseArrayString(document.getElementById('yData').value);
            const Hl = parseArrayString(document.getElementById('Hl').value);
            const Hv = parseArrayString(document.getElementById('Hv').value);
            
            const validations = [
                validateArray(xData, 'xData'), validateArray(yData, 'yData'),
                validateArray(Hl, 'Hl'), validateArray(Hv, 'Hv')
            ];
            for (let v of validations) if (!v.valid) { alert('Error: ' + v.error); return; }
            
            if (xData.length !== yData.length || xData.length !== Hl.length || xData.length !== Hv.length) {
                alert('Error: Semua array harus sama panjang!'); 
                return;
            }
            
            let q = document.getElementById('q').value;
            q = q === 'custom' ? parseFloat(document.getElementById('customQ').value) : parseFloat(q);
            
            const inputData = {
                xData, yData, Hl, Hv,
                zF: parseFloat(document.getElementById('zF').value),
                F: parseFloat(document.getElementById('F').value),
                xD: parseFloat(document.getElementById('xD').value),
                xB: parseFloat(document.getElementById('xB').value),
                q, R: parseFloat(document.getElementById('R').value)
            };
            
            btn.disabled = true;
            if (loading) loading.style.display = 'block';
            
            try {
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
                if (loading) loading.style.display = 'none';
            }
        });
    }
    
    // Export button
    const exportBtn = document.getElementById('exportBtn');
    if (exportBtn) {
        exportBtn.addEventListener('click', function() {
            if (!currentResults) return;
            
            const systemType = document.getElementById('systemType')?.value || 'user-defined';
            const isBritish = systemType.includes('ethanol');
            const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
            
            let csv = 'Stage,x (Liquid),y (Vapor),HL (' + enthalpyUnit + '),HV (' + enthalpyUnit + ')\n';
            if (currentResults.stage_compositions?.length) {
                currentResults.stage_compositions.forEach((stage, i) => {
                    const idxLiq = Math.round(stage.x * 199);
                    const idxVap = Math.round(stage.y * 199);
                    const HL = currentResults.HL_curve?.[idxLiq]?.toFixed(2) || 'N/A';
                    const HV = currentResults.HV_curve?.[idxVap]?.toFixed(2) || 'N/A';
                    
                    csv += `${i+1},${stage.x.toFixed(4)},${stage.y.toFixed(4)},${HL},${HV}\n`;
                });
            }
            
            const blob = new Blob([csv], { type: 'text/csv' });
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url; 
            a.download = `ponchon_savarit_results_${systemType}_${new Date().toISOString().slice(0,10)}.csv`; 
            a.click();
            window.URL.revokeObjectURL(url);
            
            showToast('📊 Results exported to CSV', 'success');
        });
    }
    
    // Initialize
    updatePreview();
    initPyodide();
});


