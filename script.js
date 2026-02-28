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
// CREATE PLOT - PONCHON-SAVARIT (VLE STEPPING UPWARD - FIXED)
// ==========================================
function createPlot(results) {
    const stageColors = [
        '#E63946', '#2196F3', '#E9A825', '#7B2D8B', '#00ACC1',
        '#388E3C', '#FF7043', '#8D6E63', '#EC407A', '#546E7A',
        '#26A69A', '#AB47BC'
    ];

    const traces = [];
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish = systemType.includes('ethanol');
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';

    // ========== SATURATED LIQUID & VAPOR CURVES ==========
    traces.push({
        x: results.x_range, y: results.HL_curve,
        mode: 'lines', name: 'Saturated Liquid',
        line: { color: '#1565C0', width: 3 },
        xaxis: 'x', yaxis: 'y'
    });
    traces.push({
        x: results.x_range, y: results.HV_curve,
        mode: 'lines', name: 'Saturated Vapor',
        line: { color: '#880E4F', width: 3 },
        xaxis: 'x', yaxis: 'y'
    });

    // ========== VERTICAL DASHED LINES xD, xB, zF ==========
    const verticals = [
        { val: results.xD, label: 'x<sub>D</sub>', color: '#9E9E9E' },
        { val: results.xB, label: 'x<sub>B</sub>', color: '#9E9E9E' },
        { val: results.zF, label: 'z<sub>F</sub> (Feed)', color: '#2E7D32' }
    ];
    verticals.forEach(v => {
        traces.push({
            x: [v.val, v.val], y: [results.yMin, results.yMax],
            mode: 'lines', name: v.label,
            line: { color: v.color, width: 2, dash: 'dash' },
            xaxis: 'x', yaxis: 'y'
        });
        traces.push({
            x: [v.val, v.val], y: [0, 1],
            mode: 'lines', showlegend: false,
            line: { color: v.color, width: 2, dash: 'dash' },
            xaxis: 'x2', yaxis: 'y2'
        });
    });

    // ========== OPERATING LINE ==========
    traces.push({
        x: [results.xDeltaR, results.zF, results.xDeltaS],
        y: [results.HDeltaR, results.HF, results.HDeltaS],
        mode: 'lines+markers', name: 'Operating Line',
        line: { color: '#00897B', width: 2.5 },
        marker: { size: 7, color: '#00897B' },
        xaxis: 'x', yaxis: 'y'
    });

    // ========== DIFFERENCE POINTS ==========
    traces.push({
        x: [results.xDeltaR, results.xDeltaS],
        y: [results.HDeltaR, results.HDeltaS],
        mode: 'markers+text', name: 'Difference Points',
        marker: { color: '#FF6D00', size: 16, symbol: 'star' },
        text: ['Δ<sub>R</sub>', 'Δ<sub>S</sub>'],
        textposition: ['top center', 'bottom center'],
        textfont: { size: 13, color: '#FF6D00' },
        xaxis: 'x', yaxis: 'y'
    });

    // ========== MINIMUM REFLUX LINE ==========
    if (results.yFMin !== undefined && results.HVyF !== undefined) {
        traces.push({
            x: [results.xB, results.zF, results.yFMin, results.xD],
            y: [results.QDoublePrimeMin, results.HF, results.HVyF, results.QPrimeMin],
            mode: 'lines+markers', name: 'Minimum Reflux Line',
            line: { color: '#F9A825', width: 2, dash: 'dash' },
            marker: { size: 7, color: '#F9A825', symbol: 'diamond' },
            xaxis: 'x', yaxis: 'y'
        });
        traces.push({
            x: [results.xD, results.xB],
            y: [results.QPrimeMin, results.QDoublePrimeMin],
            mode: 'markers+text', name: 'Min Diff Points',
            marker: { color: '#F9A825', size: 12, symbol: 'diamond' },
            text: ["Δ'<sub>R,min</sub>", "Δ'<sub>S,min</sub>"],
            textposition: ['top right', 'bottom left'],
            textfont: { size: 11, color: '#F9A825' },
            xaxis: 'x', yaxis: 'y'
        });
    }

    // ========== CONSTRUCTION LINES ==========
    if (results.construction_lines?.length) {
        results.construction_lines.forEach((cl, i) => {
            traces.push({
                x: cl.x, y: cl.y, mode: 'lines',
                name: i === 0 ? 'Construction Line 1' : undefined,
                showlegend: i === 0,
                line: { color: '#BDBDBD', width: 1.5, dash: 'dot' },
                xaxis: 'x', yaxis: 'y'
            });
        });
    }

    // ========== VLE EQUILIBRIUM CURVE (drawn first, behind stages) ==========
    traces.push({
        x: results.x_range, y: results.y_equilibrium,
        mode: 'lines', name: 'Equilibrium Curve',
        line: { color: '#1A1A1A', width: 3.5 },
        xaxis: 'x2', yaxis: 'y2'
    });

    // ========== y = x DIAGONAL ==========
    traces.push({
        x: [0, 1], y: [0, 1],
        mode: 'lines', name: 'y = x',
        line: { color: '#9E9E9E', width: 1.5, dash: 'dash' },
        xaxis: 'x2', yaxis: 'y2'
    });

    // ========== STAGE TIE LINES + PROJECTIONS + VLE STEPPING ==========
    //
    // Python calculate_stages() builds tie_lines ordered from xD toward xB.
    // tie_lines[0] is closest to xD (top of column), tie_lines[N-1] closest to xB.
    //
    // Each tie line: tie.x[0] = x_liq, tie.x[1] = y_vap
    //   x_liq < y_vap always (liquid leaner than vapor)
    //
    // On VLE diagram, correct Ponchon-Savarit construction stepping UPWARD:
    //
    //   Start from xB on y=x diagonal, step UP to equilibrium curve, then
    //   horizontal RIGHT to y=x, repeat toward xD.
    //
    // For stage i (0 = closest to xD):
    //   Equilibrium point: (x_liq, y_vap)  ← on the equilibrium curve
    //   Tie line:  horizontal from (x_liq, y_vap) RIGHT to (y_vap, y_vap) on y=x
    //   Operating step: vertical from (y_vap, y_vap) UP to next equilibrium point
    //     → next stage equilibrium point is stage[i+1] = (x_liq_next, y_vap_next)
    //     → but we go from y_vap upward to y_vap_next (the y value of next stage)
    //
    // Since tie_lines are ordered xD→xB, we reverse to step xB→xD (upward):

    const tieLines = results.tie_lines || [];
    const numStages = tieLines.length;

    // Reverse so we step from bottom (xB side) to top (xD side)
    const tieLinesAsc = [...tieLines].reverse();

    let projLegendAdded = false;

    for (let i = 0; i < numStages; i++) {
        // Use original index for color (stage 1 = closest to xD = tieLines[0])
        const origIdx = numStages - 1 - i;
        const tie = tieLinesAsc[i];
        const color = stageColors[origIdx % stageColors.length];
        const stageNum = origIdx + 1;

        const x_liq = tie.x[0];   // liquid composition on HL
        const y_vap = tie.x[1];   // vapor composition on HV
        const H_liq = tie.y[0];
        const H_vap = tie.y[1];

        // ── H-x-y: TIE LINE ──
        traces.push({
            x: [x_liq, y_vap], y: [H_liq, H_vap],
            mode: 'lines', name: `Stage ${stageNum}`,
            line: { color: color, width: 2.5 },
            legendgroup: `stage_${stageNum}`,
            xaxis: 'x', yaxis: 'y'
        });
        traces.push({
            x: [x_liq, y_vap], y: [H_liq, H_vap],
            mode: 'markers',
            marker: { color: color, size: 10, symbol: 'circle', line: { color: 'white', width: 1.5 } },
            showlegend: false, legendgroup: `stage_${stageNum}`,
            xaxis: 'x', yaxis: 'y'
        });

        // ── H-x-y: PROJECTION LINES ──
        traces.push({
            x: [x_liq, x_liq], y: [H_liq, results.yMin],
            mode: 'lines',
            name: !projLegendAdded ? 'Projection Lines' : undefined,
            showlegend: !projLegendAdded,
            line: { color: color, width: 1.2, dash: 'dot' },
            legendgroup: `stage_${stageNum}`,
            xaxis: 'x', yaxis: 'y'
        });
        projLegendAdded = true;
        traces.push({
            x: [y_vap, y_vap], y: [H_vap, results.yMin],
            mode: 'lines', showlegend: false,
            line: { color: color, width: 1.2, dash: 'dot' },
            legendgroup: `stage_${stageNum}`,
            xaxis: 'x', yaxis: 'y'
        });

        // ══════════════════════════════════════════════════════
        // VLE DIAGRAM — Ponchon-Savarit stepping UPWARD
        //
        // Step sequence for stage i (ascending from xB to xD):
        //
        //  A = (x_liq, y_vap)        ← on equilibrium curve   [ABOVE]
        //  B = (y_vap, y_vap)        ← on y=x diagonal        [RIGHT of A, same height]
        //  C = (y_vap, y_vap_next)   ← on operating line      [ABOVE B, vertical up]
        //
        // where y_vap_next = y_vap of the NEXT stage (i+1) going toward xD
        // Last stage: y_vap_next = xD (distillate)
        //
        // So:
        //   Tie line:      A → B  (horizontal RIGHT)
        //   Operating line: B → C  (vertical UP)
        // ══════════════════════════════════════════════════════

        // y_vap of next stage (going toward xD = stage i+1 in ascending order)
        const nextTie = (i + 1 < numStages) ? tieLinesAsc[i + 1] : null;
        const y_vap_next = nextTie ? nextTie.x[1] : results.xD;

        // Point A: equilibrium curve (x_liq, y_vap)
        traces.push({
            x: [x_liq], y: [y_vap],
            mode: 'markers',
            marker: { color: color, size: 10, symbol: 'circle', line: { color: 'white', width: 1.5 } },
            showlegend: false, legendgroup: `stage_${stageNum}`,
            xaxis: 'x2', yaxis: 'y2'
        });

        // TIE LINE on VLE: horizontal RIGHT from A(x_liq, y_vap) → B(y_vap, y_vap)
        traces.push({
            x: [x_liq, y_vap], y: [y_vap, y_vap],
            mode: 'lines',
            line: { color: color, width: 2.5 },
            showlegend: false, legendgroup: `stage_${stageNum}`,
            xaxis: 'x2', yaxis: 'y2'
        });

        // Point B: y=x diagonal (y_vap, y_vap)
        traces.push({
            x: [y_vap], y: [y_vap],
            mode: 'markers',
            marker: { color: color, size: 8, symbol: 'diamond', line: { color: 'white', width: 1 } },
            showlegend: false, legendgroup: `stage_${stageNum}`,
            xaxis: 'x2', yaxis: 'y2'
        });

        // OPERATING LINE STEP: vertical UP from B(y_vap, y_vap) → C(y_vap, y_vap_next)
        // y_vap_next > y_vap so this goes UPWARD
        traces.push({
            x: [y_vap, y_vap], y: [y_vap, y_vap_next],
            mode: 'lines',
            line: { color: color, width: 2.5 },
            showlegend: false, legendgroup: `stage_${stageNum}`,
            xaxis: 'x2', yaxis: 'y2'
        });
    }

    // ========== LAYOUT ==========
    const containerWidth = document.querySelector('.main-content')?.clientWidth || 1050;

    const layout = {
        title: {
            text: '<b>Ponchon–Savarit Diagram: Ethanol-Water System</b>',
            font: { size: 18, family: 'Arial, sans-serif' },
            x: 0.5, xanchor: 'center'
        },
        annotations: [{
            text: '<b>Ponchon–Savarit Diagram (H-x-y)</b>',
            xref: 'paper', yref: 'paper',
            x: 0.5, y: 0.995,
            xanchor: 'center', yanchor: 'top',
            showarrow: false,
            font: { size: 12, color: '#1565C0' }
        }],
        grid: { rows: 2, columns: 1, pattern: 'independent', roworder: 'top to bottom' },
        margin: { l: 80, r: 220, t: 70, b: 70 },

        xaxis: {
            domain: [0.0, 0.95], anchor: 'y',
            range: [0, 1], tickformat: '.2f',
            showline: true, linecolor: '#555', mirror: true,
            showgrid: true, gridcolor: '#E0E0E0',
            zeroline: true, zerolinecolor: '#333', zerolinewidth: 1.5,
            title: ''
        },
        yaxis: {
            domain: [0.52, 0.97], anchor: 'x',
            title: { text: `<b>Enthalpy (${enthalpyUnit})</b>`, font: { size: 13 } },
            range: [results.yMin, results.yMax],
            showline: true, linecolor: '#555', mirror: true,
            showgrid: true, gridcolor: '#E0E0E0',
            zeroline: true, zerolinecolor: '#333', zerolinewidth: 1
        },
        xaxis2: {
            domain: [0.0, 0.95], anchor: 'y2',
            range: [0, 1], tickformat: '.2f',
            showline: true, linecolor: '#555', mirror: true,
            showgrid: true, gridcolor: '#E0E0E0',
            title: { text: '<b>Mole Fraction (x or y)</b>', font: { size: 13 } }
        },
        yaxis2: {
            domain: [0.03, 0.48], anchor: 'x2',
            title: { text: '<b>y (Vapor Fraction)</b>', font: { size: 13 } },
            range: [0, 1.02],
            showline: true, linecolor: '#555', mirror: true,
            showgrid: true, gridcolor: '#E0E0E0',
            zeroline: true, zerolinecolor: '#333', zerolinewidth: 1
        },
        legend: {
            x: 1.01, y: 1.0, xanchor: 'left', yanchor: 'top',
            font: { size: 10, family: 'Arial, sans-serif' },
            bgcolor: 'rgba(255,255,255,0.95)',
            bordercolor: '#BDBDBD', borderwidth: 1,
            tracegroupgap: 1
        },
        plot_bgcolor: '#FAFAFA',
        paper_bgcolor: '#FFFFFF',
        height: 780,
        width: containerWidth - 40
    };

    Plotly.newPlot('plotDiv', traces, layout, { responsive: true });

    window.addEventListener('resize', () => {
        const w = document.querySelector('.main-content')?.clientWidth || 1000;
        Plotly.relayout('plotDiv', { width: w - 40 });
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








