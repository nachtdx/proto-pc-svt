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

    while stages < 100:
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
                tie_lines.append({'x': [xB, y], 'y': [HLxB, cubic_interpolate(y, data['yData'], data['Hv'])]})
                stage_compositions.append({'x': xB, 'y': y})
                stages += 1
            break

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
            HEnd = cubic_interpolate(yNext, data['yData'], data['Hv'])
        
        construction_lines.append({'x': [xDelta, xEnd], 'y': [HDelta, HEnd]})
        y = yNext

    if stages == 0:
        error = "Failed to calculate stages. Check input data."
    
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
        return {'error': "Interpolation failed."}
    
    D = F * (zF - xB) / (xD - xB)
    W = F - D
    
    HVxD = cubic_interpolate(xD, data['yData'], data['Hv'])
    Qc = D * (HVxD - HD) * (R + 1)
    QcKW = Qc * 0.27778
    xDeltaR = xD
    HDeltaR = HD + Qc / D
    
    xDeltaS = xB
    slope = (HDeltaR - HF) / (xDeltaR - zF)
    HDeltaS = HF + slope * (xDeltaS - zF)
    
    Qr = W * (HW - HDeltaS)
    QrKW = Qr * 0.27778
    
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
        'QcKW': round(QcKW, 2), 'QrKW': round(QrKW, 2),
        'QPrimeMin': round(QPrimeMin, 2), 'QDoublePrimeMin': round(QDoublePrimeMin, 2),
        'RMin': round(RMin, 2), 'stages': stage_results['stages'],
        'feed_stage': stage_results['feed_stage'],
        'stage_compositions': stage_results['stage_compositions'],
        'x_range': x_range, 'HL_curve': HL_curve, 'HV_curve': HV_curve,
        'y_equilibrium': y_equilibrium, 'yMin': yMin, 'yMax': yMax,
        'tie_lines': stage_results['tie_lines'],
        'construction_lines': stage_results['construction_lines'],
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
// CREATE PLOT - PONCHON-SAVARIT (SUBPOLOT MEPET + PROYEKSI NYAMBUNG)
// ==========================================
function createPlot(results) {
    const stageColors = [
        '#FF6B6B', '#4ECDC4', '#FF9F1C', '#6A4C93', '#2E86AB', 
        '#A23B72', '#F18F01', '#2D6A4F', '#9E2A2B', '#540D6E'
    ];
    
    const traces = [];
    
    // ========== KURVA DASAR H-x-y ==========
    traces.push({
        x: results.x_range,
        y: results.HL_curve,
        mode: 'lines',
        name: 'Saturated Liquid',
        line: {color: '#2E86AB', width: 4},
        legendgroup: 'liquid',
        xaxis: 'x',
        yaxis: 'y'
    });
    
    traces.push({
        x: results.x_range,
        y: results.HV_curve,
        mode: 'lines',
        name: 'Saturated Vapor',
        line: {color: '#A23B72', width: 4},
        legendgroup: 'vapor',
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
        marker: {color: '#F97316', size: 14, symbol: 'star', line: {color: 'white', width: 1}},
        text: ['Δ<sub>R</sub>', 'Δ<sub>S</sub>'],
        textposition: ['top center', 'bottom center'],
        textfont: {size: 14, color: '#F97316'},
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
        marker: {size: 8, color: '#0A9396', line: {color: 'white', width: 1}},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // ========== MINIMUM REFLUX LINE ==========
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
    
    // ========== MINIMUM DIFFERENCE POINTS ==========
    traces.push({
        x: [results.xD, results.xB],
        y: [results.QPrimeMin, results.QDoublePrimeMin],
        mode: 'markers+text',
        name: 'Min Diff Points',
        marker: {color: '#E9C46A', size: 12, symbol: 'star-diamond', line: {color: 'white', width: 1}},
        text: ['Δ<sub>R,min</sub>', 'Δ<sub>S,min</sub>'],
        textposition: ['middle left', 'middle right'],
        textfont: {size: 12, color: '#E9C46A'},
        xaxis: 'x',
        yaxis: 'y'
    });
    
    // ========== CONSTRUCTION LINES ==========
    if (results.construction_lines?.length) {
        results.construction_lines.forEach((line, i) => {
            traces.push({
                x: line.x,
                y: line.y,
                mode: 'lines',
                name: i === 0 ? 'Construction Lines' : undefined,
                line: {color: '#B0B0B0', width: 1.5, dash: 'dot'},
                showlegend: i === 0,
                xaxis: 'x',
                yaxis: 'y'
            });
        });
    }
    
    // ========== STAGE TIE LINES ==========
    results.tie_lines.forEach((tie, i) => {
        const color = stageColors[i % stageColors.length];
        traces.push({
            x: tie.x,
            y: tie.y,
            mode: 'lines',
            name: `Stage ${i+1}`,
            line: {color, width: 3},
            legendgroup: `stage_${i+1}`,
            xaxis: 'x',
            yaxis: 'y'
        });
        traces.push({
            x: tie.x,
            y: tie.y,
            mode: 'markers',
            showlegend: false,
            marker: {color, size: 10, symbol: ['circle', 'diamond'], line: {color: 'white', width: 1}},
            legendgroup: `stage_${i+1}`,
            xaxis: 'x',
            yaxis: 'y'
        });
    });
    
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
    
    // ========== VLE TRACING ==========
    results.stage_compositions.forEach((stage, i) => {
        const color = stageColors[i % stageColors.length];
        traces.push({
            x: [stage.x, stage.y],
            y: [stage.y, stage.y],
            mode: 'lines',
            line: {color, width: 2.5},
            showlegend: false,
            legendgroup: `stage_${i+1}`,
            xaxis: 'x2',
            yaxis: 'y2'
        });
        traces.push({
            x: [stage.x, stage.y],
            y: [stage.y, stage.y],
            mode: 'markers',
            showlegend: false,
            marker: {color, size: 12, symbol: ['circle', 'diamond'], line: {color: 'white', width: 1.5}},
            hovertemplate: `<b>Stage ${i+1}</b><br>Liquid: x = %{x[0]:.3f}<br>Vapor: y = %{y[1]:.3f}<extra></extra>`,
            legendgroup: `stage_${i+1}`,
            xaxis: 'x2',
            yaxis: 'y2'
        });
    });
    
    // ========== GARIS PROYEKSI ==========
   // ========== GARIS PROYEKSI VERTIKAL (MULAI DARI KURVA) ==========
// ========== GARIS PROYEKSI VERTIKAL (MULAI DARI KURVA) ==========
results.stage_compositions.forEach((stage, i) => {
    const color = stageColors[i % stageColors.length];
    const x_liq = stage.x;
    const y_liq = stage.y;
    
    // Cari enthalpy yang sesuai di KURVA
    const idxLiq = Math.round(x_liq * 199);
    const idxVap = Math.round(y_liq * 199);
    const H_liq_stage = results.HL_curve[idxLiq];
    const H_vap_stage = results.HV_curve[idxVap];
    
    if (H_liq_stage && H_vap_stage) {
        // Proyeksi liquid
        traces.push({
            x: [x_liq, x_liq],
            y: [H_liq_stage, results.yMin],
            mode: 'lines',
            showlegend: false,
            line: {color, width: 1.8, dash: 'dot'},
            xaxis: 'x',
            yaxis: 'y'
        });
        traces.push({
            x: [x_liq, x_liq],
            y: [1.0, y_liq],
            mode: 'lines',
            showlegend: false,
            line: {color, width: 1.8, dash: 'dot'},
            xaxis: 'x2',
            yaxis: 'y2'
        });
        
        // Proyeksi vapor
        traces.push({
            x: [y_liq, y_liq],
            y: [H_vap_stage, results.yMin],
            mode: 'lines',
            showlegend: false,
            line: {color, width: 1.8, dash: 'dot'},
            xaxis: 'x',
            yaxis: 'y'
        });
        traces.push({
            x: [y_liq, y_liq],
            y: [1.0, y_liq],
            mode: 'lines',
            showlegend: false,
            line: {color, width: 1.8, dash: 'dot'},
            xaxis: 'x2',
            yaxis: 'y2'
        });
    }
});
    
    // ========== LAYOUT ==========
// ========== LAYOUT - SUBPLOT MEPET + TITLE X ATAS DIHAPUS ==========
// ========== LAYOUT - SUBPLOT CONCATENATED ==========
const layout = {
    title: {
        text: '<b>Ponchon–Savarit Diagram: Binary Distillation Analysis</b>',
        font: {size: 18, family: 'Arial', color: '#1E1E1E'},
        x: 0.5,
        y: 0.98  // Title agak naik
    },
    grid: {
        rows: 2,
        columns: 1,
        pattern: 'independent',
        roworder: 'top to bottom'
    },
    margin: {
        l: 70,   // kiri
        r: 130,  // kanan (buat legend)
        t: 50,   // atas dikit
        b: 70,   // bawah agak besar biar ga nutup
        pad: 0
    },
    xaxis: {
        domain: [0.1, 0.9],
        title: '',  // No title
        range: [0, 1],
        tickformat: '.2f',
        tickfont: {size: 10},
        showline: true,
        linecolor: '#1E1E1E',
        mirror: true
    },
    yaxis: {
        domain: [0.5, 0.95],  // Subplot atas (50% ke atas)
        title: '<b>Enthalpy (MJ/kmol)</b>',
        range: [results.yMin, results.yMax],
        tickfont: {size: 10},
        titlefont: {size: 12},
        showline: true,
        linecolor: '#1E1E1E',
        mirror: true
    },
    xaxis2: {
        domain: [0.1, 0.9],
        title: '<b>Mole Fraction (x or y)</b>',
        range: [0, 1],
        tickformat: '.2f',
        tickfont: {size: 10},
        titlefont: {size: 12},
        showline: true,
        linecolor: '#1E1E1E',
        mirror: true
    },
    yaxis2: {
        domain: [0.05, 0.45],  // Subplot bawah (5% ke 45%)
        title: '<b>y (Vapor Fraction)</b>',
        range: [0, 1],
        tickfont: {size: 10},
        titlefont: {size: 12},
        showline: true,
        linecolor: '#1E1E1E',
        mirror: true
    },
    // Legend di kanan
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
    // Ukuran plot
    height: 800,
    width: document.querySelector('.main-content')?.clientWidth - 40 || 1000
};
    Plotly.newPlot('plotDiv', traces, layout, {responsive: true});
    
    window.addEventListener('resize', () => {
        const containerWidth = document.querySelector('.main-content')?.clientWidth || 900;
        Plotly.relayout('plotDiv', {width: containerWidth - 40});
    });
}

// ==========================================
// DISPLAY RESULTS
// ==========================================
function displayResults(results) {
    createPlot(results);
    
    const summaryHtml = `
        <tr><td>Distillate Flow Rate (D)</td><td>${results.D} kmol/hr</td></tr>
        <tr><td>Bottoms Flow Rate (W)</td><td>${results.W} kmol/hr</td></tr>
        <tr><td>Δ_R</td><td>(${results.xDeltaR}, ${results.HDeltaR} MJ/kmol)</td></tr>
        <tr><td>Δ_S</td><td>(${results.xDeltaS}, ${results.HDeltaS} MJ/kmol)</td></tr>
        <tr><td>Condenser Duty (Qc)</td><td>${results.QcKW} kW</td></tr>
        <tr><td>Reboiler Duty (Qr)</td><td>${results.QrKW} kW</td></tr>
        <tr><td>Δ_R min</td><td>(${results.xD}, ${results.QPrimeMin} MJ/kmol)</td></tr>
        <tr><td>Δ_S min</td><td>(${results.xB}, ${results.QDoublePrimeMin} MJ/kmol)</td></tr>
        <tr><td>Minimum Reflux Ratio</td><td>${results.RMin}</td></tr>
        <tr><td>Number of Stages</td><td>${results.stages}</td></tr>
        <tr><td>Feed Stage</td><td>${results.feed_stage}</td></tr>
    `;
    document.getElementById('summaryBody').innerHTML = summaryHtml;
    
    let stagesRows = '';
    results.stage_compositions.forEach((stage, i) => {
        stagesRows += `<tr><td>Stage ${i+1}</td><td>${stage.x.toFixed(4)}</td><td>${stage.y.toFixed(4)}</td></tr>`;
    });
    document.getElementById('stagesBody').innerHTML = stagesRows;
    
    document.getElementById('exportBtn').disabled = false;
    currentResults = results;
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
    
    loadingDiv.style.display = 'flex';
    
    try {
        statusText.innerText = 'Loading Pyodide core...';
        progressBar.style.width = '20%'; progressBar.innerText = '20%';
        
        pyodide = await loadPyodide({indexURL: "https://cdn.jsdelivr.net/pyodide/v0.23.4/full/"});
        
        statusText.innerText = 'Loading numpy...';
        progressBar.style.width = '40%'; progressBar.innerText = '40%';
        await pyodide.loadPackage('numpy');
        
        statusText.innerText = 'Loading scipy...';
        progressBar.style.width = '60%'; progressBar.innerText = '60%';
        await pyodide.loadPackage('scipy');
        
        statusText.innerText = 'Loading pandas...';
        progressBar.style.width = '80%'; progressBar.innerText = '80%';
        await pyodide.loadPackage('pandas');
        
        statusText.innerText = 'Initializing calculator...';
        progressBar.style.width = '90%'; progressBar.innerText = '90%';
        await loadCalculatorCode();
        
        progressBar.style.width = '100%'; progressBar.innerText = '100%';
        statusText.innerText = 'Ready!';
        statusDiv.innerHTML = '<i class="fas fa-check-circle text-success"></i> Python ready!';
        pyodideReady = true;
        calculateBtn.disabled = false;
        
        setTimeout(() => loadingDiv.style.display = 'none', 500);
        
    } catch (error) {
        console.error(error);
        statusDiv.innerHTML = `<i class="fas fa-exclamation-circle text-danger"></i> Error: ${error.message}`;
        loadingDiv.style.display = 'none';
    }
}

// ==========================================
// EVENT LISTENERS
// ==========================================
['xData', 'yData', 'Hl', 'Hv'].forEach(id => {
    document.getElementById(id)?.addEventListener('input', updatePreview);
});

document.getElementById('q').addEventListener('change', function() {
    document.getElementById('customQDiv').style.display = 
        this.value === 'custom' ? 'block' : 'none';
});

document.getElementById('calculateBtn').addEventListener('click', async function() {
    const loading = document.getElementById('loading');
    const btn = this;
    
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
        alert('Error: Semua array harus sama panjang!'); return;
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
    loading.style.display = 'block';
    
    try {
        const results = await runCalculation(inputData);
        results.error ? alert('Error: ' + results.error) : displayResults(results);
    } catch (error) {
        alert('Error: ' + error.message);
    } finally {
        btn.disabled = false;
        loading.style.display = 'none';
    }
});

document.getElementById('exportBtn').addEventListener('click', function() {
    if (!currentResults) return;
    let csv = 'Stage,x (Liquid),y (Vapor)\n';
    currentResults.stage_compositions.forEach((stage, i) => {
        csv += `${i+1},${stage.x.toFixed(4)},${stage.y.toFixed(4)}\n`;
    });
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'ponchon_savarit_results.csv'; a.click();
    window.URL.revokeObjectURL(url);
});

updatePreview();
initPyodide();







