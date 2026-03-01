// ==========================================
// STATE & VARIABLES GLOBAL
// ==========================================
let pyodide = null;
let pyodideReady = false;
let currentResults = null;
let currentSystem = 'user-defined';

// ==========================================
// DATASETS (Ethanol-Water + Ammonia-Water)
// ==========================================
const DATASETS = {
    'ethanol-1atm': {
        name: 'Ethanol-Water at 1 atm (Faust)',
        // Hl (BTU/lbmol): liquid enthalpy at bubble point, ref = liquid at 32°F
        //   HL(x=0.05)=3200, HL(x=0.8)=3611 consistent with Faust Fig 3.4 / handwritten soln
        Hl: [3240, 3200,  3280,  3350,  3400,  3440,  3480,  3530,  3570,  3611,  3750,  4215],
        // Hv (BTU/lbmol): vapor enthalpy at bubble point, indexed by x (liquid composition)
        //   i.e. Hv[i] = enthalpy of saturated vapor in equilibrium with liquid x[i]
        //   Hv is plotted vs y (vapor composition) in the H-x-y diagram
        //   Key: Hv(y=0.654)=22227, Hv(y=0.8)=20939 → Rmin=1.1 for zF=0.5, xD=0.8, xB=0.05
        Hv: [21240,23500, 23200, 22900, 22600, 22450, 22227, 21700, 21500, 20726, 19000, 13200],
        T:  [212.0,196.2, 190.2, 187.0, 184.5, 182.8, 181.3, 180.1, 178.5, 176.5, 173.1, 172.6],
        units: { enthalpy: 'BTU/lbmol', temperature: '°F', duty: 'BTU/hr' },
        azeotrope: { x: 0.8943, T: 172.6 }
    },
    'ethanol-76mmHg': {
        name: 'Ethanol-Water at 76 mmHg (Faust)',
        x:  [0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 1.0],
        y:  [0, 0.192, 0.377, 0.527, 0.713, 0.746, 0.771, 0.794, 0.822, 0.912, 0.942, 0.959, 0.978, 0.990, 1.0],
        Hl: [1037, 1019, 1002, 974, 928, 893, 865, 842, 819, 802, 791, 779, 773, 768, 762],
        Hv: [2778, 2745, 2704, 2635, 2544, 2478, 2426, 2380, 2338, 2302, 2274, 2248, 2234, 2225, 2221],
        T:  [212.0, 203.4, 197.2, 189.2, 184.5, 181.7, 179.6, 177.8, 176.2, 174.3, 173.0, 172.8, 172.7, 172.8, 173.0],
        units: { enthalpy: 'BTU/lbmol', temperature: '°F', duty: 'BTU/hr' },
        azeotrope: { x: 0.86, T: 172.8 }
    },
    'ammonia-1atm': {
        name: 'Ammonia-Water at 100 psia (Faust Fig 3.5)',
        // ── VLE data (Perry / Faust Fig 3.5, 100 psia) ──
        // Ammonia is highly volatile: y >> x
        x:  [0,      0.02,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.40,  0.50,  0.60,  0.70,  0.80,  0.90,  1.0 ],
        y:  [0,      0.250, 0.440, 0.590, 0.660, 0.710, 0.752, 0.790, 0.851, 0.898, 0.932, 0.960, 0.979, 0.992, 1.0 ],
        // ── Hl: saturated LIQUID enthalpy (BTU/lbmol) ──
        // Ref: liquid water at 32°F, liquid ammonia at -40°F
        // Goes NEGATIVE in midrange due to large negative enthalpy of NH3 dissolution
        // HL(x=0) ~ bubble-pt water at 327°F sensible heat from 32°F
        // HL minimum ~ -5800 at x≈0.4
        Hl: [1000,  -300,  -1800, -3800, -5000, -5500, -5750, -5800, -5500, -4600, -3400, -2000, -600,   300,   800 ],
        // ── Hv: saturated VAPOR enthalpy (BTU/lbmol), indexed by x (liquid composition) ──
        // Hv[i] = enthalpy of saturated vapor in equilibrium with liquid x[i]
        // Plotted vs y (vapor composition) in H-x-y diagram
        // HV(x=0) ~ pure water steam at 327°F, 100 psia ≈ 20200 BTU/lbmol
        // HV(x=1) ~ pure NH3 vapor at 82°F, 100 psia ≈ 4800 BTU/lbmol
        Hv: [20200, 19800, 19000, 17800, 16800, 15800, 14900, 14000, 12300, 10700, 9200,  7900,  6700,  5600,  4800],
        // ── Bubble point temperatures (°F) at 100 psia ──
        T:  [327,   315,   295,   270,   252,   237,   223,   211,   190,   172,   157,   143,   130,   107,    82  ],
        units: { enthalpy: 'BTU/lbmol', temperature: '°F', duty: 'BTU/hr' },
        azeotrope: null
    }
};

// Backward-compat alias
const ETHANOL_WATER_DATA = {
    '1atm':   DATASETS['ethanol-1atm'],
    '76mmHg': DATASETS['ethanol-76mmHg']
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
    } catch (e) { return []; }
}

// ==========================================
// VALIDATE ARRAY
// ==========================================
function validateArray(arr, name) {
    if (!Array.isArray(arr) || arr.length === 0) return { valid: false, error: `${name} tidak boleh kosong` };
    if (arr.some(isNaN)) return { valid: false, error: `${name} harus berupa angka` };
    if (name.includes('x') || name.includes('y')) {
        if (arr.some(v => v < 0 || v > 1)) return { valid: false, error: `${name} harus antara 0 dan 1` };
    }
    return { valid: true };
}

// ==========================================
// SHOW TOAST
// ==========================================
function showToast(message, type = 'info') {
    const toastContainer = document.getElementById('toastContainer');
    if (!toastContainer) return;
    const toast = document.createElement('div');
    toast.className = `toast align-items-center text-white bg-${type} border-0`;
    toast.setAttribute('role', 'alert');
    toast.setAttribute('aria-atomic', 'true');
    toast.innerHTML = `<div class="d-flex"><div class="toast-body">${message}</div><button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button></div>`;
    toastContainer.appendChild(toast);
    const bsToast = new bootstrap.Toast(toast, { delay: 3000 });
    bsToast.show();
    toast.addEventListener('hidden.bs.toast', () => toast.remove());
}

// ==========================================
// UPDATE UNIT DISPLAY
// ==========================================
function updateUnitDisplay(systemType) {
    const isBritish = systemType !== 'user-defined';
    const hlUnit    = document.getElementById('hlUnit');
    const hvUnit    = document.getElementById('hvUnit');
    const previewUnitBadge = document.getElementById('previewUnitBadge');
    const unitNote  = document.getElementById('unitNote');
    const flowUnit  = document.getElementById('flowUnit');
    const flowNote  = document.getElementById('flowNote');
    const dutyUnit  = document.getElementById('dutyUnit');
    const unitDisplay = document.getElementById('unitDisplay');
    const footerUnit  = document.getElementById('footerUnit');

    // Determine label
    const isAmmonia  = systemType.includes('ammonia');
    const isEthanol  = systemType.includes('ethanol');
    const systemLabel = isAmmonia ? 'Ammonia-Water 100 psia (Faust Fig 3.5)' : isEthanol ? 'Ethanol-Water 1 atm (Faust Fig 3.4)' : 'User Defined';

    if (unitDisplay) unitDisplay.textContent = systemLabel;
    if (footerUnit)  footerUnit.textContent  = systemLabel;

    if (isBritish) {
        if (hlUnit) { hlUnit.textContent = 'BTU/lbmole'; hlUnit.className = 'unit-badge bg-primary text-white'; }
        if (hvUnit) { hvUnit.textContent = 'BTU/lbmole'; hvUnit.className = 'unit-badge bg-primary text-white'; }
        if (previewUnitBadge) { previewUnitBadge.textContent = 'BTU/lbmole'; previewUnitBadge.className = 'unit-badge bg-primary text-white'; }
        if (unitNote)  unitNote.innerHTML  = '<b>British Units:</b> Enthalpy in BTU/lbmole, Flow in lbmol/hr, Duty in BTU/hr';
        if (flowUnit)  flowUnit.textContent = 'lbmol/hr';
        if (flowNote)  flowNote.innerHTML  = 'Flow rate in lbmol/hr (British units)';
        if (dutyUnit)  dutyUnit.textContent = 'BTU/hr';

        document.getElementById('xData').readOnly = true;
        document.getElementById('yData').readOnly = true;
        document.getElementById('Hl').readOnly = true;
        document.getElementById('Hv').readOnly = true;
    } else {
        if (hlUnit) { hlUnit.textContent = 'any units'; hlUnit.className = 'unit-badge'; }
        if (hvUnit) { hvUnit.textContent = 'any units'; hvUnit.className = 'unit-badge'; }
        if (previewUnitBadge) { previewUnitBadge.textContent = 'any units'; previewUnitBadge.className = 'unit-badge'; }
        if (unitNote)  unitNote.innerHTML  = 'User Defined: bebas menggunakan satuan apapun';
        if (flowUnit)  flowUnit.textContent = 'kmol/hr';
        if (flowNote)  flowNote.innerHTML  = 'Flow rate in kmol/hr (SI units)';
        if (dutyUnit && dutyUnit.textContent) dutyUnit.textContent = 'kW';

        document.getElementById('xData').readOnly = false;
        document.getElementById('yData').readOnly = false;
        document.getElementById('Hl').readOnly = false;
        document.getElementById('Hv').readOnly = false;
    }

    // Update enthalpy table header
    const hlHeader = document.getElementById('hlHeader');
    const hvHeader = document.getElementById('hvHeader');
    const unit = isBritish ? 'BTU/lbmole' : 'any units';
    if (hlHeader) hlHeader.textContent = `HL (${unit})`;
    if (hvHeader) hvHeader.textContent = `HV (${unit})`;
}

// ==========================================
// LOAD DATASET (generic)
// ==========================================
function loadDataset(key) {
    const dataset = DATASETS[key];
    if (!dataset) return;

    console.log(`Loading dataset: ${dataset.name}`);

    document.getElementById('xData').value = JSON.stringify(dataset.x);
    document.getElementById('yData').value = JSON.stringify(dataset.y);
    document.getElementById('Hl').value    = JSON.stringify(dataset.Hl);
    document.getElementById('Hv').value    = JSON.stringify(dataset.Hv);
    document.getElementById('TData').value = JSON.stringify(dataset.T);

    // Info panel
    const datasetInfoPanel = document.getElementById('datasetInfoPanel');
    const datasetName      = document.getElementById('datasetName');
    const datasetDetails   = document.getElementById('datasetDetails');

    if (datasetInfoPanel) {
        datasetInfoPanel.style.display = 'block';
        datasetInfoPanel.className = 'dataset-info ' + (key.includes('ammonia') ? 'system-ammonia' : 'system-ethanol');
    }
    if (datasetName) datasetName.textContent = dataset.name;

    const azeoInfo = dataset.azeotrope
        ? ` | Azeotrope: x = ${dataset.azeotrope.x.toFixed(3)} at ${dataset.azeotrope.T.toFixed(1)}°F`
        : ' | No azeotrope';

    if (datasetDetails) {
        datasetDetails.innerHTML = `
            <i class="fas fa-database me-1"></i> ${dataset.x.length} points
            <i class="fas fa-fire ms-2 me-1"></i> Hl: ${Math.min(...dataset.Hl)}–${Math.max(...dataset.Hl)} BTU/lbmole
            <i class="fas fa-temperature-high ms-2 me-1"></i> T: ${Math.min(...dataset.T).toFixed(1)}–${Math.max(...dataset.T).toFixed(1)}°F
            ${azeoInfo}
        `;
    }

    const systemInfo     = document.getElementById('systemInfo');
    const systemInfoText = document.getElementById('systemInfoText');
    if (systemInfo) systemInfo.style.display = 'block';
    if (systemInfoText) {
        const emoji = key.includes('ammonia') ? '🧪' : '🍸';
        systemInfoText.innerHTML = `
            ${emoji} <b>${dataset.name}</b><br>
            📊 ${dataset.x.length} data points<br>
            🔥 Enthalpy: <b>${dataset.units.enthalpy}</b><br>
            🌡️ Temperature: <b>${dataset.units.temperature}</b>
        `;
    }

    updateUnitDisplay(key);
    updatePreview();
    showToast(`✅ Loaded: ${dataset.name}`, 'success');
}

// Backward compat
function loadEthanolWaterDataset(pressure) {
    loadDataset(`ethanol-${pressure}`);
}

// ==========================================
// UPDATE PREVIEW TABLE
// ==========================================
function updatePreview() {
    const xData = parseArrayString(document.getElementById('xData').value);
    const yData = parseArrayString(document.getElementById('yData').value);
    const Hl    = parseArrayString(document.getElementById('Hl').value);
    const Hv    = parseArrayString(document.getElementById('Hv').value);

    const lengths   = [xData.length, yData.length, Hl.length, Hv.length];
    const maxLength = Math.max(...lengths);
    const systemType = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish  = systemType !== 'user-defined';

    let html = '';
    for (let i = 0; i < maxLength; i++) {
        html += '<tr>';
        html += `<td>${i + 1}</td>`;
        html += `<td>${i < xData.length ? xData[i].toFixed(3) : '-'}</td>`;
        html += `<td>${i < yData.length ? yData[i].toFixed(3) : '-'}</td>`;
        html += `<td>${i < Hl.length ? (isBritish ? Hl[i].toFixed(0) : Hl[i].toFixed(2)) : '-'}</td>`;
        html += `<td>${i < Hv.length ? (isBritish ? Hv[i].toFixed(0) : Hv[i].toFixed(2)) : '-'}</td>`;
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
// LOAD PYTHON CALCULATOR CODE
// ==========================================
async function loadCalculatorCode() {
    const pythonCode = `
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import root_scalar
import json

# ════════════════════════════════════════════════════════════════════
#  INTERPOLATION HELPERS
# ════════════════════════════════════════════════════════════════════

def linear_interpolate(x, x_val, y_val):
    x_val = list(x_val); y_val = list(y_val)
    if len(x_val) != len(y_val) or len(x_val) < 2:
        return float('nan')
    if x <= x_val[0]:  return float(y_val[0])
    if x >= x_val[-1]: return float(y_val[-1])
    f = interp1d(x_val, y_val, kind='linear', fill_value='extrapolate')
    return float(f(x))

def cubic_interpolate(x, x_val, y_val):
    """Cubic spline; falls back to linear on any failure."""
    x_val = list(x_val); y_val = list(y_val)
    if len(x_val) != len(y_val) or len(x_val) < 2:
        return float('nan')
    xc = float(np.clip(x, x_val[0], x_val[-1]))
    try:
        cs = CubicSpline(x_val, y_val, extrapolate=False)
        v  = float(cs(xc))
        return v if np.isfinite(v) else linear_interpolate(x, x_val, y_val)
    except Exception:
        return linear_interpolate(x, x_val, y_val)

def safe_hv(y, yData, Hv):
    """HV at vapor composition y.  Hv is indexed by yData (vapor mole fraction).
    Linear interpolation for stability."""
    return linear_interpolate(y, yData, Hv)

def equil_y(x, xData, yData):
    """Equilibrium y at liquid x (cubic -> linear fallback)."""
    v = cubic_interpolate(x, xData, yData)
    if not np.isfinite(v):
        v = linear_interpolate(x, xData, yData)
    return v

# ════════════════════════════════════════════════════════════════════
#  MINIMUM REFLUX  (Ponchon-Savarit, handles all VLE shapes)
# ════════════════════════════════════════════════════════════════════

def calc_rmin(xData, yData, Hl, Hv, zF, xD, xB, q):
    """
    Correct Ponchon-Savarit Rmin.

    Method: scan all equilibrium tie-lines from zF to xD.
    For each tie-line at liquid composition xi:
      - Tie-line connects (xi, HL_xi) on liquid curve to (yi, HV_yi) on vapor curve
      - Slope: sl = (HV_yi - HL_xi) / (yi - xi)
      - Extend to x=xD: Q_prime = HL_xi + sl*(xD - xi)
      - Rmin candidate: Ri = (Q_prime - HV_xD) / (HV_xD - HL_xD)
    Take the MAXIMUM Ri (most constraining tie-line).
    For normal convex VLE: pinch is at feed tie-line (xi=zF).
    For concave/non-ideal VLE: tangent pinch may occur above feed.
    """
    yF   = equil_y(zF, xData, yData)
    HLzF = linear_interpolate(zF, xData, Hl)
    HVzF = safe_hv(yF, yData, Hv)
    HF   = q * HLzF + (1.0 - q) * HVzF

    HD   = linear_interpolate(xD, xData, Hl)
    HVxD = safe_hv(xD, yData, Hv)
    HW   = linear_interpolate(xB, xData, Hl)

    denom = HVxD - HD
    if abs(denom) < 1e-8:
        return 0.0, HVxD, HW, False

    # Scan tie-lines from zF to xD (only rectifying section)
    scan_xs = sorted(set(
        [float(x) for x in xData if zF - 1e-6 <= float(x) < xD - 1e-6] +
        list(np.linspace(zF, xD - 1e-5, 100)) +
        [float(zF)]
    ))

    best_R   = -1e10
    best_QP  = HVxD
    best_QDP = HW
    feed_is_pinch = True

    for xi in scan_xs:
        yi  = equil_y(xi, xData, yData)
        if yi <= xi + 1e-9 or yi > xD + 1e-6:
            continue
        HLi = linear_interpolate(xi, xData, Hl)
        HVi = safe_hv(yi, yData, Hv)
        if not (np.isfinite(HLi) and np.isfinite(HVi)):
            continue
        if abs(yi - xi) < 1e-9:
            continue
        sl   = (HVi - HLi) / (yi - xi)
        QP   = HLi + sl * (xD - xi)
        QDP  = HLi + sl * (xB - xi)
        Ri   = (QP - HVxD) / denom
        if Ri > best_R:
            best_R   = Ri
            best_QP  = QP
            best_QDP = QDP
            feed_is_pinch = (abs(xi - zF) < 0.02)

    if best_R < -1e-4:
        return 0.0, HVxD, HW, False

    tangent_pinch = not feed_is_pinch
    return max(0.0, best_R), best_QP, best_QDP, tangent_pinch

# ════════════════════════════════════════════════════════════════════
#  STAGE STEPPING
# ════════════════════════════════════════════════════════════════════

def calculate_stages(xD, xB, zF, HF, q, R, D, W,
                     xDeltaR, HDeltaR, xDeltaS, HDeltaS, data):
    xData = list(data['xData'])
    yData = list(data['yData'])
    Hl    = list(data['Hl'])
    Hv    = list(data['Hv'])

    y      = float(xD)
    stages = 0
    hv_xD  = safe_hv(xD, yData, Hv)
    stage_points       = [{'x': float(xD), 'y': float(hv_xD)}]
    tie_lines          = []
    construction_lines = []
    stage_compositions = []
    in_rectifying = True
    feed_stage    = 0
    error         = ""

    while stages < 100:
        # 1. Find x_n from y_n via equilibrium curve
        def f_equil(xi):
            return cubic_interpolate(xi, xData, yData) - y

        try:
            sol = root_scalar(f_equil, bracket=[1e-9, 1.0 - 1e-9], method='bisect')
            if not sol.converged:
                error = f"Stage {stages+1}: could not find x from y={y:.4f}"
                break
            x_n = float(sol.root)
        except Exception as e:
            error = f"Stage {stages+1}: {e}"
            break

        if not np.isfinite(x_n):
            break

        # 2. Check bottoms termination
        if x_n <= xB + 1e-6:
            HLxB = linear_interpolate(xB, xData, Hl)
            HVy  = safe_hv(y, yData, Hv)
            stage_points.append({'x': float(xB), 'y': float(HLxB)})
            tie_lines.append({'x': [float(xB), float(y)],
                              'y': [float(HLxB), float(HVy)]})
            stage_compositions.append({'x': float(xB), 'y': float(y)})
            stages += 1
            break

        # 3. Enthalpies at this stage
        HLx_n = linear_interpolate(x_n, xData, Hl)
        HVy_n = safe_hv(y, yData, Hv)

        if not (np.isfinite(HLx_n) and np.isfinite(HVy_n)):
            error = f"Stage {stages+1}: enthalpy interpolation failed"
            break

        stage_points.append({'x': float(x_n), 'y': float(HLx_n)})
        tie_lines.append({'x': [float(x_n), float(y)],
                          'y': [float(HLx_n), float(HVy_n)]})
        stage_compositions.append({'x': float(x_n), 'y': float(y)})
        stages += 1

        # 4. Switch sections at feed stage (BEFORE choosing difference point)
        if in_rectifying and x_n <= zF:
            in_rectifying = False
            feed_stage    = stages

        xDelta = xDeltaR if in_rectifying else xDeltaS
        HDelta = HDeltaR if in_rectifying else HDeltaS

        if abs(x_n - xDelta) < 1e-6:
            break

        slope = (HDelta - HLx_n) / (xDelta - x_n)
        if not np.isfinite(slope):
            break

        # 5. Find y_next: HV(y) = HDelta + slope*(y - xDelta)
        #    Scan full [0,1] range to find sign change, then bisect
        def f_y(yi):
            return safe_hv(yi, yData, Hv) - (HDelta + slope * (yi - xDelta))

        n_scan = 500
        ys_scan = np.linspace(1e-6, 1.0 - 1e-6, n_scan)
        fs_scan = [f_y(yi) for yi in ys_scan]

        bracket = None
        for k in range(len(fs_scan) - 1):
            if (np.isfinite(fs_scan[k]) and np.isfinite(fs_scan[k+1])
                    and fs_scan[k] * fs_scan[k+1] < 0
                    and ys_scan[k] < y - 1e-6):
                bracket = [float(ys_scan[k]), float(ys_scan[k+1])]
                break
        if bracket is None:
            for k in range(len(fs_scan) - 1):
                if (np.isfinite(fs_scan[k]) and np.isfinite(fs_scan[k+1])
                        and fs_scan[k] * fs_scan[k+1] < 0):
                    bracket = [float(ys_scan[k]), float(ys_scan[k+1])]
                    break

        if bracket is None:
            error = f"Stage {stages+1}: could not bracket y_next (y={y:.4f})"
            break

        try:
            sol = root_scalar(f_y, bracket=bracket, method='bisect')
            if not sol.converged:
                error = f"Stage {stages+1}: could not find y_next"
                break
            yNext = float(sol.root)
        except Exception as e:
            error = f"Stage {stages+1}: {e}"
            break

        HVyNext = safe_hv(yNext, yData, Hv)
        if not np.isfinite(HVyNext):
            break

        stage_points.append({'x': float(yNext), 'y': float(HVyNext)})

        # 6. Construction lines
        if in_rectifying:
            def f_cl(xi):
                return linear_interpolate(xi, xData, Hl)                        - (HDelta + slope * (xi - xDelta))
            try:
                sol = root_scalar(f_cl, bracket=[0.0, 1.0], method='bisect')
                if sol.converged:
                    xEnd = float(sol.root)
                    HEnd = linear_interpolate(xEnd, xData, Hl)
                    if np.isfinite(HEnd):
                        construction_lines.append({
                            'x': [float(xDelta), float(xEnd)],
                            'y': [float(HDelta),  float(HEnd)]
                        })
            except Exception:
                pass
        else:
            if np.isfinite(HVyNext):
                construction_lines.append({
                    'x': [float(xDelta), float(yNext)],
                    'y': [float(HDelta),  float(HVyNext)]
                })

        prev_y = y
        y = yNext

        if yNext >= prev_y - 1e-6 and stages > 1:
            error = f"Stage {stages}: y not decreasing ({prev_y:.4f} -> {yNext:.4f})"
            break

    if stages == 0 or len(stage_points) < 2:
        error = error or "Failed to calculate stages."

    return {
        'stages': stages, 'feed_stage': feed_stage,
        'stage_points': stage_points, 'tie_lines': tie_lines,
        'construction_lines': construction_lines,
        'stage_compositions': stage_compositions,
        'error': error
    }

# ════════════════════════════════════════════════════════════════════
#  MAIN CALCULATE
# ════════════════════════════════════════════════════════════════════

def calculate(data, params):
    zF = float(params['zF']); F  = float(params['F'])
    xD = float(params['xD']); xB = float(params['xB'])
    q  = float(params['q']);  R  = float(params['R'])

    xData = list(data['xData'])
    yData = list(data['yData'])

    xD = float(np.clip(xD, xData[0], xData[-1]))
    xB = float(np.clip(xB, xData[0], xData[-1]))
    zF = float(np.clip(zF, xData[0], xData[-1]))

    # Feed enthalpies
    yF   = equil_y(zF, xData, yData)
    HLzF = linear_interpolate(zF, xData, data['Hl'])
    HVzF = safe_hv(yF, yData, data['Hv'])
    HF   = q * HLzF + (1.0 - q) * HVzF
    HD   = linear_interpolate(xD, xData, data['Hl'])
    HW   = linear_interpolate(xB, xData, data['Hl'])

    bad = {k: v for k, v in
           [('yF',yF),('HLzF',HLzF),('HVzF',HVzF),('HF',HF),('HD',HD),('HW',HW)]
           if not np.isfinite(v)}
    if bad:
        return {'error': f'Interpolation failed for: {list(bad.keys())}'}

    # Material balance
    D = F * (zF - xB) / (xD - xB)
    W = F - D
    if not (np.isfinite(D) and np.isfinite(W) and D > 1e-9 and W > 1e-9):
        return {'error': f'Material balance error: D={D:.3f}, W={W:.3f}'}

    # Condenser & rectifying difference point
    HVxD = safe_hv(xD, yData, data['Hv'])
    if not np.isfinite(HVxD):
        HVxD = float(data['Hv'][-1])

    Qc      = D * (HVxD - HD) * (R + 1.0)
    xDeltaR = float(xD)
    HDeltaR = float(HD + Qc / D)
    if not (np.isfinite(Qc) and np.isfinite(HDeltaR)):
        return {'error': 'Condenser duty / Delta_R calculation failed'}

    # Stripping difference point & reboiler duty
    xDeltaS  = float(xB)
    slope_op = (HDeltaR - HF) / (xDeltaR - zF)
    HDeltaS  = float(HF + slope_op * (xDeltaS - zF))
    if not np.isfinite(HDeltaS):
        return {'error': 'Delta_S calculation failed'}
    Qr = W * (HW - HDeltaS)

    # Minimum reflux
    RMin, QPrimeMin, QDoublePrimeMin, tangent_pinch = calc_rmin(
        xData, yData, data['Hl'], data['Hv'], zF, xD, xB, q)

    if RMin > 1e-4 and R < RMin - 1e-4:
        return {'error': f'R = {R:.3f} is below Rmin = {RMin:.3f}. Increase R.'}

    # Stage stepping
    sr = calculate_stages(xD, xB, zF, HF, q, R, D, W,
                          xDeltaR, HDeltaR, xDeltaS, HDeltaS, data)
    if sr['error']:
        return {'error': sr['error']}

    # Curves for plot
    x_range  = np.linspace(0, 1, 500).tolist()
    HL_curve = [linear_interpolate(xi, xData, data['Hl']) for xi in x_range]
    HV_curve = [safe_hv(xi, yData, data['Hv'])            for xi in x_range]

    try:
        cs_eq = CubicSpline(xData, yData, extrapolate=False)
        y_eq  = [float(np.clip(float(cs_eq(xi)), 0.0, 1.0))
                 if np.isfinite(float(cs_eq(xi)))
                 else equil_y(xi, xData, yData)
                 for xi in x_range]
    except Exception:
        y_eq  = [equil_y(xi, xData, yData) for xi in x_range]

    all_H = (list(data['Hl']) + list(data['Hv']) +
             [HF, HD, HW, HDeltaR, HDeltaS, QPrimeMin, QDoublePrimeMin])
    all_H = [v for v in all_H if np.isfinite(v)]
    span  = max(all_H) - min(all_H)
    pad   = span * 0.10 + 1.0
    yMin  = min(all_H) - pad
    yMax  = max(all_H) + pad

    return {
        'D': round(D, 4),      'W': round(W, 4),
        'xDeltaR': round(xDeltaR, 4), 'HDeltaR': round(HDeltaR, 4),
        'xDeltaS': round(xDeltaS, 4), 'HDeltaS': round(HDeltaS, 4),
        'QcDuty':  round(Qc, 4),      'QrDuty':  round(Qr, 4),
        'QPrimeMin':       round(QPrimeMin, 4),
        'QDoublePrimeMin': round(QDoublePrimeMin, 4),
        'RMin':            round(RMin, 4),
        'tangent_pinch':   tangent_pinch,
        'stages':          sr['stages'],
        'feed_stage':      sr['feed_stage'],
        'stage_compositions': sr['stage_compositions'],
        'tie_lines':          sr['tie_lines'],
        'construction_lines': sr['construction_lines'],
        'x_range':    x_range,
        'HL_curve':   HL_curve,
        'HV_curve':   HV_curve,
        'y_equilibrium': y_eq,
        'yMin': round(yMin, 4), 'yMax': round(yMax, 4),
        'HF':  round(HF, 4),   'zF': zF, 'xD': xD, 'xB': xB,
        'yFMin': round(yF, 4), 'HVyF': round(HVzF, 4)
    }

def calculate_from_js(xData, yData, Hl, Hv, zF, F, xD, xB, q, R):
    try:
        data = {
            'xData': [float(v) for v in xData],
            'yData': [float(v) for v in yData],
            'Hl':    [float(v) for v in Hl],
            'Hv':    [float(v) for v in Hv],
        }
        params = {'zF': zF, 'F': F, 'xD': xD, 'xB': xB, 'q': q, 'R': R}
        return json.dumps(calculate(data, params))
    except Exception as e:
        import traceback
        return json.dumps({'error': str(e) + ' | ' + traceback.format_exc()})
`;
    pyodide.runPython(pythonCode);
    console.log('✅ Ponchon-Savarit calculator loaded');
}


// ==========================================
// RUN CALCULATION
// ==========================================
async function runCalculation(inputData) {
    if (!pyodideReady) throw new Error('Pyodide belum siap.');
    return JSON.parse(pyodide.runPython(`
        calculate_from_js(
            ${JSON.stringify(inputData.xData)}, ${JSON.stringify(inputData.yData)},
            ${JSON.stringify(inputData.Hl)},    ${JSON.stringify(inputData.Hv)},
            ${inputData.zF}, ${inputData.F}, ${inputData.xD},
            ${inputData.xB}, ${inputData.q}, ${inputData.R}
        )
    `));
}

// ==========================================
// CREATE PLOT - PONCHON-SAVARIT
// ==========================================
function createPlot(results) {
    const stageColors = [
        '#E63946', '#2196F3', '#E9A825', '#7B2D8B', '#00ACC1',
        '#388E3C', '#FF7043', '#8D6E63', '#EC407A', '#546E7A',
        '#26A69A', '#AB47BC'
    ];

    const traces = [];
    const systemType    = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish     = systemType !== 'user-defined';
    const enthalpyUnit  = isBritish ? 'BTU/lbmole' : 'MJ/kmol';

    // ── Saturated curves ──
    traces.push({ x: results.x_range, y: results.HL_curve, mode: 'lines', name: 'Saturated Liquid', line: { color: '#1565C0', width: 3 }, xaxis: 'x', yaxis: 'y' });
    traces.push({ x: results.x_range, y: results.HV_curve, mode: 'lines', name: 'Saturated Vapor',  line: { color: '#880E4F', width: 3 }, xaxis: 'x', yaxis: 'y' });

    // ── Vertical dashed lines xD, xB, zF ──
    const verticals = [
        { val: results.xD, label: 'x<sub>D</sub>',      color: '#9E9E9E' },
        { val: results.xB, label: 'x<sub>B</sub>',      color: '#9E9E9E' },
        { val: results.zF, label: 'z<sub>F</sub> (Feed)', color: '#2E7D32' }
    ];
    verticals.forEach(v => {
        traces.push({ x: [v.val, v.val], y: [results.yMin, results.yMax], mode: 'lines', name: v.label, line: { color: v.color, width: 2, dash: 'dash' }, xaxis: 'x', yaxis: 'y' });
        traces.push({ x: [v.val, v.val], y: [0, 1], mode: 'lines', showlegend: false, line: { color: v.color, width: 2, dash: 'dash' }, xaxis: 'x2', yaxis: 'y2' });
    });

    // ── Operating line ──
    traces.push({ x: [results.xDeltaR, results.zF, results.xDeltaS], y: [results.HDeltaR, results.HF, results.HDeltaS], mode: 'lines+markers', name: 'Operating Line', line: { color: '#00897B', width: 2.5 }, marker: { size: 7, color: '#00897B' }, xaxis: 'x', yaxis: 'y' });

    // ── Difference points ──
    traces.push({ x: [results.xDeltaR, results.xDeltaS], y: [results.HDeltaR, results.HDeltaS], mode: 'markers+text', name: 'Difference Points', marker: { color: '#FF6D00', size: 16, symbol: 'star' }, text: ['Δ<sub>R</sub>', 'Δ<sub>S</sub>'], textposition: ['top center', 'bottom center'], textfont: { size: 13, color: '#FF6D00' }, xaxis: 'x', yaxis: 'y' });

    // ── Minimum reflux line ──
    if (results.yFMin !== undefined && results.HVyF !== undefined) {
        traces.push({ x: [results.xB, results.zF, results.yFMin, results.xD], y: [results.QDoublePrimeMin, results.HF, results.HVyF, results.QPrimeMin], mode: 'lines+markers', name: 'Minimum Reflux Line', line: { color: '#F9A825', width: 2, dash: 'dash' }, marker: { size: 7, color: '#F9A825', symbol: 'diamond' }, xaxis: 'x', yaxis: 'y' });
        traces.push({ x: [results.xD, results.xB], y: [results.QPrimeMin, results.QDoublePrimeMin], mode: 'markers+text', name: 'Min Diff Points', marker: { color: '#F9A825', size: 12, symbol: 'diamond' }, text: ["Δ'<sub>R,min</sub>", "Δ'<sub>S,min</sub>"], textposition: ['top right', 'bottom left'], textfont: { size: 11, color: '#F9A825' }, xaxis: 'x', yaxis: 'y' });
    }

    // ── Construction lines ──
    if (results.construction_lines?.length) {
        results.construction_lines.forEach((cl, i) => {
            traces.push({ x: cl.x, y: cl.y, mode: 'lines', name: i === 0 ? 'Construction Lines' : undefined, showlegend: i === 0, line: { color: '#BDBDBD', width: 1.5, dash: 'dot' }, xaxis: 'x', yaxis: 'y' });
        });
    }

    // ── VLE Equilibrium Curve ──
    traces.push({ x: results.x_range, y: results.y_equilibrium, mode: 'lines', name: 'Equilibrium Curve', line: { color: '#1A1A1A', width: 3.5 }, xaxis: 'x2', yaxis: 'y2' });

    // ── y = x diagonal ──
    traces.push({ x: [0, 1], y: [0, 1], mode: 'lines', name: 'y = x', line: { color: '#9E9E9E', width: 1.5, dash: 'dash' }, xaxis: 'x2', yaxis: 'y2' });

    // ── Stage construction (tie lines + connectors + VLE stepping) ──
    const tieLines    = results.tie_lines || [];
    const numStages   = tieLines.length;
    const tieLinesAsc = [...tieLines].reverse(); // xB → xD order

    let connectorLegendAdded = false;

    for (let i = 0; i < numStages; i++) {
        const origIdx  = numStages - 1 - i;
        const tie      = tieLinesAsc[i];
        const color    = stageColors[origIdx % stageColors.length];
        const stageNum = origIdx + 1;

        const x_liq = tie.x[0]; // liquid composition
        const y_vap = tie.x[1]; // vapor composition
        const H_liq = tie.y[0];
        const H_vap = tie.y[1];

        // H-x-y: tie line
        traces.push({ x: [x_liq, y_vap], y: [H_liq, H_vap], mode: 'lines', name: `Stage ${stageNum}`, line: { color, width: 2.5 }, legendgroup: `stage_${stageNum}`, xaxis: 'x', yaxis: 'y' });
        traces.push({ x: [x_liq, y_vap], y: [H_liq, H_vap], mode: 'markers', marker: { color, size: 10, symbol: 'circle', line: { color: 'white', width: 1.5 } }, showlegend: false, legendgroup: `stage_${stageNum}`, xaxis: 'x', yaxis: 'y' });

        // VLE: equilibrium point A = (x_liq, y_vap)
        traces.push({ x: [x_liq], y: [y_vap], mode: 'markers', marker: { color, size: 10, symbol: 'circle', line: { color: 'white', width: 1.5 } }, showlegend: false, legendgroup: `stage_${stageNum}`, xaxis: 'x2', yaxis: 'y2' });

        // VLE: tie line horizontal RIGHT from (x_liq, y_vap) → (y_vap, y_vap)
        traces.push({ x: [x_liq, y_vap], y: [y_vap, y_vap], mode: 'lines', line: { color, width: 2.5 }, showlegend: false, legendgroup: `stage_${stageNum}`, xaxis: 'x2', yaxis: 'y2' });

        // Point B on y=x: (y_vap, y_vap)
        traces.push({ x: [y_vap], y: [y_vap], mode: 'markers', marker: { color, size: 8, symbol: 'diamond', line: { color: 'white', width: 1 } }, showlegend: false, legendgroup: `stage_${stageNum}`, xaxis: 'x2', yaxis: 'y2' });

        // VLE: operating step vertical UP from (y_vap, y_vap) → (y_vap, y_vap_next)
        const nextTie   = (i + 1 < numStages) ? tieLinesAsc[i + 1] : null;
        const y_vap_next = nextTie ? nextTie.x[1] : results.xD;
        traces.push({ x: [y_vap, y_vap], y: [y_vap, y_vap_next], mode: 'lines', line: { color, width: 2.5 }, showlegend: false, legendgroup: `stage_${stageNum}`, xaxis: 'x2', yaxis: 'y2' });

        // ── Vertical connector lines: VLE → Enthalpy ──
        // Bottom plot: from equilibrium point upward to top edge
        traces.push({
            x: [x_liq, x_liq], y: [y_vap, 1.05],
            mode: 'lines',
            name: !connectorLegendAdded ? 'VLE↔Enthalpy Links' : undefined,
            showlegend: !connectorLegendAdded,
            line: { color, width: 1.2, dash: 'dot' },
            legendgroup: `stage_${stageNum}`,
            xaxis: 'x2', yaxis: 'y2'
        });
        connectorLegendAdded = true;

        // Top plot: from bottom edge up to HL curve
        traces.push({ x: [x_liq, x_liq], y: [results.yMin, H_liq], mode: 'lines', showlegend: false, line: { color, width: 1.2, dash: 'dot' }, legendgroup: `stage_${stageNum}`, xaxis: 'x', yaxis: 'y' });

        // Bottom plot: y_vap upward
        traces.push({ x: [y_vap, y_vap], y: [y_vap, 1.05], mode: 'lines', showlegend: false, line: { color, width: 1.2, dash: 'dot' }, legendgroup: `stage_${stageNum}`, xaxis: 'x2', yaxis: 'y2' });

        // Top plot: from bottom edge up to HV curve
        traces.push({ x: [y_vap, y_vap], y: [results.yMin, H_vap], mode: 'lines', showlegend: false, line: { color, width: 1.2, dash: 'dot' }, legendgroup: `stage_${stageNum}`, xaxis: 'x', yaxis: 'y' });
    }

    // ── Layout ──
    const containerWidth = document.querySelector('.main-content')?.clientWidth || 1050;
    const dataset = DATASETS[systemType];
    const titleText = dataset
        ? `<b>Ponchon–Savarit: ${dataset.name}</b>`
        : '<b>Ponchon–Savarit Diagram</b>';

    const layout = {
        title: { text: titleText, font: { size: 17, family: 'Arial, sans-serif' }, x: 0.5, xanchor: 'center' },
        annotations: [{ text: '<b>Ponchon–Savarit Diagram (H-x-y)</b>', xref: 'paper', yref: 'paper', x: 0.5, y: 0.995, xanchor: 'center', yanchor: 'top', showarrow: false, font: { size: 12, color: '#1565C0' } }],
        grid: { rows: 2, columns: 1, pattern: 'independent', roworder: 'top to bottom' },
        margin: { l: 80, r: 220, t: 60, b: 70 },
        xaxis:  { domain: [0.0, 0.95], anchor: 'y',  range: [0, 1], tickformat: '.2f', showline: true, linecolor: '#555', mirror: true, showgrid: true, gridcolor: '#E0E0E0', zeroline: true, zerolinecolor: '#333', zerolinewidth: 1.5, title: '' },
        yaxis:  { domain: [0.52, 0.97], anchor: 'x', title: { text: `<b>Enthalpy (${enthalpyUnit})</b>`, font: { size: 13 } }, range: [results.yMin, results.yMax], showline: true, linecolor: '#555', mirror: true, showgrid: true, gridcolor: '#E0E0E0', zeroline: true, zerolinecolor: '#333', zerolinewidth: 1 },
        xaxis2: { domain: [0.0, 0.95], anchor: 'y2', range: [0, 1], tickformat: '.2f', showline: true, linecolor: '#555', mirror: true, showgrid: true, gridcolor: '#E0E0E0', title: { text: '<b>Mole Fraction (x or y)</b>', font: { size: 13 } } },
        yaxis2: { domain: [0.03, 0.48], anchor: 'x2', title: { text: '<b>y (Vapor Fraction)</b>', font: { size: 13 } }, range: [0, 1.02], showline: true, linecolor: '#555', mirror: true, showgrid: true, gridcolor: '#E0E0E0', zeroline: true, zerolinecolor: '#333', zerolinewidth: 1 },
        legend: { x: 1.01, y: 1.0, xanchor: 'left', yanchor: 'top', font: { size: 10, family: 'Arial, sans-serif' }, bgcolor: 'rgba(255,255,255,0.95)', bordercolor: '#BDBDBD', borderwidth: 1, tracegroupgap: 1 },
        plot_bgcolor: '#FAFAFA', paper_bgcolor: '#FFFFFF',
        height: 780, width: containerWidth - 40
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
    const isBritish  = systemType !== 'user-defined';
    const flowUnit   = isBritish ? 'lbmol/hr' : 'kmol/hr';
    const dutyUnit   = isBritish ? 'BTU/hr' : 'kW';

    const rminLabel = results.RMin === 0
        ? (results.tangent_pinch ? '≈ 0.00 (tangent)' : '≈ 0.00 (strip-ltd)')
        : results.RMin.toFixed(4);
    quickSummaryContent.innerHTML = `
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>Stages:</span><strong>${results.stages}</strong></div></div>
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>Feed Stage:</span><strong>${results.feed_stage}</strong></div></div>
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>R<sub>min</sub>:</span><strong>${rminLabel}</strong></div></div>
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>Qc:</span><strong>${results.QcDuty?.toLocaleString() || 0} ${dutyUnit}</strong></div></div>
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>D:</span><strong>${results.D} ${flowUnit}</strong></div></div>
        <div class="col-md-2 col-sm-4"><div class="summary-item"><span>W:</span><strong>${results.W} ${flowUnit}</strong></div></div>
    `;
    quickSummary.style.display = 'block';
}

// ==========================================
// DISPLAY RESULTS
// ==========================================
function displayResults(results) {
    createPlot(results);
    updateQuickSummary(results);

    const systemType   = document.getElementById('systemType')?.value || 'user-defined';
    const isBritish    = systemType !== 'user-defined';
    const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';
    const flowUnit     = isBritish ? 'lbmol/hr' : 'kmol/hr';
    const dutyUnit     = isBritish ? 'BTU/hr' : 'kW';

    // Summary table
    const summaryBody = document.getElementById('summaryBody');
    if (summaryBody) {
        summaryBody.innerHTML = `
            <tr><td>Distillate Flow (D)</td><td>${results.D} ${flowUnit}</td></tr>
            <tr><td>Bottoms Flow (W)</td><td>${results.W} ${flowUnit}</td></tr>
            <tr><td>Δ_R (xR, HR)</td><td>(${results.xDeltaR}, ${results.HDeltaR.toLocaleString()} ${enthalpyUnit})</td></tr>
            <tr><td>Δ_S (xS, HS)</td><td>(${results.xDeltaS}, ${results.HDeltaS.toLocaleString()} ${enthalpyUnit})</td></tr>
            <tr><td>Condenser Duty (Qc)</td><td>${results.QcDuty?.toLocaleString() || 0} ${dutyUnit}</td></tr>
            <tr><td>Reboiler Duty (Qr)</td><td>${results.QrDuty?.toLocaleString() || 0} ${dutyUnit}</td></tr>
            <tr><td>Δ_R,min</td><td>(${results.xD}, ${results.QPrimeMin?.toLocaleString()} ${enthalpyUnit})</td></tr>
            <tr><td>Δ_S,min</td><td>(${results.xB}, ${results.QDoublePrimeMin?.toLocaleString()} ${enthalpyUnit})</td></tr>
            <tr><td>Minimum Reflux Ratio (R<sub>min</sub>)</td><td>${
                results.RMin === 0
                    ? (results.tangent_pinch
                        ? '≈ 0 <small class="text-muted">(tangent pinch — concave VLE)</small>'
                        : '≈ 0 <small class="text-muted">(stripping-limited — xD ≤ y<sub>F</sub>)</small>')
                    : results.RMin.toFixed(4) + (results.tangent_pinch ? ' <small class="text-muted">(tangent pinch)</small>' : '')
            }</td></tr>
            <tr><td>Number of Stages</td><td>${results.stages}</td></tr>
            <tr><td>Feed Stage</td><td>${results.feed_stage}</td></tr>
        `;
    }

    // Stages table
    const stagesBody = document.getElementById('stagesBody');
    if (stagesBody) {
        stagesBody.innerHTML = results.stage_compositions?.length
            ? results.stage_compositions.map((s, i) => `<tr><td>Stage ${i + 1}</td><td>${s.x.toFixed(4)}</td><td>${s.y.toFixed(4)}</td></tr>`).join('')
            : '<tr><td colspan="3" class="text-center text-muted">No stage data</td></tr>';
    }

    // Enthalpy table
    const enthalpyBody = document.getElementById('enthalpyBody');
    if (enthalpyBody) {
        enthalpyBody.innerHTML = results.stage_compositions?.length
            ? results.stage_compositions.map((stage, i) => {
                const idxLiq = Math.min(Math.round(stage.x * 199), 199);
                const idxVap = Math.min(Math.round(stage.y * 199), 199);
                const HL = results.HL_curve?.[idxLiq]?.toFixed(0) || 'N/A';
                const HV = results.HV_curve?.[idxVap]?.toFixed(0) || 'N/A';
                return `<tr><td>Stage ${i + 1}</td><td>${HL}</td><td>${HV}</td></tr>`;
            }).join('')
            : '<tr><td colspan="3" class="text-center text-muted">No enthalpy data</td></tr>';
    }

    document.getElementById('exportBtn').disabled = false;
    currentResults = results;
    showToast('✅ Calculation completed!', 'success');
}

// ==========================================
// INITIALIZE PYODIDE
// ==========================================
async function initPyodide() {
    const loadingDiv   = document.getElementById('pyodide-loading');
    const progressBar  = document.getElementById('loading-progress');
    const statusText   = document.getElementById('loading-status');
    const statusDiv    = document.getElementById('pyodide-status');
    const calculateBtn = document.getElementById('calculateBtn');

    if (!loadingDiv) return;
    loadingDiv.style.display = 'flex';

    try {
        const steps = [
            { msg: 'Loading Pyodide core...', pct: '20%',  fn: async () => { pyodide = await loadPyodide({ indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.23.4/full/' }); } },
            { msg: 'Loading numpy...',         pct: '45%',  fn: async () => { await pyodide.loadPackage('numpy'); } },
            { msg: 'Loading scipy...',         pct: '70%',  fn: async () => { await pyodide.loadPackage('scipy'); } },
            { msg: 'Loading pandas...',        pct: '85%',  fn: async () => { await pyodide.loadPackage('pandas'); } },
            { msg: 'Initializing calc...',     pct: '95%',  fn: async () => { await loadCalculatorCode(); } },
        ];

        for (const step of steps) {
            statusText.innerText = step.msg;
            progressBar.style.width = step.pct;
            progressBar.innerText = step.pct;
            await step.fn();
        }

        progressBar.style.width = '100%';
        progressBar.innerText = '100%';
        statusDiv.innerHTML = '<i class="fas fa-check-circle text-success"></i> Python ready!';
        statusDiv.className = 'alert alert-success py-2';
        pyodideReady = true;
        calculateBtn.disabled = false;
        showToast('✅ Pyodide loaded!', 'success');
        setTimeout(() => loadingDiv.style.display = 'none', 500);

    } catch (error) {
        console.error('Pyodide error:', error);
        statusDiv.innerHTML = `<i class="fas fa-exclamation-circle text-danger"></i> Error: ${error.message}`;
        statusDiv.className = 'alert alert-danger py-2';
        loadingDiv.style.display = 'none';
    }
}

// ==========================================
// EVENT LISTENERS
// ==========================================
document.addEventListener('DOMContentLoaded', function () {

    // Preview on input change
    ['xData', 'yData', 'Hl', 'Hv'].forEach(id => {
        document.getElementById(id)?.addEventListener('input', updatePreview);
    });

    // System type change
    document.getElementById('systemType')?.addEventListener('change', function () {
        const type = this.value;
        currentSystem = type;
        document.getElementById('quickSummary').style.display = 'none';

        if (type === 'user-defined') {
            document.getElementById('systemInfo').style.display = 'none';
            document.getElementById('datasetInfoPanel').style.display = 'none';
            document.getElementById('xData').value = '[0, 0.08, 0.18, 0.25, 0.49, 0.65, 0.79, 0.91, 1.0]';
            document.getElementById('yData').value = '[0, 0.28, 0.43, 0.51, 0.73, 0.83, 0.90, 0.96, 1.0]';
            document.getElementById('Hl').value    = '[24.3, 24.1, 23.2, 22.8, 22.05, 21.75, 21.7, 21.6, 21.4]';
            document.getElementById('Hv').value    = '[61.2, 59.6, 58.5, 58.1, 56.5, 55.2, 54.4, 53.8, 53.3]';
            document.getElementById('TData').value = '';
            updateUnitDisplay('user-defined');
            showToast('📝 User Defined mode', 'info');
        } else {
            loadDataset(type);
        }
        updatePreview();
    });

    // q selector
    document.getElementById('q')?.addEventListener('change', function () {
        document.getElementById('customQDiv').style.display = this.value === 'custom' ? 'block' : 'none';
    });

    // Calculate button
    document.getElementById('calculateBtn')?.addEventListener('click', async function () {
        if (!pyodideReady) { alert('Pyodide belum siap.'); return; }

        const xData = parseArrayString(document.getElementById('xData').value);
        const yData = parseArrayString(document.getElementById('yData').value);
        const Hl    = parseArrayString(document.getElementById('Hl').value);
        const Hv    = parseArrayString(document.getElementById('Hv').value);

        for (const v of [validateArray(xData, 'xData'), validateArray(yData, 'yData'), validateArray(Hl, 'Hl'), validateArray(Hv, 'Hv')]) {
            if (!v.valid) { alert('Error: ' + v.error); return; }
        }
        if (new Set([xData.length, yData.length, Hl.length, Hv.length]).size > 1) {
            alert('Error: Semua array harus sama panjang!'); return;
        }

        let q = document.getElementById('q').value;
        q = q === 'custom' ? parseFloat(document.getElementById('customQ').value) : parseFloat(q);

        const inputData = {
            xData, yData, Hl, Hv,
            zF: parseFloat(document.getElementById('zF').value),
            F:  parseFloat(document.getElementById('F').value),
            xD: parseFloat(document.getElementById('xD').value),
            xB: parseFloat(document.getElementById('xB').value),
            q, R: parseFloat(document.getElementById('R').value)
        };

        this.disabled = true;
        document.getElementById('loading').style.display = 'flex';

        try {
            const results = await runCalculation(inputData);
            if (results.error) { alert('Error: ' + results.error); }
            else { displayResults(results); }
        } catch (error) {
            alert('Error: ' + error.message);
        } finally {
            this.disabled = false;
            document.getElementById('loading').style.display = 'none';
        }
    });

    // Export button
    document.getElementById('exportBtn')?.addEventListener('click', function () {
        if (!currentResults) return;
        const systemType   = document.getElementById('systemType')?.value || 'user-defined';
        const isBritish    = systemType !== 'user-defined';
        const enthalpyUnit = isBritish ? 'BTU/lbmole' : 'MJ/kmol';

        let csv = `Stage,x (Liquid),y (Vapor),HL (${enthalpyUnit}),HV (${enthalpyUnit})\n`;
        currentResults.stage_compositions?.forEach((stage, i) => {
            const idxLiq = Math.min(Math.round(stage.x * 199), 199);
            const idxVap = Math.min(Math.round(stage.y * 199), 199);
            const HL = currentResults.HL_curve?.[idxLiq]?.toFixed(2) || 'N/A';
            const HV = currentResults.HV_curve?.[idxVap]?.toFixed(2) || 'N/A';
            csv += `${i + 1},${stage.x.toFixed(4)},${stage.y.toFixed(4)},${HL},${HV}\n`;
        });

        const blob = new Blob([csv], { type: 'text/csv' });
        const url  = URL.createObjectURL(blob);
        const a    = document.createElement('a');
        a.href = url;
        a.download = `ponchon_savarit_${systemType}_${new Date().toISOString().slice(0, 10)}.csv`;
        a.click();
        URL.revokeObjectURL(url);
        showToast('📊 Results exported to CSV', 'success');
    });

    // Init
    updatePreview();
    initPyodide();
});

