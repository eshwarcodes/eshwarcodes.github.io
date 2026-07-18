/* ==========================================================================
   Photosynthesis Explorer — app.js
   Interactive teaching tool for photosynthesis concepts
   ========================================================================== */

// ── Constants ─────────────────────────────────────────────────────────────
const R_GAS = 8.314;          // J mol⁻¹ K⁻¹
const O2_AMBIENT = 210;       // mmol/mol (21%)

// ── Plotly defaults ───────────────────────────────────────────────────────
const PLOT_CFG = { responsive: true, displayModeBar: false };
const LAYOUT_BASE = {
    font: { family: 'Segoe UI, system-ui, sans-serif', size: 13 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    margin: { l: 65, r: 25, t: 40, b: 55 },
    xaxis: { gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
    yaxis: { gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
    legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(255,255,255,0.85)',
              bordercolor: '#ddd', borderwidth: 1, font: { size: 11.5 } },
    hovermode: 'x unified',
};

function mergeLayout(extra) {
    return Object.assign({}, JSON.parse(JSON.stringify(LAYOUT_BASE)), extra);
}

// ==========================================================================
// PHOTOSYNTHESIS MODELS
// ==========================================================================

// ── Arrhenius with high-T deactivation (Kattge & Knorr 2007) ─────────────
function arrh(val25, Ha, Hd, Ds, T_C) {
    const Tk = T_C + 273.15, Tr = 298.15;
    const f = Math.exp(Ha * (Tk - Tr) / (Tr * R_GAS * Tk));
    const numD = 1 + Math.exp((Tr * Ds - Hd) / (Tr * R_GAS));
    const denD = 1 + Math.exp((Tk * Ds - Hd) / (Tk * R_GAS));
    return val25 * f * numD / denD;
}

// ── Simple Arrhenius (no deactivation) ────────────────────────────────────
function arrh_s(val25, Ha, T_C) {
    const Tk = T_C + 273.15, Tr = 298.15;
    return val25 * Math.exp(Ha * (Tk - Tr) / (Tr * R_GAS * Tk));
}

// ── Crop parameter sets ───────────────────────────────────────────────────
const CROPS = {
    lettuce: {
        name: 'Lettuce', Vcmax25: 60, Jmax25: 110,
        Ha_V: 65330, Hd_V: 149250, Ds_V: 485,
        Ha_J: 43540, Hd_J: 152040, Ds_J: 495,
        Rd25: 1.2, Ha_Rd: 46390,
        Gamma25: 42.75, Ha_Gamma: 37830,
        Kc25: 404, Ha_Kc: 79430, Ko25: 278, Ha_Ko: 36380,
        alpha_j: 0.3, theta: 0.7,
    },
    tomato: {
        name: 'Tomato', Vcmax25: 100, Jmax25: 180,
        Ha_V: 65330, Hd_V: 220000, Ds_V: 640,
        Ha_J: 43540, Hd_J: 200000, Ds_J: 610,
        Rd25: 1.8, Ha_Rd: 46390,
        Gamma25: 42.75, Ha_Gamma: 37830,
        Kc25: 404, Ha_Kc: 79430, Ko25: 278, Ha_Ko: 36380,
        alpha_j: 0.3, theta: 0.7,
    },
    cucumber: {
        name: 'Cucumber', Vcmax25: 90, Jmax25: 165,
        Ha_V: 65330, Hd_V: 210000, Ds_V: 625,
        Ha_J: 43540, Hd_J: 195000, Ds_J: 600,
        Rd25: 1.5, Ha_Rd: 46390,
        Gamma25: 42.75, Ha_Gamma: 37830,
        Kc25: 404, Ha_Kc: 79430, Ko25: 278, Ha_Ko: 36380,
        alpha_j: 0.3, theta: 0.7,
    },
};

// ── FvCB model with variable O2 ──────────────────────────────────────────
function fvcb(Ci, T_C, PPFD, p, O2_mmol) {
    if (O2_mmol === undefined) O2_mmol = O2_AMBIENT;
    const Vcmax = arrh(p.Vcmax25, p.Ha_V, p.Hd_V, p.Ds_V, T_C);
    const Jmax  = arrh(p.Jmax25,  p.Ha_J, p.Hd_J, p.Ds_J, T_C);
    const Rd    = arrh_s(p.Rd25, p.Ha_Rd, T_C);
    const Kc    = arrh_s(p.Kc25, p.Ha_Kc, T_C);
    const Ko    = arrh_s(p.Ko25, p.Ha_Ko, T_C);
    // Scale Gamma* by O2 ratio
    const Gamma_base = arrh_s(p.Gamma25, p.Ha_Gamma, T_C);
    const Gamma = Gamma_base * (O2_mmol / O2_AMBIENT);
    const Km = Kc * (1 + O2_mmol / Ko);

    // Rubisco-limited
    const Ac = Vcmax * (Ci - Gamma) / (Ci + Km);

    // Electron transport rate (quadratic)
    const aQ = p.alpha_j * PPFD;
    const a = p.theta;
    const b = -(aQ + Jmax);
    const c = aQ * Jmax;
    const disc = Math.max(b * b - 4 * a * c, 0);
    const J = (-b - Math.sqrt(disc)) / (2 * a);

    // RuBP-limited
    const Aj = (J / 4) * (Ci - Gamma) / (Ci + 2 * Gamma);

    const An = Math.min(Ac, Aj) - Rd;
    return { An, Ac: Ac - Rd, Aj: Aj - Rd, Rd, Gamma, Vcmax, Jmax, J };
}

// ── Light response curve (rectangular hyperbola) ─────────────────────────
function lightResponse(PPFD, Amax, alpha, Rd) {
    return (alpha * PPFD * Amax) / (alpha * PPFD + Amax) - Rd;
}

// ── Medlyn stomatal conductance ──────────────────────────────────────────
function medlyn_gs(An, Ca, VPD_kPa, g0, g1) {
    const D = Math.max(VPD_kPa, 0.05);
    if (An <= 0) return g0;
    return g0 + 1.6 * (1 + g1 / Math.sqrt(D)) * An / Ca;
}

// ── Saturation vapor pressure (Tetens) ───────────────────────────────────
function esat(T_C) {
    return 0.6108 * Math.exp(17.27 * T_C / (T_C + 237.3)); // kPa
}

// ==========================================================================
// TAB SWITCHING
// ==========================================================================
const tabs = document.querySelectorAll('.tab');
const panels = document.querySelectorAll('.panel');
let currentTab = 'light-quality';
const rendered = new Set();

tabs.forEach(t => {
    t.addEventListener('click', () => {
        tabs.forEach(x => x.classList.remove('active'));
        panels.forEach(x => x.classList.remove('active'));
        t.classList.add('active');
        const id = t.dataset.tab;
        document.getElementById(id).classList.add('active');
        currentTab = id;
        if (!rendered.has(id)) {
            rendered.add(id);
            renderTab(id);
        }
    });
});

function renderTab(id) {
    switch (id) {
        case 'light-quality':    renderLightQuality(); break;
        case 'light-intensity':  renderLightIntensity(); break;
        case 'co2':              renderCO2(); break;
        case 'temperature':      renderTemperature(); break;
        case 'water':            renderWater(); break;
        case 'photorespiration': renderPhotorespiration(); break;
        case 'mechanisms':       renderMechanisms(); break;
        case 'source-sink':      renderSourceSink(); break;
    }
}

// ==========================================================================
// HELPER: bind sliders
// ==========================================================================
function bindSlider(id, valId, cb, fmt) {
    const sl = document.getElementById(id);
    const vl = document.getElementById(valId);
    const update = () => {
        const v = parseFloat(sl.value);
        vl.textContent = fmt ? fmt(v) : v;
        cb(v);
    };
    sl.addEventListener('input', update);
}

function bindCheckbox(id, cb) {
    document.getElementById(id).addEventListener('change', cb);
}

// ==========================================================================
// TAB 1: LIGHT QUALITY
// ==========================================================================
function buildAbsorptionSpectra() {
    const wl = [];
    for (let w = 300; w <= 780; w++) wl.push(w);

    function gaussian(center, sigma, peak) {
        return wl.map(w => peak * Math.exp(-0.5 * ((w - center) / sigma) ** 2));
    }

    // Chlorophyll a: peaks at 430 & 662 nm
    const chla_b = gaussian(430, 15, 0.90);
    const chla_r = gaussian(662, 14, 0.85);
    const chla = wl.map((_, i) => Math.min(chla_b[i] + chla_r[i], 1.0));

    // Chlorophyll b: peaks at 453 & 642 nm
    const chlb_b = gaussian(453, 16, 0.75);
    const chlb_r = gaussian(642, 15, 0.60);
    const chlb = wl.map((_, i) => Math.min(chlb_b[i] + chlb_r[i], 1.0));

    // Carotenoids: broad 400-500
    const car1 = gaussian(450, 30, 0.70);
    const car2 = gaussian(475, 25, 0.55);
    const carot = wl.map((_, i) => Math.min(car1[i] + car2[i], 1.0));

    // Action spectrum (McCree 1972 approximation)
    const act_b = gaussian(440, 22, 0.90);
    const act_r = gaussian(670, 20, 1.0);
    const act_g = gaussian(550, 35, 0.35);
    const action = wl.map((_, i) => Math.min(act_b[i] + act_r[i] + act_g[i], 1.0));

    return { wl, chla, chlb, carot, action };
}

function renderLightQuality() {
    const spec = buildAbsorptionSpectra();

    function draw() {
        const traces = [];
        const showPAR = document.getElementById('lq-par').checked;

        if (showPAR) {
            traces.push({
                x: [400, 400, 700, 700],
                y: [0, 1.05, 1.05, 0],
                fill: 'toself',
                fillcolor: 'rgba(255,215,0,0.12)',
                line: { color: 'rgba(200,170,0,0.4)', width: 1, dash: 'dot' },
                name: 'PAR band',
                hoverinfo: 'skip',
            });
        }
        if (document.getElementById('lq-chla').checked) {
            traces.push({ x: spec.wl, y: spec.chla, name: 'Chlorophyll a',
                line: { color: '#1b9e77', width: 2.5 } });
        }
        if (document.getElementById('lq-chlb').checked) {
            traces.push({ x: spec.wl, y: spec.chlb, name: 'Chlorophyll b',
                line: { color: '#7570b3', width: 2.5 } });
        }
        if (document.getElementById('lq-carot').checked) {
            traces.push({ x: spec.wl, y: spec.carot, name: 'Carotenoids',
                line: { color: '#d95f02', width: 2.5 } });
        }
        if (document.getElementById('lq-action').checked) {
            traces.push({ x: spec.wl, y: spec.action, name: 'Action Spectrum',
                line: { color: '#e7298a', width: 3, dash: 'dash' } });
        }

        // Spectral rainbow background
        const shapes = [];
        const rainbow = [
            [380, 420, 'rgba(108,0,180,0.08)'],  // violet
            [420, 480, 'rgba(0,70,255,0.08)'],    // blue
            [480, 510, 'rgba(0,180,180,0.06)'],   // cyan
            [510, 565, 'rgba(0,180,0,0.06)'],     // green
            [565, 590, 'rgba(255,210,0,0.06)'],   // yellow
            [590, 625, 'rgba(255,140,0,0.06)'],   // orange
            [625, 750, 'rgba(255,0,0,0.06)'],     // red
        ];
        rainbow.forEach(([x0, x1, col]) => {
            shapes.push({ type: 'rect', xref: 'x', yref: 'paper',
                x0, x1, y0: 0, y1: 1, fillcolor: col, line: { width: 0 }, layer: 'below' });
        });

        Plotly.react('plot-lq', traces, mergeLayout({
            title: { text: 'Pigment Absorption & Photosynthetic Action Spectra', font: { size: 15 } },
            xaxis: { title: 'Wavelength (nm)', range: [300, 780], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'Relative Absorption / Action', range: [0, 1.08], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            shapes,
        }), PLOT_CFG);
    }

    draw();
    ['lq-chla', 'lq-chlb', 'lq-carot', 'lq-action', 'lq-par'].forEach(id =>
        bindCheckbox(id, draw));
}

// ==========================================================================
// TAB 2: LIGHT INTENSITY
// ==========================================================================
function renderLightIntensity() {
    function draw() {
        const Amax = parseFloat(document.getElementById('li-amax').value);
        const alpha = parseFloat(document.getElementById('li-alpha').value);
        const Rd = parseFloat(document.getElementById('li-rd').value);

        const ppfd = [], an = [];
        for (let p = 0; p <= 2500; p += 10) {
            ppfd.push(p);
            an.push(lightResponse(p, Amax, alpha, Rd));
        }

        // Find compensation point (An = 0)
        const compPPFD = Rd / (alpha * (1 - Rd / Amax));
        // Approximate saturation point (90% of Amax - Rd)
        const target90 = 0.9 * (Amax - Rd);
        let satPPFD = 0;
        for (let p = 0; p <= 2500; p += 5) {
            if (lightResponse(p, Amax, alpha, Rd) >= target90) { satPPFD = p; break; }
        }

        const traces = [
            { x: ppfd, y: an, name: 'Net photosynthesis A<sub>n</sub>',
              line: { color: '#2d6a4f', width: 3 } },
            { x: [0, 2500], y: [0, 0], name: 'Zero line',
              line: { color: '#aaa', width: 1, dash: 'dot' }, showlegend: false },
        ];

        const annotations = [
            { x: compPPFD, y: 0, text: `Compensation<br>${Math.round(compPPFD)} µmol`,
              showarrow: true, arrowhead: 2, ax: 50, ay: -40,
              font: { size: 11, color: '#c0392b' }, arrowcolor: '#c0392b' },
        ];
        if (satPPFD > 0 && satPPFD < 2400) {
            annotations.push({
                x: satPPFD, y: target90, text: `~Saturation<br>${Math.round(satPPFD)} µmol`,
                showarrow: true, arrowhead: 2, ax: 50, ay: -30,
                font: { size: 11, color: '#2980b9' }, arrowcolor: '#2980b9',
            });
        }

        // Shade the saturated region
        const shapes = satPPFD > 0 ? [{
            type: 'rect', xref: 'x', yref: 'paper',
            x0: satPPFD, x1: 2500, y0: 0, y1: 1,
            fillcolor: 'rgba(52,152,219,0.06)', line: { width: 0 }, layer: 'below',
        }] : [];

        Plotly.react('plot-li', traces, mergeLayout({
            title: { text: 'Light Response Curve (Rectangular Hyperbola)', font: { size: 15 } },
            xaxis: { title: 'PPFD (µmol photons / m² / s)', range: [0, 2500], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            annotations, shapes,
        }), PLOT_CFG);
    }

    draw();
    bindSlider('li-amax', 'li-amax-val', draw);
    bindSlider('li-alpha', 'li-alpha-val', draw, v => v.toFixed(3));
    bindSlider('li-rd', 'li-rd-val', draw);
}

// ==========================================================================
// TAB 3: CO2 CONCENTRATION (A–Ci curve)
// ==========================================================================
function renderCO2() {
    function draw() {
        const Vcmax25 = parseFloat(document.getElementById('co2-vcmax').value);
        const Jmax25  = parseFloat(document.getElementById('co2-jmax').value);
        const PPFD    = parseFloat(document.getElementById('co2-ppfd').value);
        const T_C     = parseFloat(document.getElementById('co2-temp').value);
        const showComp = document.getElementById('co2-show-components').checked;

        const p = Object.assign({}, CROPS.tomato, { Vcmax25, Jmax25 });
        const ci = [], an = [], ac = [], aj = [];
        for (let c = 0; c <= 1500; c += 5) {
            ci.push(c);
            const r = fvcb(c, T_C, PPFD, p);
            an.push(r.An);
            ac.push(r.Ac);
            aj.push(r.Aj);
        }

        const traces = [
            { x: ci, y: an, name: 'Net A<sub>n</sub> (min of limitations)',
              line: { color: '#2d6a4f', width: 3 } },
        ];
        if (showComp) {
            traces.push({ x: ci, y: ac, name: 'Rubisco-limited (A<sub>c</sub>) — controlled by V<sub>cmax</sub>',
                line: { color: '#e76f51', width: 2, dash: 'dash' } });
            traces.push({ x: ci, y: aj, name: 'RuBP-limited (A<sub>j</sub>) — controlled by J<sub>max</sub>',
                line: { color: '#457b9d', width: 2, dash: 'dash' } });
        }
        traces.push({ x: [0, 1500], y: [0, 0], showlegend: false,
            line: { color: '#aaa', width: 1, dash: 'dot' } });

        // Find crossover point where Ac and Aj intersect
        let crossCi = 300;  // default
        for (let c = 5; c <= 1500; c += 5) {
            const r = fvcb(c, T_C, PPFD, p);
            if (r.Aj <= r.Ac) { crossCi = c; break; }
        }
        const crossAn = fvcb(crossCi, T_C, PPFD, p).An;

        // Ambient Ci ≈ 0.7 × Ca
        const ambientCi = 0.7 * 420;
        const annotations = [{
            x: ambientCi, y: fvcb(ambientCi, T_C, PPFD, p).An,
            text: 'Ambient<br>(~420 ppm)',
            showarrow: true, arrowhead: 2, ax: -60, ay: -30,
            font: { size: 11, color: '#8e44ad' }, arrowcolor: '#8e44ad',
        }];

        // Shaded regions + labels showing which parameter controls which zone
        const shapes = [];
        if (showComp) {
            // Rubisco-limited zone (left of crossover)
            shapes.push({ type: 'rect', xref: 'x', yref: 'paper',
                x0: 0, x1: crossCi, y0: 0, y1: 1,
                fillcolor: 'rgba(231,111,81,0.06)', line: { width: 0 }, layer: 'below' });
            // RuBP-limited zone (right of crossover)
            shapes.push({ type: 'rect', xref: 'x', yref: 'paper',
                x0: crossCi, x1: 1500, y0: 0, y1: 1,
                fillcolor: 'rgba(69,123,157,0.06)', line: { width: 0 }, layer: 'below' });
            // Crossover line
            shapes.push({ type: 'line', x0: crossCi, x1: crossCi, y0: 0, y1: 1, yref: 'paper',
                line: { color: '#999', width: 1.5, dash: 'dashdot' } });

            annotations.push({
                x: crossCi * 0.4, y: 1.0, yref: 'paper', yanchor: 'top',
                text: '<b>V<sub>cmax</sub> controls<br>this region</b><br>(Rubisco-limited)',
                showarrow: false, font: { size: 11, color: '#c0392b' },
            });
            annotations.push({
                x: crossCi + (1500 - crossCi) * 0.5, y: 1.0, yref: 'paper', yanchor: 'top',
                text: '<b>J<sub>max</sub> controls<br>this region</b><br>(RuBP-limited)',
                showarrow: false, font: { size: 11, color: '#2471a3' },
            });
            annotations.push({
                x: crossCi, y: crossAn,
                text: 'Transition',
                showarrow: true, arrowhead: 2, ax: 45, ay: 25,
                font: { size: 10, color: '#666' }, arrowcolor: '#999',
            });
        }

        Plotly.react('plot-co2', traces, mergeLayout({
            title: { text: 'A–C<sub>i</sub> Curve (FvCB Model)', font: { size: 15 } },
            xaxis: { title: 'Intercellular CO₂, C<sub>i</sub> (µmol/mol)', range: [0, 1500], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            annotations, shapes,
        }), PLOT_CFG);
    }

    draw();
    bindSlider('co2-vcmax', 'co2-vcmax-val', draw);
    bindSlider('co2-jmax', 'co2-jmax-val', draw);
    bindSlider('co2-ppfd', 'co2-ppfd-val', draw);
    bindSlider('co2-temp', 'co2-temp-val', draw);
    bindCheckbox('co2-show-components', draw);
}

// ==========================================================================
// TAB 4: TEMPERATURE
// ==========================================================================
function renderTemperature() {
    function draw() {
        const PPFD = parseFloat(document.getElementById('temp-ppfd').value);
        const Ca   = parseFloat(document.getElementById('temp-co2').value);
        const Ci   = 0.7 * Ca;
        const showVJ = document.getElementById('temp-show-vcmax').checked;

        const temps = [];
        for (let t = 0; t <= 50; t += 0.5) temps.push(t);

        const traces = [];
        const cropList = [
            { id: 'temp-lettuce', key: 'lettuce', color: '#2d6a4f' },
            { id: 'temp-tomato',  key: 'tomato',  color: '#e63946' },
            { id: 'temp-cucumber',key: 'cucumber', color: '#457b9d' },
        ];

        cropList.forEach(({ id, key, color }) => {
            if (!document.getElementById(id).checked) return;
            const p = CROPS[key];
            const an = temps.map(t => fvcb(Ci, t, PPFD, p).An);
            traces.push({ x: temps, y: an, name: p.name + ' A<sub>n</sub>',
                line: { color, width: 2.5 } });

            if (showVJ) {
                const vc = temps.map(t => arrh(p.Vcmax25, p.Ha_V, p.Hd_V, p.Ds_V, t));
                const jm = temps.map(t => arrh(p.Jmax25, p.Ha_J, p.Hd_J, p.Ds_J, t));
                traces.push({ x: temps, y: vc, name: p.name + ' V<sub>cmax</sub>',
                    line: { color, width: 1.5, dash: 'dot' }, yaxis: 'y2' });
                traces.push({ x: temps, y: jm, name: p.name + ' J<sub>max</sub>',
                    line: { color, width: 1.5, dash: 'dashdot' }, yaxis: 'y2' });
            }
        });

        traces.push({ x: [0, 50], y: [0, 0], showlegend: false,
            line: { color: '#aaa', width: 1, dash: 'dot' } });

        const layout = mergeLayout({
            title: { text: 'Temperature Response of Photosynthesis', font: { size: 15 } },
            xaxis: { title: 'Leaf Temperature (°C)', range: [0, 50], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
        });

        if (showVJ) {
            layout.yaxis2 = {
                title: 'V<sub>cmax</sub>, J<sub>max</sub> (µmol/m²/s)',
                overlaying: 'y', side: 'right', showgrid: false,
            };
        }

        Plotly.react('plot-temp', traces, layout, PLOT_CFG);
    }

    draw();
    bindSlider('temp-ppfd', 'temp-ppfd-val', draw);
    bindSlider('temp-co2', 'temp-co2-val', draw);
    ['temp-lettuce', 'temp-tomato', 'temp-cucumber', 'temp-show-vcmax'].forEach(id =>
        bindCheckbox(id, draw));
}

// ==========================================================================
// TAB 5: WATER & STOMATA
// ==========================================================================
function renderWater() {
    function draw() {
        const g1   = parseFloat(document.getElementById('water-g1').value);
        const g0   = parseFloat(document.getElementById('water-g0').value);
        const PPFD = parseFloat(document.getElementById('water-ppfd').value);
        const T_C  = parseFloat(document.getElementById('water-temp').value);
        const RH   = parseFloat(document.getElementById('water-rh').value);

        const p = CROPS.tomato;
        const Ca = 400;
        const Ci_approx = 0.7 * Ca;
        const An_pot = Math.max(fvcb(Ci_approx, T_C, PPFD, p).An, 0.5);

        // Air VPD from temperature and RH
        const e_air = esat(T_C) * RH / 100;
        const VPD_air = esat(T_C) - e_air;

        // ── Leaf energy balance constants ────────────────────────────
        const P_atm = 101.3;
        const lambda = 44000;   // J/mol latent heat
        const cp_mol = 29.2;    // J/mol/K
        const g_bh = 1.4;       // mol/m²/s boundary layer conductance
        const R_abs = PPFD * 0.22;  // W/m² absorbed radiation on leaf

        // ── Coupled solver: gs-An-Tleaf all interact ────────────────
        // Uses damped iteration + warm-start from previous SWF to avoid
        // bistability jumps in the leaf energy balance.
        const DAMP = 0.4;  // relaxation: new = damp*computed + (1-damp)*old

        function solveCoupled(swf_val, init_T, init_An, init_gs) {
            let T_leaf = init_T;
            let An = init_An;
            let gs = init_gs * swf_val;

            for (let iter = 0; iter < 30; iter++) {
                const vpd_leaf = Math.max(esat(T_leaf) - e_air, 0.05);

                const gs_new = Math.max(medlyn_gs(An, Ca, vpd_leaf, g0, g1) * swf_val, g0);
                gs = DAMP * gs_new + (1 - DAMP) * gs;

                const E_mol = gs * vpd_leaf / P_atm;
                const LE = E_mol * lambda;

                const T_new = T_C + (R_abs - LE) / (g_bh * cp_mol);
                T_leaf = DAMP * Math.max(T_C - 2, Math.min(T_C + 18, T_new)) + (1 - DAMP) * T_leaf;

                const Ci = Math.max(Ca - An / Math.max(gs * 1.6, 0.001), 50);
                const An_new = Math.max(fvcb(Ci, T_leaf, PPFD, p).An, 0);
                An = DAMP * An_new + (1 - DAMP) * An;
            }
            return { An, gs, T_leaf };
        }

        // ── Isothermal solver: stomatal effect only, no heating ─────
        function solveIsothermal(swf_val, init_An, init_gs) {
            let An = init_An;
            let gs = init_gs * swf_val;
            for (let iter = 0; iter < 15; iter++) {
                const gs_new = Math.max(medlyn_gs(An, Ca, VPD_air, g0, g1) * swf_val, g0);
                gs = DAMP * gs_new + (1 - DAMP) * gs;
                const Ci = Math.max(Ca - An / Math.max(gs * 1.6, 0.001), 50);
                const An_new = Math.max(fvcb(Ci, T_C, PPFD, p).An, 0);
                An = DAMP * An_new + (1 - DAMP) * An;
            }
            return { An, gs, T_leaf: T_C };
        }

        // ── Sweep SWF from 1.0 (wet) → 0.1 (dry), warm-starting ───
        const swfs = [];
        for (let s = 1.0; s >= 0.08; s -= 0.02) swfs.push(s);

        const coupled = [], isothermal = [];
        // Start at SWF=1.0 with well-watered initial conditions
        let prev_T = T_C, prev_An = An_pot;
        let prev_gs = medlyn_gs(An_pot, Ca, VPD_air, g0, g1);
        let prev_An_iso = An_pot, prev_gs_iso = prev_gs;

        for (let i = 0; i < swfs.length; i++) {
            const rc = solveCoupled(swfs[i], prev_T, prev_An, prev_gs / Math.max(swfs[i], 0.01));
            coupled.push(rc);
            prev_T = rc.T_leaf;
            prev_An = rc.An;
            prev_gs = rc.gs;

            const ri = solveIsothermal(swfs[i], prev_An_iso, prev_gs_iso / Math.max(swfs[i], 0.01));
            isothermal.push(ri);
            prev_An_iso = ri.An;
            prev_gs_iso = ri.gs;
        }

        // ── Plot 1: Leaf temperature vs soil water ──────────────────
        const traces1 = [];
        const tl_arr = coupled.map(r => r.T_leaf);

        // Shaded fill: T_leaf above T_air
        traces1.push({
            x: swfs.concat([...swfs].reverse()),
            y: tl_arr.concat(Array(swfs.length).fill(T_C)),
            fill: 'toself', fillcolor: 'rgba(231,76,60,0.20)',
            line: { width: 0 }, name: 'Lost evaporative cooling',
            hoverinfo: 'skip',
        });
        traces1.push({ x: swfs, y: tl_arr,
            name: 'T<sub>leaf</sub>',
            line: { color: '#e74c3c', width: 3 } });
        traces1.push({ x: swfs, y: Array(swfs.length).fill(T_C),
            name: `T<sub>air</sub> = ${T_C}°C`,
            line: { color: '#888', width: 2, dash: 'dash' } });

        const maxTl = Math.max(...tl_arr);
        const minSwfIdx = tl_arr.indexOf(maxTl);

        Plotly.react('plot-water-gs', traces1, mergeLayout({
            title: { text: `Leaf Heats Up as Soil Dries  (air ${T_C}°C, ${RH}% RH, VPD ${VPD_air.toFixed(1)} kPa)`, font: { size: 13 } },
            xaxis: { title: 'Soil Water Fraction (1 = wet, 0 = dry)', range: [1.02, 0.08],
                     gridcolor: '#e8e8e8', zerolinecolor: '#ccc', autorange: false },
            yaxis: { title: 'Temperature (°C)', range: [T_C - 2, Math.max(maxTl + 2, T_C + 3)],
                     gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            margin: { l: 55, r: 25, t: 42, b: 50 },
            legend: { x: 0.60, y: 0.30 },
            annotations: [{
                x: swfs[minSwfIdx], y: maxTl,
                text: `<b>+${(maxTl - T_C).toFixed(1)}°C</b>`,
                showarrow: true, arrowhead: 2, ax: 40, ay: -20,
                font: { size: 13, color: '#c0392b' }, arrowcolor: '#c0392b',
            }],
        }), PLOT_CFG);

        // ── Plot 2: An with vs without thermal coupling ─────────────
        const traces2 = [];
        const an_coupled = coupled.map(r => r.An);
        const an_iso     = isothermal.map(r => r.An);

        // Shaded gap = thermal penalty
        traces2.push({
            x: swfs.concat([...swfs].reverse()),
            y: an_iso.concat([...an_coupled].reverse()),
            fill: 'toself', fillcolor: 'rgba(231,76,60,0.18)',
            line: { width: 0 }, name: 'Thermal penalty',
            hoverinfo: 'skip',
        });
        traces2.push({ x: swfs, y: an_iso,
            name: 'A<sub>n</sub> stomata only (leaf = air temp)',
            line: { color: '#888', width: 2.5, dash: 'dash' } });
        traces2.push({ x: swfs, y: an_coupled,
            name: 'A<sub>n</sub> actual (leaf heats up)',
            line: { color: '#2d6a4f', width: 3 } });

        // Show the penalty at driest point
        const lastIdx = an_iso.length - 1;
        const penalty = an_iso[lastIdx] - an_coupled[lastIdx];
        const penaltyPct = (penalty / Math.max(an_iso[0], 1) * 100).toFixed(0);
        const totalDrop = ((an_coupled[lastIdx] - an_coupled[0]) / Math.max(an_coupled[0], 1) * 100).toFixed(0);
        const stomataOnlyDrop = ((an_iso[lastIdx] - an_iso[0]) / Math.max(an_iso[0], 1) * 100).toFixed(0);

        // Find mid-SWF for label placement
        const midIdx = Math.round(swfs.length * 0.65);

        Plotly.react('plot-water-an', traces2, mergeLayout({
            title: { text: 'Photosynthesis Drops More Than Stomatal Closure Alone Predicts', font: { size: 13 } },
            xaxis: { title: 'Soil Water Fraction (1 = wet, 0 = dry)', range: [1.02, 0.08],
                     gridcolor: '#e8e8e8', zerolinecolor: '#ccc', autorange: false },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', rangemode: 'tozero',
                     gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            margin: { l: 55, r: 25, t: 55, b: 50 },
            legend: { x: 0.55, y: 0.98, font: { size: 11 } },
            annotations: [
                { x: swfs[midIdx],
                  y: (an_iso[midIdx] + an_coupled[midIdx]) / 2,
                  text: '<b>Thermal<br>penalty</b>',
                  showarrow: false, font: { size: 12, color: '#c0392b' } },
                { x: 0.5, y: 1.08, xref: 'paper', yref: 'paper',
                  text: `At SWF=0.1:  stomata alone ${stomataOnlyDrop}%  |  with heating ${totalDrop}%  |  thermal penalty = ${penaltyPct}% extra loss`,
                  showarrow: false, font: { size: 11, color: '#555' } },
            ],
        }), PLOT_CFG);
    }

    draw();
    ['water-g1', 'water-g0', 'water-ppfd', 'water-temp', 'water-rh'].forEach(id => {
        const valId = id + '-val';
        const fmt = id === 'water-g0' ? v => v.toFixed(3) : null;
        bindSlider(id, valId, draw, fmt);
    });
}

// ==========================================================================
// TAB 6: PHOTORESPIRATION
// ==========================================================================
function renderPhotorespiration() {
    function draw() {
        const O2_pct  = parseFloat(document.getElementById('pr-o2').value);
        const T_C     = parseFloat(document.getElementById('pr-temp').value);
        const PPFD    = parseFloat(document.getElementById('pr-ppfd').value);
        const showRef = document.getElementById('pr-show-ref').checked;

        const O2_mmol = O2_pct * 10;  // % → mmol/mol
        const p = CROPS.tomato;

        const ci = [], an = [], an_ref = [];
        for (let c = 0; c <= 1200; c += 5) {
            ci.push(c);
            an.push(fvcb(c, T_C, PPFD, p, O2_mmol).An);
            if (showRef) an_ref.push(fvcb(c, T_C, PPFD, p, O2_AMBIENT).An);
        }

        const gamma_curr = fvcb(0, T_C, PPFD, p, O2_mmol).Gamma;
        const gamma_ref  = fvcb(0, T_C, PPFD, p, O2_AMBIENT).Gamma;

        const traces = [
            { x: ci, y: an, name: `O₂ = ${O2_pct}%`,
              line: { color: '#2d6a4f', width: 3 } },
        ];
        if (showRef && O2_pct !== 21) {
            traces.push({ x: ci, y: an_ref, name: 'O₂ = 21% (reference)',
                line: { color: '#aaa', width: 2, dash: 'dash' } });
        }
        traces.push({ x: [0, 1200], y: [0, 0], showlegend: false,
            line: { color: '#aaa', width: 1, dash: 'dot' } });

        const annotations = [{
            x: gamma_curr, y: 0,
            text: `Γ* = ${gamma_curr.toFixed(1)}`,
            showarrow: true, arrowhead: 2, ax: 40, ay: -35,
            font: { size: 11, color: '#c0392b' }, arrowcolor: '#c0392b',
        }];
        if (showRef && O2_pct !== 21) {
            annotations.push({
                x: gamma_ref, y: 0,
                text: `Γ*<sub>21%</sub> = ${gamma_ref.toFixed(1)}`,
                showarrow: true, arrowhead: 2, ax: -40, ay: -35,
                font: { size: 11, color: '#888' }, arrowcolor: '#888',
            });
        }

        Plotly.react('plot-pr', traces, mergeLayout({
            title: { text: 'Photorespiration: Effect of O₂ on A–C<sub>i</sub> Curve', font: { size: 15 } },
            xaxis: { title: 'C<sub>i</sub> (µmol/mol)', range: [0, 1200], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            annotations,
        }), PLOT_CFG);
    }

    draw();
    bindSlider('pr-o2', 'pr-o2-val', draw);
    bindSlider('pr-temp', 'pr-temp-val', draw);
    bindSlider('pr-ppfd', 'pr-ppfd-val', draw);
    bindCheckbox('pr-show-ref', draw);
}

// ==========================================================================
// TAB 7: C3 / C4 / CAM
// ==========================================================================
let mechView = 'mech-light';

function renderMechanisms() {
    const pills = document.querySelectorAll('#mechanisms .pill');
    pills.forEach(pill => {
        pill.addEventListener('click', () => {
            pills.forEach(p => p.classList.remove('active'));
            pill.classList.add('active');
            mechView = pill.dataset.view;
            drawMech();
        });
    });

    ['mech-c3', 'mech-c4', 'mech-cam'].forEach(id => bindCheckbox(id, drawMech));

    drawMech();
}

// Simplified C4 model
function c4_An(Ci, T_C, PPFD) {
    // PEPc concentrates CO2 ~10x at Rubisco site
    const Vcmax25 = 40;
    const Vpmax25 = 120; // PEPc capacity
    const Kp = 80;       // PEPc Km for CO2

    // PEPc
    const Vp = Vpmax25 * Ci / (Ci + Kp);

    // Rubisco in bundle sheath sees concentrated CO2
    const Ci_bs = Math.min(Ci * 10, 3000);  // CO2 concentrating
    const p = Object.assign({}, CROPS.tomato, { Vcmax25: Vcmax25 });
    const r = fvcb(Ci_bs, T_C, PPFD, p);

    // C4 limited by min(PEPc, Rubisco at elevated Ci)
    const Ac = r.An;
    const Ap = Vp - 1.5;  // PEPc limited with cost
    const An = Math.min(Ac, Ap);

    // C4 light response (higher saturation)
    const alpha_c4 = 0.06;
    const Amax_c4 = 45;
    const Rd_c4 = 1.0;
    const An_light = (alpha_c4 * PPFD * Amax_c4) / (alpha_c4 * PPFD + Amax_c4) - Rd_c4;

    return Math.min(An, An_light);
}

// Simplified C4 temperature response
function c4_An_temp(T_C, PPFD, Ca) {
    const Topt = 35;
    const Amax = 45;
    const sigma = 12;
    const scale = (PPFD / 2000) * (Ca / 400);
    return Amax * scale * Math.exp(-0.5 * ((T_C - Topt) / sigma) ** 2) - 1.0;
}

function drawMech() {
    const showC3  = document.getElementById('mech-c3').checked;
    const showC4  = document.getElementById('mech-c4').checked;
    const showCAM = document.getElementById('mech-cam').checked;

    const traces = [];

    if (mechView === 'mech-light') {
        // Light response comparison
        const ppfd = [];
        for (let p = 0; p <= 2500; p += 10) ppfd.push(p);

        if (showC3) {
            const an = ppfd.map(p => lightResponse(p, 25, 0.05, 1.5));
            traces.push({ x: ppfd, y: an, name: 'C3 (e.g., wheat)',
                line: { color: '#2d6a4f', width: 2.5 } });
        }
        if (showC4) {
            const an = ppfd.map(p => lightResponse(p, 45, 0.06, 1.0));
            traces.push({ x: ppfd, y: an, name: 'C4 (e.g., corn)',
                line: { color: '#e76f51', width: 2.5 } });
        }
        if (showCAM) {
            const an = ppfd.map(p => lightResponse(p, 12, 0.04, 0.8));
            traces.push({ x: ppfd, y: an, name: 'CAM (e.g., pineapple)',
                line: { color: '#264653', width: 2.5 } });
        }

        Plotly.react('plot-mech', traces, mergeLayout({
            title: { text: 'Light Response: C3 vs C4 vs CAM', font: { size: 15 } },
            xaxis: { title: 'PPFD (µmol / m² / s)', range: [0, 2500], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
        }), PLOT_CFG);

    } else if (mechView === 'mech-co2') {
        // CO2 response comparison
        const ca = [];
        for (let c = 0; c <= 1200; c += 5) ca.push(c);

        if (showC3) {
            const an = ca.map(c => fvcb(0.7 * c, 25, 1000, CROPS.tomato).An);
            traces.push({ x: ca, y: an, name: 'C3',
                line: { color: '#2d6a4f', width: 2.5 } });
        }
        if (showC4) {
            const an = ca.map(c => c4_An(0.4 * c, 25, 1000));
            traces.push({ x: ca, y: an, name: 'C4',
                line: { color: '#e76f51', width: 2.5 } });
        }
        if (showCAM) {
            // CAM similar to C3 but lower Amax
            const cam_p = Object.assign({}, CROPS.lettuce, { Vcmax25: 25, Jmax25: 50 });
            const an = ca.map(c => fvcb(0.7 * c, 25, 1000, cam_p).An);
            traces.push({ x: ca, y: an, name: 'CAM',
                line: { color: '#264653', width: 2.5 } });
        }

        // Mark ambient
        traces.push({ x: [0, 1200], y: [0, 0], showlegend: false,
            line: { color: '#aaa', width: 1, dash: 'dot' } });

        Plotly.react('plot-mech', traces, mergeLayout({
            title: { text: 'CO₂ Response: C3 vs C4 vs CAM', font: { size: 15 } },
            xaxis: { title: 'Ambient CO₂ (ppm)', range: [0, 1200], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            shapes: [{ type: 'line', x0: 420, x1: 420, y0: 0, y1: 1, yref: 'paper',
                line: { color: '#8e44ad', width: 1.5, dash: 'dot' } }],
            annotations: [{ x: 420, y: 1, yref: 'paper', text: 'Ambient (420 ppm)',
                showarrow: false, yanchor: 'bottom', font: { size: 10, color: '#8e44ad' } }],
        }), PLOT_CFG);

    } else if (mechView === 'mech-temp') {
        // Temperature response
        const temps = [];
        for (let t = 0; t <= 50; t += 0.5) temps.push(t);

        if (showC3) {
            const an = temps.map(t => fvcb(0.7 * 400, t, 1000, CROPS.tomato).An);
            traces.push({ x: temps, y: an, name: 'C3 (tomato)',
                line: { color: '#2d6a4f', width: 2.5 } });
        }
        if (showC4) {
            const an = temps.map(t => c4_An_temp(t, 1000, 400));
            traces.push({ x: temps, y: an, name: 'C4 (corn)',
                line: { color: '#e76f51', width: 2.5 } });
        }
        if (showCAM) {
            const cam_p = Object.assign({}, CROPS.lettuce, { Vcmax25: 20, Jmax25: 40 });
            const an = temps.map(t => fvcb(0.7 * 400, t, 500, cam_p).An);
            traces.push({ x: temps, y: an, name: 'CAM (low light adapted)',
                line: { color: '#264653', width: 2.5 } });
        }

        traces.push({ x: [0, 50], y: [0, 0], showlegend: false,
            line: { color: '#aaa', width: 1, dash: 'dot' } });

        Plotly.react('plot-mech', traces, mergeLayout({
            title: { text: 'Temperature Response: C3 vs C4 vs CAM', font: { size: 15 } },
            xaxis: { title: 'Temperature (°C)', range: [0, 50], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'A<sub>n</sub> (µmol CO₂ / m² / s)', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
        }), PLOT_CFG);

    } else if (mechView === 'mech-diurnal') {
        // CAM diurnal pattern
        const hours = [];
        for (let h = 0; h < 24; h += 0.25) hours.push(h);

        // Sunlight pattern (bell curve centered at 12)
        const light = hours.map(h => {
            if (h < 6 || h > 18) return 0;
            return 1500 * Math.exp(-0.5 * ((h - 12) / 3) ** 2);
        });

        if (showC3 || showC4) {
            // C3/C4: stomata open during day
            const gs_c3 = hours.map((h, i) => light[i] > 10 ? 0.3 * (light[i] / 1500) : 0.01);
            if (showC3) traces.push({ x: hours, y: gs_c3, name: 'C3/C4 stomata (g<sub>s</sub>)',
                line: { color: '#2d6a4f', width: 2.5 }, yaxis: 'y' });
        }
        if (showCAM) {
            // CAM: stomata open at night, closed during day
            const gs_cam = hours.map((h, i) => {
                if (h >= 20 || h <= 6) return 0.15 + 0.1 * Math.sin((h < 12 ? h + 4 : h - 20) / 10 * Math.PI);
                if (h > 6 && h < 8) return 0.15 * (1 - (h - 6) / 2);   // closing
                if (h > 18 && h < 20) return 0.15 * ((h - 18) / 2);     // opening
                return 0.02;  // closed during day
            });
            traces.push({ x: hours, y: gs_cam, name: 'CAM stomata (g<sub>s</sub>)',
                line: { color: '#264653', width: 2.5 }, yaxis: 'y' });
        }

        // Add light as secondary axis
        traces.push({ x: hours, y: light, name: 'PPFD',
            line: { color: '#f4a261', width: 1.5, dash: 'dot' }, yaxis: 'y2' });

        // Night shading
        const shapes = [
            { type: 'rect', xref: 'x', yref: 'paper', x0: 0, x1: 6, y0: 0, y1: 1,
              fillcolor: 'rgba(30,30,60,0.08)', line: { width: 0 }, layer: 'below' },
            { type: 'rect', xref: 'x', yref: 'paper', x0: 18, x1: 24, y0: 0, y1: 1,
              fillcolor: 'rgba(30,30,60,0.08)', line: { width: 0 }, layer: 'below' },
        ];

        Plotly.react('plot-mech', traces, mergeLayout({
            title: { text: 'Diurnal Stomatal Pattern: C3/C4 vs CAM', font: { size: 15 } },
            xaxis: { title: 'Hour of Day', range: [0, 24], dtick: 3, gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'g<sub>s</sub> (mol / m² / s)', rangemode: 'tozero', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis2: { title: 'PPFD (µmol/m²/s)', overlaying: 'y', side: 'right', showgrid: false, rangemode: 'tozero' },
            shapes,
            annotations: [
                { x: 3, y: 1, yref: 'paper', text: 'NIGHT', showarrow: false,
                  font: { size: 11, color: '#555' } },
                { x: 21, y: 1, yref: 'paper', text: 'NIGHT', showarrow: false,
                  font: { size: 11, color: '#555' } },
            ],
        }), PLOT_CFG);
    }
}

// ==========================================================================
// TAB 8: SOURCE-SINK
// ==========================================================================
const SS_CROPS = {
    lettuce: {
        name: 'Lettuce', cycle: 60, harvestable: 'leaves',
        stages: [
            { name: 'Seedling',    end: 0.15, root: 0.40, stem: 0.10, leaf: 0.50, fruit: 0.00 },
            { name: 'Vegetative',  end: 0.60, root: 0.20, stem: 0.10, leaf: 0.70, fruit: 0.00 },
            { name: 'Head fill',   end: 1.00, root: 0.10, stem: 0.05, leaf: 0.85, fruit: 0.00 },
        ],
    },
    tomato: {
        name: 'Tomato', cycle: 180, harvestable: 'fruit',
        stages: [
            { name: 'Seedling',    end: 0.10, root: 0.35, stem: 0.25, leaf: 0.40, fruit: 0.00 },
            { name: 'Vegetative',  end: 0.30, root: 0.20, stem: 0.25, leaf: 0.55, fruit: 0.00 },
            { name: 'Flowering',   end: 0.45, root: 0.12, stem: 0.18, leaf: 0.35, fruit: 0.35 },
            { name: 'Fruiting',    end: 0.80, root: 0.08, stem: 0.12, leaf: 0.20, fruit: 0.60 },
            { name: 'Late fruit',  end: 1.00, root: 0.05, stem: 0.10, leaf: 0.15, fruit: 0.70 },
        ],
    },
    strawberry: {
        name: 'Strawberry', cycle: 150, harvestable: 'fruit',
        stages: [
            { name: 'Establishment', end: 0.15, root: 0.35, stem: 0.20, leaf: 0.45, fruit: 0.00 },
            { name: 'Vegetative',    end: 0.35, root: 0.20, stem: 0.15, leaf: 0.65, fruit: 0.00 },
            { name: 'Flowering',     end: 0.50, root: 0.12, stem: 0.13, leaf: 0.35, fruit: 0.40 },
            { name: 'Fruiting',      end: 0.85, root: 0.08, stem: 0.10, leaf: 0.25, fruit: 0.57 },
            { name: 'Late',          end: 1.00, root: 0.06, stem: 0.09, leaf: 0.25, fruit: 0.60 },
        ],
    },
};

function getPartitioning(crop, dayFrac, photo) {
    const stages = crop.stages;
    // Interpolate base partitioning from stage table
    let prev = stages[0];
    let base;
    for (const s of stages) {
        if (dayFrac <= s.end) {
            const startFrac = prev === s ? 0 : prev.end;
            const t = (dayFrac - startFrac) / (s.end - startFrac);
            base = {
                stage: s.name,
                root: prev.root + t * (s.root - prev.root),
                stem: prev.stem + t * (s.stem - prev.stem),
                leaf: prev.leaf + t * (s.leaf - prev.leaf),
                fruit: prev.fruit + t * (s.fruit - prev.fruit),
            };
            break;
        }
        prev = s;
    }
    if (!base) {
        const last = stages[stages.length - 1];
        base = { stage: last.name, root: last.root, stem: last.stem, leaf: last.leaf, fruit: last.fruit };
    }

    // --- Carbon-availability modulation ---
    // Reference photosynthesis where the stage table applies as-is
    const photoRef = 30;  // g CH2O/m²/d
    // Minimum carbon for fruit set (below this, no fruit partitioning)
    const photoMinFruit = 10;

    if (base.fruit > 0 && photo < photoRef) {
        // Fruit sink strength declines under carbon limitation
        // Below photoMinFruit: fruit fraction → 0 (abortion / no fruit set)
        // Between min and ref: linear scaling
        const fruitScale = Math.max(0, (photo - photoMinFruit) / (photoRef - photoMinFruit));
        const fruitLost = base.fruit * (1 - fruitScale);
        base.fruit *= fruitScale;

        // Redistribute lost fruit carbon to vegetative organs
        // Under stress: roots get priority (foraging for resources)
        const toRoot = fruitLost * 0.45;
        const toLeaf = fruitLost * 0.40;
        const toStem = fruitLost * 0.15;
        base.root += toRoot;
        base.leaf += toLeaf;
        base.stem += toStem;
    }

    // Under very low carbon (< 15), boost root allocation further at expense of stem
    // (survival priority: root foraging)
    if (photo < 15) {
        const stressLevel = 1 - photo / 15;  // 0 at 15, 1 at 0
        const shift = Math.min(base.stem * 0.3 * stressLevel, 0.05);
        base.root += shift;
        base.stem -= shift;
    }

    // Normalize to ensure fractions sum to 1.0
    const total = base.root + base.stem + base.leaf + base.fruit;
    if (total > 0) {
        base.root /= total;
        base.stem /= total;
        base.leaf /= total;
        base.fruit /= total;
    }

    return base;
}

function renderSourceSink() {
    function draw() {
        const cropKey = document.getElementById('ss-crop').value;
        const photo   = parseFloat(document.getElementById('ss-photo').value);
        const day     = parseInt(document.getElementById('ss-day').value);

        const crop = SS_CROPS[cropKey];
        const cycle = crop.cycle;

        // Time series of cumulative partitioning
        const days = [], roots = [], stems = [], leaves = [], fruits = [];
        let cumRoot = 0, cumStem = 0, cumLeaf = 0, cumFruit = 0;

        for (let d = 1; d <= Math.min(day, 200); d++) {
            const frac = Math.min(d / cycle, 1.0);
            const part = getPartitioning(crop, frac, photo);
            // Simple growth: daily net = photo * growth_eff * (1 - respiration_frac)
            const net = photo * 0.7 * 0.75;  // growth efficiency × available
            cumRoot  += net * part.root;
            cumStem  += net * part.stem;
            cumLeaf  += net * part.leaf;
            cumFruit += net * part.fruit;
            days.push(d);
            roots.push(cumRoot);
            stems.push(cumStem);
            leaves.push(cumLeaf);
            fruits.push(cumFruit);
        }

        const traces1 = [
            { x: days, y: fruits, name: 'Fruit / Harvestable', fill: 'tonexty',
              fillcolor: 'rgba(231,111,81,0.5)', line: { color: '#e76f51', width: 0 },
              stackgroup: 'one' },
            { x: days, y: leaves, name: 'Leaves', fill: 'tonexty',
              fillcolor: 'rgba(45,106,79,0.5)', line: { color: '#2d6a4f', width: 0 },
              stackgroup: 'one' },
            { x: days, y: stems, name: 'Stems', fill: 'tonexty',
              fillcolor: 'rgba(168,162,120,0.5)', line: { color: '#a8a278', width: 0 },
              stackgroup: 'one' },
            { x: days, y: roots, name: 'Roots', fill: 'tonexty',
              fillcolor: 'rgba(139,90,43,0.5)', line: { color: '#8b5a2b', width: 0 },
              stackgroup: 'one' },
        ];

        // Add stage boundaries
        const shapes = [];
        let prevEnd = 0;
        const stageColors = ['rgba(200,200,200,0.05)', 'rgba(180,220,180,0.05)', 'rgba(255,240,200,0.05)', 'rgba(255,220,200,0.05)', 'rgba(255,200,200,0.05)'];
        crop.stages.forEach((s, i) => {
            const d0 = prevEnd * cycle;
            const d1 = s.end * cycle;
            if (d1 <= day) {
                shapes.push({ type: 'line', x0: d1, x1: d1, y0: 0, y1: 1, yref: 'paper',
                    line: { color: '#bbb', width: 1, dash: 'dot' } });
            }
            prevEnd = s.end;
        });

        Plotly.react('plot-ss-time', traces1, mergeLayout({
            title: { text: `${crop.name}: Cumulative Biomass Partitioning`, font: { size: 14 } },
            xaxis: { title: 'Days After Transplant', range: [1, Math.max(day, 10)], gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            yaxis: { title: 'Cumulative Dry Matter (g/m²)', rangemode: 'tozero', gridcolor: '#e8e8e8', zerolinecolor: '#ccc' },
            margin: { l: 60, r: 20, t: 40, b: 50 },
            legend: { x: 0.02, y: 0.98 },
            shapes,
        }), PLOT_CFG);

        // Pie chart for current day
        const frac = Math.min(day / cycle, 1.0);
        const part = getPartitioning(crop, frac, photo);

        const pieLabels = ['Roots', 'Stems', 'Leaves'];
        const pieValues = [part.root, part.stem, part.leaf];
        const pieColors = ['#8b5a2b', '#a8a278', '#2d6a4f'];

        if (crop.harvestable === 'fruit') {
            pieLabels.push('Fruit');
            pieValues.push(part.fruit);
            pieColors.push('#e76f51');
        }

        // Carbon status label
        const net = photo * 0.7 * 0.75;
        let carbonStatus = 'Surplus carbon';
        let statusColor = '#2d6a4f';
        if (photo < 10) { carbonStatus = 'Severe carbon deficit — no fruit set'; statusColor = '#c0392b'; }
        else if (photo < 15) { carbonStatus = 'Carbon limited — reduced fruiting'; statusColor = '#e67e22'; }
        else if (photo < 25) { carbonStatus = 'Moderate carbon — some fruit suppression'; statusColor = '#f39c12'; }

        Plotly.react('plot-ss-pie', [{
            labels: pieLabels,
            values: pieValues,
            type: 'pie',
            marker: { colors: pieColors },
            textinfo: 'label+percent',
            textfont: { size: 12 },
            hole: 0.40,
            sort: false,
        }], mergeLayout({
            title: { text: `Day ${day}: ${part.stage}`, font: { size: 14 } },
            margin: { l: 20, r: 20, t: 45, b: 30 },
            showlegend: false,
            annotations: [
                {
                    text: `${net.toFixed(1)}<br>g/m²/d`,
                    showarrow: false, font: { size: 13, color: '#333' },
                    x: 0.5, y: 0.52, xref: 'paper', yref: 'paper',
                },
                {
                    text: carbonStatus,
                    showarrow: false, font: { size: 11, color: statusColor },
                    x: 0.5, y: -0.02, xref: 'paper', yref: 'paper',
                },
            ],
        }), PLOT_CFG);
    }

    draw();
    document.getElementById('ss-crop').addEventListener('change', () => {
        // Update day slider max based on crop
        const crop = SS_CROPS[document.getElementById('ss-crop').value];
        document.getElementById('ss-day').max = crop.cycle + 20;
        draw();
    });
    bindSlider('ss-photo', 'ss-photo-val', draw);
    bindSlider('ss-day', 'ss-day-val', draw);
}

// ==========================================================================
// INIT
// ==========================================================================
document.addEventListener('DOMContentLoaded', () => {
    rendered.add('light-quality');
    renderLightQuality();
});
