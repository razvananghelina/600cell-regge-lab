/**
 * UI panel for Particle Physics Visualizer.
 *
 * Sections:
 * 1. Display Mode (coset / gauge / soliton / cooling)
 * 2. 4D View Rotation (6 sliders + auto-rotate)
 * 3. Coset Toggles (5 checkboxes)
 * 4. Gauge Toggles (U1 / SU2 / SU3)
 * 5. Soliton Controls (click to flip, pair mode, reset)
 * 6. Potts Cooling (temperature, speed, cool/reset/randomize)
 * 7. Info Readout (counts, energy, defects)
 */

import { PLANES } from './rotation4d.js';
import { COSET_NAMES, COSET_COLORS } from './colormap.js';
import { TOUR_STEPS } from './tour.js';

export class UI {
    constructor(app) {
        this.app = app;
        this.panel = document.getElementById('panel');
        this._build();
    }

    _build() {
        this.panel.innerHTML = '';

        // --- Display Mode ---
        this._addSection('Display Mode');
        this.modeSelect = this._addSelect('display-mode', [
            { value: 'tour',    label: '>>> Theory Tour (start here)' },
            { value: 'coset',   label: 'Cosets (5 x 24-cell)' },
            { value: 'gauge',   label: 'Gauge Channels' },
            { value: 'soliton', label: 'Soliton Creation' },
            { value: 'cooling', label: 'Kibble-Zurek Cooling' },
        ]);
        this.modeSelect.addEventListener('change', () => {
            this.app.setMode(this.modeSelect.value);
            this._updateVisibility();
        });

        // Mode description
        this.modeDesc = document.createElement('div');
        this.modeDesc.className = 'info-readout';
        this.modeDesc.style.marginTop = '4px';
        this.panel.appendChild(this.modeDesc);

        // --- 4D View Rotation ---
        this._addSection('4D View Rotation');
        this.sliders = {};
        for (const plane of PLANES) {
            const slider = this._addSlider(plane.toUpperCase(), -Math.PI, Math.PI, 0, 0.01);
            this.sliders[plane] = slider;
            slider.addEventListener('input', () => {
                this.app.rotation4d.setAngle(plane, parseFloat(slider.value));
            });
        }

        const autoDiv = document.createElement('div');
        autoDiv.className = 'control-row';
        const autoLabel = document.createElement('label');
        autoLabel.textContent = 'Auto-rot';
        const autoCheck = document.createElement('input');
        autoCheck.type = 'checkbox';
        autoCheck.checked = false;
        autoCheck.addEventListener('change', () => {
            this.app.autoRotate = autoCheck.checked;
        });
        autoDiv.appendChild(autoLabel);
        autoDiv.appendChild(autoCheck);
        this.panel.appendChild(autoDiv);

        this.speedSlider = this._addSlider('Speed', 0, 0.05, 0.005, 0.001);
        this.speedSlider.addEventListener('input', () => {
            this.app.autoRotateSpeed = parseFloat(this.speedSlider.value);
        });

        // --- Coset Toggles ---
        this.cosetHeader = this._addSection('Cosets');
        this.cosetSection = document.createElement('div');
        this.cosetSection.id = 'coset-section';
        this.cosetToggles = [];
        for (let c = 0; c < 5; c++) {
            const row = document.createElement('div');
            row.className = 'control-row';
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.checked = true;
            cb.dataset.coset = c;
            cb.addEventListener('change', () => {
                this.app.toggleCoset(c, cb.checked);
            });
            const lbl = document.createElement('label');
            lbl.textContent = `Coset ${c}`;
            const swatch = document.createElement('span');
            swatch.className = 'color-swatch';
            const col = COSET_COLORS[c];
            swatch.style.backgroundColor = `rgb(${col[0] * 255 | 0}, ${col[1] * 255 | 0}, ${col[2] * 255 | 0})`;
            const count = document.createElement('span');
            count.className = 'slider-value';
            count.textContent = `${COSET_NAMES[c]} (24)`;

            row.appendChild(cb);
            row.appendChild(swatch);
            row.appendChild(lbl);
            row.appendChild(count);
            this.cosetSection.appendChild(row);
            this.cosetToggles.push(cb);
        }
        this.panel.appendChild(this.cosetSection);

        // --- Gauge Toggles ---
        this.gaugeHeader = this._addSection('Gauge Channels');
        this.gaugeSection = document.createElement('div');
        this.gaugeSection.id = 'gauge-section';

        const gaugeItems = [
            { key: 'U1',  label: 'U(1)',  color: '#ffe633', count: '48' },
            { key: 'SU2', label: 'SU(2)', color: '#33ccff', count: '288' },
            { key: 'SU3', label: 'SU(3)', color: '#ff4d33', count: '384' },
        ];
        this.gaugeToggles = {};
        for (const g of gaugeItems) {
            const row = document.createElement('div');
            row.className = 'control-row';
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.checked = true;
            cb.addEventListener('change', () => {
                this.app.toggleGauge(g.key, cb.checked);
            });
            const swatch = document.createElement('span');
            swatch.className = 'color-swatch';
            swatch.style.backgroundColor = g.color;
            const lbl = document.createElement('label');
            lbl.textContent = g.label;
            const info = document.createElement('span');
            info.className = 'slider-value';
            info.textContent = `${g.count}`;
            info.id = `gauge-count-${g.key}`;

            row.appendChild(cb);
            row.appendChild(swatch);
            row.appendChild(lbl);
            row.appendChild(info);
            this.gaugeSection.appendChild(row);
            this.gaugeToggles[g.key] = cb;
        }
        this.panel.appendChild(this.gaugeSection);

        // --- Soliton Controls ---
        this.solitonHeader = this._addSection('Soliton Physics');
        this.solitonSection = document.createElement('div');
        this.solitonSection.id = 'soliton-section';

        // Explanation box
        const solitonExplain = document.createElement('div');
        solitonExplain.className = 'info-readout';
        solitonExplain.innerHTML = `
            <div style="color:#88aaff; margin-bottom:4px"><b>What is a soliton?</b></div>
            <div>In the ground state, every vertex has a "coset" color, and no two neighbors share the same color (proper 5-coloring).</div>
            <div style="margin-top:3px">A <b>soliton</b> = you change one vertex's color. Now it matches 3 of its 12 neighbors &rarr; 3 "broken" edges appear (shown in white).</div>
            <div style="margin-top:3px">This costs <b>energy E = 3</b>, always, for every vertex. The soliton is a topological defect - like a particle!</div>
            <div style="margin-top:3px; color:#6688aa">Look at the spectrum panel (top right): each eigenvalue level splits into 3 sub-levels when a soliton exists.</div>
        `;
        this.solitonSection.appendChild(solitonExplain);

        // Demo button - creates a soliton automatically
        const demoBtn = document.createElement('button');
        demoBtn.textContent = 'Create Demo Soliton';
        demoBtn.className = 'btn';
        demoBtn.style.background = '#2a3a5a';
        demoBtn.style.borderColor = '#4a6a9a';
        demoBtn.addEventListener('click', () => {
            this.app.createDemoSoliton();
        });
        this.solitonSection.appendChild(demoBtn);
        this.demoBtn = demoBtn;

        // Demo pair button
        const demoPairBtn = document.createElement('button');
        demoPairBtn.textContent = 'Create Demo Pair (2 adj. solitons)';
        demoPairBtn.className = 'btn';
        demoPairBtn.style.background = '#2a3a5a';
        demoPairBtn.style.borderColor = '#4a6a9a';
        demoPairBtn.addEventListener('click', () => {
            this.app.createDemoPair();
        });
        this.solitonSection.appendChild(demoPairBtn);

        // Manual mode section
        const manualLabel = document.createElement('div');
        manualLabel.style.cssText = 'color:#6666aa; font-size:10px; margin-top:10px; text-transform:uppercase; letter-spacing:1px;';
        manualLabel.textContent = 'Manual mode';
        this.solitonSection.appendChild(manualLabel);

        this.flipBtn = document.createElement('button');
        this.flipBtn.textContent = 'Click vertex to flip';
        this.flipBtn.className = 'btn';
        this.flipBtn.addEventListener('click', () => {
            this.app.toggleFlipMode();
            this.flipBtn.classList.toggle('active', this.app.flipMode);
            this.flipBtn.textContent = this.app.flipMode
                ? 'Now click any vertex on the 3D model...'
                : 'Click vertex to flip';
        });
        this.solitonSection.appendChild(this.flipBtn);

        this.pairBtn = document.createElement('button');
        this.pairBtn.textContent = 'Pair mode (click 2 neighbors)';
        this.pairBtn.className = 'btn';
        this.pairBtn.addEventListener('click', () => {
            this.app.togglePairMode();
            this.pairBtn.classList.toggle('active', this.app.pairMode);
            this.pairBtn.textContent = this.app.pairMode
                ? 'Click first vertex...'
                : 'Pair mode (click 2 neighbors)';
        });
        this.solitonSection.appendChild(this.pairBtn);

        const resetSolBtn = document.createElement('button');
        resetSolBtn.textContent = 'Reset (back to ground state)';
        resetSolBtn.className = 'btn';
        resetSolBtn.addEventListener('click', () => {
            this.app.resetSolitons();
            this.flipBtn.classList.remove('active');
            this.flipBtn.textContent = 'Click vertex to flip';
            this.pairBtn.classList.remove('active');
            this.pairBtn.textContent = 'Pair mode (click 2 neighbors)';
        });
        this.solitonSection.appendChild(resetSolBtn);

        // Soliton info (dynamic, updated each frame)
        this.solitonInfo = document.createElement('div');
        this.solitonInfo.className = 'info-readout';
        this.solitonInfo.id = 'soliton-info';
        this.solitonSection.appendChild(this.solitonInfo);

        this.panel.appendChild(this.solitonSection);

        // --- Potts Cooling Controls ---
        this.coolingHeader = this._addSection('Kibble-Zurek Cooling');
        this.coolingSection = document.createElement('div');
        this.coolingSection.id = 'cooling-section';

        this.tempSlider = this._addSliderTo(this.coolingSection, 'Temp', 0.01, 5.0, 5.0, 0.01);
        this.tempSlider.addEventListener('input', () => {
            this.app.setTemperature(parseFloat(this.tempSlider.value));
        });

        this.coolSpeedSlider = this._addSliderTo(this.coolingSection, 'Speed', 1, 20, 5, 1);

        const coolBtn = document.createElement('button');
        coolBtn.textContent = 'Start Cooling';
        coolBtn.className = 'btn';
        coolBtn.addEventListener('click', () => this.app.startCooling());
        this.coolingSection.appendChild(coolBtn);
        this.coolBtn = coolBtn;

        const randomBtn = document.createElement('button');
        randomBtn.textContent = 'Randomize (T=inf)';
        randomBtn.className = 'btn';
        randomBtn.addEventListener('click', () => this.app.randomizePotts());
        this.coolingSection.appendChild(randomBtn);

        const gsBtn = document.createElement('button');
        gsBtn.textContent = 'Ground State (T=0)';
        gsBtn.className = 'btn';
        gsBtn.addEventListener('click', () => this.app.setGroundState());
        this.coolingSection.appendChild(gsBtn);

        this.coolingInfo = document.createElement('div');
        this.coolingInfo.className = 'info-readout';
        this.coolingInfo.id = 'cooling-info';
        this.coolingSection.appendChild(this.coolingInfo);

        this.panel.appendChild(this.coolingSection);

        // --- Theory Tour ---
        this.tourHeader = this._addSection('Theory Tour');
        this.tourSection = document.createElement('div');
        this.tourSection.id = 'tour-section';

        // Step indicator
        this.tourStepIndicator = document.createElement('div');
        this.tourStepIndicator.style.cssText = 'text-align:center; color:#6688aa; font-size:11px; margin:4px 0;';
        this.tourSection.appendChild(this.tourStepIndicator);

        // Navigation buttons
        const navRow = document.createElement('div');
        navRow.style.cssText = 'display:flex; gap:6px; margin:6px 0;';

        this.tourPrevBtn = document.createElement('button');
        this.tourPrevBtn.textContent = 'Previous';
        this.tourPrevBtn.className = 'btn';
        this.tourPrevBtn.style.flex = '1';
        this.tourPrevBtn.addEventListener('click', () => this.app.prevTourStep());

        this.tourNextBtn = document.createElement('button');
        this.tourNextBtn.textContent = 'Next';
        this.tourNextBtn.className = 'btn';
        this.tourNextBtn.style.flex = '1';
        this.tourNextBtn.style.background = '#2a3a5a';
        this.tourNextBtn.style.borderColor = '#4a6a9a';
        this.tourNextBtn.addEventListener('click', () => this.app.nextTourStep());

        navRow.appendChild(this.tourPrevBtn);
        navRow.appendChild(this.tourNextBtn);
        this.tourSection.appendChild(navRow);

        // Tour title
        this.tourTitle = document.createElement('div');
        this.tourTitle.style.cssText = 'color:#88aaff; font-size:14px; font-weight:bold; margin:10px 0 6px 0; border-bottom:1px solid #2a2a4a; padding-bottom:4px;';
        this.tourSection.appendChild(this.tourTitle);

        // Tour text
        this.tourText = document.createElement('div');
        this.tourText.className = 'info-readout tour-text';
        this.tourSection.appendChild(this.tourText);

        this.panel.appendChild(this.tourSection);

        // --- Info Readout ---
        this._addSection('Info');
        this.infoDiv = document.createElement('div');
        this.infoDiv.className = 'info-readout';
        this.infoDiv.id = 'info-readout';
        this.panel.appendChild(this.infoDiv);

        // Initial visibility
        this._updateVisibility();
    }

    _updateVisibility() {
        const mode = this.modeSelect.value;

        // Show/hide section headers + content divs together
        const showCosets = (mode === 'coset' || mode === 'soliton');
        const showGauge = (mode === 'gauge');
        const showSoliton = (mode === 'soliton');
        const showCooling = (mode === 'cooling');
        const showTour = (mode === 'tour');

        this.tourHeader.style.display = showTour ? '' : 'none';
        this.tourSection.style.display = showTour ? '' : 'none';

        this.cosetHeader.style.display = showCosets ? '' : 'none';
        this.cosetSection.style.display = showCosets ? '' : 'none';

        this.gaugeHeader.style.display = showGauge ? '' : 'none';
        this.gaugeSection.style.display = showGauge ? '' : 'none';

        this.solitonHeader.style.display = showSoliton ? '' : 'none';
        this.solitonSection.style.display = showSoliton ? '' : 'none';

        this.coolingHeader.style.display = showCooling ? '' : 'none';
        this.coolingSection.style.display = showCooling ? '' : 'none';

        // Update mode description
        const descriptions = {
            tour: 'Step-by-step: how does the Standard Model emerge from 600-cell geometry? Use Next/Previous below.',
            coset: 'Vertices colored by coset (5 inscribed 24-cells). Edges by gauge type. Toggle cosets below.',
            gauge: 'Vertices by type: A (gauge, white), B (Higgs, gray), C (fermion, dim). Edges prominent by gauge channel.',
            soliton: 'Create topological defects on the 600-cell. Each defect = a "particle" with energy E=3. See the explanation below.',
            cooling: 'Antiferromagnetic Potts model. Randomize then cool to see ~88 frozen defects (Kibble-Zurek).',
        };
        this.modeDesc.innerHTML = `<small style="color:#6688aa">${descriptions[mode] || ''}</small>`;
    }

    /**
     * Returns the h3 element so we can toggle its visibility.
     */
    _addSection(title) {
        const h = document.createElement('h3');
        h.textContent = title;
        this.panel.appendChild(h);
        return h;
    }

    _addSlider(label, min, max, value, step) {
        return this._addSliderTo(this.panel, label, min, max, value, step);
    }

    _addSliderTo(container, label, min, max, value, step) {
        const row = document.createElement('div');
        row.className = 'control-row';
        const lbl = document.createElement('label');
        lbl.textContent = label;
        const input = document.createElement('input');
        input.type = 'range';
        input.min = min;
        input.max = max;
        input.value = value;
        input.step = step;
        const val = document.createElement('span');
        val.className = 'slider-value';
        val.textContent = parseFloat(value).toFixed(2);
        input.addEventListener('input', () => {
            val.textContent = parseFloat(input.value).toFixed(2);
        });
        row.appendChild(lbl);
        row.appendChild(input);
        row.appendChild(val);
        container.appendChild(row);
        return input;
    }

    _addSelect(id, options) {
        const sel = document.createElement('select');
        sel.id = id;
        for (const opt of options) {
            const o = document.createElement('option');
            o.value = opt.value;
            o.textContent = opt.label;
            sel.appendChild(o);
        }
        this.panel.appendChild(sel);
        return sel;
    }

    updateInfo(info) {
        this.infoDiv.innerHTML = `
            <div><b>600-cell</b></div>
            <div>V: ${info.vertices} | E: ${info.edges} | deg: 12</div>
            <div>FPS: ${info.fps}</div>
            <div>Mode: ${info.mode}</div>
            <hr>
            <div>Cosets: ${info.cosetValid ? '5x24 OK' : 'INVALID'}</div>
            <div>Edges: ${info.gaugeU1} U(1) + ${info.gaugeSU2} SU(2) + ${info.gaugeSU3} SU(3)</div>
            <div>Type A: ${info.typeA} | B: ${info.typeB} | C: ${info.typeC}</div>
        `;
    }

    updateSolitonInfo(info) {
        if (info.count === 0) {
            this.solitonInfo.innerHTML = `
                <hr>
                <div style="color:#6688aa">No solitons yet. Press "Create Demo Soliton" or use manual mode above.</div>
            `;
            return;
        }

        const lastFlipHtml = info.lastFlip ? `
            <hr>
            <div style="color:#88aaff"><b>Last action:</b></div>
            <div>Vertex #${info.lastFlip.vertex}: changed color</div>
            <div>Broken edges created: <b>${info.lastFlip.energy}</b></div>
            <div style="color:#6688aa">(${info.lastFlip.energy} neighbors now share same color = defect)</div>
        ` : '';

        this.solitonInfo.innerHTML = `
            <hr>
            <div style="color:#88aaff"><b>Status</b></div>
            <div>Active solitons: <b style="color:#ff8866">${info.count}</b> (pulsing vertices)</div>
            <div>Total energy: <b>${info.totalEnergy}</b> (broken edges in whole graph)</div>
            <div>White edges: ${info.brokenEdges} (same-color neighbors)</div>
            ${info.count === 1 && info.totalEnergy === 3 ? '<div style="color:#44aa66">E = 3 exactly. Always! This is topological.</div>' : ''}
            ${info.count === 2 && info.totalEnergy === 4 ? '<div style="color:#44aa66">E_pair = 4 < 2x3 = 6. Binding energy = 2!</div>' : ''}
            ${lastFlipHtml}
            <hr>
            <div style="color:#555577; font-size:10px">Spectrum panel (top-right) shows eigenvalue splitting: each level &rarr; 3 sub-levels (bulk/shell/core).</div>
        `;
    }

    updateCoolingInfo(info) {
        this.coolingInfo.innerHTML = `
            <hr>
            <div>Temperature: <b>${info.temperature.toFixed(2)}</b></div>
            <div>Energy: <b>${info.energy}</b></div>
            <div>Defects: <b>${info.defects}</b></div>
            <div>Sweeps: ${info.sweeps}</div>
        `;
        // Sync temp slider
        this.tempSlider.value = info.temperature;
        const valSpan = this.tempSlider.nextElementSibling;
        if (valSpan) valSpan.textContent = info.temperature.toFixed(2);
    }

    getCoolSpeed() {
        return parseInt(this.coolSpeedSlider.value);
    }

    updateTourStep(stepIndex, totalSteps, title, text) {
        this.tourStepIndicator.textContent = `Step ${stepIndex + 1} of ${totalSteps}`;
        this.tourTitle.textContent = title;
        this.tourText.innerHTML = text;

        // Disable prev on first step, next on last step
        this.tourPrevBtn.disabled = (stepIndex === 0);
        this.tourPrevBtn.style.opacity = (stepIndex === 0) ? '0.3' : '1';
        this.tourNextBtn.disabled = (stepIndex === totalSteps - 1);
        this.tourNextBtn.style.opacity = (stepIndex === totalSteps - 1) ? '0.3' : '1';
    }

    syncSliders(angles) {
        for (const plane of PLANES) {
            if (this.sliders[plane]) {
                this.sliders[plane].value = angles[plane];
                const valSpan = this.sliders[plane].nextElementSibling;
                if (valSpan) valSpan.textContent = angles[plane].toFixed(3);
            }
        }
    }
}
