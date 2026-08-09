/**
 * UI panel: rotation sliders, polytope selector, STG controls, info readouts.
 */

import { PLANES } from './rotation4d.js';
import { PHI } from './constants.js';

export class UI {
    constructor(app) {
        this.app = app;
        this.panel = document.getElementById('panel');
        this._build();
    }

    _build() {
        this.panel.innerHTML = '';

        // --- Polytope selector ---
        this._addSection('Polytope');
        const select = this._addSelect('polytope-select', [
            { value: '600cell', label: '600-cell (120V, 720E)' },
            { value: '120cell', label: '120-cell (600V, 1200E)' },
        ]);
        select.addEventListener('change', () => this.app.loadPolytope(select.value));

        // --- 4D View Rotation (visualization only) ---
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
        autoLabel.textContent = 'Auto-rotate';
        const autoCheck = document.createElement('input');
        autoCheck.type = 'checkbox';
        autoCheck.id = 'auto-rotate';
        autoCheck.checked = true;
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

        // --- STG Controls ---
        this._addSection('STG (no free params)');

        // Mass slider (in Planck masses) - NOT alpha_G
        this.massSlider = this._addSlider('Mass (M_P)', 0.1, 5, 1.0, 0.1);

        const massBtn = document.createElement('button');
        massBtn.textContent = 'Click vertex to place mass';
        massBtn.className = 'btn';
        massBtn.addEventListener('click', () => {
            this.app.placeMassMode = !this.app.placeMassMode;
            massBtn.classList.toggle('active', this.app.placeMassMode);
            massBtn.textContent = this.app.placeMassMode
                ? 'Click vertex now...'
                : 'Click vertex to place mass';
        });
        this.panel.appendChild(massBtn);
        this.massModeBtn = massBtn;

        const clearBtn = document.createElement('button');
        clearBtn.textContent = 'Clear masses';
        clearBtn.className = 'btn';
        clearBtn.addEventListener('click', () => {
            this.app.stg.clearMasses();
            this.app.particle.active = false;
        });
        this.panel.appendChild(clearBtn);

        // Note about physics
        const note = document.createElement('div');
        note.className = 'info-readout';
        note.innerHTML = '<small>C(d) = 1 + M/d<br>a = nabla(C)/C^3<br>No free parameters (DERIVED)</small>';
        note.style.marginTop = '4px';
        note.style.color = '#6688aa';
        this.panel.appendChild(note);

        // --- Particle Controls ---
        this._addSection('Geodesic on S^3');

        const particleBtn = document.createElement('button');
        particleBtn.textContent = 'Launch particle';
        particleBtn.className = 'btn';
        particleBtn.addEventListener('click', () => {
            this.app.launchParticle();
        });
        this.panel.appendChild(particleBtn);

        this.dtSlider = this._addSlider('dt', 0.001, 0.02, 0.005, 0.001);
        this.dtSlider.addEventListener('input', () => {
            this.app.particle.dt = parseFloat(this.dtSlider.value);
        });

        // --- Info readouts ---
        this._addSection('Info');
        this.infoDiv = document.createElement('div');
        this.infoDiv.className = 'info-readout';
        this.infoDiv.id = 'info-readout';
        this.panel.appendChild(this.infoDiv);
    }

    /**
     * Get mass value from slider (in Planck units).
     */
    getMassValue() {
        return parseFloat(this.massSlider.value);
    }

    _addSection(title) {
        const h = document.createElement('h3');
        h.textContent = title;
        this.panel.appendChild(h);
    }

    _addSlider(label, min, max, value, step) {
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
        val.textContent = parseFloat(value).toFixed(3);
        input.addEventListener('input', () => {
            val.textContent = parseFloat(input.value).toFixed(3);
        });

        row.appendChild(lbl);
        row.appendChild(input);
        row.appendChild(val);
        this.panel.appendChild(row);

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
            <div><b>${info.name || ''}</b></div>
            <div>V: ${info.vertices} | E: ${info.edges} | deg: ${info.degree}</div>
            <div>Diameter: ${info.diameter}</div>
            <div>FPS: ${info.fps}</div>
            <hr>
            <div>RCS (chi): ${info.chi} ${info.chiValid ? 'OK' : 'FAIL'}</div>
            <div>Masses: ${info.massCount}</div>
            <div>C range: [${info.cMin}, ${info.cMax}]</div>
            <div>TDI drift: ${info.drift} / ${info.driftThreshold}</div>
            ${info.driftAlert ? '<div class="alert">TDI ALERT: drift > ln(phi)</div>' : ''}
            <hr>
            <div>Particle: ${info.particleActive ? 'ACTIVE' : 'idle'}</div>
            ${info.particleActive ? `<div>v = ${info.particleSpeed}</div>` : ''}
            ${info.particleActive ? `<div>Nearest vertex: ${info.particleNearest}</div>` : ''}
        `;
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
