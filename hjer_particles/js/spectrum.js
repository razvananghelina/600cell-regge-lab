/**
 * Spectral structure of the 600-cell adjacency matrix.
 *
 * 9 eigenvalues, ALL of the form a + b*phi:
 *   12, 6phi, 4phi, 3, 0, -2, 4-4phi, -3, 6-6phi
 *
 * Multiplicities ALL perfect squares:
 *   1, 4, 9, 16, 25, 36, 9, 16, 4
 *
 * Sum of multiplicities = 120 (number of vertices).
 *
 * Key physics:
 * - Phi-sector: dimensions 4+9+9+4 = 26 (sin^2(tW) denominator)
 * - Sign change at d=3 (confinement-like)
 * - Each soliton splits every eigenvalue group into (m-3, 2, 1) sub-levels
 */

import { PHI } from './constants.js';

/**
 * Pre-computed eigenvalue data.
 */
export const EIGENVALUES = [
    { value: 12,             mult: 1,  label: '12',       isPhiSector: false },
    { value: 6 * PHI,        mult: 4,  label: '6phi',     isPhiSector: true  },
    { value: 4 * PHI,        mult: 9,  label: '4phi',     isPhiSector: true  },
    { value: 3,              mult: 16, label: '3',        isPhiSector: false },
    { value: 0,              mult: 25, label: '0',        isPhiSector: false },
    { value: -2,             mult: 36, label: '-2',       isPhiSector: false },
    { value: 4 - 4 * PHI,    mult: 9,  label: '4-4phi',   isPhiSector: true  },
    { value: -3,             mult: 16, label: '-3',       isPhiSector: false },
    { value: 6 - 6 * PHI,    mult: 4,  label: '6-6phi',   isPhiSector: true  },
];

/**
 * Total dimensions.
 */
export const TOTAL_DIM = 120;
export const PHI_SECTOR_DIM = 26; // 4 + 9 + 9 + 4

/**
 * When a soliton is created, each eigenvalue group with multiplicity m
 * splits into sub-levels: (m-3, 2, 1).
 * Exception: m=1 doesn't split (stays as 1).
 *
 * @param {number} mult - original multiplicity
 * @returns {number[]} sub-level multiplicities
 */
export function computeSplitting(mult) {
    if (mult <= 1) return [mult];
    if (mult <= 3) return [mult - 2, 1, 1];
    return [mult - 3, 2, 1];
}

/**
 * Get splitting colors for sub-levels.
 * bulk = main (stays in place), shell = doublet, core = singlet
 */
export const SPLIT_COLORS = {
    bulk:  [0.3, 0.5, 0.9],  // blue
    shell: [0.2, 0.8, 0.8],  // cyan
    core:  [0.9, 0.3, 0.2],  // red
};

/**
 * Compute the spectral data for canvas rendering.
 * Returns arrays ready for bar chart display.
 *
 * @param {boolean} solitonActive - whether to show splitting
 * @param {number} solitonCount - number of solitons (affects splitting depth)
 * @returns {Array<{value: number, mult: number, label: string, subLevels: number[], isPhiSector: boolean}>}
 */
export function getSpectralData(solitonActive, solitonCount = 0) {
    return EIGENVALUES.map(ev => ({
        ...ev,
        subLevels: solitonActive && solitonCount > 0
            ? computeSplitting(ev.mult)
            : [ev.mult],
    }));
}

/**
 * Draw the spectral panel on a canvas element.
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} width
 * @param {number} height
 * @param {boolean} solitonActive
 * @param {number} solitonCount
 * @param {number} animProgress - 0 to 1 for splitting animation
 */
export function drawSpectrum(ctx, width, height, solitonActive, solitonCount, animProgress = 1) {
    ctx.clearRect(0, 0, width, height);

    // Background
    ctx.fillStyle = '#0a0a1a';
    ctx.fillRect(0, 0, width, height);

    const data = getSpectralData(solitonActive, solitonCount);
    const padding = { top: 30, bottom: 30, left: 60, right: 15 };
    const plotW = width - padding.left - padding.right;
    const plotH = height - padding.top - padding.bottom;

    // Value range: from 6-6phi (~-3.7) to 12
    const minVal = 6 - 6 * PHI - 1;
    const maxVal = 13;
    const valRange = maxVal - minVal;

    // Title
    ctx.fillStyle = '#88aaff';
    ctx.font = '11px Courier New';
    ctx.textAlign = 'center';
    ctx.fillText('EIGENVALUE SPECTRUM', width / 2, 14);

    // Axis line
    ctx.strokeStyle = '#2a2a4a';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, height - padding.bottom);
    ctx.stroke();

    // Zero line
    const zeroY = padding.top + plotH * (1 - (0 - minVal) / valRange);
    ctx.strokeStyle = '#3a3a5a';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(padding.left, zeroY);
    ctx.lineTo(width - padding.right, zeroY);
    ctx.stroke();
    ctx.setLineDash([]);

    // Draw each eigenvalue
    const maxMult = 36;
    const barHeight = 4;

    for (const ev of data) {
        const y = padding.top + plotH * (1 - (ev.value - minVal) / valRange);

        if (ev.subLevels.length === 1 || animProgress < 0.01) {
            // Single bar
            const barW = (ev.mult / maxMult) * plotW;
            ctx.fillStyle = ev.isPhiSector
                ? 'rgba(180, 140, 255, 0.8)'
                : 'rgba(100, 140, 200, 0.7)';
            ctx.fillRect(padding.left + 2, y - barHeight / 2, barW, barHeight);
        } else {
            // Split bars
            const splitOffset = animProgress * 6; // pixels between sub-levels
            const colors = [
                SPLIT_COLORS.bulk,
                SPLIT_COLORS.shell,
                SPLIT_COLORS.core,
            ];

            for (let s = 0; s < ev.subLevels.length; s++) {
                const subMult = ev.subLevels[s];
                if (subMult === 0) continue;
                const barW = (subMult / maxMult) * plotW;
                const offsetY = (s - (ev.subLevels.length - 1) / 2) * splitOffset;
                const c = colors[s] || colors[0];
                ctx.fillStyle = `rgba(${c[0] * 255 | 0}, ${c[1] * 255 | 0}, ${c[2] * 255 | 0}, 0.8)`;
                ctx.fillRect(padding.left + 2, y - barHeight / 2 + offsetY, barW, barHeight);
            }
        }

        // Label
        ctx.fillStyle = '#8888aa';
        ctx.font = '9px Courier New';
        ctx.textAlign = 'right';
        ctx.fillText(ev.label, padding.left - 4, y + 3);

        // Multiplicity
        ctx.textAlign = 'left';
        const totalBarW = (ev.mult / maxMult) * plotW;
        ctx.fillStyle = '#6666aa';
        ctx.fillText(`x${ev.mult}`, padding.left + totalBarW + 6, y + 3);
    }

    // Legend
    const legendY = height - 14;
    ctx.font = '9px Courier New';
    ctx.textAlign = 'left';

    if (solitonActive && solitonCount > 0) {
        const items = [
            { color: SPLIT_COLORS.bulk, label: 'bulk' },
            { color: SPLIT_COLORS.shell, label: 'shell' },
            { color: SPLIT_COLORS.core, label: 'core' },
        ];
        let x = 10;
        for (const item of items) {
            ctx.fillStyle = `rgb(${item.color[0] * 255 | 0}, ${item.color[1] * 255 | 0}, ${item.color[2] * 255 | 0})`;
            ctx.fillRect(x, legendY - 6, 8, 8);
            ctx.fillStyle = '#8888aa';
            ctx.fillText(item.label, x + 12, legendY + 1);
            x += 55;
        }
    } else {
        ctx.fillStyle = '#6666aa';
        ctx.fillText('phi-sector: purple | integer: blue', 10, legendY + 1);
    }
}
