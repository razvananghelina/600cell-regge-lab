/**
 * Coset decomposition of 600-cell into 5 inscribed 24-cells.
 *
 * The 120 vertices of the 600-cell = binary icosahedral group 2I.
 * The 24 vertices of a 24-cell = binary tetrahedral group 2T (subgroup of 2I).
 * Index [2I : 2T] = 5, giving 5 cosets = 5 inscribed 24-cells.
 *
 * Vertex types:
 *   A-type (8):  permutations of (+-1, 0, 0, 0)
 *   B-type (16): all (+-1/2, +-1/2, +-1/2, +-1/2)
 *   C-type (96): even permutations of (+-1/2, +-phi/2, +-1/(2*phi), 0)
 *
 * Key result: ALL 720 edges are inter-coset (5-coloring is proper).
 */

import { PHI } from './constants.js';

const PHI_INV = 1 / PHI;
const EPS = 1e-6;

/**
 * Classify a vertex as A, B, or C type based on coordinates.
 * @param {number[]} v - [x,y,z,w]
 * @returns {'A'|'B'|'C'}
 */
function classifyVertex(v) {
    const abs = v.map(c => Math.abs(c));
    abs.sort((a, b) => b - a);  // descending

    // A-type: one coord ~1, rest ~0
    if (Math.abs(abs[0] - 1.0) < EPS &&
        abs[1] < EPS && abs[2] < EPS && abs[3] < EPS) {
        return 'A';
    }

    // B-type: all coords ~0.5
    if (Math.abs(abs[0] - 0.5) < EPS &&
        Math.abs(abs[1] - 0.5) < EPS &&
        Math.abs(abs[2] - 0.5) < EPS &&
        Math.abs(abs[3] - 0.5) < EPS) {
        return 'B';
    }

    // C-type: everything else (golden ratio coords)
    return 'C';
}

/**
 * Quaternion multiplication: q1 * q2
 * Quaternion as [x, y, z, w] where w is real part.
 * Convention: q = w + x*i + y*j + z*k
 */
function qmul(q1, q2) {
    const [x1, y1, z1, w1] = q1;
    const [x2, y2, z2, w2] = q2;
    return [
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
    ];
}

/**
 * Find the closest vertex to a given 4D point.
 */
function findClosestVertex(vertices, target) {
    let bestIdx = -1;
    let bestDot = -Infinity;
    for (let i = 0; i < vertices.length; i++) {
        const v = vertices[i];
        const dot = v[0] * target[0] + v[1] * target[1] +
                    v[2] * target[2] + v[3] * target[3];
        if (dot > bestDot) {
            bestDot = dot;
            bestIdx = i;
        }
    }
    return bestIdx;
}

/**
 * Compute the 5 cosets of 2T in 2I.
 * Returns coset assignment for each vertex (0-4).
 *
 * Method: Start with coset 0 = {identity vertex + its 2T orbit}.
 * The binary tetrahedral group 2T has 24 elements.
 * We pick the 24-cell containing vertex 0, then use quaternion
 * multiplication by coset representatives to find the other 4 cosets.
 *
 * @param {Array<number[]>} vertices - 120 vertices of 600-cell
 * @param {Array<number[]>} adjacency - adjacency lists
 * @returns {number[]} coset[i] = 0..4 for each vertex
 */
export function computeCosets(vertices, adjacency) {
    const n = vertices.length;
    const coset = new Array(n).fill(-1);

    // The 24 vertices of a 24-cell inscribed in the 600-cell form
    // the binary tetrahedral group 2T. We identify coset 0 by finding
    // vertices that form a 24-cell (degree 8 subgraph, distance-regular).

    // Strategy: Pick vertex 0 as seed. Find which 24-cell it belongs to
    // by checking the induced subgraph on vertices at even graph distances.
    // Actually, we use the quaternion method: 2T consists of the 8 A-type
    // vertices and 16 B-type vertices. So coset 0 = A + B.

    // Find A and B type vertices
    const types = vertices.map(v => classifyVertex(v));
    const aVertices = [];
    const bVertices = [];
    const cVertices = [];

    for (let i = 0; i < n; i++) {
        if (types[i] === 'A') aVertices.push(i);
        else if (types[i] === 'B') bVertices.push(i);
        else cVertices.push(i);
    }

    // Coset 0 = the 24-cell formed by A + B vertices (24 total)
    // This is 2T itself.
    const coset0 = [...aVertices, ...bVertices];
    for (const idx of coset0) {
        coset[idx] = 0;
    }

    // The remaining 96 C-type vertices form 4 cosets of 24 each.
    // Use quaternion multiplication: pick a C-type vertex as coset representative,
    // multiply all coset0 vertices by it to get the next coset.

    // Find coset representatives: pick unassigned C-vertices
    const representatives = [];
    const assigned = new Set(coset0);

    for (let i = 0; i < n; i++) {
        if (assigned.has(i)) continue;

        // This vertex starts a new coset
        representatives.push(i);

        // Generate the full coset: {rep * t : t in 2T}
        // In quaternion terms: left-multiply all 2T elements by rep
        const rep = vertices[i];
        const cosetMembers = new Set();

        for (const t of coset0) {
            const tVert = vertices[t];
            const product = qmul(rep, tVert);
            // Find closest vertex to this product
            const closest = findClosestVertex(vertices, product);
            cosetMembers.add(closest);
        }

        // Also try right multiplication and take the union approach
        // Actually, for coset decomposition we need left cosets: g * 2T
        // But let's also verify with right: 2T * g
        if (cosetMembers.size < 24) {
            // Try right multiplication instead
            cosetMembers.clear();
            for (const t of coset0) {
                const tVert = vertices[t];
                const product = qmul(tVert, rep);
                const closest = findClosestVertex(vertices, product);
                cosetMembers.add(closest);
            }
        }

        const cosetIdx = representatives.length; // 1, 2, 3, 4

        for (const idx of cosetMembers) {
            if (!assigned.has(idx)) {
                coset[idx] = cosetIdx;
                assigned.add(idx);
            }
        }

        if (assigned.size >= n) break;
    }

    // Fallback: if any vertices still unassigned, use BFS-based coloring
    // This shouldn't happen with correct quaternion multiplication
    let nextCoset = representatives.length + 1;
    for (let i = 0; i < n; i++) {
        if (coset[i] === -1) {
            // Find nearest assigned vertex in same would-be coset
            // Simple fallback: assign remaining by greedy proper coloring
            const neighborCosets = new Set();
            for (const nb of adjacency[i]) {
                if (coset[nb] >= 0) neighborCosets.add(coset[nb]);
            }
            // Since all edges are inter-coset, pick any coset not used by neighbors
            for (let c = 0; c < 5; c++) {
                if (!neighborCosets.has(c)) {
                    coset[i] = c;
                    break;
                }
            }
        }
    }

    return { coset, types };
}

/**
 * Verify coset properties.
 * @returns {{ valid: boolean, counts: number[], allInterCoset: boolean }}
 */
export function verifyCosets(coset, types, edges) {
    // Count per coset
    const counts = [0, 0, 0, 0, 0];
    for (const c of coset) counts[c]++;

    // Check all edges are inter-coset (proper 5-coloring)
    let allInterCoset = true;
    for (const [i, j] of edges) {
        if (coset[i] === coset[j]) {
            allInterCoset = false;
            break;
        }
    }

    // Check each coset has 24 vertices
    const valid = counts.every(c => c === 24) && allInterCoset;

    return { valid, counts, allInterCoset };
}
