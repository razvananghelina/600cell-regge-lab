/**
 * Gauge edge classification for 600-cell.
 *
 * Each of 720 edges classified by shared-neighbor profile:
 *   U(1):  48 CC edges with (shared_A=1, shared_B=0)
 *   SU(2): 288 edges = AC(96) + BC(192)
 *   SU(3): 384 CC edges (remaining after U(1))
 *
 * Key results:
 * - 1+3+8 = 12 = dim(SM gauge group) = vertex degree
 * - U(1) edges form perfect matching (degree 1 per C-vertex)
 * - SU(3) has 288 pure plaquettes (only gauge type with plaquettes)
 * - ALL 720 edges are inter-coset
 *
 * CC shared-neighbor profiles:
 *   (1,0): 48  -> U(1)
 *   (0,1): 96  -> SU(3) T1
 *   (0,2): 96  -> SU(3) T2
 *   (1,1): 192 -> SU(3) T3
 */

/**
 * For each edge (i,j), count shared neighbors by vertex type.
 */
function countSharedNeighbors(i, j, adjacency, types) {
    const nbI = new Set(adjacency[i]);
    let sharedA = 0, sharedB = 0, sharedC = 0;

    for (const nb of adjacency[j]) {
        if (nbI.has(nb)) {
            if (types[nb] === 'A') sharedA++;
            else if (types[nb] === 'B') sharedB++;
            else sharedC++;
        }
    }

    return { sharedA, sharedB, sharedC };
}

/**
 * Classify all 720 edges into gauge types.
 *
 * Classification:
 *   U(1):  CC edges with shared_A=1, shared_B=0 (48 edges)
 *   SU(2): All AC edges (96) + all BC edges (192) = 288
 *   SU(3): All remaining CC edges (384 = 96+96+192 by profile)
 *
 * @param {Array<number[]>} edges
 * @param {Array<number[]>} adjacency
 * @param {string[]} types - vertex types
 * @returns {{ edgeType: string[], edgeSubType: string[], counts: object, subCounts: object }}
 */
export function classifyEdges(edges, adjacency, types) {
    const edgeType = [];
    const edgeSubType = [];
    const counts = { U1: 0, SU2: 0, SU3: 0 };
    const subCounts = {};

    for (let e = 0; e < edges.length; e++) {
        const [i, j] = edges[e];
        const ti = types[i];
        const tj = types[j];
        const pairType = [ti, tj].sort().join('');
        const shared = countSharedNeighbors(i, j, adjacency, types);

        let gauge, sub;

        // No AA, AB, or BB edges exist (these are intra-coset in 24-cell)
        if (pairType === 'CC' && shared.sharedA === 1 && shared.sharedB === 0) {
            // U(1): CC edges with exactly 1 shared A neighbor, 0 shared B
            gauge = 'U1';
            sub = 'U1_CC';
        } else if (pairType === 'AC') {
            // SU(2): A-C edges
            gauge = 'SU2';
            sub = 'SU2_AC';
        } else if (pairType === 'BC') {
            // SU(2): B-C edges
            gauge = 'SU2';
            sub = 'SU2_BC';
        } else {
            // SU(3): all remaining CC edges (profiles: (0,1), (0,2), (1,1))
            gauge = 'SU3';
            const key = `${shared.sharedA}${shared.sharedB}`;
            if (key === '01') sub = 'SU3_T1';       // (0,1): 96 edges
            else if (key === '02') sub = 'SU3_T2';   // (0,2): 96 edges
            else if (key === '11') sub = 'SU3_T3';   // (1,1): 192 edges
            else sub = `SU3_${key}`;                  // unexpected profile
        }

        edgeType.push(gauge);
        edgeSubType.push(sub);
        counts[gauge]++;
        subCounts[sub] = (subCounts[sub] || 0) + 1;
    }

    return { edgeType, edgeSubType, counts, subCounts };
}

/**
 * Get color for each gauge type.
 * @param {string} gaugeType - 'U1', 'SU2', 'SU3'
 * @returns {[number, number, number]} RGB in [0,1]
 */
export function gaugeColor(gaugeType) {
    switch (gaugeType) {
        case 'U1':  return [1.0, 0.9, 0.2];   // yellow
        case 'SU2': return [0.2, 0.8, 1.0];   // cyan
        case 'SU3': return [1.0, 0.3, 0.2];   // red
        default:    return [0.5, 0.5, 0.5];   // gray
    }
}

/**
 * Get the set of edges for a specific gauge type.
 * @param {Array<number[]>} edges
 * @param {string[]} edgeType
 * @param {string} type - 'U1', 'SU2', 'SU3'
 * @returns {number[]} indices into edges array
 */
export function getEdgesByType(edges, edgeType, type) {
    const result = [];
    for (let e = 0; e < edges.length; e++) {
        if (edgeType[e] === type) result.push(e);
    }
    return result;
}
