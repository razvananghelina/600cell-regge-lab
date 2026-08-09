/**
 * PolytopeGraph: loads JSON polytope data and provides graph operations.
 */
export class PolytopeGraph {
    constructor() {
        this.vertices4D = [];   // Float64Array-like: [[x,y,z,w], ...]
        this.edges = [];        // [[i,j], ...]
        this.adjacency = [];    // [[], [], ...] neighbor lists
        this.numVertices = 0;
        this.numEdges = 0;
        this.degree = 0;
        this.graphDiameter = 0;
        this.name = '';
    }

    /**
     * Load polytope from JSON file.
     * @param {string} url - Path to JSON file
     * @returns {Promise<PolytopeGraph>}
     */
    async load(url) {
        const resp = await fetch(url);
        const data = await resp.json();

        this.name = data.name;
        this.vertices4D = data.vertices4D;
        this.edges = data.edges;
        this.adjacency = data.adjacency;
        this.numVertices = data.numVertices;
        this.numEdges = data.numEdges;
        this.degree = data.degree;
        this.graphDiameter = data.graphDiameter;

        return this;
    }

    /**
     * Get neighbors of vertex i.
     * @param {number} i
     * @returns {number[]}
     */
    getNeighbors(i) {
        return this.adjacency[i];
    }

    /**
     * BFS from source vertex, returns distance array.
     * @param {number} source
     * @returns {Int32Array} distances (-1 = unreachable)
     */
    bfs(source) {
        const dist = new Int32Array(this.numVertices).fill(-1);
        dist[source] = 0;
        const queue = [source];
        let head = 0;

        while (head < queue.length) {
            const u = queue[head++];
            for (const v of this.adjacency[u]) {
                if (dist[v] === -1) {
                    dist[v] = dist[u] + 1;
                    queue.push(v);
                }
            }
        }
        return dist;
    }

    /**
     * Get all edges as flat array [i0,j0, i1,j1, ...] for rendering.
     * @returns {Uint16Array}
     */
    getEdgeIndices() {
        const arr = new Uint16Array(this.numEdges * 2);
        for (let e = 0; e < this.numEdges; e++) {
            arr[e * 2] = this.edges[e][0];
            arr[e * 2 + 1] = this.edges[e][1];
        }
        return arr;
    }
}
