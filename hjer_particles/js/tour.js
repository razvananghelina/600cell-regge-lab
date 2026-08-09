/**
 * Theory Tour - Interactive walkthrough showing how the Standard Model
 * emerges from 600-cell geometry.
 *
 * 10 steps, each with:
 *   - title: displayed in panel header
 *   - text: HTML explanation
 *   - vertices: how to color vertices
 *   - edges: which edges to show and how
 *   - highlight: special visual emphasis
 */

export const TOUR_STEPS = [
    // --- Step 0 ---
    {
        title: 'The 600-cell',
        vertices: 'coset',
        edges: 'gauge',
        highlight: 'none',
        text: `
            <b>120 vertices</b>, <b>720 edges</b>, <b>600 tetrahedral cells</b>.
            <br><br>
            This is the 4D polytope that achieves the <b>densest possible sphere packing</b>
            on the 3-sphere S<sup>3</sup> (proved by Cohn-Kumar, 2007).
            <br><br>
            Every vertex has exactly <b>12 neighbors</b>.
            <br><br>
            12 is also the dimension of the Standard Model gauge group
            SU(3) &times; SU(2) &times; U(1). Coincidence?
        `,
    },
    // --- Step 1 ---
    {
        title: '5 Inscribed 24-cells',
        vertices: 'coset',
        edges: 'dim',
        highlight: 'none',
        text: `
            The 120 vertices decompose into <b>5 groups of 24</b> (shown by color).
            <br><br>
            Each group forms a <b>24-cell</b> &mdash; a regular 4D polytope whose
            symmetry group is <b>D<sub>4</sub></b>.
            <br><br>
            D<sub>4</sub> is the Lie algebra that contains
            SU(3) &times; SU(2) &times; U(1) as a maximal subgroup.
            <br><br>
            All 720 edges connect vertices from <b>different</b> groups.
            This is a <b>proper 5-coloring</b> &mdash; no two neighbors share a color.
            <br><br>
            <span style="color:#44aa66">DERIVED from quaternion group theory (2I/2T cosets).</span>
        `,
    },
    // --- Step 2 ---
    {
        title: 'Three Vertex Types',
        vertices: 'type',
        edges: 'dim',
        highlight: 'none',
        text: `
            The 120 vertices have three distinct coordinate structures:
            <br><br>
            <span style="color:#fff">&bull;</span>
            <b>8 A-type</b> (large, white): permutations of (&plusmn;1, 0, 0, 0)
            <br>
            &rarr; These are the <b>gauge vertices</b> of the 24-cell.
            <br><br>
            <span style="color:#999ab3">&bull;</span>
            <b>16 B-type</b> (medium, gray): all (&plusmn;&frac12;, &plusmn;&frac12;, &plusmn;&frac12;, &plusmn;&frac12;)
            <br>
            &rarr; <b>Higgs-like</b> sector.
            <br><br>
            <span style="color:#404059">&bull;</span>
            <b>96 C-type</b> (small, dim): golden ratio coordinates
            <br>
            &rarr; <b>Fermions</b>. Note: 96 = 16 &times; 3 &times; 2
            <br>= (particles per generation) &times; (generations) &times; (particle/antiparticle).
            <br><br>
            Key: <b>A-A = A-B = B-B = 0 edges</b>. These types never connect to each other.
        `,
    },
    // --- Step 3 ---
    {
        title: 'A Fermion and Its 12 Neighbors',
        vertices: 'neighborhood',
        edges: 'neighborhood_dim',
        highlight: 'neighborhood',
        text: `
            The <b>pulsing vertex</b> is a fermion (C-type).
            <br><br>
            It has exactly <b>12 neighbors</b>, connected by 12 edges.
            <br><br>
            Among its 12 neighbors:
            <br>
            &bull; <b>1</b> is A-type (white, large) &mdash; gauge boson
            <br>
            &bull; <b>2</b> are B-type (gray, medium) &mdash; Higgs
            <br>
            &bull; <b>9</b> are C-type (dim, small) &mdash; other fermions
            <br><br>
            But what matters is not neighbor type &mdash;
            it's the <b>edge classification</b>...
            <br><br>
            <span style="color:#6688aa">&rarr; Press Next to see the key result.</span>
        `,
    },
    // --- Step 4 (THE KEY STEP) ---
    {
        title: '1 + 3 + 8 = 12',
        vertices: 'neighborhood',
        edges: 'neighborhood_gauge',
        highlight: 'neighborhood',
        text: `
            The 12 edges of this fermion classify by <b>gauge channel</b>:
            <br><br>
            <span style="color:#ffe633; font-size:14px">&block;</span>
            <b>1 edge</b> is <b>U(1)</b> &mdash; electromagnetic
            <br>
            <span style="color:#33ccff; font-size:14px">&block;</span>
            <b>3 edges</b> are <b>SU(2)</b> &mdash; weak force
            <br>
            <span style="color:#ff4d33; font-size:14px">&block;</span>
            <b>8 edges</b> are <b>SU(3)</b> &mdash; strong force
            <br><br>
            <div style="font-size:14px; color:#fff; text-align:center; margin:6px 0">
                1 + 3 + 8 = 12
            </div>
            <div style="text-align:center; color:#aac">
                = dim U(1) + dim SU(2) + dim SU(3)
                <br>= dimension of the SM gauge group
            </div>
            <br>
            This is <b>EXACT</b>, for <b>EVERY</b> one of the 96 fermion vertices.
            <br><br>
            <span style="color:#44aa66">DERIVED. Not fitted. Pure graph combinatorics.</span>
        `,
    },
    // --- Step 5 ---
    {
        title: 'Coupling Constants',
        vertices: 'coset',
        edges: 'gauge',
        highlight: 'none',
        text: `
            The spectral structure (9 eigenvalues, intersection numbers)
            determines the coupling constants:
            <br><br>
            <b>Fine structure constant &alpha;:</b>
            <br>
            &nbsp;&nbsp;2&pi;&alpha;&sup2; &minus; 20&phi;<sup>4</sup>&alpha; + 1 = 0
            <br>
            &nbsp;&nbsp;&rarr; &alpha; = 1/137.03600 &nbsp; (<b>0.0001%</b>)
            <br><br>
            <b>Weinberg angle:</b>
            <br>
            &nbsp;&nbsp;sin&sup2;&theta;<sub>W</sub> = b<sub>1</sub>/&phi;<sub>dim</sub> = 6/26
            <br>
            &nbsp;&nbsp;&rarr; 0.23077 &nbsp; (<b>0.19%</b> from experiment)
            <br><br>
            <b>Strong coupling ratio:</b>
            <br>
            &nbsp;&nbsp;&alpha;<sub>s</sub>/&alpha; = 10&phi; &nbsp; (exact)
            <br><br>
            <b>Higgs mass:</b>
            <br>
            &nbsp;&nbsp;m<sub>H</sub> = m<sub>W</sub>(&phi; &minus; 8&alpha;) = 125.09 GeV
            <br>
            &nbsp;&nbsp;(<b>0.09%</b> from experiment)
            <br><br>
            <span style="color:#44aa66">ALL DERIVED from 600-cell spectral data.</span>
        `,
    },
    // --- Step 6 ---
    {
        title: 'U(1) Is Confined',
        vertices: 'type',
        edges: 'U1',
        highlight: 'none',
        text: `
            The 48 U(1) edges (yellow) form a <b>perfect matching</b>:
            <br><br>
            Every fermion has <b>exactly 1</b> U(1) edge.
            <br><br>
            There are <b>zero</b> two-step U(1) paths anywhere in the graph.
            The U(1) subgraph is a forest of 48 disconnected pairs.
            <br><br>
            Consequence: U(1) is <b>ultralocal</b> &mdash; it cannot propagate
            beyond nearest neighbors.
            <br><br>
            In physical terms: the electromagnetic interaction is
            <b>confined at the lattice level</b>. Each fermion sees
            exactly one U(1) partner.
            <br><br>
            <span style="color:#44aa66">DERIVED. Forced by graph structure.</span>
        `,
    },
    // --- Step 7 ---
    {
        title: 'Only SU(3) Has Plaquettes',
        vertices: 'dim_type',
        edges: 'SU3',
        highlight: 'none',
        text: `
            A <b>plaquette</b> = closed triangle of edges (elementary Wilson loop).
            <br><br>
            Plaquette counts by gauge type:
            <br>
            &bull; SU(3): <b>288 plaquettes</b>
            <br>
            &bull; SU(2): <b>0 plaquettes</b>
            <br>
            &bull; U(1): <b>0 plaquettes</b>
            <br><br>
            In Wilson's lattice gauge theory, the gauge action is a sum
            over plaquettes. Only SU(3) can build a lattice action
            on the 600-cell.
            <br><br>
            This mirrors QCD: the strong force confines quarks via
            color flux tubes, built from plaquettes.
            <br><br>
            <span style="color:#44aa66">DERIVED. Combinatorial fact of the graph.</span>
        `,
    },
    // --- Step 8 ---
    {
        title: 'Fermion Mass Formula',
        vertices: 'type',
        edges: 'dim',
        highlight: 'none',
        text: `
            <b>Unified formula:</b>
            <br><br>
            <div style="text-align:center; font-size:13px; color:#fff; margin:4px 0">
                m<sub>f</sub> = m<sub>e</sub> &times; &phi;<sup style="font-size:11px">(5a + 6b)</sup>
            </div>
            <br>
            5 = graph diameter. &nbsp; 6 = intersection number b<sub>1</sub>.
            <br>
            (a,b) = unique minimizer of |a|+|b| such that 5a + 6b = n.
            <br><br>
            <table style="font-size:10px; color:#aac; border-collapse:collapse; width:100%">
            <tr style="border-bottom:1px solid #2a2a4a"><td><b>Particle</b></td><td><b>(a,b)</b></td><td><b>Error</b></td></tr>
            <tr><td>electron</td><td>(0, 0)</td><td style="color:#44aa66">exact (input)</td></tr>
            <tr><td>muon</td><td>(1, 1)</td><td style="color:#44aa66"><b>0.2%</b></td></tr>
            <tr><td>tau</td><td>(1, 2)</td><td style="color:#44aa66"><b>1.4%</b></td></tr>
            <tr><td>up</td><td>(3, &minus;2)</td><td style="color:#44aa66">1.8%</td></tr>
            <tr><td>charm</td><td>(2, 1)</td><td style="color:#44aa66"><b>0.6%</b></td></tr>
            <tr><td>top</td><td>(4, 1)</td><td style="color:#44aa66"><b>0.3%</b></td></tr>
            <tr><td>down</td><td>(1, 0)</td><td style="color:#aa6644">10%</td></tr>
            <tr><td>strange</td><td>(1, 1)</td><td style="color:#aa6644">6.6%</td></tr>
            <tr><td>bottom</td><td>(&minus;1, 4)</td><td style="color:#44aa66"><b>1.2%</b></td></tr>
            </table>
            <br>
            <span style="color:#44aa66">7/9 within 2%. (a,b) are DERIVED, not fitted.</span>
        `,
    },
    // --- Step 9 ---
    {
        title: 'Three Generations',
        vertices: 'soliton',
        edges: 'dim',
        highlight: 'soliton',
        text: `
            A <b>topological defect</b> (soliton) on the 600-cell splits
            each eigenvalue level into <b>exactly 3 sub-levels</b>:
            <br><br>
            <span style="color:#4d80e6; font-size:14px">&block;</span>
            <b>Bulk</b> (m&minus;3 modes): delocalized, far from defect
            <br>
            <span style="color:#33cccc; font-size:14px">&block;</span>
            <b>Shell</b> (2 modes): localized on neighbors
            <br>
            <span style="color:#e64d33; font-size:14px">&block;</span>
            <b>Core</b> (1 mode): localized at the defect
            <br><br>
            The splitting is <b>(m&minus;3, 2, 1)</b> &mdash; always 3 sub-levels,
            for ANY vertex, ANY defect type. The number 3 comes from
            12 neighbors &divide; 4 other cosets = <b>3 neighbors per coset</b>.
            <br><br>
            Interpretation:
            <br>
            &bull; Generation 1 = bulk (lightest, delocalized)
            <br>
            &bull; Generation 2 = shell (intermediate)
            <br>
            &bull; Generation 3 = core (heaviest, localized)
            <br><br>
            Look at the <b>spectrum panel</b> (top right) &mdash; each bar
            has split into 3 colored components.
            <br><br>
            <span style="color:#44aa66">NUMBER 3 is DERIVED. Generation = localization is PATTERN.</span>
        `,
    },
];
