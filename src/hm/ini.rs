use prng_mt::MT19937;
use crate::hm::HeisenbergModel;

impl HeisenbergModel {
    pub fn new(para_lx: usize, para_ly: usize, para_beta: f64, para_seed: u32) -> Self {
        Self {
            // ----------------------------------------------------------------
            //  Basic params
            // ----------------------------------------------------------------
            lx: para_lx,
            ly: para_ly,
            beta: para_beta,
            n: 0,
            m: 20,

            // --------------------------------------------------------
            //  Lattice
            // --------------------------------------------------------
            n_sites: para_lx * para_ly,
            n_bonds: 0,                 // default
            b_sites: Vec::new(),         // default

            // --------------------------------------------------------
            //  Two frequently-used factors
            // --------------------------------------------------------
            prob_add_factor: 0.0,       // default
            prob_remove_factor: 0.0,    // default

            // ----------------------------------------
            //  Random number generator
            // ----------------------------------------
            rng: MT19937::new(para_seed),

            // ----------------------------------------
            //  Data structures for configuration
            // ----------------------------------------
            spins: Vec::new(),
            op_string: Vec::new(),
            v_first: Vec::new(),
            v_last: Vec::new(),
            vertex_list: Vec::new(),

            // ----------------------------------------
            //  For measurements
            // ----------------------------------------
            energy: 0.0,
            energy2: 0.0,
            heat_capacity: 0.0,
            zz_correlation: vec![0.0; para_lx * para_ly],
        }
    }

    pub fn initialization(&mut self){
        self.n_bonds = if self.ly == 1 { self.lx } else { 2 * self.lx * self.ly };

        self.b_sites = {
            let mut lattice: Vec<Vec<usize>> = Vec::new();

            //  make the 1D lattice (chain)
            if self.ly == 1 {
                for i in 0..self.n_sites {
                    lattice.push(vec![i, (i + 1) % self.n_sites]);
                }
            }

            //  make the 2D square lattice
            else {
                for j in 0..self.ly {
                    for i in 0..self.lx {
                        lattice.push(vec![
                            self.lx * j + i,
                            self.lx * j + (i + 1) % self.lx,
                        ]);
                    }
                }

                for j in 0..self.ly {
                    for i in 0..self.lx {
                        lattice.push(vec![
                            self.lx * j + i,
                            (self.lx * j + self.lx + i) % (self.lx * self.ly),
                        ])
                    }
                }
            }
            lattice.clone()
        };

        self.prob_add_factor = 0.5 * self.beta * self.n_bonds as f64;
        self.prob_remove_factor = 1.0 / self.prob_add_factor;

        // Initialize other data structures
        self.spins = (0..self.n_sites)
            .map(|_| if self.rand_prob() > 0.5 { 1 } else { -1 })
            .collect();

        self.op_string = vec![-1; self.m];
        self.v_first = vec![-1; self.n_sites];
        self.v_last = vec![-1; self.n_sites];
        self.vertex_list = vec![-1; 4 * self.m];
    }
}