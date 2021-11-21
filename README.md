# ZeemanSpectra.jl

ZeemanSpectra.jl is a julia package for computing hyperfine structure of alkali atoms (and iodine atom) in weak ambient magnetic field.

## Roadmap
1. Atom data
    - alkali atoms
        - [x] 1. $^6\mathrm{Li}$
        - [ ] 2. $^{39}\mathrm{K}$
        - [ ] ...
    - iodine atoms
        1. [x] 1. $^{127}\mathrm{I}$
        2. [ ] 2. $^{129}\mathrm{I}$
2. Angular momentum operators
    - [x] ladder operator: $J_+$ and $J_-$ ($I$ as well)
    - [x] $J^2$, $J_z$ ($I$ as well)
    - [x] any combination of above oprators, such as $J_zJ_+I_zI_-$
3. Perturbed hamiltonian matrix of zeeman effect and the hyperfine interaction
    - [x] $\hat{H'} = \hat{H}_{Zeeman} + \hat{H}_{hfs}$
4. Degenerate perturbation theory (hamiltonian diagonalization)
    - [x] level splits (eigenvalues)
    - [x] quantum states (eigenstates)
5. Irreducible tensor operators
    - [x] Wigner-Eckart theorem
    - [x] ...
6. Transition line strength
    - [ ] E1 Electric diople radiation
    - [ ] M1 Magnetic diople radiation
7. Spectral line broadenings
    - [x] doppler
    - [x] pressure
    - [x] voigt
