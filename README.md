# PoMiN
PoMiN: a Post-Minkowskian N-Body Solver

## Description

N-body code which solves Hamilton's equations for a system of particles using
the N-body Hamiltonian in the following paper:

https://arxiv.org/abs/0807.0214 [Ledvinka, T., Schafer, G., Bicak, J., Phys. Rev. Lett. 100, 251101 (2008)]

The Hamiltonian presented in the paper above models a system of weakly
gravitating point particles. It is fully special relativistic and 
includes first-order effects in (Newton's constant) G from general 
relativity.

This version is a minimal Julia language refactoring of PoMiN for use in
modeling trajectories for laser-propelled spacecraft missions to nearby 
stellar destinations. A more complete version of the Julia language 
refactoring of PoMiN will be released at a later date.

The original C language code is described in the paper:

https://arxiv.org/abs/1805.00813 [Feng, J., Baumann, M., Hall, B., Doss, J., Spencer, L., Matzner, R., Ap. J. 859, 130 (2018)]

## Usage

Example scripts are provided in the src/examples directory. To run them,
first install Julia, then execute a command of the form:

    > julia MXIC_massive.jl

to run the MXIC_massive.jl script.
    
## License

This program is licensed under the MIT License. See LICENSE file.

