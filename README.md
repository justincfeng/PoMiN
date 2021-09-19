# PoMiN
PoMiN: a Post-Minkowskian N-Body Solver

## Description

N-body code which solves Hamilton's equations for a system of particles using
the N-body Hamiltonian in the following paper:

https://arxiv.org/abs/0807.0214 [Ledvinka, T., Schafer, G., Bicak, J., Phys. Rev. Lett. 100, 251101 (2008)]

The Hamiltonian presented in the paper above models a system of weakly
gravitating point particles. It is fully special relativistic, but includes
only first-order effects in (Newton's constant) G from General Relativity.

The code uses a fourth-order Runge-Kutta integrator, and has a computational
complexity that scales as N^2. The code implements a basic, global adaptive
timestepping scheme, based on the distances and velocities of the two closest
particles. Support for parallel computing has not yet been implemented.

For more, see our PoMiN paper:

https://arxiv.org/abs/1805.00813 [Feng, J., Baumann, M., Hall, B., Doss, J., Spencer, L., Matzner, R., Ap. J. 859, 130 (2018)]


## Input

Input must be read from a CSV file with one line in the following format:

    description,start-time,end-time,timestep,iterations,courant-number,gravitational-constant,number-of-particles,mass_1,qx_1,qy_1,qz_1,px_1,py_1,pz_1,mass_2,qx_2,...


    INPUT,                  TYPE,    DESCRIPTION  
    description,            string,  ignored by program--except it cannot be empty  
    start-time,             float,   time to start at (physics time)  
    end-time,               float,   time to end at (physics time)  
    timestep,               float,   (maximum) length of each computational iteration  
    iterations,             integer, upper bound of iterations (0 = unbounded)  
    courant-number,         float,   adaptive-timestep parameter (0 = off)  
    gravitational-constant, float,   you know: G  
    number-of-particles     integer, number of interacting particles
    mass_1,                 float,   mass of particle 1  
    qx_1,                   float,   x-position of particle 1  
    ...,                    ...,     ...  
    py_3,                   float,   y-momentum of particle 3  
    ...,                    ...,     ...  

NOTE: the speed of light c is assumed to have a value of 1.

NOTE: the 'number-of-particles' is the number of fully interacting particles in the 
system. Any additional particles in the input line after the first 'number-of-particles' 
particles will be treated as test particles. Set 'number-of-particles' to zero if all
particles in the input line are fully interacting particles (no test particles).

NOTE: the 'courant-number' is the ratio of the change in distance between the
two closest particles (in a single timestep) to the distance between those
particles. A positive 'courant-number' turns on adaptive timestepping and computes
the timestep using the distances between all particles, including test particles. 
Setting 'courant-number' to zero turns off adaptive timestepping. A negative 
'courant-number' turns on adaptive timestepping, but computes the timestep (with 
the absolute value of 'courant-number') using only the distances between fully 
interacting particles.


## Output

Output is written to stdout (except for errors, which are written to stderr) in
this form:

    iteration_#,time,qx_1,qy_1,qz_1,px_1,py_1,pz_1,qx_2,qy_2,...

## License

This program is licensed under the MIT License. See LICENSE file.

