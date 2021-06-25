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

For more, see link:doc/pomin.pdf[pomin.pdf].

pomin takes input from the standard input (stdin) and outputs to the
standard output (stdout); both of which are in the comma-separated values (CSV)
format (see the *INPUT* and *OUTPUT* sections).

For more information about using stdin and stdout, see the *NOTES* section below.
Also see the *EXAMPLES* section for some useful uses.

## Options

Print debug (e.g., auxiliary values and metadata) information to stdout.  

    -d  

Print verbose (e.g., wall-time) information to stderr.  

    -v  

Which parts of the Hamiltonian to use:  

    -p <parts>  
    
    where <parts> is one of the following (see paper for definitions of H1, H2, H3):  

    1    to use H1 only  
    12   to use H1+H2  
    13   to use H1+H3 or  
    123  to use H1+H2+H3 (full Hamiltonian, default)  


## Input

Input must be read from stdin as one line in the following format:

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

WARNING: the 'description' field cannot be empty due to the strtok(1) function
parsing; however, using a unique value here is helpful when selecting input to
use with, e.g., grep(1).

To generate a header line for n particles, run the “header” make(1) target:

    make header N=n

The pomin code will automatically re-run itself if it detects another line of input
(i.e., recursion).


## Output

Output is written to stdout (except for errors, which are written to stderr) in
this form:

    iteration_#,time,qx_1,qy_1,qz_1,px_1,py_1,pz_1,qx_2,qy_2,...

## Validation

## Debug

## Notes

### Obtaining

The Git repository can be obtained through the HTTPS protocol with

    git clone https://github.com/justincfeng/PoMiN.git


### Compiling

To compile run

    make [double | quad]

where double (the default) or quad is the desired precision. Then, to install
the program on an operating system that abides by the Filesystem Heirarchy
Standard (FHS; e.g., a GNU/Linux distribution) run make install and the
executible will be put into /usr/local/bin/ and this man page will be put into
/usr/local/man/man1/ so that you can run the program from any directory and
view its manual at any time with man pomin


### Scripts

## Examples

## License

This program is licensed under the MIT License. See LICENSE file.

