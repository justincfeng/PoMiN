#-----------------------------------------------------------------------
#
#   PoMiN: A Hamiltonian Post-Minkowski N-Body Code
#
#       Version 2.0
#
#-----------------------------------------------------------------------

module pomin

using LinearAlgebra
using CSV
using ForwardDiff
using OrdinaryDiffEq
using Logging
using DoubleFloats

# Core functionality
include("core/types.jl")                 # Type definitions
include("core/hamiltonian/HamPM.jl")    # Post-Minkowski Hamiltonian
include("core/hamiltonian/HamPM_massless.jl") # Massless particle support

# Physics modules
include("physics/gravitational_waves/gwsc.jl")  # GW strain calculator
include("physics/post_minkowski/HO.jl")        # Harmonic oscillator

# Numerical methods
include("integrators/epsi.jl")          # Symplectic Integrator (Tao 2016)
include("integrators/rk4i.jl")          # 4th order Runge-Kutta
include("integrators/intjul.jl")        # Julia integrator functions

# Utilities
include("utils/idxer.jl")               # Index management
include("utils/idgen.jl")               # ID generation
include("utils/broyden.jl")             # Broyden method solver
include("utils/misc.jl")                # Miscellaneous utilities

# Input/Output
include("io/io.jl")                     # I/O routines

# Export public interface

# Types
export RealVec, Particles, soln

# Core functionality
export H, dH, Jsympl

# Integration methods
export jlintegrator, jlintegratorfull

end
