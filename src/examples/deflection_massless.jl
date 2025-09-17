#-----------------------------------------------------------------------
#
#   Light deflection calculation
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, DoubleFloats, Printf
using Plots

# Import Z2q function from HamTools
using .pomin: Z2q
tpfl  =  Double64

ν0 = zero(tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR CENTRAL MASS
#-----------------------------------------------------------------------

# Mass, position, momentum
M    = one(tpfl)
qC   = zeros(tpfl,3)
pC   = zeros(tpfl,3)

# Particle object for central mass
PC   = pomin.setup_single_particle(M, qC, pC, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR LIGHT
#-----------------------------------------------------------------------

# Scattering parameters
b    = tpfl(1e6)
dist = tpfl(1e12)

# Position
qL   = tpfl.([-dist, b, ν0])

# Momentum
pL   = tpfl.([tpfl(1e-12), ν0, ν0])

# Particle object for light (as test particle)
PL   = pomin.setup_single_particle(zero(tpfl), qL, pL, tpfl)

#-----------------------------------------------------------------------
#   INTEGRATION
#-----------------------------------------------------------------------

tspan = (zero(tpfl), 2*tpfl(1e12))

tols  = tpfl(1e-16)

params       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)                                

sol          = pomin.solve(PC, params; testparticles=PL)

ΔZ = sol.u[end]-sol.u[1]

# Extract momentum components for the light particle (test particle)
# Test particle momentum is at indices 10:12 (after central mass q,p and test particle q)
Pi = sol.u[1][10:12]    # Initial momentum
Pf = sol.u[end][10:12]  # Final momentum

# Calculate deflection angle
deflection_angle = abs(asin(Pf[2]/norm(Pf)))

# Theoretical prediction for small angle approximation
theoretical_deflection = 4*M/b

# Results
println("Light Deflection Calculation Results:")
println("=====================================")
println("Impact parameter b = ", b, " length units")
println("Central mass M = ", M, " mass units")
println("")
println("Numerical deflection angle: ", deflection_angle, " radians")
println("Theoretical deflection:     ", theoretical_deflection, " radians") 
println("Relative error:             ", abs(deflection_angle - theoretical_deflection)/theoretical_deflection * 100, " %")