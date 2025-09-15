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
ν1 = one(tpfl)

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
#   INITIAL DATA SETUP FOR MASSIVE PARTICLE
#-----------------------------------------------------------------------

# Scattering parameters
b    = tpfl(1e6)
dist = tpfl(1e12)

# Position
q   = tpfl.([-dist, b, ν0])

# Mass
m   = tpfl(1e-12)

# Speed
v   = 0.999

γ   = ν1/sqrt(ν1-v^2)

# Momentum
p   = tpfl.([γ*v, ν0, ν0])

# Particle object for massive particle
PM   = pomin.setup_single_particle(m, q, p, tpfl)

#-----------------------------------------------------------------------
#   INTEGRATION
#-----------------------------------------------------------------------

tspan = (zero(tpfl), 2*tpfl(1e12))

tols  = tpfl(1e-16)

params       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)                                

sol          = pomin.solveT(PC, params, PM)

ΔZ = sol.u[end]-sol.u[1]

# Extract momentum components for the massive particle
Pi = sol.u[1][10:12]    # Initial momentum
Pf = sol.u[end][10:12]  # Final momentum

# Calculate deflection angle
deflection_angle = abs(asin(Pf[2]/norm(Pf)))

# Theoretical prediction for small angle approximation (ultrarelativistic)
theoretical_deflection = 4*M/b

# Results
println("ultrarelativistic Deflection Calculation Results:")
println("=====================================")
println("Impact parameter b = ", b, " length units")
println("Central mass M = ", M, " mass units")
println("Massive particle mass m = ", m, " mass units")
println("")
println("Numerical deflection angle: ", deflection_angle, " radians")
println("Theoretical deflection:     ", theoretical_deflection, " radians") 
println("Relative error:             ", abs(deflection_angle - theoretical_deflection)/theoretical_deflection * 100, " %")
