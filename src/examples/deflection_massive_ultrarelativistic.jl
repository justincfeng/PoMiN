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
p   = m .* tpfl.([γ*v, ν0, ν0])
pN  = m .* tpfl.([v, ν0, ν0])

# Particle object for massive particle
PM   = pomin.setup_single_particle(m, q, p, tpfl)
PN   = pomin.setup_single_particle(m, q, pN, tpfl)

PCombined = pomin.merge_particle_systems(PC, PN)

#-----------------------------------------------------------------------
#   INTEGRATION
#-----------------------------------------------------------------------

tspan = (zero(tpfl), 2*tpfl(1e12))

tols  = tpfl(1e-16)

params       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)                                

sol          = pomin.solve(PC, params; testparticles=PM)
solN         = pomin.solve(PCombined, params; Newtonian=true)

# Extract momentum components for the massive particle (Test Particle approach)
Pi = sol.u[1][10:12]    # Initial momentum
Pf = sol.u[end][10:12]  # Final momentum

# Calculate deflection angle for test particle
deflection_angle_T = abs(asin(Pf[2]/norm(Pf)))

# Extract momentum components for the massive particle (Newtonian approach)
# For Newtonian solution: [qC(1:3), qM(4:6), pC(7:9), pM(10:12)]
Pi_N = solN.u[1][10:12]    # Initial momentum (second particle)
Pf_N = solN.u[end][10:12]  # Final momentum (second particle)

# Calculate deflection angle for Newtonian using momentum directly
deflection_angle_N = abs(asin(Pf_N[2]/norm(Pf_N)))

# Theoretical predictions
theoretical_deflection_relativistic = 4*M/b  # Ultrarelativistic limit
theoretical_deflection_newtonian = 2*M/b     # Classical Newtonian scattering

# Results
println("Massive Particle Deflection Calculation Results:")
println("===============================================")
println("Impact parameter b = ", b, " length units")
println("Central mass M = ", M, " mass units")
println("Massive particle mass m = ", m, " mass units")
println("Particle speed v = ", v, " (relativistic)")
println("")
println("Test Particle (Relativistic) Results:")
println("  Numerical deflection angle: ", deflection_angle_T, " radians")
println("  Theoretical deflection:     ", theoretical_deflection_relativistic, " radians") 
println("  Relative error:             ", abs(deflection_angle_T - theoretical_deflection_relativistic)/theoretical_deflection_relativistic * 100, " %")
println("")
println("Newtonian Results:")
println("  Numerical deflection angle: ", deflection_angle_N, " radians")
println("  Theoretical deflection:     ", theoretical_deflection_newtonian, " radians") 
println("  Relative error:             ", abs(deflection_angle_N - theoretical_deflection_newtonian)/theoretical_deflection_newtonian * 100, " %")
println("")
println("Comparison:")
println("  Relativistic/Newtonian ratio: ", deflection_angle_T/deflection_angle_N)
println("  Theoretical ratio (4M/b)/(2M/b): ", theoretical_deflection_relativistic/theoretical_deflection_newtonian)
