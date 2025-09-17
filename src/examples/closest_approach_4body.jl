#-----------------------------------------------------------------------
#
#   4-BODY CLOSEST APPROACH CALCULATION
#   Spacecraft, Proxima Centauri, Sun, and Jupiter
#   Modernized version compatible with current PoMiN framework
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots

# Precision options tested:
# Float64    - Standard precision, works well
# Double64   - Extended precision (DoubleFloats.jl), works well  
# ArbFloat   - Arbitrary precision (ArbNumerics.jl), cache issues with ODE integrators
# BigFloat   - Julia built-in arbitrary precision, should work better

using DoubleFloats
tpfl  = Double64

# using ArbNumerics
# setprecision(ArbFloat, 256)
# tpfl  = ArbFloat

#setprecision(BigFloat, 256)  # Set precision to 256 bits
#tpfl  = BigFloat

# Import Z2q function from HamTools
using .pomin: Z2q
using .pomin: Z2Part

println("=== 4-Body Closest Approach Calculation ===")
println("Setting up particle system...")

AU = tpfl(99731913.0);

# Spacecraft (chip) parameters - increased from 1e-33 to 1e-8 for better numerical scaling
# while still remaining a test particle (negligible compared to other masses)
mchip = tpfl(1.0E-30)
qchip = tpfl.([-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07])
# Momentum scaled by same factor as mass to preserve velocity
pchip = mchip .* tpfl.([-7.443174824381408e-2, -5.690948861061366e-2, -1.8135019058747826e-1])

Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)
println("✓ Spacecraft particle created with mass: ", mchip)

# Proxima Centauri parameters
mProx = tpfl(0.1221)
qProx = tpfl.([-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]) 
pProx = mProx .* tpfl.([-3.34779933990201E-5, 7.24198145771899E-5, 6.82553747232694E-5])

PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)
println("✓ Proxima Centauri particle created with mass: ", mProx)

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)
println("✓ Sun particle created with mass: ", msol)

# Jupiter parameters (realistic orbital position)
# Jupiter mass relative to Sun: ~0.000954
# Average distance from Sun: ~7.78 × 10^8 km = 7.78 × 10^11 m
# Orbital velocity: ~13.1 km/s = 1.31 × 10^4 m/s
mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([7.78E+11, 0.0, 0.0])  # Jupiter position (x-axis for simplicity)
pjup = tpfl.([0.0, 1.31E+4 * mjup, 0.0])  # Jupiter momentum (y-direction for circular orbit)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)
println("✓ Jupiter particle created with mass: ", mjup)

# Merge heavy particles
P3body = pomin.merge_particle_systems(PProx, Psol, Pjup)
println("✓ 3-body system merged successfully (Proxima + Sun + Jupiter)")

Ptest  = Pchip

P2body = pomin.merge_particle_systems(Pchip, PProx)
println("✓ 2-body system merged successfully (Spacecraft + Proxima)")

# Integration parameters
tspan = (tpfl(0), tpfl(1.410E+14)) # 22 years
δ = tpfl(7.31E+8) # about one hour

params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tpfl(1e-16), rtol=tpfl(1e-16), Nrec=100)
println("✓ Integration parameters set: ", tspan[2]/1e14, " × 10^14 time units")

println("\nStarting integration...")
@time sol = pomin.solve(P3body, params; testparticles=Ptest)

#@time sol = pomin.solve(Psol, params; testparticles=Ptest)

# @time sol = pomin.solve(P2body, params)
println("✓ Integration completed with ", length(sol.t), " time steps")

# Find closest approach
println("\nAnalyzing closest approach...")
mindist = tpfl(1E100)
minvec = zeros(tpfl, 3)
min_time = tpfl(0.0)
numsteps = length(sol.t)


qtestFunc = t->sol(t)[19:21]
qproxFunc = t->sol(t)[1:3]

dq = t->norm(qtestFunc(t)-qproxFunc(t))/(AU)

tend = sol.t[end]

dq(0.9635391124 * tend)

# dq = 1.28113381156829800727317956434636241e-01 AU @ mchip=1e-8, tol=1e-26
# dq = 1.28113381156829809438081060116200672e-01 AU @ mchip=1e-6, tol=1e-20
# dq = 1.28113381157121865267484762416170654e-01 AU @ mchip=1e-4, tol=1e-16

# ----------------------------------------------------------------------

plot(dq, zero(tpfl) , tend , label="Separation")

plot(dq, 0.9635391122 * tend , 0.9635391126 * tend , label="Separation")
