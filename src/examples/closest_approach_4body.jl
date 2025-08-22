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
using LinearAlgebra, Printf

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

println("=== 4-Body Closest Approach Calculation ===")
println("Setting up particle system...")

# Spacecraft (chip) parameters - increased from 1e-33 to 1e-8 for better numerical scaling
# while still remaining a test particle (negligible compared to other masses)
mchip = tpfl(1.0E-12)
qchip = tpfl.([-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07])
# Momentum scaled by same factor as mass to preserve velocity
pchip = mchip .* tpfl.([-7.443174824381408e-2, -5.690948861061366e-2, -1.8135019058747826e-1])

Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)
println("✓ Spacecraft particle created with mass: ", mchip)

# Proxima Centauri parameters
mProx = tpfl(0.1221*1e-14)
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

# Merge all particles into 4-body system
P4body = pomin.merge_particle_systems(Pchip, PProx, Psol, Pjup)
println("✓ 4-body system merged successfully (Spacecraft + Proxima + Sun + Jupiter)")

P2body = pomin.merge_particle_systems(Pchip, PProx)
println("✓ 2-body system merged successfully (Spacecraft + Proxima)")

# Integration parameters
tspan = (tpfl(0), tpfl(1.410E+14)) # 22 years
δ = tpfl(7.31E+8) # about one hour

params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tpfl(1e-22), rtol=tpfl(1e-22), Nrec=100)
println("✓ Integration parameters set: ", tspan[2]/1e14, " × 10^14 time units")

println("\nStarting integration...")
@time sol = pomin.solve(P4body, params)
# @time sol = pomin.solve(P2body, params)
println("✓ Integration completed with ", length(sol.t), " time steps")

# Find closest approach
println("\nAnalyzing closest approach...")
mindist = tpfl(1E100)
minvec = zeros(tpfl, 3)
min_time = tpfl(0.0)
numsteps = length(sol.t)

Npar = 4
# Npar = 2

for tn in 1:numsteps  # tn is timestep number
    # get q vector for particle 1 (spacecraft) at timestep tn
    q1 = Z2q(Npar, 3, 1, sol.u[tn])  # Npar particles, 3 dimensions, particle 1
    # get q vector for particle 2 (Proxima Centauri) at timestep tn
    q2 = Z2q(Npar, 3, 2, sol.u[tn])  # Npar particles, 3 dimensions, particle 2
    # get distance and save it if it's the minimum distance
    dist = norm(q1-q2)
    if dist < mindist
        global mindist = dist
        global minvec = q1-q2
        global min_time = sol.t[tn]
    end
end

# Print results
println("\n=== RESULTS ===")
println("Minimum distance: ", mindist)
println("Separation vector at closest approach: ", minvec)
println("Time of closest approach: ", min_time)
println("Closest approach in scientific notation: ", @sprintf("%.3e", mindist))
println("=== END RESULTS ===")

