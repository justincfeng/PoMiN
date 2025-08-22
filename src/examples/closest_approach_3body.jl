#-----------------------------------------------------------------------
#
#   3-BODY CLOSEST APPROACH CALCULATION
#   Modernized version compatible with current PoMiN framework
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, DoubleFloats, Printf

# Import Z2q function from HamTools
using .pomin: Z2q

println("=== 3-Body Closest Approach Calculation ===")
println("Setting up particle system...")

tpfl  =  Float64  # Use Float64 for better compatibility

# Spacecraft (chip) parameters
mchip = tpfl(1.0057832537E-33)
qchip = tpfl.([-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07])
pchip = tpfl.([-7.486220592724260E-35, -5.723861062118610E-35, -1.823989847481890E-34])

Pchip = pomin.setup_single_particle(mchip, qchip, pchip)
println("✓ Spacecraft particle created with mass: ", mchip)

# Proxima Centauri parameters
mProx = tpfl(0.1221)
qProx = tpfl.([-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13])
pProx = tpfl.([-3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45])

PProx = pomin.setup_single_particle(mProx, qProx, pProx)
println("✓ Proxima Centauri particle created with mass: ", mProx)

# Solar system (third body) parameters
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol)
println("✓ Solar system particle created with mass: ", msol)

P3body = pomin.merge_particle_systems(Pchip, PProx, Psol)
println("✓ 3-body system merged successfully")

# Integration parameters
tspan = (tpfl(0), tpfl(1.410E+14)) # 22 years
δ = tpfl(7.31E+8) # about one hour

params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tpfl(1e-12), rtol=tpfl(1e-12), Nrec=100)
println("✓ Integration parameters set: ", tspan[2]/1e14, " × 10^14 time units")

println("\nStarting integration...")
@time sol = pomin.solve(P3body, params)
println("✓ Integration completed with ", length(sol.t), " time steps")

# Find closest approach
println("\nAnalyzing closest approach...")
mindist = tpfl(1E100)
minvec = zeros(tpfl, 3)
min_time = tpfl(0.0)
numsteps = length(sol.t)

for tn in 1:numsteps  # tn is timestep number
    # get q vector for particle 1 (spacecraft) at timestep tn
    q1 = Z2q(3, 3, 1, sol.u[tn])  # 3 particles, 3 dimensions, particle 1
    # get q vector for particle 2 (Proxima Centauri) at timestep tn
    q2 = Z2q(3, 3, 2, sol.u[tn])  # 3 particles, 3 dimensions, particle 2
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

