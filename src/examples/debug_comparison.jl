#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf
using DoubleFloats
tpfl = Double64

include("target_proxima.jl")

println("=== DEBUGGING ALPHA CENTAURI MISS DISTANCE DIFFERENCE ===")
println()

# Extract the exact same initial conditions used in both scripts
println("1. ALPHA CENTAURI SETUP:")
malpha = tpfl(2.0429)
qalpha = tpfl.([-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13])
valpha = tpfl.([-3.11496332572849e-05, 7.38228261899770e-05, 7.22023391262230e-05])
γalpha = one(tpfl) / sqrt(one(tpfl) - norm(valpha)^2)
palpha = tpfl.(malpha * γalpha .* valpha)
Palpha = pomin.setup_single_particle(malpha, qalpha, palpha, tpfl)

println("Alpha Centauri mass: ", malpha)
println("Alpha Centauri position: ", qalpha)
println("Alpha Centauri velocity: ", valpha)
println()

println("2. SPACECRAFT SETUP (Original velocity from target_proxima.jl):")
mchip = tpfl(1.0E-30)
qchip = Xst0
println("Spacecraft mass: ", mchip)
println("Spacecraft position: ", qchip)
println("Original velocity Vst: ", Vst)

# Setup spacecraft with original velocity (like AlphaAB.jl)
vst = norm(Vst)
γst = one(tpfl) / sqrt(one(tpfl) - (vst / c)^2)
pchip_orig = (mchip * γst) .* Vst
Pchip_orig = pomin.setup_single_particle(mchip, qchip, pchip_orig, tpfl)

# Setup spacecraft with pchip_rel0 (like All_FTID.jl)
γst0 = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel0 = (mchip*γst0) .* Vst
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)

println("Original momentum (AlphaAB style): ", pchip_orig)
println("pchip_rel0 momentum (FTID style): ", pchip_rel0)
println("Momentum difference: ", norm(pchip_orig - pchip_rel0))
println()

println("3. PROXIMA SETUP:")
mProx = tpfl(0.1221)
qProx = Xpx0
pProx = mProx .* γVpx
PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)

println("Proxima mass: ", mProx)
println("Proxima position: ", qProx)
println("Proxima momentum: ", pProx)
println()

println("4. INTEGRATION PARAMETERS:")
tspan = (tpfl(0), tcl * tpfl(1.1))
tols = tpfl(1e-16)
params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tols, rtol=tols)

println("Integration time span: ", tspan)
println("tcl (closest approach time): ", tcl)
# time_closest_approach from AlphaAB.jl
time_closest_approach = pfs[4]  # Same as tcl
println("time_closest_approach: ", time_closest_approach)
println("Time difference: ", abs(tcl - time_closest_approach))
println()

println("5. RUNNING SIMULATIONS:")

# AlphaAB.jl style simulation
println("AlphaAB.jl style:")
main_particle_system = pomin.merge_particle_systems(Palpha)
test_particle_system = pomin.merge_particle_systems(Pchip_orig, PProx)
sol_AB = pomin.solve(main_particle_system, params; testparticles=test_particle_system)

Zend_AB = sol_AB(time_closest_approach)
q_starchip_AB = Zend_AB[7:9]
q_proxima_AB = Zend_AB[10:12]
q_target_AB = q_proxima_AB + bv
miss_dist_AB = norm(q_starchip_AB - q_target_AB) / AU

println("Spacecraft final position (AlphaAB): ", q_starchip_AB)
println("Proxima final position (AlphaAB): ", q_proxima_AB)
println("Target position (AlphaAB): ", q_target_AB)
println("Miss distance (AlphaAB): ", miss_dist_AB, " AU")
println()

# All_FTID.jl style simulation
println("All_FTID.jl style:")
sol_FTID = pomin.solve(Palpha, params; testparticles=Pchip_rel0)
zend_FTID = sol_FTID(tcl)
q_starchip_FTID = zend_FTID[7:9]

# Get flat space spacecraft final position
XstF,XpxF,XbF,tcl_from_target = pfs
flat_space_spacecraft_final = XstF(tcl)

qmissFlat = flat_space_spacecraft_final - q_starchip_FTID
dmissFlat = norm(qmissFlat) / AU

println("Spacecraft final position (FTID): ", q_starchip_FTID)
println("Flat space target position: ", flat_space_spacecraft_final)
println("Miss distance (FTID): ", dmissFlat, " AU")
println()

println("6. COMPARISON:")
println("AlphaAB miss distance: ", miss_dist_AB, " AU")
println("FTID miss distance: ", dmissFlat, " AU")
println("Ratio (AB/FTID): ", miss_dist_AB / dmissFlat)
println("Difference: ", abs(miss_dist_AB - dmissFlat), " AU")
println()

println("7. TARGET COMPARISON:")
println("AlphaAB target: ", q_target_AB)
println("FTID target: ", flat_space_spacecraft_final)
println("Target difference: ", norm(q_target_AB - flat_space_spacecraft_final) / AU, " AU")
println()

println("8. SPACECRAFT POSITION COMPARISON:")
println("AlphaAB spacecraft: ", q_starchip_AB)
println("FTID spacecraft: ", q_starchip_FTID)
println("Spacecraft difference: ", norm(q_starchip_AB - q_starchip_FTID) / AU, " AU")
