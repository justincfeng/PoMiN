#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf
using DoubleFloats
tpfl = Double64

include("target_proxima.jl")

println("=== PROXIMA TREATMENT COMPARISON ===")
println()

println("1. ALL_FTID.jl ALPHA CENTAURI SECTION:")
println("   sol0 = pomin.solve( Palpha, params; testparticles=Pchip_rel0 )")
println("   - Main gravitational source: Alpha Centauri (Palpha)")
println("   - Test particles: Spacecraft only (Pchip_rel0)")
println("   - Proxima: NOT included in simulation")
println()

println("2. AlphaAB.jl:")
println("   main_particle_system = pomin.merge_particle_systems(Palpha)")
println("   test_particle_system = pomin.merge_particle_systems(Pchip, PProx)")
println("   sol = pomin.solve(main_particle_system, params; testparticles=test_particle_system)")
println("   - Main gravitational source: Alpha Centauri (Palpha)")
println("   - Test particles: Spacecraft (Pchip) AND Proxima (PProx)")
println("   - Proxima: IS included as test particle")
println()

println("3. KEY DIFFERENCE:")
println("   - FTID: Proxima is NOT simulated - uses targeting problem position")
println("   - AlphaAB: Proxima IS simulated under Alpha Centauri gravity")
println()

# Setup the same initial conditions
malpha = tpfl(2.0429)
qalpha = tpfl.([-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13])
valpha = tpfl.([-3.11496332572849e-05, 7.38228261899770e-05, 7.22023391262230e-05])
γalpha = one(tpfl) / sqrt(one(tpfl) - norm(valpha)^2)
palpha = tpfl.(malpha * γalpha .* valpha)
Palpha = pomin.setup_single_particle(malpha, qalpha, palpha, tpfl)

mchip = tpfl(1.0E-30)
qchip = Xst0
vst = norm(Vst)
γst = one(tpfl) / sqrt(one(tpfl) - (vst / c)^2)
pchip_orig = (mchip * γst) .* Vst
Pchip_orig = pomin.setup_single_particle(mchip, qchip, pchip_orig, tpfl)

mProx = tpfl(0.1221)
qProx = Xpx0
pProx = mProx .* γVpx
PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)

tspan = (tpfl(0), tcl * tpfl(1.1))
tols = tpfl(1e-16)
params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tols, rtol=tols)

time_closest_approach = pfs[4]

println("4. SIMULATION COMPARISON:")

# FTID style - no Proxima simulation
println("FTID style (no Proxima simulation):")
sol_FTID = pomin.solve(Palpha, params; testparticles=Pchip_orig)
zend_FTID = sol_FTID(tcl)
q_starchip_FTID = zend_FTID[7:9]

# Target from targeting problem (where Proxima "should" be)
XstF,XpxF,XbF,tcl_from_target = pfs
proxima_targeting = XpxF(tcl)
target_FTID = proxima_targeting + bv

println("  Spacecraft final: ", q_starchip_FTID)
println("  Proxima (targeting): ", proxima_targeting)
println("  Target (targeting + bv): ", target_FTID)

miss_FTID = norm(q_starchip_FTID - target_FTID) / AU
println("  Miss distance: ", miss_FTID, " AU")
println()

# AlphaAB style - with Proxima simulation
println("AlphaAB style (with Proxima simulation):")
main_particle_system = pomin.merge_particle_systems(Palpha)
test_particle_system = pomin.merge_particle_systems(Pchip_orig, PProx)
sol_AB = pomin.solve(main_particle_system, params; testparticles=test_particle_system)

Zend_AB = sol_AB(time_closest_approach)
q_starchip_AB = Zend_AB[7:9]
q_proxima_AB = Zend_AB[10:12]
target_AB = q_proxima_AB + bv

println("  Spacecraft final: ", q_starchip_AB)
println("  Proxima (simulated): ", q_proxima_AB)
println("  Target (simulated + bv): ", target_AB)

miss_AB = norm(q_starchip_AB - target_AB) / AU
println("  Miss distance: ", miss_AB, " AU")
println()

println("5. PROXIMA POSITION COMPARISON:")
println("  Proxima (targeting): ", proxima_targeting)
println("  Proxima (simulated): ", q_proxima_AB)
proxima_diff = norm(q_proxima_AB - proxima_targeting) / AU
println("  Proxima difference: ", proxima_diff, " AU")
println()

println("6. TARGET POSITION COMPARISON:")
println("  Target (FTID): ", target_FTID)
println("  Target (AlphaAB): ", target_AB)
target_diff = norm(target_AB - target_FTID) / AU
println("  Target difference: ", target_diff, " AU")
println()

println("7. CONCLUSION:")
println("  Miss distance difference: ", abs(miss_AB - miss_FTID), " AU")
println("  This difference comes from Proxima's gravitational evolution")
println("  under Alpha Centauri's influence vs. the original targeting assumption")
