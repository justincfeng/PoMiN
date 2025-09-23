#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf
using DoubleFloats
tpfl = Double64

include("target_proxima.jl")

println("=== INVESTIGATING TARGET DEFINITIONS ===")
println()

# From target_proxima.jl, we have:
# pfs = parfuncs(ics)
# XstF,XpxF,XbF,tcl = pfs

XstF,XpxF,XbF,tcl_from_target = pfs

println("1. TARGET PROBLEM FUNCTIONS:")
println("XstF: spacecraft position function")
println("XpxF: Proxima position function") 
println("XbF: target position function")
println("tcl: closest approach time = ", tcl_from_target)
println()

println("2. POSITIONS AT CLOSEST APPROACH TIME:")
flat_space_spacecraft_final = XstF(tcl)
proxima_final_from_func = XpxF(tcl)
target_final_from_func = XbF(tcl)

println("flat_space_spacecraft_final = XstF(tcl): ", flat_space_spacecraft_final)
println("Proxima final = XpxF(tcl): ", proxima_final_from_func)
println("Target final = XbF(tcl): ", target_final_from_func)
println()

println("3. COMPARISON WITH bv CONSTRUCTION:")
println("bv (target offset): ", bv)
println("Xpx0 (initial Proxima): ", Xpx0)
println("Proxima + bv: ", Xpx0 + bv)
println()

println("4. CHECKING IF XbF(tcl) = XpxF(tcl) + bv:")
calculated_target = proxima_final_from_func + bv
println("XpxF(tcl) + bv: ", calculated_target)
println("XbF(tcl): ", target_final_from_func)
println("Difference: ", norm(calculated_target - target_final_from_func))
println()

println("5. CHECKING WHAT flat_space_spacecraft_final REPRESENTS:")
println("Is flat_space_spacecraft_final = XbF(tcl)?")
println("XstF(tcl): ", flat_space_spacecraft_final)
println("XbF(tcl): ", target_final_from_func)
println("Difference: ", norm(flat_space_spacecraft_final - target_final_from_func))
println()

println("6. UNDERSTANDING THE TARGETING PROBLEM:")
# From target_proxima.jl, the targeting problem sets up initial conditions
# such that the spacecraft reaches the target position XbF(tcl) = Xpx0 + bv
# at time tcl in flat space (no gravity)

println("The targeting problem is designed so that:")
println("- Spacecraft starts at Xst0 with velocity Vst")
println("- In flat space, spacecraft reaches target Xpx0 + bv at time tcl")
println("- Therefore XstF(tcl) should equal the target position")
println()

println("7. VERIFICATION:")
println("Target from targeting: ", target_final_from_func)
println("Spacecraft final (flat): ", flat_space_spacecraft_final)
println("Are they equal? ", norm(target_final_from_func - flat_space_spacecraft_final) < 1e-10)
