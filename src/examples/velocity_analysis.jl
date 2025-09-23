#!/usr/bin/env julia

using DoubleFloats, LinearAlgebra
tpfl = Double64

# Velocity from All.jl (fine-tuned for all bodies)
Vst_All = [tpfl("-7.29280133630346695175238301027084351e-02"), 
           tpfl("-5.57596267171519947482140178590898884e-02"), 
           tpfl("-1.77686212323486749509362995970259857e-01")]

# Velocity from All_FTID.jl (fine-tuned specifically for Alpha Centauri)
Vst_FTID = [tpfl("-7.29279538728389700446721236547889991e-02"), 
            tpfl("-5.57596913770760242903838491049976388e-02"), 
            tpfl("-1.77686156920466641922704142460267104e-01")]

# Calculate differences
dV = Vst_All - Vst_FTID
dV_magnitude = norm(dV)

# Calculate relative differences
rel_diff = dV ./ Vst_All
rel_diff_magnitude = norm(rel_diff)

# Calculate speeds
speed_All = norm(Vst_All)
speed_FTID = norm(Vst_FTID)
speed_diff = speed_All - speed_FTID

println("=== VELOCITY COMPARISON ===")
println("All.jl velocity:     ", Vst_All)
println("All_FTID.jl velocity:", Vst_FTID)
println()
println("Absolute differences: ", dV)
println("Difference magnitude: ", dV_magnitude)
println("Relative differences: ", rel_diff)
println("Relative diff magnitude: ", rel_diff_magnitude)
println()
println("Speed All.jl:     ", speed_All, " c")
println("Speed All_FTID.jl: ", speed_FTID, " c")
println("Speed difference:  ", speed_diff, " c")
println("Relative speed diff: ", speed_diff/speed_All)

# Integration time and distance calculation
AU = tpfl(1.495978707E+11)  # meters
c = tpfl(2.99792458E+8)     # m/s
tcl = tpfl(1.35858929869399063971729181215165883e+14)  # geometric time units
Msol2m = tpfl(1476.67)      # geometric unit conversion

# Convert integration time to seconds
tcl_seconds = tcl * Msol2m / c
tcl_years = tcl_seconds / (365.25 * 24 * 3600)

println()
println("=== TRAJECTORY ANALYSIS ===")
println("Integration time: ", tcl_seconds, " seconds")
println("Integration time: ", tcl_years, " years")

# Calculate trajectory difference over integration time
# Convert velocities from units of c to m/s
Vst_All_ms = Vst_All * c
Vst_FTID_ms = Vst_FTID * c
dV_ms = dV * c

# Position difference after integration time
position_diff = dV_ms * tcl_seconds
position_diff_AU = norm(position_diff) / AU

println("Velocity difference: ", norm(dV_ms), " m/s")
println("Position difference after ", tcl_years, " years:")
println("  ", norm(position_diff), " meters")
println("  ", position_diff_AU, " AU")

# Compare with observed miss distance differences
miss_All = 0.135419  # AU
miss_FTID = 1.986e-07  # AU
miss_ratio = miss_All / miss_FTID

println()
println("=== MISS DISTANCE COMPARISON ===")
println("Miss distance All.jl:     ", miss_All, " AU")
println("Miss distance All_FTID.jl: ", miss_FTID, " AU")
println("Ratio (All/FTID):         ", miss_ratio)
println("Expected from velocity:    ", position_diff_AU, " AU")
