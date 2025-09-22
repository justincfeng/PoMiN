#!/usr/bin/env julia

# Convert Alpha Centauri position from parsecs to geometric units
using Printf

# Constants
cMKS = 299792458.0  # m/s
AUMKS = 149597870700.0  # m
Msol2m = 1476.67  # m (solar mass to meters conversion)

# Position in parsecs
qalphapc = [0.95845, -0.93402, -0.01601]  # pc

println("="^60)
println("ALPHA CENTAURI POSITION CONVERSION")
println("="^60)

println("Position in parsecs:")
println("  qalphapc = [$(qalphapc[1]), $(qalphapc[2]), $(qalphapc[3])] pc")

# Convert parsecs to meters
# 1 parsec = 3.0857e16 m
pc_to_m = 3.0857e16  # m/pc

qalpha_m = qalphapc .* pc_to_m
println("\nPosition in meters:")
println("  qalpha_m = [$(Printf.@sprintf("%.6e", qalpha_m[1])), $(Printf.@sprintf("%.6e", qalpha_m[2])), $(Printf.@sprintf("%.6e", qalpha_m[3]))] m")

# Convert to geometric units
# Geometric units: length in units of GM/c^2
# For solar mass: GM/c^2 = Msol2m = 1476.67 m
qalpha_geom = qalpha_m ./ Msol2m

println("\nPosition in geometric units:")
println("  qalpha_geom = [$(Printf.@sprintf("%.12e", qalpha_geom[1])), $(Printf.@sprintf("%.12e", qalpha_geom[2])), $(Printf.@sprintf("%.12e", qalpha_geom[3]))]")

println("\nFor Julia code:")
println("qalpha = tpfl.([$(Printf.@sprintf("%.12E", qalpha_geom[1])),")
println("                $(Printf.@sprintf("%.12E", qalpha_geom[2])),")
println("                $(Printf.@sprintf("%.12E", qalpha_geom[3]))])         # in geometric units")

# Verify by comparing with existing values
qalpha_existing = [-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13]
println("\nComparison with existing values:")
println("Existing: [$(Printf.@sprintf("%.12E", qalpha_existing[1])), $(Printf.@sprintf("%.12E", qalpha_existing[2])), $(Printf.@sprintf("%.12E", qalpha_existing[3]))]")
println("New:      [$(Printf.@sprintf("%.12E", qalpha_geom[1])), $(Printf.@sprintf("%.12E", qalpha_geom[2])), $(Printf.@sprintf("%.12E", qalpha_geom[3]))]")

diff = qalpha_geom .- qalpha_existing
println("Difference: [$(Printf.@sprintf("%.6e", diff[1])), $(Printf.@sprintf("%.6e", diff[2])), $(Printf.@sprintf("%.6e", diff[3]))]")

# Calculate relative differences
rel_diff = abs.(diff ./ qalpha_existing) .* 100
println("Relative difference: [$(Printf.@sprintf("%.3f", rel_diff[1]))%, $(Printf.@sprintf("%.3f", rel_diff[2]))%, $(Printf.@sprintf("%.3f", rel_diff[3]))%]")

println("\n" * "="^60)
