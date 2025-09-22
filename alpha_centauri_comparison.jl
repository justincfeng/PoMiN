#!/usr/bin/env julia

# Compare Alpha Centauri position data from different sources
using Printf, LinearAlgebra

println("="^80)
println("ALPHA CENTAURI POSITION COMPARISON")
println("="^80)

# 1. Existing values in the script (geometric units)
qalpha_existing = [-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13]

# 2. Values from qalphapc converted to geometric units
qalphapc_script = [0.95845, -0.93402, -0.01601]  # pc (from script)
pc_to_m = 3.0857e16  # m/pc
Msol2m = 1476.67  # m (GM_sun/c^2)
qalpha_from_script_pc = qalphapc_script .* pc_to_m ./ Msol2m

# 3. Values from astropy calculation (ICRS coordinates)
qalpha_astropy = [-1.042788e+13, -8.719694e+12, -2.435599e+13]  # geometric units
qalphapc_astropy = [-0.499029, -0.417283, -1.165563]  # pc

println("POSITION DATA COMPARISON:")
println("-"^80)

println("1. EXISTING VALUES (from script):")
println("   qalpha = [$(Printf.@sprintf("%.6e", qalpha_existing[1])), $(Printf.@sprintf("%.6e", qalpha_existing[2])), $(Printf.@sprintf("%.6e", qalpha_existing[3]))] (geom)")
println("   |qalpha| = $(Printf.@sprintf("%.6e", norm(qalpha_existing))) (geom)")

println("\n2. FROM SCRIPT qalphapc CONVERTED:")
println("   qalphapc = [$(qalphapc_script[1]), $(qalphapc_script[2]), $(qalphapc_script[3])] (pc)")
println("   qalpha = [$(Printf.@sprintf("%.6e", qalpha_from_script_pc[1])), $(Printf.@sprintf("%.6e", qalpha_from_script_pc[2])), $(Printf.@sprintf("%.6e", qalpha_from_script_pc[3]))] (geom)")
println("   |qalpha| = $(Printf.@sprintf("%.6e", norm(qalpha_from_script_pc))) (geom)")

println("\n3. FROM ASTROPY CALCULATION (ICRS):")
println("   qalphapc = [$(qalphapc_astropy[1]), $(qalphapc_astropy[2]), $(qalphapc_astropy[3])] (pc)")
println("   qalpha = [$(Printf.@sprintf("%.6e", qalpha_astropy[1])), $(Printf.@sprintf("%.6e", qalpha_astropy[2])), $(Printf.@sprintf("%.6e", qalpha_astropy[3]))] (geom)")
println("   |qalpha| = $(Printf.@sprintf("%.6e", norm(qalpha_astropy))) (geom)")

println("\n" * "="^80)
println("COMPARISON ANALYSIS:")
println("="^80)

# Compare existing vs astropy
diff_existing_astropy = qalpha_existing .- qalpha_astropy
rel_diff_existing_astropy = abs.(diff_existing_astropy ./ qalpha_existing) .* 100

println("\nEXISTING vs ASTROPY:")
println("   Difference: [$(Printf.@sprintf("%.3e", diff_existing_astropy[1])), $(Printf.@sprintf("%.3e", diff_existing_astropy[2])), $(Printf.@sprintf("%.3e", diff_existing_astropy[3]))] (geom)")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_existing_astropy[1]))%, $(Printf.@sprintf("%.1f", rel_diff_existing_astropy[2]))%, $(Printf.@sprintf("%.1f", rel_diff_existing_astropy[3]))%]")
println("   |Δq| = $(Printf.@sprintf("%.3e", norm(diff_existing_astropy))) (geom)")
println("   Relative |Δq| = $(Printf.@sprintf("%.1f", norm(diff_existing_astropy)/norm(qalpha_existing)*100))%")

# Compare script pc vs astropy
diff_script_astropy = qalpha_from_script_pc .- qalpha_astropy
rel_diff_script_astropy = abs.(diff_script_astropy ./ qalpha_from_script_pc) .* 100

println("\nSCRIPT qalphapc vs ASTROPY:")
println("   Difference: [$(Printf.@sprintf("%.3e", diff_script_astropy[1])), $(Printf.@sprintf("%.3e", diff_script_astropy[2])), $(Printf.@sprintf("%.3e", diff_script_astropy[3]))] (geom)")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_script_astropy[1]))%, $(Printf.@sprintf("%.1f", rel_diff_script_astropy[2]))%, $(Printf.@sprintf("%.1f", rel_diff_script_astropy[3]))%]")
println("   |Δq| = $(Printf.@sprintf("%.3e", norm(diff_script_astropy))) (geom)")
println("   Relative |Δq| = $(Printf.@sprintf("%.1f", norm(diff_script_astropy)/norm(qalpha_from_script_pc)*100))%")

# Compare existing vs script pc
diff_existing_script = qalpha_existing .- qalpha_from_script_pc
rel_diff_existing_script = abs.(diff_existing_script ./ qalpha_existing) .* 100

println("\nEXISTING vs SCRIPT qalphapc:")
println("   Difference: [$(Printf.@sprintf("%.3e", diff_existing_script[1])), $(Printf.@sprintf("%.3e", diff_existing_script[2])), $(Printf.@sprintf("%.3e", diff_existing_script[3]))] (geom)")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_existing_script[1]))%, $(Printf.@sprintf("%.1f", rel_diff_existing_script[2]))%, $(Printf.@sprintf("%.1f", rel_diff_existing_script[3]))%]")
println("   |Δq| = $(Printf.@sprintf("%.3e", norm(diff_existing_script))) (geom)")
println("   Relative |Δq| = $(Printf.@sprintf("%.1f", norm(diff_existing_script)/norm(qalpha_existing)*100))%")

println("\n" * "="^80)
println("SUMMARY & RECOMMENDATIONS:")
println("="^80)

println("\nBest Match Analysis:")
if norm(diff_existing_astropy) < norm(diff_existing_script)
    println("   ✓ EXISTING values match ASTROPY better than script qalphapc")
    println("   → The existing qalpha values appear to be from a similar ICRS calculation")
else
    println("   ✓ EXISTING values match SCRIPT qalphapc better than astropy")
    println("   → The existing qalpha values may be from the script's qalphapc conversion")
end

println("\nAccuracy Assessment:")
if norm(diff_existing_astropy)/norm(qalpha_existing) < 0.01
    println("   ✓ EXCELLENT: Existing vs Astropy difference < 1%")
elseif norm(diff_existing_astropy)/norm(qalpha_existing) < 0.05
    println("   ✓ GOOD: Existing vs Astropy difference < 5%")
elseif norm(diff_existing_astropy)/norm(qalpha_existing) < 0.1
    println("   ⚠ MODERATE: Existing vs Astropy difference < 10%")
else
    println("   ⚠ LARGE: Existing vs Astropy difference > 10%")
end

println("\nRecommendation:")
if norm(diff_existing_astropy)/norm(qalpha_existing) < 0.05
    println("   → Keep existing qalpha values - they're already quite accurate")
    println("   → The small differences may be due to:")
    println("     • Different epoch (proper motion effects)")
    println("     • Different astrometric catalog")
    println("     • Rounding differences")
else
    println("   → Consider updating to astropy-derived values for better accuracy")
    println("   → Use: qalpha = [$(Printf.@sprintf("%.6e", qalpha_astropy[1])), $(Printf.@sprintf("%.6e", qalpha_astropy[2])), $(Printf.@sprintf("%.6e", qalpha_astropy[3]))]")
end

println("\n" * "="^80)
