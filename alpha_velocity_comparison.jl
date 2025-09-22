#!/usr/bin/env julia

# Compare Alpha Centauri velocity data from different sources
using Printf, LinearAlgebra

println("="^80)
println("ALPHA CENTAURI VELOCITY COMPARISON")
println("="^80)

# Constants
cMKS = 299792458.0  # m/s
malpha = 2.0429  # solar masses

# 1. valpha from script
valpha_script = [-29.291, 1.710, 13.589]  # km/s

# 2. palpha from script (convert to velocity)
palpha_script = [-6.363558578130740E-05, 1.508126516235040E-04, 1.475021586009610E-04]  # geometric units
valpha_from_palpha_c = palpha_script ./ malpha  # velocity in units of c
valpha_from_palpha_km_s = valpha_from_palpha_c .* cMKS / 1000  # convert to km/s

# 3. vAlphaCen from astropy (new values)
vAlphaCen_astropy = [-21.600000, -23.271732, 3.048579]  # km/s

println("VELOCITY DATA COMPARISON:")
println("-"^80)

println("1. valpha (from script):")
println("   valpha = [$(valpha_script[1]), $(valpha_script[2]), $(valpha_script[3])] km/s")
println("   |valpha| = $(Printf.@sprintf("%.3f", norm(valpha_script))) km/s")

println("\n2. palpha converted to velocity:")
println("   palpha = [$(Printf.@sprintf("%.6e", palpha_script[1])), $(Printf.@sprintf("%.6e", palpha_script[2])), $(Printf.@sprintf("%.6e", palpha_script[3]))] (geom)")
println("   valpha = [$(Printf.@sprintf("%.3f", valpha_from_palpha_km_s[1])), $(Printf.@sprintf("%.3f", valpha_from_palpha_km_s[2])), $(Printf.@sprintf("%.3f", valpha_from_palpha_km_s[3]))] km/s")
println("   |valpha| = $(Printf.@sprintf("%.3f", norm(valpha_from_palpha_km_s))) km/s")

println("\n3. vAlphaCen (from astropy proper motion + radial velocity):")
println("   vAlphaCen = [$(vAlphaCen_astropy[1]), $(vAlphaCen_astropy[2]), $(vAlphaCen_astropy[3])] km/s")
println("   |vAlphaCen| = $(Printf.@sprintf("%.3f", norm(vAlphaCen_astropy))) km/s")

println("\n" * "="^80)
println("COMPARISON ANALYSIS:")
println("="^80)

# Compare valpha vs astropy
diff_valpha_astropy = valpha_script .- vAlphaCen_astropy
rel_diff_valpha_astropy = abs.(diff_valpha_astropy ./ valpha_script) .* 100

println("\nvalpha vs vAlphaCen (astropy):")
println("   Difference: [$(Printf.@sprintf("%.3f", diff_valpha_astropy[1])), $(Printf.@sprintf("%.3f", diff_valpha_astropy[2])), $(Printf.@sprintf("%.3f", diff_valpha_astropy[3]))] km/s")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_valpha_astropy[1]))%, $(Printf.@sprintf("%.1f", rel_diff_valpha_astropy[2]))%, $(Printf.@sprintf("%.1f", rel_diff_valpha_astropy[3]))%]")
println("   |Δv| = $(Printf.@sprintf("%.3f", norm(diff_valpha_astropy))) km/s")
println("   Relative |Δv| = $(Printf.@sprintf("%.1f", norm(diff_valpha_astropy)/norm(valpha_script)*100))%")

# Compare palpha vs astropy
diff_palpha_astropy = valpha_from_palpha_km_s .- vAlphaCen_astropy
rel_diff_palpha_astropy = abs.(diff_palpha_astropy ./ valpha_from_palpha_km_s) .* 100

println("\npalpha (converted) vs vAlphaCen (astropy):")
println("   Difference: [$(Printf.@sprintf("%.3f", diff_palpha_astropy[1])), $(Printf.@sprintf("%.3f", diff_palpha_astropy[2])), $(Printf.@sprintf("%.3f", diff_palpha_astropy[3]))] km/s")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_palpha_astropy[1]))%, $(Printf.@sprintf("%.1f", rel_diff_palpha_astropy[2]))%, $(Printf.@sprintf("%.1f", rel_diff_palpha_astropy[3]))%]")
println("   |Δv| = $(Printf.@sprintf("%.3f", norm(diff_palpha_astropy))) km/s")
println("   Relative |Δv| = $(Printf.@sprintf("%.1f", norm(diff_palpha_astropy)/norm(valpha_from_palpha_km_s)*100))%")

# Compare valpha vs palpha
diff_valpha_palpha = valpha_script .- valpha_from_palpha_km_s
rel_diff_valpha_palpha = abs.(diff_valpha_palpha ./ valpha_script) .* 100

println("\nvalpha vs palpha (converted):")
println("   Difference: [$(Printf.@sprintf("%.3f", diff_valpha_palpha[1])), $(Printf.@sprintf("%.3f", diff_valpha_palpha[2])), $(Printf.@sprintf("%.3f", diff_valpha_palpha[3]))] km/s")
println("   Relative diff: [$(Printf.@sprintf("%.1f", rel_diff_valpha_palpha[1]))%, $(Printf.@sprintf("%.1f", rel_diff_valpha_palpha[2]))%, $(Printf.@sprintf("%.1f", rel_diff_valpha_palpha[3]))%]")
println("   |Δv| = $(Printf.@sprintf("%.3f", norm(diff_valpha_palpha))) km/s")
println("   Relative |Δv| = $(Printf.@sprintf("%.1f", norm(diff_valpha_palpha)/norm(valpha_script)*100))%")

println("\n" * "="^80)
println("CONSISTENCY CHECK:")
println("="^80)

# Check if palpha is consistent with valpha
println("\nChecking if palpha = malpha * valpha_c:")
valpha_c_from_script = valpha_script .* 1000 / cMKS
palpha_expected = malpha .* valpha_c_from_script
palpha_diff = palpha_script .- palpha_expected

println("   valpha (script) in units of c: [$(Printf.@sprintf("%.8e", valpha_c_from_script[1])), $(Printf.@sprintf("%.8e", valpha_c_from_script[2])), $(Printf.@sprintf("%.8e", valpha_c_from_script[3]))]")
println("   Expected palpha: [$(Printf.@sprintf("%.8e", palpha_expected[1])), $(Printf.@sprintf("%.8e", palpha_expected[2])), $(Printf.@sprintf("%.8e", palpha_expected[3]))]")
println("   Actual palpha:   [$(Printf.@sprintf("%.8e", palpha_script[1])), $(Printf.@sprintf("%.8e", palpha_script[2])), $(Printf.@sprintf("%.8e", palpha_script[3]))]")
println("   Difference:      [$(Printf.@sprintf("%.8e", palpha_diff[1])), $(Printf.@sprintf("%.8e", palpha_diff[2])), $(Printf.@sprintf("%.8e", palpha_diff[3]))]")

if norm(palpha_diff) / norm(palpha_script) < 0.01
    println("   ✓ CONSISTENT: palpha matches valpha within 1%")
else
    println("   ⚠ INCONSISTENT: palpha doesn't match valpha ($(Printf.@sprintf("%.1f", norm(palpha_diff)/norm(palpha_script)*100))% difference)")
end

println("\n" * "="^80)
println("SUMMARY & RECOMMENDATIONS:")
println("="^80)

println("\nAccuracy Assessment:")

# Determine which is most accurate
best_match = ""
best_diff = Inf

if norm(diff_valpha_astropy) < best_diff
    best_diff = norm(diff_valpha_astropy)
    best_match = "valpha"
end

if norm(diff_palpha_astropy) < best_diff
    best_diff = norm(diff_palpha_astropy)
    best_match = "palpha"
end

println("   Best match to astropy: $(best_match) ($(Printf.@sprintf("%.1f", best_diff/norm(vAlphaCen_astropy)*100))% difference)")

println("\nRecommendations:")
if norm(diff_valpha_astropy)/norm(valpha_script) < 0.1
    println("   → valpha is reasonably accurate (< 10% difference)")
else
    println("   → valpha has significant differences from astropy data")
end

if norm(diff_palpha_astropy)/norm(valpha_from_palpha_km_s) < 0.1
    println("   → palpha is reasonably accurate (< 10% difference)")
else
    println("   → palpha has significant differences from astropy data")
end

println("\n   For best accuracy, consider using:")
println("   vAlphaCen = [$(vAlphaCen_astropy[1]), $(vAlphaCen_astropy[2]), $(vAlphaCen_astropy[3])] km/s")
println("   (from proper motion + radial velocity measurements)")

println("\n" * "="^80)
