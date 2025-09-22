#!/usr/bin/env julia

# Compare existing momenta with new velocity data
using Printf, LinearAlgebra

# Constants (from target_proxima.jl)
cMKS = 299792458.0  # m/s
tpfl = Float64

# Existing data from the script
println("="^70)
println("VELOCITY COMPARISON: EXISTING MOMENTA vs NEW ASTROPY DATA")
println("="^70)

# Alpha Centauri comparison
println("\nALPHA CENTAURI:")
println("-"^50)

malpha = 2.0429  # solar masses
palpha_existing = [-6.363558578130740E-05, 1.508126516235040E-04, 1.475021586009610E-04]  # geometric units

# Convert existing momentum to velocity (assuming non-relativistic: p = mv)
valpha_existing_c = palpha_existing ./ malpha  # velocity in units of c
valpha_existing_km_s = valpha_existing_c .* cMKS / 1000  # convert to km/s

# New astropy-derived velocity
valpha_new_km_s = [-21.600000, -23.271732, 3.048579]  # km/s
valpha_new_c = valpha_new_km_s .* 1000 / cMKS  # convert to units of c

println("Existing momentum-derived velocity:")
println("  p_alpha = $(palpha_existing)")
println("  v_alpha = $(valpha_existing_c) c")
println("  v_alpha = $(valpha_existing_km_s) km/s")
println("  |v_alpha| = $(norm(valpha_existing_km_s)) km/s")

println("\nNew astropy-derived velocity:")
println("  v_alpha = $(valpha_new_c) c")
println("  v_alpha = $(valpha_new_km_s) km/s")
println("  |v_alpha| = $(norm(valpha_new_km_s)) km/s")

println("\nDifference:")
diff_alpha_km_s = valpha_new_km_s .- valpha_existing_km_s
diff_alpha_c = valpha_new_c .- valpha_existing_c
println("  Δv = $(diff_alpha_km_s) km/s")
println("  Δv = $(diff_alpha_c) c")
println("  |Δv| = $(norm(diff_alpha_km_s)) km/s")
println("  Relative difference: $(norm(diff_alpha_km_s)/norm(valpha_existing_km_s)*100)%")

# Earth comparison
println("\n\nEARTH:")
println("-"^50)

mEarth = 3.0033693739E-06  # solar masses
pEarth_existing = [-2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11]  # geometric units

# Convert existing momentum to velocity (assuming non-relativistic: p = mv)
vEarth_existing_c = pEarth_existing ./ mEarth  # velocity in units of c
vEarth_existing_km_s = vEarth_existing_c .* cMKS / 1000  # convert to km/s

# New astropy-derived velocity
vEarth_new_km_s = [-29.784947, -5.029754, -2.180644]  # km/s
vEarth_new_c = vEarth_new_km_s .* 1000 / cMKS  # convert to units of c

println("Existing momentum-derived velocity:")
println("  p_Earth = $(pEarth_existing)")
println("  v_Earth = $(vEarth_existing_c) c")
println("  v_Earth = $(vEarth_existing_km_s) km/s")
println("  |v_Earth| = $(norm(vEarth_existing_km_s)) km/s")

println("\nNew astropy-derived velocity:")
println("  v_Earth = $(vEarth_new_c) c")
println("  v_Earth = $(vEarth_new_km_s) km/s")
println("  |v_Earth| = $(norm(vEarth_new_km_s)) km/s")

println("\nDifference:")
diff_Earth_km_s = vEarth_new_km_s .- vEarth_existing_km_s
diff_Earth_c = vEarth_new_c .- vEarth_existing_c
println("  Δv = $(diff_Earth_km_s) km/s")
println("  Δv = $(diff_Earth_c) c")
println("  |Δv| = $(norm(diff_Earth_km_s)) km/s")
println("  Relative difference: $(norm(diff_Earth_km_s)/norm(vEarth_existing_km_s)*100)%")

println("\n" * "="^70)
println("SUMMARY:")
println("="^70)
println("Alpha Centauri velocity difference: $(Printf.@sprintf("%.3f", norm(diff_alpha_km_s))) km/s ($(Printf.@sprintf("%.1f", norm(diff_alpha_km_s)/norm(valpha_existing_km_s)*100))%)")
println("Earth velocity difference: $(Printf.@sprintf("%.3f", norm(diff_Earth_km_s))) km/s ($(Printf.@sprintf("%.1f", norm(diff_Earth_km_s)/norm(vEarth_existing_km_s)*100))%)")

# Check if differences are significant
println("\nSignificance Assessment:")
if norm(diff_alpha_km_s)/norm(valpha_existing_km_s) > 0.1
    println("  Alpha Centauri: LARGE difference (>10%) - consider updating")
elseif norm(diff_alpha_km_s)/norm(valpha_existing_km_s) > 0.01
    println("  Alpha Centauri: MODERATE difference (1-10%) - may want to update")
else
    println("  Alpha Centauri: SMALL difference (<1%) - existing values are good")
end

if norm(diff_Earth_km_s)/norm(vEarth_existing_km_s) > 0.1
    println("  Earth: LARGE difference (>10%) - consider updating")
elseif norm(diff_Earth_km_s)/norm(vEarth_existing_km_s) > 0.01
    println("  Earth: MODERATE difference (1-10%) - may want to update")
else
    println("  Earth: SMALL difference (<1%) - existing values are good")
end
