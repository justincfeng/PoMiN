#-----------------------------------------------------------------------
#
#   NEWTONIAN VS RELATIVISTIC COMPARISON WITH FINE-TUNING
#   Spacecraft trajectory to Proxima Centauri
#   Includes Broyden optimization for both cases
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots
using DoubleFloats, ForwardDiff
tpfl  = Double64

include("target_proxima.jl")
include("../utils/broyden.jl")

#-----------------------------------------------------------------------
#   COMMON SETUP
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0
vst = norm(Vst)

# Proxima Centauri
mProx = tpfl(0.1221)
qProx = Xpx0
pProx = mProx .* γVpx
PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)

# Sun
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])
Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)

# Jupiter
mjup = tpfl(0.000954)
qjup = tpfl.([5.2*AU, 0.0, 0.0])
vorb = tpfl(13.1E+3) / cMKS

# Integration parameters
tspan = (tpfl(0), tcl*tpfl(1.1))
tols  = tpfl(1e-16)

# Target position
q_target_pos_Flat = XbF(tcl)

#-----------------------------------------------------------------------
#   RELATIVISTIC CASE
#-----------------------------------------------------------------------

println("="^70)
println("RELATIVISTIC CASE")
println("="^70)

# Relativistic Lorentz factors
γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
γjup = one(tpfl)/sqrt(one(tpfl)-vorb^2)

# Relativistic momenta
pchip_rel = (mchip*γst) .* Vst
pjup_rel = tpfl.([0.0, mjup * γjup * vorb, 0.0])

# Particle objects
Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pjup_rel = pomin.setup_single_particle(mjup, qjup, pjup_rel, tpfl)

# System setup
PintSPJ_rel = pomin.merge_particle_systems(Psol, PProx, Pjup_rel)
PtestSPJ_rel = Pchip_rel

paramsSPJ_rel = pomin.ParametersJulia(tspan, integrator="Vern9", 
                                     atol=tpfl(1e-24), rtol=tpfl(1e-24))

# Solve uncorrected trajectory
solSPJ_rel = pomin.solve(PintSPJ_rel, paramsSPJ_rel; testparticles=PtestSPJ_rel)

ZendSPJ_rel = solSPJ_rel(tcl)
q_spacecraft_SPJ_rel = ZendSPJ_rel[19:21]
q_proxima_SPJ_rel = ZendSPJ_rel[4:6]
q_target_pos_SPJ_rel = q_proxima_SPJ_rel + bv
miss_SPJ_rel = norm(q_spacecraft_SPJ_rel - q_target_pos_SPJ_rel) / AU

println("Uncorrected miss distance: $(miss_SPJ_rel) AU")

# Crude correction
dxe_rel = q_spacecraft_SPJ_rel - q_target_pos_SPJ_rel
δV_rel = - dxe_rel ./ tcl
vcorr_rel = (Vst .+ (δV_rel ./ γst))
γstcorr_rel = one(tpfl)/sqrt(one(tpfl)-(norm(vcorr_rel)/c)^2)
pchipcorr_rel = (mchip*γstcorr_rel) .* vcorr_rel

Pchipcorr_rel = pomin.setup_single_particle(mchip, qchip, pchipcorr_rel, tpfl)

solSPJc_rel = pomin.solve(PintSPJ_rel, paramsSPJ_rel; testparticles=Pchipcorr_rel)
ZendSPJc_rel = solSPJc_rel(tcl)
q_spacecraft_SPJc_rel = ZendSPJc_rel[19:21]
q_target_pos_SPJc_rel = ZendSPJc_rel[4:6] + bv
miss_SPJc_rel = norm(q_spacecraft_SPJc_rel - q_target_pos_SPJc_rel) / AU

println("Corrected miss distance: $(miss_SPJc_rel) AU")

# Target function for relativistic case
function targetrootconstructor_rel(Pint,Ptest0,params,bv,tcl)
    return function f(v)
        v_mag = norm(v)
        γ = one(eltype(v)) / sqrt(one(eltype(v)) - (v_mag/c)^2)
        pv = (Ptest0.m[1] * γ) .* v
        Ptest = pomin.setup_single_particle(Ptest0.m[1], Ptest0.q[1], pv, eltype(Ptest0.m))
        sol_modified = pomin.solve(Pint, params; testparticles=Ptest)
        Zend_modified = sol_modified(tcl)
        q_spacecraft = Zend_modified[19:21]
        q_proxima = Zend_modified[4:6]
        q_target = q_proxima + bv
        return q_spacecraft - q_target
    end
end

FSPJc_rel = targetrootconstructor_rel(PintSPJ_rel, Pchipcorr_rel, paramsSPJ_rel, bv, tcl)

# Broyden optimization (relativistic)
println("Computing Jacobian (relativistic)...")
J_init_rel = ForwardDiff.jacobian(FSPJc_rel, vcorr_rel)

println("Broyden iterations (relativistic)...")
v_finetuned_rel = bsolve(FSPJc_rel, J_init_rel, FSPJc_rel(vcorr_rel), vcorr_rel, 10)

# Final relativistic result
f_final_rel = FSPJc_rel(v_finetuned_rel)
miss_final_rel = norm(f_final_rel) / AU

# Compute fine-tuned trajectory for relativistic case
v_magft_rel = norm(v_finetuned_rel)
γft_rel = one(tpfl)/sqrt(one(tpfl)-(v_magft_rel/c)^2)
pvft_rel = (mchip*γft_rel) .* v_finetuned_rel

Ptestft_rel = pomin.setup_single_particle(mchip, qchip, pvft_rel, tpfl)
sol_ft_rel = pomin.solve(PintSPJ_rel, paramsSPJ_rel; testparticles=Ptestft_rel)
Zend_ft_rel = sol_ft_rel(tcl)
q_spacecraft_ft_rel = Zend_ft_rel[19:21]
q_proxima_ft_rel = Zend_ft_rel[4:6]

println("Fine-tuned miss distance (relativistic initial data): $(miss_final_rel) AU")
println("Improvement factor (relativistic initial data): $(norm(FSPJc_rel(vcorr_rel)) / norm(f_final_rel))")

#-----------------------------------------------------------------------
#   NEWTONIAN CASE
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("NEWTONIAN CASE")
println("="^70)

# Newtonian case: no Lorentz factors
γst_newt = one(tpfl)
γjup_newt = one(tpfl)

# Newtonian momenta
pchip_newt = (mchip*γst_newt) .* Vst
pjup_newt = tpfl.([0.0, mjup * γjup_newt * vorb, 0.0])

# Particle objects
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)
Pjup_newt = pomin.setup_single_particle(mjup, qjup, pjup_newt, tpfl)

# System setup
PintSPJ_newt = pomin.merge_particle_systems(Psol, PProx, Pjup_newt)
PtestSPJ_newt = Pchip_newt

paramsSPJ_newt = pomin.ParametersJulia(tspan, integrator="Vern9", 
                                      atol=tols, rtol=tols)

# Solve uncorrected trajectory
solSPJ_newt = pomin.solve(PintSPJ_newt, paramsSPJ_newt; 
                         testparticles=PtestSPJ_newt, Newtonian=true)

ZendSPJ_newt = solSPJ_newt(tcl)
q_spacecraft_SPJ_newt = ZendSPJ_newt[19:21]
q_proxima_SPJ_newt = ZendSPJ_newt[4:6]
q_target_pos_SPJ_newt = q_proxima_SPJ_newt + bv
miss_SPJ_newt = norm(q_spacecraft_SPJ_newt - q_target_pos_SPJ_newt) / AU

println("Uncorrected miss distance: $(miss_SPJ_newt) AU")

# Crude correction
dxe_newt = q_spacecraft_SPJ_newt - q_target_pos_SPJ_newt
δV_newt = - dxe_newt ./ tcl
vcorr_newt = (Vst .+ (δV_newt ./ γst_newt))
pchipcorr_newt = (mchip*γst_newt) .* vcorr_newt

Pchipcorr_newt = pomin.setup_single_particle(mchip, qchip, pchipcorr_newt, tpfl)

solSPJc_newt = pomin.solve(PintSPJ_newt, paramsSPJ_newt; 
                          testparticles=Pchipcorr_newt, Newtonian=true)
ZendSPJc_newt = solSPJc_newt(tcl)
q_spacecraft_SPJc_newt = ZendSPJc_newt[19:21]
q_target_pos_SPJc_newt = ZendSPJc_newt[4:6] + bv
miss_SPJc_newt = norm(q_spacecraft_SPJc_newt - q_target_pos_SPJc_newt) / AU

println("Corrected miss distance: $(miss_SPJc_newt) AU")

# Target function for Newtonian case
function targetrootconstructor_newt(Pint,Ptest0,params,bv,tcl)
    return function f(v)
        pv = Ptest0.m[1] .* v  # Newtonian momentum
        Ptest = pomin.setup_single_particle(Ptest0.m[1], Ptest0.q[1], pv, eltype(Ptest0.m))
        sol_modified = pomin.solve(Pint, params; testparticles=Ptest, Newtonian=true)
        Zend_modified = sol_modified(tcl)
        q_spacecraft = Zend_modified[19:21]
        q_proxima = Zend_modified[4:6]
        q_target = q_proxima + bv
        return q_spacecraft - q_target
    end
end

FSPJc_newt = targetrootconstructor_newt(PintSPJ_newt, Pchipcorr_newt, paramsSPJ_newt, bv, tcl)

# Broyden optimization (Newtonian)
println("Computing Jacobian (Newtonian)...")
J_init_newt = ForwardDiff.jacobian(FSPJc_newt, vcorr_newt)

println("Broyden iterations (Newtonian)...")
v_finetuned_newt = bsolve(FSPJc_newt, J_init_newt, FSPJc_newt(vcorr_newt), vcorr_newt, 10)

# Final Newtonian result
f_final_newt = FSPJc_newt(v_finetuned_newt)
miss_final_newt = norm(f_final_newt) / AU

# Compute fine-tuned trajectory for Newtonian case
v_magft_newt = norm(v_finetuned_newt)
pvft_newt = mchip .* v_finetuned_newt  # Newtonian momentum

Ptestft_newt = pomin.setup_single_particle(mchip, qchip, pvft_newt, tpfl)
sol_ft_newt = pomin.solve(PintSPJ_newt, paramsSPJ_newt; testparticles=Ptestft_newt, Newtonian=true)
Zend_ft_newt = sol_ft_newt(tcl)
q_spacecraft_ft_newt = Zend_ft_newt[19:21]
q_proxima_ft_newt = Zend_ft_newt[4:6]

println("Fine-tuned miss distance: $(miss_final_newt) AU")
println("Improvement factor: $(norm(FSPJc_newt(vcorr_newt)) / norm(f_final_newt))")

#-----------------------------------------------------------------------
#   RELATIVISTIC CASE WITH NEWTONIAN INITIAL VELOCITY
#-----------------------------------------------------------------------

# Relativistic Lorentz factors
γst_ftn = one(tpfl)/sqrt(one(tpfl)-(norm(v_finetuned_newt)/c)^2)

# Relativistic momenta
pchip_ftn = (mchip*γst_ftn) .* v_finetuned_newt

Pchip_ftn = pomin.setup_single_particle(mchip, qchip, pchip_ftn, tpfl)

sol_ftn_rel = pomin.solve(PintSPJ_rel, paramsSPJ_rel; testparticles=Pchip_ftn)
Zend_ftn_rel = sol_ftn_rel(tcl)
q_spacecraft_ftn_rel = Zend_ftn_rel[19:21]
q_proxima_ftn_rel = Zend_ftn_rel[4:6]
q_target_pos_ftn_rel = q_proxima_ftn_rel + bv
miss_ftn_rel = norm(q_spacecraft_ftn_rel - q_target_pos_ftn_rel) / AU

println("Fine-tuned miss distance (Newtonian initial velocity): $(miss_ftn_rel) AU")

#-----------------------------------------------------------------------
#   COMPARISON RESULTS
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("COMPARISON RESULTS")
println("="^70)
println()

println("Spacecraft Parameters:")
println("  Target distance from Proxima: $(b0/AU) AU")
println("  Time to closest approach: $(tcl) (geometric units)")
println("  Spacecraft velocity: $(vst) c")
println("  Relativistic γ factor: $(γst)")
println()

println("Miss Distances (AU):")
println("                                    Relativistic    Newtonian    Rel(Newt vel)")
println("  Uncorrected:                     $(Printf.@sprintf("%12.6e", miss_SPJ_rel))    $(Printf.@sprintf("%12.6e", miss_SPJ_newt))    N/A")
println("  Crude correction:                $(Printf.@sprintf("%12.6e", miss_SPJc_rel))    $(Printf.@sprintf("%12.6e", miss_SPJc_newt))    N/A")
println("  Fine-tuned (Broyden):            $(Printf.@sprintf("%12.6e", miss_final_rel))    $(Printf.@sprintf("%12.6e", miss_final_newt))    $(Printf.@sprintf("%12.6e", miss_ftn_rel))")
println()

println("Velocity Corrections:")
println("  Relativistic crude correction:   $(norm(vcorr_rel - Vst)) c")
println("  Relativistic fine correction:    $(norm(v_finetuned_rel - vcorr_rel)) c")
println("  Newtonian crude correction:      $(norm(vcorr_newt - Vst)) c")
println("  Newtonian fine correction:       $(norm(v_finetuned_newt - vcorr_newt)) c")
println()

println("Relativistic Effects:")
println("  Difference in uncorrected miss:  $(abs(miss_SPJ_rel - miss_SPJ_newt)) AU")
println("  Difference in corrected miss:    $(abs(miss_SPJc_rel - miss_SPJc_newt)) AU")
println("  Difference in fine-tuned miss:   $(abs(miss_final_rel - miss_final_newt)) AU")
println("  Rel(Newt vel) vs Newtonian:      $(abs(miss_ftn_rel - miss_final_newt)) AU")
println()

println("Improvement Factors:")
rel_improvement = miss_SPJ_rel / miss_final_rel
newt_improvement = miss_SPJ_newt / miss_final_newt
println("  Relativistic (uncorrected → fine-tuned): $(Printf.@sprintf("%.2e", rel_improvement))")
println("  Newtonian (uncorrected → fine-tuned):    $(Printf.@sprintf("%.2e", newt_improvement))")
println()

println("Final Velocities:")
println("  Relativistic fine-tuned speed:   $(norm(v_finetuned_rel)) c")
println("  Newtonian fine-tuned speed:      $(norm(v_finetuned_newt)) c")
println("  Speed difference:                $(abs(norm(v_finetuned_rel) - norm(v_finetuned_newt))) c")
println()

# Relativistic vs Newtonian significance
rel_significance = abs(miss_final_rel - miss_final_newt) / min(miss_final_rel, miss_final_newt)
println("Relativistic Significance:")
println("  Relative difference in final miss: $(Printf.@sprintf("%.2e", rel_significance)) ($(Printf.@sprintf("%.1f", rel_significance*100))%)")

if rel_significance > 0.01
    println("  → Relativistic effects are SIGNIFICANT (>1%)")
elseif rel_significance > 0.001
    println("  → Relativistic effects are MODERATE (0.1-1%)")
else
    println("  → Relativistic effects are SMALL (<0.1%)")
end

#-----------------------------------------------------------------------
#   POSITION DIFFERENCE ANALYSIS (FINE-TUNED TRAJECTORIES)
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("POSITION DIFFERENCE ANALYSIS (FINE-TUNED TRAJECTORIES)")
println("="^70)
println()

# Extract velocity vectors at endpoint for fine-tuned trajectories
v_ft_newton = Zend_ft_newt[16:18] / mchip  # p/m for Newtonian
v_ft_rel = Zend_ft_rel[16:18] / (mchip * γft_rel)  # p/(γm) for relativistic

# Velocity difference
Δv_ft = v_ft_rel - v_ft_newton
Δv_ft_magnitude = norm(Δv_ft)

# Use relativistic velocity as reference direction
v_ft_rel_unit = v_ft_rel / norm(v_ft_rel)

# Position difference between fine-tuned trajectories
Δq_spacecraft_ft = q_spacecraft_ft_rel - q_spacecraft_ft_newt

# Components relative to velocity direction
# Parallel component (along velocity)
Δq_parallel_ft = dot(Δq_spacecraft_ft, v_ft_rel_unit) * v_ft_rel_unit
Δq_parallel_ft_magnitude = abs(dot(Δq_spacecraft_ft, v_ft_rel_unit))

# Transverse component (perpendicular to velocity)
Δq_transverse_ft = Δq_spacecraft_ft - Δq_parallel_ft
Δq_transverse_ft_magnitude = norm(Δq_transverse_ft)

# Convert to AU and km
Δq_spacecraft_ft_AU = norm(Δq_spacecraft_ft) / AU
Δq_parallel_ft_AU = Δq_parallel_ft_magnitude / AU
Δq_transverse_ft_AU = Δq_transverse_ft_magnitude / AU
km = AU / 1.496e8
Δq_spacecraft_ft_km = norm(Δq_spacecraft_ft) / km
Δq_parallel_ft_km = Δq_parallel_ft_magnitude / km
Δq_transverse_ft_km = Δq_transverse_ft_magnitude / km

# Velocity angle difference
cos_angle_ft = dot(v_ft_newton, v_ft_rel) / (norm(v_ft_newton) * norm(v_ft_rel))
cos_angle_ft = clamp(cos_angle_ft, -1.0, 1.0)  # Ensure valid range for acos
angle_diff_ft_rad = acos(cos_angle_ft)
angle_diff_ft_deg = angle_diff_ft_rad * 180.0 / π
angle_diff_ft_arcsec = angle_diff_ft_deg * 3600.0

println("Fine-Tuned Trajectory Differences:")
println("Total position difference:")
println("  Δq_spacecraft = $(Printf.@sprintf("%.6e", Δq_spacecraft_ft_AU)) AU")
println("  Δq_spacecraft = $(Printf.@sprintf("%.6e", Δq_spacecraft_ft_km)) km")
println()
println("Components relative to velocity direction:")
println("  Parallel (along velocity):")
println("    Δq_parallel = $(Printf.@sprintf("%.6e", Δq_parallel_ft_AU)) AU")
println("    Δq_parallel = $(Printf.@sprintf("%.6e", Δq_parallel_ft_km)) km")
println("  Transverse (perpendicular to velocity):")
println("    Δq_transverse = $(Printf.@sprintf("%.6e", Δq_transverse_ft_AU)) AU")
println("    Δq_transverse = $(Printf.@sprintf("%.6e", Δq_transverse_ft_km)) km")
println()
println("Velocity comparison (fine-tuned):")
println("  |v_newton| = $(Printf.@sprintf("%.12f", norm(v_ft_newton))) c")
println("  |v_rel| = $(Printf.@sprintf("%.12f", norm(v_ft_rel))) c")
println("  Δv magnitude = $(Printf.@sprintf("%.6e", Δv_ft_magnitude)) c")
println("  Angle between velocities = $(Printf.@sprintf("%.6e", angle_diff_ft_deg)) degrees")
println("  Angle between velocities = $(Printf.@sprintf("%.6e", angle_diff_ft_arcsec)) arcseconds")
println()
println("Verification: √(parallel² + transverse²) = $(Printf.@sprintf("%.6e", sqrt(Δq_parallel_ft_AU^2 + Δq_transverse_ft_AU^2))) AU")

# Component analysis
parallel_fraction = Δq_parallel_ft_AU / Δq_spacecraft_ft_AU
transverse_fraction = Δq_transverse_ft_AU / Δq_spacecraft_ft_AU

println()
println("Component Analysis:")
println("  Parallel component: $(Printf.@sprintf("%.1f", parallel_fraction*100))% of total difference")
println("  Transverse component: $(Printf.@sprintf("%.1f", transverse_fraction*100))% of total difference")

# Physical interpretation
println()
println("Physical Interpretation:")
if Δq_parallel_ft_AU > Δq_transverse_ft_AU
    println("  → Relativistic effects are primarily LONGITUDINAL (along flight direction)")
else
    println("  → Relativistic effects are primarily TRANSVERSE (perpendicular to flight)")
end

# Scale comparison
earth_radius_AU = 6371e3 / (AU / km)  # Earth radius in AU
println("  → Total difference is $(Printf.@sprintf("%.1f", Δq_spacecraft_ft_AU / earth_radius_AU)) Earth radii")
println("  → Parallel difference is $(Printf.@sprintf("%.1f", Δq_parallel_ft_AU / earth_radius_AU)) Earth radii")
println("  → Transverse difference is $(Printf.@sprintf("%.1f", Δq_transverse_ft_AU / earth_radius_AU)) Earth radii")

#-----------------------------------------------------------------------
#   POSITION DIFFERENCE: REL(NEWT VEL) vs NEWTONIAN
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("POSITION DIFFERENCE: RELATIVISTIC(NEWT VEL) vs NEWTONIAN")
println("="^70)
println()

# Position difference between relativistic with Newtonian velocity and pure Newtonian
Δq_ftn_newt = q_spacecraft_ftn_rel - q_spacecraft_ft_newt

# Use Newtonian velocity as reference direction for this comparison
v_newt_unit = v_ft_newton / norm(v_ft_newton)

# Components relative to Newtonian velocity direction
Δq_parallel_ftn = dot(Δq_ftn_newt, v_newt_unit) * v_newt_unit
Δq_parallel_ftn_magnitude = abs(dot(Δq_ftn_newt, v_newt_unit))

Δq_transverse_ftn = Δq_ftn_newt - Δq_parallel_ftn
Δq_transverse_ftn_magnitude = norm(Δq_transverse_ftn)

# Convert to AU and km
Δq_ftn_newt_AU = norm(Δq_ftn_newt) / AU
Δq_parallel_ftn_AU = Δq_parallel_ftn_magnitude / AU
Δq_transverse_ftn_AU = Δq_transverse_ftn_magnitude / AU
Δq_ftn_newt_km = norm(Δq_ftn_newt) / km
Δq_parallel_ftn_km = Δq_parallel_ftn_magnitude / km
Δq_transverse_ftn_km = Δq_transverse_ftn_magnitude / km

println("Relativistic Physics with Newtonian Initial Velocity vs Pure Newtonian:")
println("Total position difference:")
println("  Δq_spacecraft = $(Printf.@sprintf("%.6e", Δq_ftn_newt_AU)) AU")
println("  Δq_spacecraft = $(Printf.@sprintf("%.6e", Δq_ftn_newt_km)) km")
println()
println("Components relative to Newtonian velocity direction:")
println("  Parallel (along Newtonian velocity):")
println("    Δq_parallel = $(Printf.@sprintf("%.6e", Δq_parallel_ftn_AU)) AU")
println("    Δq_parallel = $(Printf.@sprintf("%.6e", Δq_parallel_ftn_km)) km")
println("  Transverse (perpendicular to Newtonian velocity):")
println("    Δq_transverse = $(Printf.@sprintf("%.6e", Δq_transverse_ftn_AU)) AU")
println("    Δq_transverse = $(Printf.@sprintf("%.6e", Δq_transverse_ftn_km)) km")
println()
println("Verification: √(parallel² + transverse²) = $(Printf.@sprintf("%.6e", sqrt(Δq_parallel_ftn_AU^2 + Δq_transverse_ftn_AU^2))) AU")

# Component analysis for this comparison
parallel_fraction_ftn = Δq_parallel_ftn_AU / Δq_ftn_newt_AU
transverse_fraction_ftn = Δq_transverse_ftn_AU / Δq_ftn_newt_AU

println()
println("Component Analysis (Rel(Newt vel) vs Newtonian):")
println("  Parallel component: $(Printf.@sprintf("%.1f", parallel_fraction_ftn*100))% of total difference")
println("  Transverse component: $(Printf.@sprintf("%.1f", transverse_fraction_ftn*100))% of total difference")

println()
println("Physical Interpretation:")
if Δq_parallel_ftn_AU > Δq_transverse_ftn_AU
    println("  → Relativistic physics effects are primarily LONGITUDINAL")
else
    println("  → Relativistic physics effects are primarily TRANSVERSE")
end

# Compare the magnitude of this effect to the full relativistic vs Newtonian difference
ratio_to_full = Δq_ftn_newt_AU / Δq_spacecraft_ft_AU
println("  → This difference is $(Printf.@sprintf("%.1f", ratio_to_full*100))% of the full relativistic vs Newtonian difference")

# Scale comparison
println("  → Total difference is $(Printf.@sprintf("%.1f", Δq_ftn_newt_AU / earth_radius_AU)) Earth radii")

println()
println("="^70)
