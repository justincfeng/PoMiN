#-----------------------------------------------------------------------
#
#   NEWTONIAN vs RELATIVISTIC COMPARISON
#   Sun + Proxima + Jupiter scenario
#   Direct comparison of endpoints to quantify post-Minkowskian effects
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf
using DoubleFloats
tpfl  = Double64

include("target_proxima.jl")

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP
#-----------------------------------------------------------------------

# Spacecraft parameters
mchip = tpfl(1.0E-30)
qchip = Xst0
vst = norm(Vst)

# Proxima Centauri parameters
mProx = tpfl(0.1221)
qProx = Xpx0
pProx = mProx .* γVpx

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

# Jupiter parameters
mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([5.2*AU, 0.0, 0.0])  # Jupiter position in geometric units
vorb = tpfl(13.1E+3) / cMKS  # Orbital velocity in units of c
γjup = one(tpfl)/sqrt(one(tpfl)-vorb^2)  # Lorentz factor
pjup = tpfl.([0.0, mjup * γjup * vorb, 0.0])  # Jupiter momentum (y-direction)

#-----------------------------------------------------------------------
#   NEWTONIAN CALCULATION
#-----------------------------------------------------------------------

println("Running Newtonian calculation...")

# Newtonian Lorentz factor (set to 1)
γst_newton = one(tpfl)
pchip_newton = (mchip*γst_newton) .* (Vst)

# Create particle objects
Pchip_newton = pomin.setup_single_particle(mchip, qchip, pchip_newton, tpfl)
PProx_newton = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)
Psol_newton = pomin.setup_single_particle(msol, qsol, psol, tpfl)
Pjup_newton = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)

# Merge systems: Sun + Proxima + Jupiter as integrators, spacecraft as test particle
Pint_newton = pomin.merge_particle_systems(
    pomin.merge_particle_systems(Psol_newton, PProx_newton), 
    Pjup_newton
)
Ptest_newton = Pchip_newton

# Integration parameters
tspan = (tpfl(0), tcl*tpfl(1.1))
tols = tpfl(1e-16)
params_newton = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tols, rtol=tols)

# Solve with Newtonian physics
sol_newton = pomin.solve(Pint_newton, params_newton; testparticles=Ptest_newton, Newtonian=true)

# Extract endpoints
Zend_newton = sol_newton(tcl)
q_spacecraft_newton = Zend_newton[19:21]  # Spacecraft position
q_proxima_newton = Zend_newton[4:6]       # Proxima position

#-----------------------------------------------------------------------
#   RELATIVISTIC CALCULATION
#-----------------------------------------------------------------------

println("Running relativistic calculation...")

# Relativistic Lorentz factor
γst_rel = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst_rel) .* (Vst)

# Create particle objects (same masses and initial positions)
Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
PProx_rel = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)
Psol_rel = pomin.setup_single_particle(msol, qsol, psol, tpfl)
Pjup_rel = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)

# Merge systems
Pint_rel = pomin.merge_particle_systems(
    pomin.merge_particle_systems(Psol_rel, PProx_rel), 
    Pjup_rel
)
Ptest_rel = Pchip_rel

# Integration parameters (same as Newtonian)
params_rel = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tols, rtol=tols)

# Solve with post-Minkowskian relativistic physics
sol_rel = pomin.solve(Pint_rel, params_rel; testparticles=Ptest_rel)

# Extract endpoints
Zend_rel = sol_rel(tcl)
q_spacecraft_rel = Zend_rel[19:21]  # Spacecraft position
q_proxima_rel = Zend_rel[4:6]       # Proxima position

# Extract velocity vectors at endpoint
v_newton = Zend_newton[16:18] / mchip  # p/m for Newtonian
v_rel = Zend_rel[16:18] / (mchip * γst_rel)  # p/(γm) for relativistic

# Velocity difference
Δv = v_rel - v_newton
Δv_magnitude = norm(Δv)

# Use relativistic velocity as reference direction
v_rel_unit = v_rel / norm(v_rel)

# Components relative to velocity direction
# Parallel component (along velocity)
Δq_spacecraft = q_spacecraft_rel - q_spacecraft_newton
Δq_parallel = dot(Δq_spacecraft, v_rel_unit) * v_rel_unit
Δq_parallel_magnitude = abs(dot(Δq_spacecraft, v_rel_unit))

# Transverse component (perpendicular to velocity)
Δq_transverse = Δq_spacecraft - Δq_parallel
Δq_transverse_magnitude = norm(Δq_transverse)

# Convert to AU and km
Δq_spacecraft_AU = norm(Δq_spacecraft) / AU
Δq_parallel_AU = Δq_parallel_magnitude / AU
Δq_transverse_AU = Δq_transverse_magnitude / AU
km = AU / 1.496e8
Δq_spacecraft_km = norm(Δq_spacecraft) / km
Δq_parallel_km = Δq_parallel_magnitude / km
Δq_transverse_km = Δq_transverse_magnitude / km

# Velocity angle difference
cos_angle = dot(v_newton, v_rel) / (norm(v_newton) * norm(v_rel))
cos_angle = clamp(cos_angle, -1.0, 1.0)  # Ensure valid range for acos
angle_diff_rad = acos(cos_angle)
angle_diff_deg = angle_diff_rad * 180.0 / π
angle_diff_arcsec = angle_diff_deg * 3600.0

println("Computed quantities:")
println("Total position difference:")
println("  Δq_spacecraft_AU = $(Δq_spacecraft_AU) AU")
println("  Δq_spacecraft_km = $(Δq_spacecraft_km) km")
println()
println("Components relative to velocity direction:")
println("  Parallel (along velocity):")
println("    Δq_parallel_AU = $(Δq_parallel_AU) AU")
println("    Δq_parallel_km = $(Δq_parallel_km) km")
println("  Transverse (perpendicular to velocity):")
println("    Δq_transverse_AU = $(Δq_transverse_AU) AU")
println("    Δq_transverse_km = $(Δq_transverse_km) km")
println()
println("Velocity comparison:")
println("  |v_newton| = $(norm(v_newton)) c")
println("  |v_rel| = $(norm(v_rel)) c")
println("  Δv magnitude = $(Δv_magnitude) c")
println("  Angle between velocities = $(angle_diff_deg) degrees")
println("  Angle between velocities = $(angle_diff_arcsec) arcseconds")
println()
println("Verification: √(parallel² + transverse²) = $(sqrt(Δq_parallel_AU^2 + Δq_transverse_AU^2)) AU")
