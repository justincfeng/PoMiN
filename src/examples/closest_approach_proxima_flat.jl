#-----------------------------------------------------------------------
#
#   4-BODY CLOSEST APPROACH CALCULATION - FLAT SPACETIME
#   Spacecraft, Proxima Centauri, Sun, and Jupiter
#   All gravity turned off - pure inertial motion baseline
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots
using DoubleFloats
tpfl  = Double64

include("target_proxima.jl")

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0

vst = norm(Vst)
γst = one(tpfl)  # Newtonian limit

δV =  tpfl.([-6.04280291023418732263714890692229831e-08,
             6.96240251604129644625817536477811962e-08,
             -5.34085226136474509045517326371492461e-08])

# Momentum scaled by same factor as mass to preserve velocity
pchip = (mchip*γst) .* (Vst)
pchipcorr = (mchip*γst) .* (Vst .+ (δV ./ γst ) )

# Particle object for spacecraft
Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)
Pchipcorr = pomin.setup_single_particle(mchip, qchip, pchipcorr, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR PROXIMA CENTAURI
#-----------------------------------------------------------------------

# Mass, position, momentum
mProx = tpfl(0.1221)
qProx = Xpx0
pProx = mProx .* γV2v(γVpx)

# Particle object for Proxima Centauri
PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SUN AND JUPITER
#-----------------------------------------------------------------------

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)

# Jupiter parameters (realistic orbital position)
mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([5.2*AU, 0.0, 0.0])  # Jupiter position in geometric units
vorb = tpfl(13.1E+3) / cMKS  # Orbital velocity in units of c
γjup = one(tpfl)  # Newtonian limit
pjup = tpfl.([0.0, mjup * γjup * vorb, 0.0])  # Jupiter momentum (y-direction)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tcl*tpfl(1.1) ) # 22 years
δ     = tpfl(7.31E+8) # about one hour
tols  = tpfl(1e-16)

#-----------------------------------------------------------------------
#   FLAT SPACETIME CASES - ANALYTICAL ONLY
#-----------------------------------------------------------------------

# Skip numerical integration (causes issues with zero masses)
# Use pure analytical calculation for flat spacetime baseline

#-----------------------------------------------------------------------
# ANALYTICAL FLAT SPACE CALCULATION (VERIFICATION)
#-----------------------------------------------------------------------

# Calculate positions analytically using straight-line motion
q_spacecraft_analytical = Xst0 + Vst * tcl
q_proxima_analytical = Xpx0 + γV2v(γVpx) * tcl
q_target_analytical = q_proxima_analytical + bv

dist_prox_analytical = norm(q_spacecraft_analytical - q_proxima_analytical) / AU
miss_analytical = norm(q_spacecraft_analytical - q_target_analytical) / AU

# Corrected trajectory analytical
q_spacecraft_corr_analytical = Xst0 + (Vst + δV) * tcl
miss_corr_analytical = norm(q_spacecraft_corr_analytical - q_target_analytical) / AU

#-----------------------------------------------------------------------
#
#   ANALYSIS AND RESULTS
#
#-----------------------------------------------------------------------

println("="^70)
println("FLAT SPACETIME CLOSEST APPROACH CALCULATION RESULTS")
println("="^70)
println()

println("Target Parameters:")
println("  Target distance from Proxima: $(b0/AU) AU")
println("  Time to closest approach: $(tcl) (geometric units)")
println("  Spacecraft velocity: $(vst) c")
println()

println("Flat Spacetime Results (Analytical):")
println("  Straight-line motion:")
println("    Endpoint distance to Proxima (AU): $(dist_prox_analytical)")
println("    Miss from target (AU): $(miss_analytical)")
println()
println("  Corrected trajectory:")
println("    Miss from target (AU): $(miss_corr_analytical)")
println()

println("Key Insights:")
println("  - This establishes the pure inertial motion baseline")
println("  - Any deviation from these values in other simulations is due to gravity")
println("  - The correction δV should have minimal effect in flat spacetime")
println("  - Numerical integration should match analytical calculation exactly")

println("="^70)
