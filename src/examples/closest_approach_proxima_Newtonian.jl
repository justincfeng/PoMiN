#-----------------------------------------------------------------------
#
#   4-BODY CLOSEST APPROACH CALCULATION
#   Spacecraft, Proxima Centauri, Sun, and Jupiter
#   Modernized version compatible with current PoMiN framework
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots
using DoubleFloats
tpfl  = Double64

# using ArbNumerics
# setprecision(ArbFloat, 256)
# tpfl  = ArbFloat

#setprecision(BigFloat, 256)  # Set precision to 256 bits
#tpfl  = BigFloat

include("target_proxima.jl")


#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0

vst = norm(Vst)
γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)

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
pProx = mProx .* γVpx

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
# Jupiter mass relative to Sun: ~0.000954
# Average distance from Sun: ~5.2 AU in geometric units
# Orbital velocity: ~13.1 km/s converted to geometric units
mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([5.2*AU, 0.0, 0.0])  # Jupiter position in geometric units
# Convert orbital velocity to geometric units: v = 13.1 km/s / c
vorb = tpfl(13.1E+3) / cMKS  # Orbital velocity in units of c
γjup = one(tpfl)/sqrt(one(tpfl)-vorb^2)  # Lorentz factor
pjup = tpfl.([0.0, mjup * γjup * vorb, 0.0])  # Jupiter momentum (y-direction)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tcl*tpfl(1.1) ) # 22 years
δ     = tpfl(7.31E+8) # about one hour
tols  = tpfl(1e-16)

#-----------------------------------------------------------------------
#   TEST CASES
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Target position
q_target_pos_Flat = XbF(tcl)

#-----------------------------------------------------------------------
# SUN ONLY

PintSO          = Psol
PtestSO         = pomin.merge_particle_systems(Pchip,PProx)

paramsSO        = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                         atol=tols, rtol=tols)

solSO           = pomin.solve(PintSO, paramsSO; testparticles=PtestSO , 
                              Newtonian = true)

ZendSO          = solSO(tcl)
q_spacecraft_SO = ZendSO[7:9]
q_proxima_SO    = ZendSO[10:12] 
q_target_pos_SO = q_proxima_SO + bv
dist_prox_SO    = norm(q_spacecraft_SO - q_proxima_SO) / AU
miss_SO         = norm(q_spacecraft_SO - q_target_pos_SO) / AU
miss_SO_flat    = norm(q_spacecraft_SO - q_target_pos_Flat) / AU
d_proxima_SO    = norm(q_proxima_SO - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA

PintSP          = pomin.merge_particle_systems(Psol,PProx)
PtestSP         = Pchip

paramsSP        = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSP           = pomin.solve(PintSP, paramsSP; testparticles=PtestSP , 
                              Newtonian = true)

ZendSP          = solSP(tcl)
q_spacecraft_SP = ZendSP[13:15]
q_proxima_SP    = ZendSP[4:6] 

q_target_pos_SP = q_proxima_SP + bv
dist_prox_SP    = norm(q_spacecraft_SP - q_proxima_SP) / AU
miss_SP         = norm(q_spacecraft_SP - q_target_pos_SP) / AU
miss_SP_flat    = norm(q_spacecraft_SP - q_target_pos_Flat) / AU
d_proxima_SP    = norm(q_proxima_SP - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER

PintSPJ         = pomin.merge_particle_systems(PintSP,Pjup)
PtestSPJ        = Pchip

paramsSPJ       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSPJ          = pomin.solve(PintSPJ, paramsSPJ; 
                              testparticles=PtestSPJ, 
                              Newtonian = true)

ZendSPJ         = solSPJ(tcl)
q_spacecraft_SPJ = ZendSPJ[19:21]
q_proxima_SPJ    = ZendSPJ[4:6] 
q_target_pos_SPJ = q_proxima_SPJ + bv
dist_prox_SPJ    = norm(q_spacecraft_SPJ - q_proxima_SPJ) / AU
miss_SPJ         = norm(q_spacecraft_SPJ - q_target_pos_SPJ) / AU
miss_SPJ_flat    = norm(q_spacecraft_SPJ - q_target_pos_Flat) / AU
d_proxima_SPJ    = norm(q_proxima_SPJ - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER (CORRECTED)

## CORRECTION OBTAINED FROM:
#   dxe = q_spacecraft_SPJ - target_pos
#   δV  = (- dxe ./ tcl ) + Vst

# Momentum scaled by same factor as mass to preserve velocity
pchipcorr = (mchip*γst) .* (Vst .+ (δV ./ γst ) )

# Particle object for spacecraft
Pchipcorr = pomin.setup_single_particle(mchip, qchip, pchipcorr, tpfl)

PintSPJc        = pomin.merge_particle_systems(PintSP,Pjup)
PtestSPJc       = Pchipcorr

paramsSPJc      = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSPJc         = pomin.solve(PintSPJc, paramsSPJc; 
                              testparticles=PtestSPJc ,
                              Newtonian = true)

ZendSPJc        = solSPJc(tcl)
q_spacecraft_SPJc = ZendSPJc[19:21]
q_proxima_SPJc    = ZendSPJc[4:6] 
q_target_pos_SPJc = q_proxima_SPJc + bv
dist_prox_SPJc    = norm(q_spacecraft_SPJc - q_proxima_SPJc) / AU
miss_SPJc         = norm(q_spacecraft_SPJc - q_target_pos_SPJc) / AU
miss_SPJc_flat    = norm(q_spacecraft_SPJc - q_target_pos_Flat) / AU
d_proxima_SPJc    = norm(q_proxima_SPJc - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER + MILKY WAY (Sun's Rest Frame)

# Include external potential functionality
include("../core/physics/external_potentials/external.jl")

# Solar offset values in the Milky Way (in solar mass units)
origin_x = tpfl(-1.708859462494220e17)  # Solar x-offset
origin_z = tpfl(4.346342845091530e14)   # Solar z-offset  
origin_y = tpfl(0.0)                    # Solar y-offset
xo_MW = [origin_x, origin_y, origin_z]

km_s = tpfl(1000.0 / 299792458.0)

vpec    = km_s .* tpfl.([11.1, 12.24, 7.25])  # Solar peculiar motion
v_LSR   = tpfl.([0.0, 220.0 * km_s, 0.0])     # LSR circular motion
v_total = v_LSR + vpec                         # Total velocity

xo_MWSRF = xo_MW .- v_total * tcl

# Milky Way potential in the sun's rest frame
Φ_MW = ΦMilkyWay(tpfl, xo_MWSRF)

# Mass array for external potential (includes all massive bodies)
m_ext = [msol, mProx, mjup]  # Sun, Proxima, Jupiter masses

PintSPJMW       = pomin.merge_particle_systems(PintSP, Pjup)
PtestSPJMW      = Pchipcorr  # Use corrected spacecraft

paramsSPJMW     = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSPJMW        = pomin.solve(PintSPJMW, paramsSPJMW; 
                              testparticles=PtestSPJMW, Φ=Φ_MW, 
                              Newtonian = true)

ZendSPJMW       = solSPJMW(tcl)
q_spacecraft_SPJMW = ZendSPJMW[19:21]
q_proxima_SPJMW    = ZendSPJMW[4:6] 
q_target_pos_SPJMW = q_proxima_SPJMW + bv
dist_prox_SPJMW    = norm(q_spacecraft_SPJMW - q_proxima_SPJMW) / AU
miss_SPJMW         = norm(q_spacecraft_SPJMW - q_target_pos_SPJMW) / AU
miss_SPJMW_flat    = norm(q_spacecraft_SPJMW - q_target_pos_Flat) / AU
d_proxima_SPJMW    = norm(q_proxima_SPJMW - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
#
#   ANALYSIS AND RESULTS
#
#-----------------------------------------------------------------------

println("="^70)
println("4-BODY CLOSEST APPROACH CALCULATION RESULTS")
println("="^70)
println()

println("Target Parameters:")
println("  Target distance from Proxima: $(b0/AU) AU")
println("  Time to closest approach: $(tcl) (geometric units)")
println("  Spacecraft velocity: $(vst) c")
println()

println("Miss Distances (AU):")
println("  Sun Only:")
println("    Endpoint distance to Proxima (AU): $(dist_prox_SO)")
println("    Proxima position  (AU): $(dist_prox_SO)")
println("    Miss from target (AU): $(miss_SO)")
println("    Miss from flat space (AU): $(miss_SO_flat)")
println("    Proxima deviation (AU): $(d_proxima_SO)")
println()

println("  Sun + Proxima:")
println("    Endpoint distance to Proxima (AU): $(dist_prox_SP)")
println("    Miss from target (AU): $(miss_SP)")
println("    Miss from flat space (AU): $(miss_SP_flat)")
println("    Proxima deviation (AU): $(d_proxima_SP)")
println()

println("  Sun + Proxima + Jupiter:")
println("    Endpoint distance to Proxima (AU): $(dist_prox_SPJ)")
println("    Miss from target (AU): $(miss_SPJ)")
println("    Miss from flat space (AU): $(miss_SPJ_flat)")
println("    Proxima deviation (AU): $(d_proxima_SPJ)")
println()

println("  Sun + Proxima + Jupiter (Corrected):")
println("    Endpoint distance to Proxima (AU): $(dist_prox_SPJc)")
println("    Miss from target (AU): $(miss_SPJc)")
println("    Miss from flat space (AU): $(miss_SPJc_flat)")
println("    Proxima deviation (AU): $(d_proxima_SPJc)")
println()

println("  Sun + Proxima + Jupiter + Milky Way (Sun's Rest Frame):")
println("    Endpoint distance to Proxima (AU): $(dist_prox_SPJMW)")
println("    Miss from target (AU): $(miss_SPJMW)")
println("    Miss from flat space (AU): $(miss_SPJMW_flat)")
println("    Proxima deviation (AU): $(d_proxima_SPJMW)")
println()

println("Gravitational Effects:")
println("  Sun-only vs Flat space miss difference: $(abs(miss_SO_flat - miss_SO)) AU")
println("  Sun+Proxima vs Sun-only miss difference: $(abs(miss_SP - miss_SO)) AU")
println("  Sun+Proxima+Jupiter vs Sun+Proxima miss difference: $(abs(miss_SPJ - miss_SP)) AU")
println("  Correction effectiveness: $(abs(miss_SPJc - miss_SPJ)) AU improvement")
println("  Milky Way vs Sun+Proxima+Jupiter(Corrected) difference: $(abs(miss_SPJMW - miss_SPJc)) AU")
println()

println("="^70)
