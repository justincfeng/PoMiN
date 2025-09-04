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

# Import Z2q function from HamTools
using .pomin: Z2q
using .pomin: Z2Part

#-----------------------------------------------------------------------
#
#   FUNCTIONS FOR TARGETING PROBLEM
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function γV2v(γV)
    #   γ  = sqrt(1+γv^2)
    return γV ./ sqrt(one(typeof(γV[1]))+dot(γV,γV))
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function vunit(v)
    return v/norm(v)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function uvproj(u,v)
    U = vunit(u)
    return v .- dot(U,v) .* U
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function basisconstructor(vp,vs)
    epar = vunit(vp)
    e1   = vunit(uvproj(vp,vs))
    e2   = vunit(cross(e1,epar))
    return (epar,e1,e2)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function bconstructor(vecs,b,Θ)
    Xpx0,Xst0,Vpx = vecs

    tpfl = typeof(b)

    ONE,TWO  = one(tpfl),tpfl(2)

    ΔXp0  = Xst0 .- Xpx0

    (epar,e1,e2) = basisconstructor(-ΔXp0,Vpx)

    n = cos(Θ)*e1 + sin(Θ)*e2

    φ = acos(b/norm(ΔXp0))

    return b .* (cos(φ) .* epar + sin(φ) .* n)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function targprox(vecs,vst,bv)
    Xpx0,Xst0,Vpx = vecs

    tpfl = typeof(vst)

    ONE,TWO  = one(tpfl),tpfl(2)

    ΔX0  = Xst0 .- (Xpx0 + bv)
    b    = norm(bv)
    rpx  = norm(ΔX0)
    vb   = norm(Vpx)
    
    # From LaTeX: cos(ψ) = -dot(Vpx,ΔX0)/(vb*rpx)
    cosψ = -dot(Vpx,ΔX0) / (vb*rpx)
    sinψ = sqrt(ONE - cosψ^2)
    
    # From LaTeX: sin(ϑ) = (vb/vst) * sin(ψ)
    sinϑ = (vb/vst) * sinψ
    cosϑ = sqrt(ONE - sinϑ^2)
    
    # Construct basis with ΔX0 (not -ΔX0)
    epar = vunit(ΔX0)
    e1   = vunit(uvproj(ΔX0, Vpx))
    
    # Spacecraft velocity: pointing toward target with correct magnitude
    Vst  = vst * (-cosϑ .* epar .+ sinϑ .* e1)
    
    # Velocity difference
    ΔV   = Vst .- Vpx
    
    # Time to closest approach from LaTeX constraint
    tcl  = -dot(ΔX0, ΔV) / dot(ΔV, ΔV)

    return (Xpx0,Xpx0 + bv,Xst0,Vpx,Vst,tcl)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function parfuncs(ics)
    Xpx0,Xpxb,Xst0,Vpx,Vst,tcl = ics
    return ( t->(Xst0 .+ (Vst .* t)) , t->(Xpx0 .+ (Vpx .* t)) , 
             t->(Xpxb .+ (Vpx .* t)) , tcl )
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   PARAMETER DEFINITIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CONSTANTS
#-----------------------------------------------------------------------
cMKS   = tpfl(299792458.0)
AUMKS  = tpfl(149597870700.0)
lyMKS  = tpfl(9460730472580800.0)
Msol2m = tpfl(1476.67)

c      = one(tpfl)
AU     = AUMKS/Msol2m
ly     = lyMKS/Msol2m 

#-----------------------------------------------------------------------
#   FLAT SPACE TARGETING PROBLEM
#-----------------------------------------------------------------------

# Proxima Centauri
Xpx0 = tpfl.([-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13])
γVpx = tpfl.([-3.34779933990201E-5, 7.24198145771899E-5, 6.82553747232694E-5])
Vpx  = γV2v(γVpx)

# Starchip
Xst0   = tpfl.([-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07])
vst    = tpfl(0.2)*c         # Starchip velocity magnitude

# Target
b0     = 0.05*AU               # Target distance
bv     = bconstructor((Xpx0,Xst0,Vpx),b0,pi/3)

ics    = targprox( (Xpx0,Xst0,Vpx) , vst , bv )
pfs    = parfuncs(ics)

XstF,XpxF,XbF,tcl = pfs

Xpx0,Xpxb,Xst0,Vpx,Vst,tcl = ics

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

solSO           = pomin.solveT(PintSO, paramsSO, PtestSO)

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

solSP           = pomin.solveT(PintSP, paramsSP, PtestSP)

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

solSPJ          = pomin.solveT(PintSPJ, paramsSPJ, PtestSPJ)

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

solSPJc         = pomin.solveT(PintSPJc, paramsSPJc, PtestSPJc)

ZendSPJc        = solSPJc(tcl)
q_spacecraft_SPJc = ZendSPJc[19:21]
q_proxima_SPJc    = ZendSPJc[4:6] 
q_target_pos_SPJc = q_proxima_SPJc + bv
dist_prox_SPJc    = norm(q_spacecraft_SPJc - q_proxima_SPJc) / AU
miss_SPJc         = norm(q_spacecraft_SPJc - q_target_pos_SPJc) / AU
miss_SPJc_flat    = norm(q_spacecraft_SPJc - q_target_pos_Flat) / AU
d_proxima_SPJc    = norm(q_proxima_SPJc - XpxF(tcl)) / AU

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

println("Gravitational Effects:")
println("  Sun-only vs Flat space miss difference: $(abs(miss_SO_flat - miss_SO)) AU")
println("  Sun+Proxima vs Sun-only miss difference: $(abs(miss_SP - miss_SO)) AU")
println("  Sun+Proxima+Jupiter vs Sun+Proxima miss difference: $(abs(miss_SPJ - miss_SP)) AU")
println("  Correction effectiveness: $(abs(miss_SPJc - miss_SPJ)) AU improvement")
println()

println("="^70)
