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
# SUN ONLY

PintSO      = Psol
PtestSO     = pomin.merge_particle_systems(Pchip,PProx)

paramsSO    = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSO       = pomin.solveT(PintSO, paramsSO, PtestSO)

ZendSO      = solSO(tcl)

# Debug Z2Part dimensions for Sun-only case
println("=== DEBUG Z2Part for Sun-only case ===")
println("ZendSO length: ", length(ZendSO))
println("PintSO.m length: ", length(PintSO.m))
println("PtestSO.m length: ", length(PtestSO.m))
println("Expected main particle phase space: ", 2 * length(PintSO.m) * 3)
println("Expected test particle phase space: ", 2 * length(PtestSO.m) * 3)
println("Total expected: ", 2 * length(PintSO.m) * 3 + 2 * length(PtestSO.m) * 3)

# Try to extract particles with corrected indexing
try
    PintendSO = Z2Part(ZendSO[1:6*length(PintSO.m)], PintSO.m, length(PintSO.m), 3)
    PtestendSO = Z2Part(ZendSO[6*length(PintSO.m)+1:end], PtestSO.m, length(PtestSO.m), 3)
    println("✓ Z2Part extraction successful for Sun-only case")
catch e
    println("✗ Z2Part extraction failed: ", e)
end

#-----------------------------------------------------------------------
# SUN + PROXIMA

PintSP      = pomin.merge_particle_systems(Psol,PProx)
PtestSP     = Pchip

paramsSP    = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSP       = pomin.solveT(PintSP, paramsSP, PtestSP)

ZendSP      = solSP(tcl)

# Debug Z2Part dimensions for Sun+Proxima case
println("\n=== DEBUG Z2Part for Sun+Proxima case ===")
println("ZendSP length: ", length(ZendSP))
println("PintSP.m length: ", length(PintSP.m))
println("PtestSP.m length: ", length(PtestSP.m))
println("Expected main particle phase space: ", 2 * length(PintSP.m) * 3)
println("Expected test particle phase space: ", 2 * length(PtestSP.m) * 3)
println("Total expected: ", 2 * length(PintSP.m) * 3 + 2 * length(PtestSP.m) * 3)

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER

PintSPJ     = pomin.merge_particle_systems(PintSP,Pjup)
PtestSPJ    = Pchip

paramsSPJ   = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSPJ      = pomin.solveT(PintSPJ, paramsSPJ, PtestSPJ)

ZendSPJ     = solSPJ(tcl)

# Debug Z2Part dimensions for Sun+Proxima+Jupiter case
println("\n=== DEBUG Z2Part for Sun+Proxima+Jupiter case ===")
println("ZendSPJ length: ", length(ZendSPJ))
println("PintSPJ.m length: ", length(PintSPJ.m))
println("PtestSPJ.m length: ", length(PtestSPJ.m))
println("Expected main particle phase space: ", 2 * length(PintSPJ.m) * 3)
println("Expected test particle phase space: ", 2 * length(PtestSPJ.m) * 3)
println("Total expected: ", 2 * length(PintSPJ.m) * 3 + 2 * length(PtestSPJ.m) * 3)

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER (CORRECTED)

PintSPJc     = pomin.merge_particle_systems(PintSP,Pjup)
PtestSPJc    = Pchipcorr

paramsSPJc   = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

solSPJc      = pomin.solveT(PintSPJc, paramsSPJc, PtestSPJc)

ZendSPJc     = solSPJc(tcl)

# Debug Z2Part dimensions for Sun+Proxima+Jupiter case
println("\n=== DEBUG Z2Part for Sun+Proxima+Jupiter case ===")
println("ZendSPJc length: ", length(ZendSPJc))
println("PintSPJc.m length: ", length(PintSPJc.m))
println("PtestSPJc.m length: ", length(PtestSPJc.m))
println("Expected main particle phase space: ", 2 * length(PintSPJc.m) * 3)
println("Expected test particle phase space: ", 2 * length(PtestSPJc.m) * 3)
println("Total expected: ", 2 * length(PintSPJc.m) * 3 + 2 * length(PtestSPJc.m) * 3)

#-----------------------------------------------------------------------
#
#   ANALYSIS AND RESULTS
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("CLOSEST APPROACH ANALYSIS")
println("="^70)

# Function to calculate distance between spacecraft and Proxima
function calc_distance(sol, t_eval, scenario)
    state = sol(t_eval)
    
    if scenario == "SO"  # Sun only - both spacecraft and Proxima are test particles
        q_spacecraft = state[7:9]   # First test particle (spacecraft)
        q_proxima = state[10:12]    # Second test particle (Proxima)
    elseif scenario == "SP"  # Sun + Proxima - spacecraft is test particle
        q_spacecraft = state[13:15]  # Test particle (spacecraft)
        q_proxima = state[4:6]       # Main particle (Proxima)
    elseif scenario == "SPJ"  # Sun + Proxima + Jupiter - spacecraft is test particle
        q_spacecraft = state[19:21]  # Test particle (spacecraft)
        q_proxima = state[4:6]       # Main particle (Proxima)
    end
    
    return norm(q_spacecraft - q_proxima) / AU
end

# Analyze each scenario
println("\n1. SUN ONLY scenario:")
dist_SO = calc_distance(solSO, tcl, "SO")
println("   Closest approach: ", @sprintf("%.6f", dist_SO), " AU")

println("\n2. SUN + PROXIMA scenario:")
dist_SP = calc_distance(solSP, tcl, "SP")
println("   Closest approach: ", @sprintf("%.6f", dist_SP), " AU")

println("\n3. SUN + PROXIMA + JUPITER scenario:")
dist_SPJ = calc_distance(solSPJ, tcl, "SPJ")
println("   Closest approach: ", @sprintf("%.6f", dist_SPJ), " AU")

println("\n4. SUN + PROXIMA + JUPITER (CORRECTED) scenario:")
dist_SPJc = calc_distance(solSPJc, tcl, "SPJ")  # Same indexing as SPJ
println("   Closest approach: ", @sprintf("%.6f", dist_SPJc), " AU")

# Compare with flat-space targeting solution
target_pos = XbF(tcl)
target_dist = norm(target_pos - XpxF(tcl)) / AU


# Extract spacecraft positions from solutions
state_SO = solSO(tcl)
state_SP = solSP(tcl)
state_SPJ = solSPJ(tcl)
state_SPJc = solSPJc(tcl)

q_spacecraft_SO = state_SO[7:9]
q_spacecraft_SP = state_SP[13:15] 
q_spacecraft_SPJ = state_SPJ[19:21]
q_spacecraft_SPJc = state_SPJc[19:21]

# Extract Proxima positions first
q_proxima_SO = state_SO[10:12]  # Proxima as test particle in Sun-only case
q_proxima_SP = state_SP[4:6]    # Proxima as main particle in Sun+Proxima case
q_proxima_SPJ = state_SPJ[4:6]  # Proxima as main particle in Sun+Proxima+Jupiter case
q_proxima_SPJc = state_SPJc[4:6]  # Proxima as main particle in Sun+Proxima+Jupiter (corrected) case

# Compute actual target positions by adding bv to Proxima's end positions
target_pos_SO = q_proxima_SO + bv
target_pos_SP = q_proxima_SP + bv
target_pos_SPJ = q_proxima_SPJ + bv
target_pos_SPJc = q_proxima_SPJc + bv

# Calculate miss distances - spacecraft end position vs flat-space target position
miss_SO_flat = norm(q_spacecraft_SO - target_pos) / AU
miss_SP_flat = norm(q_spacecraft_SP - target_pos) / AU
miss_SPJ_flat = norm(q_spacecraft_SPJ - target_pos) / AU
miss_SPJc_flat = norm(q_spacecraft_SPJc - target_pos) / AU

# Calculate miss distances - spacecraft end position vs actual target positions (Proxima + bv)
miss_SO_actual = norm(q_spacecraft_SO - target_pos_SO) / AU
miss_SP_actual = norm(q_spacecraft_SP - target_pos_SP) / AU
miss_SPJ_actual = norm(q_spacecraft_SPJ - target_pos_SPJ) / AU
miss_SPJc_actual = norm(q_spacecraft_SPJc - target_pos_SPJc) / AU

# Summary comparison
println("\n" * "-"^50)
println("COMPARISON SUMMARY:")
println("-"^50)
println("Target distance:     ", @sprintf("%.6f", b0/AU), " AU")
println("Flat-space target:   ", @sprintf("%.6f", target_dist), " AU")
println("Sun only:           ", @sprintf("%.6f", dist_SO), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", dist_SP), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", dist_SPJ), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", dist_SPJc), " AU")

println("\nMiss distances (spacecraft vs flat-space target):")
println("Sun only:           ", @sprintf("%.6f", miss_SO_flat), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", miss_SP_flat), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", miss_SPJ_flat), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", miss_SPJc_flat), " AU")

println("\nMiss distances (spacecraft vs actual target = Proxima + bv):")
println("Sun only:           ", @sprintf("%.6f", miss_SO_actual), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", miss_SP_actual), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", miss_SPJ_actual), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", miss_SPJc_actual), " AU")

# Additional validation
println("\nFlat-space validation:")
println("Spacecraft at tcl:   ", XstF(tcl) / AU, " AU")
println("Proxima at tcl:      ", XpxF(tcl) / AU, " AU") 
println("Target point at tcl: ", target_pos / AU, " AU")

# Compare actual spacecraft positions with flat-space prediction
println("\nSpacecraft position comparison:")

println("Flat-space prediction: ", XstF(tcl) / AU, " AU")
println("Sun only result:       ", q_spacecraft_SO / AU, " AU")
println("Sun+Proxima result:    ", q_spacecraft_SP / AU, " AU")
println("Sun+Proxima+Jup result:", q_spacecraft_SPJ / AU, " AU")
println("Sun+Prox+Jup (C) result:", q_spacecraft_SPJc / AU, " AU")

# Distance from flat-space prediction
dev_SO = norm(q_spacecraft_SO - XstF(tcl)) / AU
dev_SP = norm(q_spacecraft_SP - XstF(tcl)) / AU
dev_SPJ = norm(q_spacecraft_SPJ - XstF(tcl)) / AU
dev_SPJc = norm(q_spacecraft_SPJc - XstF(tcl)) / AU

println("\nDeviation from flat-space prediction:")
println("Sun only:           ", @sprintf("%.6f", dev_SO), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", dev_SP), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", dev_SPJ), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", dev_SPJc), " AU")

# Compare Proxima's final position with XpxF(tcl)
println("\nProxima position comparison:")

println("Flat-space prediction: ", XpxF(tcl) / AU, " AU")
println("Sun only result:       ", q_proxima_SO / AU, " AU")
println("Sun+Proxima result:    ", q_proxima_SP / AU, " AU")
println("Sun+Proxima+Jup result:", q_proxima_SPJ / AU, " AU")
println("Sun+Prox+Jup (C) result:", q_proxima_SPJc / AU, " AU")

# Proxima deviation from flat-space prediction
prox_dev_SO = norm(q_proxima_SO - XpxF(tcl)) / AU
prox_dev_SP = norm(q_proxima_SP - XpxF(tcl)) / AU
prox_dev_SPJ = norm(q_proxima_SPJ - XpxF(tcl)) / AU
prox_dev_SPJc = norm(q_proxima_SPJc - XpxF(tcl)) / AU

println("\nProxima deviation from flat-space prediction:")
println("Sun only:           ", @sprintf("%.6f", prox_dev_SO), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", prox_dev_SP), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", prox_dev_SPJ), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", prox_dev_SPJc), " AU")

# Additional analysis: Target position comparison
println("\nTarget position comparison:")
println("Flat-space target:      ", target_pos / AU, " AU")
println("Sun only target:        ", target_pos_SO / AU, " AU")
println("Sun+Proxima target:     ", target_pos_SP / AU, " AU")
println("Sun+Proxima+Jup target: ", target_pos_SPJ / AU, " AU")
println("Sun+Prox+Jup (C) target:", target_pos_SPJc / AU, " AU")

# Target position deviations from flat-space
target_dev_SO = norm(target_pos_SO - target_pos) / AU
target_dev_SP = norm(target_pos_SP - target_pos) / AU
target_dev_SPJ = norm(target_pos_SPJ - target_pos) / AU
target_dev_SPJc = norm(target_pos_SPJc - target_pos) / AU

println("\nTarget position deviation from flat-space:")
println("Sun only:           ", @sprintf("%.6f", target_dev_SO), " AU")
println("Sun + Proxima:      ", @sprintf("%.6f", target_dev_SP), " AU")
println("Sun + Proxima + Jup:", @sprintf("%.6f", target_dev_SPJ), " AU")
println("Sun + Prox + Jup (C):", @sprintf("%.6f", target_dev_SPJc), " AU")

println("\n" * "="^70)

## CORRECTION OBTAINED FROM:
#   dxe = q_spacecraft_SPJ - target_pos
#   δV  = (- dxe ./ tcl ) + Vst
