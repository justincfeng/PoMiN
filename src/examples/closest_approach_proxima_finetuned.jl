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
using DoubleFloats, ForwardDiff
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

# Momentum scaled by same factor as mass to preserve velocity
pchip = (mchip*γst) .* (Vst)

# Particle object for spacecraft
Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)

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

solSO           = pomin.solve(PintSO, paramsSO; testparticles=PtestSO)

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

solSP           = pomin.solve(PintSP, paramsSP; testparticles=PtestSP)

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

solSPJ          = pomin.solve(PintSPJ, paramsSPJ; testparticles=PtestSPJ)

ZendSPJ         = solSPJ(tcl)
q_spacecraft_SPJ = ZendSPJ[19:21]
q_proxima_SPJ    = ZendSPJ[4:6] 
q_target_pos_SPJ = q_proxima_SPJ + bv
dist_prox_SPJ    = norm(q_spacecraft_SPJ - q_proxima_SPJ) / AU
miss_SPJ         = norm(q_spacecraft_SPJ - q_target_pos_SPJ) / AU
miss_SPJ_flat    = norm(q_spacecraft_SPJ - q_target_pos_Flat) / AU
d_proxima_SPJ    = norm(q_proxima_SPJ - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER (CRUDE CORRECTED)

## Corrected velocity
dxe     = q_spacecraft_SPJ - q_target_pos_SPJ
δV      = - dxe ./ tcl
vcorr   = (Vst .+ (δV ./ γst ) )
γstcorr = one(tpfl)/sqrt(one(tpfl)-(norm(vcorr)/c)^2)

# Corrected momentum
pchipcorr = (mchip*γstcorr) .* (vcorr )

# Particle object for spacecraft
Pchipcorr = pomin.setup_single_particle(mchip, qchip, pchipcorr, tpfl)

PintSPJc        = pomin.merge_particle_systems(PintSP,Pjup)
PtestSPJc       = Pchipcorr

paramsSPJc      = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                    atol=tpfl(1e-24), rtol=tpfl(1e-24))

solSPJc         = pomin.solve(PintSPJc, paramsSPJc; testparticles=PtestSPJc)

ZendSPJc        = solSPJc(tcl)
q_spacecraft_SPJc = ZendSPJc[19:21]
q_proxima_SPJc    = ZendSPJc[4:6] 
q_target_pos_SPJc = q_proxima_SPJc + bv
dist_prox_SPJc    = norm(q_spacecraft_SPJc - q_proxima_SPJc) / AU
miss_SPJc         = norm(q_spacecraft_SPJc - q_target_pos_SPJc) / AU
miss_SPJc_flat    = norm(q_spacecraft_SPJc - q_target_pos_Flat) / AU
d_proxima_SPJc    = norm(q_proxima_SPJc - XpxF(tcl)) / AU

#-----------------------------------------------------------------------
# SUN + PROXIMA + JUPITER (AUTODIFF CORRECTED)

function targetrootconstructor(Pint,Ptest0,params,bv,tcl)
    return function f(v)
        # Create a modified particle with new velocity
        # Calculate Lorentz factor for the new velocity
        v_mag = norm(v)
        γ = one(eltype(v)) / sqrt(one(eltype(v)) - (v_mag/c)^2)
        
        # Calculate new momentum: p = γmv
        pv = (Ptest0.m[1] * γ) .* v
        
        # Create new particle with modified momentum
        Ptest = pomin.setup_single_particle(Ptest0.m[1], Ptest0.q[1], pv, eltype(Ptest0.m))

        # Solve the system
        sol_modified = pomin.solve(Pint, params; testparticles=Ptest)
        
        # Calculate miss distance at final time
        Zend_modified = sol_modified(tcl)
        # Extract spacecraft and target positions (indices depend on system configuration)
        # For SPJc: PintSPJc has Sun(1-6) + Proxima(7-12) + Jupiter(13-18), testparticle is at (19-24)
        q_spacecraft = Zend_modified[19:21]
        q_proxima = Zend_modified[4:6]  # Proxima position (same as original)
        q_target = q_proxima + bv  # Target position
        
        return q_spacecraft - q_target
    end
end

FSPJc = targetrootconstructor(PintSPJc,PtestSPJc,paramsSPJc,bv,tcl)

Fcheck = FSPJc(vcorr) - (q_spacecraft_SPJc - q_target_pos_SPJc)
println("Fcheck = $Fcheck")
println("Check miss distance: $(norm(Fcheck))")

#-----------------------------------------------------------------------
# BROYDEN OPTIMIZATION
#-----------------------------------------------------------------------

# Include Broyden algorithm
include("../utils/broyden.jl")

println("\n=== BROYDEN OPTIMIZATION ===")

println("Computing Jacobian")
J_init = ForwardDiff.jacobian(FSPJc, vcorr)

println("Broyden iterations")
v_finetuned = bsolve(FSPJc, J_init, FSPJc(vcorr), vcorr, 10)  # Maximum 10 iterations

# Check final result
f_final = FSPJc(v_finetuned)
println("\n=== BROYDEN RESULTS ===")
println("Optimized velocity = $v_finetuned")
println("Final function value = $f_final")
println("Final norm |f| = $(norm(f_final))")
println("Velocity correction = $(v_finetuned - vcorr)")
println("Improvement factor = $(norm(FSPJc(vcorr)) / norm(f_final))")
println("=== END BROYDEN ===\n")

v_magft = norm(v_finetuned)
γft = one(eltype(v_finetuned)) / sqrt(one(eltype(v_finetuned)) - (v_magft/c)^2)
pvft = (PtestSPJc.m[1] * γft) .* v_finetuned

# Create new particle with modified momentum
Ptestft = pomin.setup_single_particle(PtestSPJc.m[1], PtestSPJc.q[1], pvft, eltype(PtestSPJc.m))

# Solve the system
sol_ft = pomin.solve(PintSPJc, paramsSPJc; testparticles=Ptestft)

# Calculate miss distance at final time
Zend_ft = sol_ft(tcl)
# Extract spacecraft and target positions (indices depend on system configuration)
# For SPJc: PintSPJc has Sun(1-6) + Proxima(7-12) + Jupiter(13-18), testparticle is at (19-24)
q_spacecraft_ft = Zend_ft[19:21]
q_proxima_ft = Zend_ft[4:6]  # Proxima position (same as original)
q_target_ft = q_proxima_ft + bv  # Target position

println("Fine tuned speed = $(v_magft)")
println("Fine tuned final spacecraft position = $(q_spacecraft_ft)")
println("Fine tuned final target position = $(q_target_ft)")