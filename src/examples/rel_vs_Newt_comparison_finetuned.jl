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
using DoubleFloats, ForwardDiff, Dates
tpfl  = Double64

include("target_proxima.jl")
include("../utils/broyden.jl")

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-14)
qchip = Xst0
vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR PROXIMA CENTAURI
#-----------------------------------------------------------------------

# Mass, position, momentum
mProx   = tpfl(0.1221)
qProx   = Xpx0
pProx   = mProx .* γVpx
pProxN  = mProx .* Vpx

# Particle object for Proxima Centauri
PProx   = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)
PProxN  = pomin.setup_single_particle(mProx, qProx, pProxN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SUN
#-----------------------------------------------------------------------

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR ALPHA CENTAURI A + B
#-----------------------------------------------------------------------

malpha = tpfl(2.0429)

# Alpha A+B barycenter position and velocity given by Kervella et al 2017
qalpha = tpfl.([-1.045245216607860E+13, 
                -8.74487747435090E+12, 
                -2.442178248634550E+13])         # in geometric units

valpha = tpfl.([-3.11496332572849e-05,
                 7.38228261899770e-05,
                 7.22023391262230e-05])          # in geometric units

γalpha = one(tpfl) / sqrt(one(tpfl) - norm(valpha)^2)  # Lorentz factor
palpha = tpfl.(malpha * γalpha .* valpha) # Alpha momentum (relativ.)
palphaN = malpha .* valpha                # Alpha momentum (Newtonian)

Palpha = pomin.setup_single_particle(malpha, qalpha, palpha, tpfl)
PalphaN = pomin.setup_single_particle(malpha, qalpha, palphaN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR EARTH
#-----------------------------------------------------------------------

mEarth = tpfl(3.0033693739E-06)

# Earth coordinates generated using astropy with epoch J2000
qEarth = tpfl.([-1.8667826140E+07, 
                 8.9633743993E+07, 
                 3.8883291714E+07])         # geometric units
vEarth = tpfl.([-9.9351889046e-05,
                -1.6777453395e-05,
                -7.2738469450e-06])

γEarth = one(tpfl) / sqrt(one(tpfl) - norm(vEarth)^2)  # Lorentz factor
pEarth = tpfl.(mEarth * γEarth .* vEarth) # Earth momentum (relativ.)
pEarthN = mEarth .* vEarth                # Earth momentum (Newtonian)

PEarth = pomin.setup_single_particle(mEarth, qEarth, pEarth, tpfl)
PEarthN = pomin.setup_single_particle(mEarth, qEarth, pEarthN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR JUPITER (astropy)
#-----------------------------------------------------------------------

mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([3.99442362 * AU, 
              2.73345644 * AU, 
              1.07451704 * AU])  # Jupiter position in geometric units (from astropy, J2000)

vjup = tpfl.([-7.8875477321829800, 
              1.0175858402135900E+01, 
              4.5538873116399600])  # in km/s, from astropy, J2000

vjup *= 1000 / cMKS  # convert from km/s to units of c

γjup = one(tpfl) / sqrt(one(tpfl) - norm(vjup)^2)  # Lorentz factor
pjup = tpfl.(mjup * γjup .* vjup)  # Jupiter momentum 
pjupN = mjup .* vjup                # Jupiter momentum (Newtonian)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)
PjupN = pomin.setup_single_particle(mjup, qjup, pjupN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MOON
#-----------------------------------------------------------------------

mMoon = 7.34767309E22 / 1.988416E30     # in units of solar masses

# Moon coordinates generated using astropy with epoch J2000
qMoon = tpfl.([-1.886529874442160E+07, 
                8.94531273942814E+07, 
                3.883175807801640E+07])    # geometric units
vMoon = tpfl.([-2.9141440121218000E+01, 
                -5.6958177245812600, 
                -2.4819741171379700]) 

vMoon *= 1000 / cMKS  # convert from km/s to units of c

γMoon = one(tpfl) / sqrt(one(tpfl) - norm(vMoon)^2)  # Lorentz factor
pMoon = tpfl.(mMoon * γMoon * vMoon)  # Moon momentum 
pMoonN = mMoon .* vMoon                # Moon momentum (Newtonian)

PMoon = pomin.setup_single_particle(mMoon, qMoon, pMoon, tpfl)
PMoonN = pomin.setup_single_particle(mMoon, qMoon, pMoonN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MARS
#-----------------------------------------------------------------------

mMars = tpfl(3.2271058587E-07)  # Mars mass in solar masses
qMars = tpfl.([1.38356874 * AU, 
               -0.00120916 * AU, 
               -0.03786079 * AU])  # Mars position in geometric units (from astropy, J2000)

vMars = tpfl.([1.17347755650600, 
               2.390740193498500E+01, 
               1.0934201867474400E+01])  # in km/s, from astropy, J2000

vMars *= 1000 / cMKS  # convert from km/s to units of c

γMars = one(tpfl) / sqrt(one(tpfl) - norm(vMars)^2)  # Lorentz factor
pMars = tpfl.(mMars * γMars * vMars)  # Mars momentum 
pMarsN = mMars .* vMars                # Mars momentum (Newtonian)

PMars = pomin.setup_single_particle(mMars, qMars, pMars, tpfl)
PMarsN = pomin.setup_single_particle(mMars, qMars, pMarsN, tpfl)

#-----------------------------------------------------------------------
#
#   INITIAL RUN
#
#-----------------------------------------------------------------------

# Integration parameters
tspan = (tpfl(0), tcl*tpfl(1.1))
tols  = tpfl(1e-17)

params = pomin.ParametersJulia(tspan, integrator="Vern7", 
                                     atol=tols, rtol=tols)

Ptest = Pchip_rel
PtestN = Pchip_newt

Pint = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup)
PintN = pomin.merge_particle_systems(PProxN, Psol, PalphaN, PjupN)
#Pint = pomin.merge_particle_systems(PProx, Psol, Pjup, Palpha, PEarth, PMoon, PMars)
#PintN = pomin.merge_particle_systems(PProxN, Psol, PjupN, PalphaN, PEarthN, PMoonN, PMarsN)

nInt = length(Pint.m)

#-----------------------------------------------------------------------
#   NEWTONIAN CASE
#-----------------------------------------------------------------------

# Solve
solN        = pomin.solve(PintN, params; testparticles=PtestN, Newtonian=true)

# Final positions and momenta
zendN       = solN(tcl)

# Masses
mProxN      = PintN.m[1]
mScN        = PtestN.m[1]

# Separate interacting and test particle phase space vectors
zendNint    = zendN[1:6*nInt]
zendNtest   = zendN[6*nInt+1:end]

# Proxima
qProxN      = zendNint[1:3]
pProxN      = zendNint[3*nInt+1:3*nInt+3]
vProxN      = pProxN ./ mProxN

# Target
qtarN       = qProxN + bv

# Starchip
qScN        = zendNtest[1:3]
pScN        = zendNtest[4:6]
vScN        = pScN ./ mScN

# Miss distance
qmissN      = qtarN - qScN
dmissN      = norm(qmissN)

# Velocity correction
δVN = qmissN ./ tcl
vcorrN = Vst .+ δVN 

#-----------------------------------------------------------------------
#   RELATIVISTIC CASE
#-----------------------------------------------------------------------

# Solve
solR = pomin.solve(Pint, params; testparticles=Ptest)

# Final positions and momenta
zendR = solR(tcl)

# Masses
mProxR      = Pint.m[1]
mScR        = Ptest.m[1]

# Separate interacting and test particle phase space vectors
zendRint    = zendR[1:6*nInt]
zendRtest   = zendR[6*nInt+1:end]

# Proxima
qProxR      = zendRint[1:3]
pProxR      = zendRint[3*nInt+1:3*nInt+3]
vProxR      = pProxR ./ mProxR

# Target
qtarR       = qProxR + bv

# Starchip
qScR        = zendRtest[1:3]
pScR        = zendRtest[4:6]
vScR        = pScR ./ mScR

# Miss distance
qmissR      = qtarR - qScR
dmissR      = norm(qmissR)

# Velocity correction
δVR = qmissR ./ tcl
vcorrR = Vst .+ δVR

#-----------------------------------------------------------------------
#
#   ROOT CONSTRUCTOR
#
#-----------------------------------------------------------------------

function missfuncconstructor(Pint,Ptest0,params,bv,tcl;Newt=false)
    nInt = length(Pint.m)
    if Newt
        return function(v)
            v_mag = norm(v)
            pv = (Ptest0.m[1]) .* v
            Ptest = pomin.setup_single_particle(Ptest0.m[1],
                                                Ptest0.q[1],pv,
                                                eltype(Ptest0.m))
            sol_modified  = pomin.solve(Pint, params; 
                                        testparticles=Ptest,
                                        Newtonian=true)
            Zend_modified = sol_modified(tcl)
            q_proxima     = Zend_modified[1:3]
            q_spacecraft  = Zend_modified[6*nInt+1:6*nInt+3]
            q_target      = q_proxima + bv
            return q_target - q_spacecraft
        end
    else    
        return function(v)
            v_mag = norm(v)
            γ = one(eltype(v)) / sqrt(one(eltype(v)) - (v_mag/c)^2)
            pv = (Ptest0.m[1] * γ) .* v
            Ptest = pomin.setup_single_particle(Ptest0.m[1], 
                                                Ptest0.q[1], pv, 
                                                eltype(Ptest0.m))
            sol_modified  = pomin.solve(Pint, params;
                                        testparticles=Ptest)
            Zend_modified = sol_modified(tcl)
            q_proxima     = Zend_modified[1:3]
            q_spacecraft  = Zend_modified[6*nInt+1:6*nInt+3]
            q_target      = q_proxima + bv
            return q_target - q_spacecraft
        end
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   TEST MISS FUNCTION CONSTRUCTOR
#-----------------------------------------------------------------------

fmN = missfuncconstructor(PintN,PtestN,params,bv,tcl;Newt=true)
fmR = missfuncconstructor(Pint,Ptest,params,bv,tcl;Newt=false)

fmNVst = fmN(Vst)
fmRVst = fmR(Vst)

fmNcheck = fmNVst - qmissN
fmRcheck = fmRVst - qmissR

println("Testing miss function consistency:")
println("fmN(Vst) - qmissN = ", fmNcheck)
println("fmN(Vst) = ", fmNVst)
println("qmissN = ", qmissN)
println()
println("fmR(Vst) - qmissR = ", fmRcheck)
println("fmR(Vst) = ", fmRVst)
println("qmissR = ", qmissR)

#-----------------------------------------------------------------------
#
#   FINE TUNING
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   NEWTONIAN FINE TUNING
#-----------------------------------------------------------------------

println("="^70)
println("NEWTONIAN FINE TUNING")
println("="^70)

println("Initial miss distance: $(dmissN/AU) AU")
println("vcorrN = ", vcorrN)
println("norm(vcorrN) = ", norm(vcorrN))
println("Vst = ", Vst) 
println("norm(Vst) = ", norm(Vst))
println("Computing Jacobian (Newtonian)...")
println("Using Vst as starting point instead of vcorrN")
J_init_N = ForwardDiff.jacobian(fmN, Vst)

println("Jacobian matrix:")
println(J_init_N)
println("Jacobian condition number: ", cond(J_init_N))
println("Jacobian determinant: ", det(J_init_N))

if abs(det(J_init_N)) < 1e-10
    println("WARNING: Jacobian is nearly singular, using pseudoinverse")
    J_init_N_inv = pinv(J_init_N)
    println("Using modified Broyden with pseudoinverse...")
    # Simple Newton step instead of full Broyden
    v_finetuned_N = Vst - J_init_N_inv * fmN(Vst)
else
    println("Broyden iterations (Newtonian)...")
    v_finetuned_N = bsolve(fmN, J_init_N, fmN(Vst), Vst, 10)
end

# Final Newtonian result
f_final_N = fmN(v_finetuned_N)
miss_final_N = norm(f_final_N) / AU

println("Fine-tuned miss distance (Newtonian): $(miss_final_N) AU")
println("Improvement factor (Newtonian): $(dmissN / norm(f_final_N))")

# Compute fine-tuned trajectory for Newtonian case
pvft_N = mchip .* v_finetuned_N  # Newtonian momentum
Ptestft_N = pomin.setup_single_particle(mchip, qchip, pvft_N, tpfl)
sol_ft_N = pomin.solve(PintN, params; testparticles=Ptestft_N, Newtonian=true)
Zend_ft_N = sol_ft_N(tcl)

# Extract final positions
zendNint_ft = Zend_ft_N[1:6*nInt]
zendNtest_ft = Zend_ft_N[6*nInt+1:end]
q_proxima_ft_N = zendNint_ft[1:3]
q_spacecraft_ft_N = zendNtest_ft[1:3]

#-----------------------------------------------------------------------
#   RELATIVISTIC FINE TUNING
#-----------------------------------------------------------------------

println("\n" * "="^70)
println("RELATIVISTIC FINE TUNING")
println("="^70)

println("Initial miss distance: $(dmissR/AU) AU")
println("Computing Jacobian (Relativistic)...")
println("Using Vst as starting point instead of vcorrR")
J_init_R = ForwardDiff.jacobian(fmR, Vst)

println("Jacobian matrix:")
println(J_init_R)
println("Jacobian condition number: ", cond(J_init_R))
println("Jacobian determinant: ", det(J_init_R))

if abs(det(J_init_R)) < 1e-10
    println("WARNING: Jacobian is nearly singular, using pseudoinverse")
    J_init_R_inv = pinv(J_init_R)
    println("Using modified Broyden with pseudoinverse...")
    # Simple Newton step instead of full Broyden
    v_finetuned_R = Vst - J_init_R_inv * fmR(Vst)
else
    println("Broyden iterations (Relativistic)...")
    v_finetuned_R = bsolve(fmR, J_init_R, fmR(Vst), Vst, 10)
end

# Final relativistic result
f_final_R = fmR(v_finetuned_R)
miss_final_R = norm(f_final_R) / AU

println("Fine-tuned miss distance (Relativistic): $(miss_final_R) AU")
println("Improvement factor (Relativistic): $(dmissR / norm(f_final_R))")

# Compute fine-tuned trajectory for relativistic case
v_magft_R = norm(v_finetuned_R)
γft_R = one(tpfl)/sqrt(one(tpfl)-(v_magft_R/c)^2)
pvft_R = (mchip*γft_R) .* v_finetuned_R
Ptestft_R = pomin.setup_single_particle(mchip, qchip, pvft_R, tpfl)
sol_ft_R = pomin.solve(Pint, params; testparticles=Ptestft_R)
Zend_ft_R = sol_ft_R(tcl)

# Extract final positions
zendRint_ft = Zend_ft_R[1:6*nInt]
zendRtest_ft = Zend_ft_R[6*nInt+1:end]
q_proxima_ft_R = zendRint_ft[1:3]
q_spacecraft_ft_R = zendRtest_ft[1:3]

println("Fine-Tuned Trajectory Differences:")

# Position differences
println("Position Differences:")
Δq_total = q_spacecraft_ft_R - q_spacecraft_ft_N
Δq_total_mag = norm(Δq_total)
println("  Total:       $(Printf.@sprintf("%.6e", Δq_total_mag/AU)) AU")
Δq_parallel = Δq_total .* (v_finetuned_R ./ norm(v_finetuned_R))
Δq_parallel_mag = norm(Δq_parallel)
println("  Parallel:    $(Printf.@sprintf("%.6e", Δq_parallel_mag/AU)) AU")
Δq_transverse = Δq_total - Δq_parallel
Δq_transverse_mag = norm(Δq_transverse)
println("  Transverse:  $(Printf.@sprintf("%.6e", Δq_transverse_mag/AU)) AU")

# Velocity differences
println("\nVelocity Comparison:")
v_ft_newton = v_finetuned_N
v_ft_rel = v_finetuned_R
println("  Newtonian:   $(Printf.@sprintf("%.10f", norm(v_ft_newton))) c")
println("  Relativistic: $(Printf.@sprintf("%.10f", norm(v_ft_rel))) c")
Δv_ft = v_ft_rel - v_ft_newton
println("  Difference:  $(Printf.@sprintf("%.6e", norm(Δv_ft))) c")
angle_diff_ft_deg = rad2deg(acos(dot(v_ft_newton, v_ft_rel) / (norm(v_ft_newton) * norm(v_ft_rel))))
angle_diff_ft_arcsec = angle_diff_ft_deg * 3600
println("  Angle diff:  $(Printf.@sprintf("%.6e", angle_diff_ft_deg))° ($(Printf.@sprintf("%.6e", angle_diff_ft_arcsec))\")")

# Verification
println("\nVerification: √(∥² + ⊥²) = $(Printf.@sprintf("%.6e", sqrt((Δq_parallel_mag/AU)^2 + (Δq_transverse_mag/AU)^2))) AU")

#-----------------------------------------------------------------------
#   SAVE RESULTS TO TEXT FILE
#-----------------------------------------------------------------------

println("\nSaving results to test_particle_data.txt...")

open("test_particle_data.txt", "w") do file
    println(file, "="^70)
    println(file, "TEST PARTICLE INITIAL DATA AND OPTIMIZATION RESULTS")
    println(file, "Generated: $(now())")
    println(file, "="^70)
    
    println(file, "\nINITIAL TEST PARTICLE DATA:")
    println(file, "Mass (solar masses): $(mchip)")
    println(file, "Initial position (geometric units): $(qchip)")
    println(file, "Initial velocity magnitude: $(vst) c")
    println(file, "Initial velocity vector: $(Vst)")
    println(file, "Lorentz factor: $(γst)")
    
    println(file, "\nINITIAL MOMENTA:")
    println(file, "Relativistic momentum: $(pchip_rel)")
    println(file, "Newtonian momentum: $(pchip_newt)")
    
    println(file, "\nINITIAL MISS DISTANCES:")
    println(file, "Newtonian miss distance: $(dmissN/AU) AU")
    println(file, "Relativistic miss distance: $(dmissR/AU) AU")
    println(file, "Difference: $((dmissR-dmissN)/AU) AU")
    
    println(file, "\nFINE-TUNED VELOCITIES:")
    println(file, "Newtonian optimized velocity: $(v_finetuned_N)")
    println(file, "Relativistic optimized velocity: $(v_finetuned_R)")
    println(file, "Newtonian speed: $(Printf.@sprintf("%.10f", norm(v_finetuned_N))) c")
    println(file, "Relativistic speed: $(Printf.@sprintf("%.10f", norm(v_finetuned_R))) c")
    
    println(file, "\nOPTIMIZATION RESULTS:")
    println(file, "Newtonian final miss distance: $(miss_final_N) AU")
    println(file, "Relativistic final miss distance: $(miss_final_R) AU")
    println(file, "Newtonian improvement factor: $(dmissN / norm(f_final_N))")
    println(file, "Relativistic improvement factor: $(dmissR / norm(f_final_R))")
    
    println(file, "\nFINAL TRAJECTORY DIFFERENCES:")
    println(file, "Total position difference: $(Printf.@sprintf("%.6e", Δq_total_mag/AU)) AU")
    println(file, "Parallel component: $(Printf.@sprintf("%.6e", Δq_parallel_mag/AU)) AU")
    println(file, "Transverse component: $(Printf.@sprintf("%.6e", Δq_transverse_mag/AU)) AU")
    println(file, "Velocity difference: $(Printf.@sprintf("%.6e", norm(Δv_ft))) c")
    println(file, "Angular difference: $(Printf.@sprintf("%.6e", angle_diff_ft_deg))° ($(Printf.@sprintf("%.6e", angle_diff_ft_arcsec))\")")
    
    println(file, "\nSYSTEM PARAMETERS:")
    println(file, "Integration time: $(tcl) (geometric units)")
    println(file, "Integration tolerances: $(tols)")
    println(file, "Target offset: $(bv)")
    println(file, "Number of main particles: $(nInt)")
    
    println(file, "\n" * "="^70)
end

println("Results saved to test_particle_data.txt")

# Save corrected initial data for the chip
println("Saving corrected chip initial data to chip_corrected_data.txt...")

open("chip_corrected_data.txt", "w") do file
    println(file, "="^70)
    println(file, "CORRECTED TEST PARTICLE (CHIP) INITIAL DATA")
    println(file, "Generated: $(now())")
    println(file, "="^70)
    
    println(file, "\nORIGINAL INITIAL DATA:")
    println(file, "Mass (solar masses): $(mchip)")
    println(file, "Position (geometric units): $(qchip)")
    println(file, "Original velocity: $(Vst)")
    println(file, "Original speed: $(norm(Vst)) c")
    
    println(file, "\nCORRECTED NEWTONIAN DATA:")
    println(file, "Mass (solar masses): $(mchip)")
    println(file, "Position (geometric units): $(qchip)")
    println(file, "Corrected velocity: $(v_finetuned_N)")
    println(file, "Corrected speed: $(Printf.@sprintf("%.15f", norm(v_finetuned_N))) c")
    println(file, "Corrected momentum: $(mchip .* v_finetuned_N)")
    
    println(file, "\nCORRECTED RELATIVISTIC DATA:")
    println(file, "Mass (solar masses): $(mchip)")
    println(file, "Position (geometric units): $(qchip)")
    println(file, "Corrected velocity: $(v_finetuned_R)")
    println(file, "Corrected speed: $(Printf.@sprintf("%.15f", norm(v_finetuned_R))) c")
    v_magft_R = norm(v_finetuned_R)
    γft_R = one(tpfl)/sqrt(one(tpfl)-(v_magft_R/c)^2)
    println(file, "Lorentz factor: $(Printf.@sprintf("%.15f", γft_R))")
    println(file, "Corrected momentum: $((mchip*γft_R) .* v_finetuned_R)")
    
    println(file, "\nVELOCITY CORRECTIONS:")
    δv_N = v_finetuned_N - Vst
    δv_R = v_finetuned_R - Vst
    println(file, "Newtonian velocity correction: $(δv_N)")
    println(file, "Relativistic velocity correction: $(δv_R)")
    println(file, "Newtonian correction magnitude: $(Printf.@sprintf("%.6e", norm(δv_N))) c")
    println(file, "Relativistic correction magnitude: $(Printf.@sprintf("%.6e", norm(δv_R))) c")
    
    println(file, "\nPERFORMANCE METRICS:")
    println(file, "Newtonian final miss: $(Printf.@sprintf("%.6e", miss_final_N)) AU")
    println(file, "Relativistic final miss: $(Printf.@sprintf("%.6e", miss_final_R)) AU")
    println(file, "Newtonian improvement: $(Printf.@sprintf("%.6e", dmissN / norm(f_final_N)))")
    println(file, "Relativistic improvement: $(Printf.@sprintf("%.6e", dmissR / norm(f_final_R)))")
    
    println(file, "\n" * "="^70)
end

println("Corrected chip data saved to chip_corrected_data.txt")
