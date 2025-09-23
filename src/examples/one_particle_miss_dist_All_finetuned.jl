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
using LinearAlgebra, Printf, Plots, Statistics
using DoubleFloats, ForwardDiff, Dates
tpfl  = Double64

include("target_proxima.jl")
include("../utils/broyden.jl")

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0
vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_newt0 = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)


Vchip_rel_FT  = [tpfl("-7.29280133630346695175238301027084351e-02"), 
                 tpfl("-5.57596267171519947482140178590898884e-02"), 
                 tpfl("-1.77686212323486749509362995970259857e-01")]
Vchip_newt_FT = [tpfl("-7.2928011848008478283629726326358021e-02"), 
                 tpfl("-5.57596297969013443215818900874537475e-02"), 
                 tpfl("-1.77686212073603716257005595684594549e-01")]

γchip_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_rel_FT))^2)
pchip_rel_FT = (mchip*γchip_rel_FT) .* Vchip_rel_FT
pchip_newt_FT = (mchip) .* Vchip_newt_FT

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel_FT, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt_FT, tpfl)

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
#   INITIAL DATA SETUP FOR EARTH (astropy, J2000)
#-----------------------------------------------------------------------

mEarth = tpfl(3.0033693739E-06)

# Earth coordinates generated using astropy with epoch J2000
qEarth = tpfl.([-1.8667826140E+07, 
                 8.9633743993E+07, 
                 3.8883291714E+07])         # geometric units
vEarth = tpfl.([-9.9351889046e-05,
                -1.6777453395e-05,
                -7.2738469450e-06])         # geometric units

γEarth = one(tpfl) / sqrt(one(tpfl) - norm(vEarth)^2)  # Lorentz factor
pEarth = tpfl.(mEarth * γEarth .* vEarth) # Earth momentum (relativ.)
pEarthN = mEarth .* vEarth                # Earth momentum (Newtonian)

PEarth = pomin.setup_single_particle(mEarth, qEarth, pEarth, tpfl)
PEarthN = pomin.setup_single_particle(mEarth, qEarth, pEarthN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR JUPITER (astropy, J2000)
#-----------------------------------------------------------------------

mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([3.99442362 * AU, 
              2.73345644 * AU, 
              1.07451704 * AU])             # geometric units

vjup = tpfl.([-7.8875477321829800, 
              1.0175858402135900E+01, 
              4.5538873116399600])  # in km/s (astropy, J2000)

vjup *= tpfl(1000 / cMKS)  # convert from km/s to units of c

γjup = one(tpfl) / sqrt(one(tpfl) - norm(vjup)^2)  # Lorentz factor
pjup = tpfl.(mjup * γjup .* vjup)  # Jupiter momentum 
pjupN = mjup .* vjup                # Jupiter momentum (Newtonian)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)
PjupN = pomin.setup_single_particle(mjup, qjup, pjupN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MOON (astropy, J2000)
#-----------------------------------------------------------------------

mMoon = tpfl(7.34767309E22 / 1.988416E30)     # in units of solar masses

# Moon coordinates generated using astropy with epoch J2000
qMoon = tpfl.([-1.886529874442160E+07, 
                8.94531273942814E+07, 
                3.883175807801640E+07])    # geometric units
vMoon = tpfl.([-2.9141440121218000E+01, 
                -5.6958177245812600, 
                -2.4819741171379700])  # in km/s (astropy, J2000)

vMoon *= tpfl(1000 / cMKS)  # convert from km/s to units of c

γMoon = one(tpfl) / sqrt(one(tpfl) - norm(vMoon)^2)  # Lorentz factor
pMoon = tpfl.(mMoon * γMoon * vMoon)  # Moon momentum 
pMoonN = mMoon .* vMoon                # Moon momentum (Newtonian)

PMoon = pomin.setup_single_particle(mMoon, qMoon, pMoon, tpfl)
PMoonN = pomin.setup_single_particle(mMoon, qMoon, pMoonN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MARS (astropy, J2000)
#-----------------------------------------------------------------------

mMars = tpfl(3.2271058587E-07)  # Mars mass in solar masses
qMars = tpfl.([1.38356874 * AU, 
               -0.00120916 * AU, 
               -0.03786079 * AU])  # Mars position in geometric units (from astropy, J2000)

vMars = tpfl.([1.17347755650600, 
               2.390740193498500E+01, 
               1.0934201867474400E+01])  # in km/s (astropy, J2000)

vMars *= tpfl(1000 / cMKS)  # convert from km/s to units of c

γMars = one(tpfl) / sqrt(one(tpfl) - norm(vMars)^2)  # Lorentz factor
pMars = tpfl.(mMars * γMars * vMars)  # Mars momentum 
pMarsN = mMars .* vMars                # Mars momentum (Newtonian)

PMars = pomin.setup_single_particle(mMars, qMars, pMars, tpfl)
PMarsN = pomin.setup_single_particle(mMars, qMars, pMarsN, tpfl)

#-----------------------------------------------------------------------
#   BROYDEN OPTIMIZATION FUNCTIONS FOR FINE-TUNING
#-----------------------------------------------------------------------

# Miss function constructor - following your existing pattern
function missfuncconstructor_single(body_particle, Ptest0, params, target_pos, tcl; Newt=false)
    if Newt
        return function(v)
            pv = (Ptest0.m[1]) .* v
            Ptest = pomin.setup_single_particle(Ptest0.m[1], Ptest0.q[1], pv, eltype(Ptest0.m))
            sol_modified = pomin.solve(body_particle, params; testparticles=Ptest, Newtonian=true)
            Zend_modified = sol_modified(tcl)
            q_spacecraft = Zend_modified[7:9]  # Test particle position
            return target_pos - q_spacecraft
        end
    else    
        return function(v)
            v_mag = norm(v)
            γ = one(eltype(v)) / sqrt(one(eltype(v)) - (v_mag/c)^2)
            pv = (Ptest0.m[1] * γ) .* v
            Ptest = pomin.setup_single_particle(Ptest0.m[1], Ptest0.q[1], pv, eltype(Ptest0.m))
            sol_modified = pomin.solve(body_particle, params; testparticles=Ptest)
            Zend_modified = sol_modified(tcl)
            q_spacecraft = Zend_modified[7:9]  # Test particle position
            return target_pos - q_spacecraft
        end
    end
end

# Fine-tuning function using your existing Broyden implementation
function fine_tune_velocity(initial_velocity, body_particle, is_newtonian::Bool, body_name::String)
    println("\n" * "="^60)
    println("FINE-TUNING $body_name ($(is_newtonian ? "Newtonian" : "Relativistic"))")
    println("="^60)
    
    # Create test particle template
    if is_newtonian
        pchip_template = mchip .* initial_velocity
    else
        γ_template = one(tpfl)/sqrt(one(tpfl) - norm(initial_velocity)^2)
        pchip_template = (mchip * γ_template) .* initial_velocity
    end
    Ptest_template = pomin.setup_single_particle(mchip, qchip, pchip_template, tpfl)
    
    # Define miss function using your pattern
    miss_func = missfuncconstructor_single(body_particle, Ptest_template, params, flat_space_spacecraft_final, tcl; Newt=is_newtonian)
    
    # Initial guess and function evaluation
    v0 = tpfl.(initial_velocity)
    f0 = miss_func(v0)
    initial_miss_distance = norm(f0)
    
    println("Initial miss distance: ", initial_miss_distance/AU, " AU")
    println("Initial velocity: ", v0)
    println("Computing Jacobian using ForwardDiff...")
    
    try
        # Compute Jacobian using ForwardDiff (following your pattern)
        J_init = ForwardDiff.jacobian(miss_func, v0)
        
        println("Jacobian condition number: ", cond(J_init))
        println("Jacobian determinant: ", det(J_init))
        
        # Check for singularity and apply Broyden
        if abs(det(J_init)) < tpfl(1e-10)
            println("WARNING: Jacobian is nearly singular, using pseudoinverse")
            J_init_inv = pinv(J_init)
            println("Using Newton step with pseudoinverse...")
            v_finetuned = v0 - J_init_inv * f0
        else
            println("Broyden iterations...")
            v_finetuned = bsolve(miss_func, J_init, f0, v0, 10)
        end
        
        # Final result
        f_final = miss_func(v_finetuned)
        final_miss_distance = norm(f_final)
        
        println("Optimization completed!")
        println("Final velocity: ", v_finetuned)
        println("Final miss distance: ", final_miss_distance/AU, " AU")
        println("Improvement factor: ", initial_miss_distance / final_miss_distance)
        
        return v_finetuned, final_miss_distance
        
    catch e
        println("Optimization failed: ", e)
        println("Using initial velocity")
        return v0, initial_miss_distance
    end
end

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

#Pint = pomin.merge_particle_systems(PProx, Psol)
#PintN = pomin.merge_particle_systems(PProxN, Psol)
#Pint = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup)
#PintN = pomin.merge_particle_systems(PProxN, Psol, PalphaN, PjupN)
#Pint = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup, PEarth)
#PintN = pomin.merge_particle_systems(PProxN, Psol, PalphaN, PjupN, PEarthN)
#Pint = pomin.merge_particle_systems(PProx, Psol, Pjup, Palpha, PEarth, PMoon, PMars)
#PintN = pomin.merge_particle_systems(PProxN, Psol, PjupN, PalphaN, PEarthN, PMoonN, PMarsN)

#nInt = length(Pint.m)

# Get flat space spacecraft final position
XstF,XpxF,XbF,tcl_from_target = pfs
flat_space_spacecraft_final = XstF(tcl)

#-----------------------------------------------------------------------
#   SUN - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Sun case
Vchip_sun_newt_FT, dmiss_sun_newt = fine_tune_velocity(Vst, Psol, true, "SUN")
Vchip_sun_rel_FT, dmiss_sun_rel = fine_tune_velocity(Vst, Psol, false, "SUN")

# Create fine-tuned particles
γchip_sun_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_sun_rel_FT))^2)
pchip_sun_rel_FT = (mchip*γchip_sun_rel_FT) .* Vchip_sun_rel_FT
pchip_sun_newt_FT = (mchip) .* Vchip_sun_newt_FT

Pchip_sun_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_sun_rel_FT, tpfl)
Pchip_sun_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_sun_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(Psol, params; testparticles=Pchip_sun_newt_FT, Newtonian=true)
sol = pomin.solve(Psol, params; testparticles=Pchip_sun_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("SUN - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table (using proper precision type)
miss_distances_N = tpfl.([dmissN/AU])
miss_distances_R = tpfl.([dmiss/AU])
body_names_completed = ["Sun"]

#-----------------------------------------------------------------------
#   ALPHA CENTAURI - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Alpha Centauri case
Vchip_alpha_newt_FT, dmiss_alpha_newt = fine_tune_velocity(Vst, Palpha, true, "ALPHA CENTAURI")
Vchip_alpha_rel_FT, dmiss_alpha_rel = fine_tune_velocity(Vst, Palpha, false, "ALPHA CENTAURI")

# Create fine-tuned particles
γchip_alpha_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_alpha_rel_FT))^2)
pchip_alpha_rel_FT = (mchip*γchip_alpha_rel_FT) .* Vchip_alpha_rel_FT
pchip_alpha_newt_FT = (mchip) .* Vchip_alpha_newt_FT

Pchip_alpha_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_alpha_rel_FT, tpfl)
Pchip_alpha_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_alpha_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(Palpha, params; testparticles=Pchip_alpha_newt_FT, Newtonian=true)
sol = pomin.solve(Palpha, params; testparticles=Pchip_alpha_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("ALPHA CENTAURI - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Alpha Centauri")

#-----------------------------------------------------------------------
#   JUPITER - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Jupiter case
Vchip_jup_newt_FT, dmiss_jup_newt = fine_tune_velocity(Vst, Pjup, true, "JUPITER")
Vchip_jup_rel_FT, dmiss_jup_rel = fine_tune_velocity(Vst, Pjup, false, "JUPITER")

# Create fine-tuned particles
γchip_jup_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_jup_rel_FT))^2)
pchip_jup_rel_FT = (mchip*γchip_jup_rel_FT) .* Vchip_jup_rel_FT
pchip_jup_newt_FT = (mchip) .* Vchip_jup_newt_FT

Pchip_jup_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_jup_rel_FT, tpfl)
Pchip_jup_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_jup_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(Pjup, params; testparticles=Pchip_jup_newt_FT, Newtonian=true)
sol = pomin.solve(Pjup, params; testparticles=Pchip_jup_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("JUPITER - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Jupiter")

#-----------------------------------------------------------------------
#   EARTH - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Earth case
Vchip_earth_newt_FT, dmiss_earth_newt = fine_tune_velocity(Vst, PEarth, true, "EARTH")
Vchip_earth_rel_FT, dmiss_earth_rel = fine_tune_velocity(Vst, PEarth, false, "EARTH")

# Create fine-tuned particles
γchip_earth_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_earth_rel_FT))^2)
pchip_earth_rel_FT = (mchip*γchip_earth_rel_FT) .* Vchip_earth_rel_FT
pchip_earth_newt_FT = (mchip) .* Vchip_earth_newt_FT

Pchip_earth_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_earth_rel_FT, tpfl)
Pchip_earth_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_earth_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(PEarth, params; testparticles=Pchip_earth_newt_FT, Newtonian=true)
sol = pomin.solve(PEarth, params; testparticles=Pchip_earth_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("EARTH - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Earth")

#-----------------------------------------------------------------------
#   PROXIMA - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Proxima case
Vchip_prox_newt_FT, dmiss_prox_newt = fine_tune_velocity(Vst, PProx, true, "PROXIMA")
Vchip_prox_rel_FT, dmiss_prox_rel = fine_tune_velocity(Vst, PProx, false, "PROXIMA")

# Create fine-tuned particles
γchip_prox_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_prox_rel_FT))^2)
pchip_prox_rel_FT = (mchip*γchip_prox_rel_FT) .* Vchip_prox_rel_FT
pchip_prox_newt_FT = (mchip) .* Vchip_prox_newt_FT

Pchip_prox_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_prox_rel_FT, tpfl)
Pchip_prox_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_prox_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(PProx, params; testparticles=Pchip_prox_newt_FT, Newtonian=true)
sol = pomin.solve(PProx, params; testparticles=Pchip_prox_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("PROXIMA - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Proxima")

#-----------------------------------------------------------------------
#   MOON - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Moon case
Vchip_moon_newt_FT, dmiss_moon_newt = fine_tune_velocity(Vst, PMoon, true, "MOON")
Vchip_moon_rel_FT, dmiss_moon_rel = fine_tune_velocity(Vst, PMoon, false, "MOON")

# Create fine-tuned particles
γchip_moon_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_moon_rel_FT))^2)
pchip_moon_rel_FT = (mchip*γchip_moon_rel_FT) .* Vchip_moon_rel_FT
pchip_moon_newt_FT = (mchip) .* Vchip_moon_newt_FT

Pchip_moon_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_moon_rel_FT, tpfl)
Pchip_moon_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_moon_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(PMoon, params; testparticles=Pchip_moon_newt_FT, Newtonian=true)
sol = pomin.solve(PMoon, params; testparticles=Pchip_moon_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("MOON - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Moon")

#-----------------------------------------------------------------------
#   MARS - WITH FINE-TUNING
#-----------------------------------------------------------------------

# Fine-tune velocities for Mars case
Vchip_mars_newt_FT, dmiss_mars_newt = fine_tune_velocity(Vst, PMars, true, "MARS")
Vchip_mars_rel_FT, dmiss_mars_rel = fine_tune_velocity(Vst, PMars, false, "MARS")

# Create fine-tuned particles
γchip_mars_rel_FT = one(tpfl)/sqrt(one(tpfl)-(norm(Vchip_mars_rel_FT))^2)
pchip_mars_rel_FT = (mchip*γchip_mars_rel_FT) .* Vchip_mars_rel_FT
pchip_mars_newt_FT = (mchip) .* Vchip_mars_newt_FT

Pchip_mars_rel_FT = pomin.setup_single_particle(mchip, qchip, pchip_mars_rel_FT, tpfl)
Pchip_mars_newt_FT = pomin.setup_single_particle(mchip, qchip, pchip_mars_newt_FT, tpfl)

# Verify final trajectories
solN = pomin.solve(PMars, params; testparticles=Pchip_mars_newt_FT, Newtonian=true)
sol = pomin.solve(PMars, params; testparticles=Pchip_mars_rel_FT)

zendN = solN(tcl)
zend = sol(tcl)

qmissN = zendN[7:9] - flat_space_spacecraft_final
qmiss = zend[7:9] - flat_space_spacecraft_final
dmissN = norm(qmissN)
dmiss = norm(qmiss)

println("\n" * "="^60)
println("MARS - FINAL RESULTS AFTER FINE-TUNING")
println("="^60)
println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("Relativistic miss: ", qmiss)
println("Relativistic miss distance: ", dmiss, " (", dmiss/AU, " AU)")

# Store results for summary table
push!(miss_distances_N, tpfl(dmissN/AU))
push!(miss_distances_R, tpfl(dmiss/AU))
push!(body_names_completed, "Mars")

#-----------------------------------------------------------------------
#   SAVE OPTIMIZED VELOCITIES FOR FUTURE USE
#-----------------------------------------------------------------------

println("\n" * "="^80)
println("OPTIMIZED VELOCITIES - FOR FUTURE REFERENCE")
println("="^80)

println("# Fine-tuned Newtonian velocities:")
println("Vchip_sun_newt_FT    = ", Vchip_sun_newt_FT)
println("Vchip_alpha_newt_FT  = ", Vchip_alpha_newt_FT)
println("Vchip_jup_newt_FT    = ", Vchip_jup_newt_FT)
println("Vchip_earth_newt_FT  = ", Vchip_earth_newt_FT)
println("Vchip_prox_newt_FT   = ", Vchip_prox_newt_FT)
println("Vchip_moon_newt_FT   = ", Vchip_moon_newt_FT)
println("Vchip_mars_newt_FT   = ", Vchip_mars_newt_FT)

println("\n# Fine-tuned Relativistic velocities:")
println("Vchip_sun_rel_FT     = ", Vchip_sun_rel_FT)
println("Vchip_alpha_rel_FT   = ", Vchip_alpha_rel_FT)
println("Vchip_jup_rel_FT     = ", Vchip_jup_rel_FT)
println("Vchip_earth_rel_FT   = ", Vchip_earth_rel_FT)
println("Vchip_prox_rel_FT    = ", Vchip_prox_rel_FT)
println("Vchip_moon_rel_FT    = ", Vchip_moon_rel_FT)
println("Vchip_mars_rel_FT    = ", Vchip_mars_rel_FT)

println("="^80)

#-----------------------------------------------------------------------
#   COMPREHENSIVE SUMMARY TABLE WITH FINE-TUNED RESULTS
#-----------------------------------------------------------------------

println("\n" * "="^100)
println("COMPREHENSIVE MISS DISTANCE SUMMARY - AFTER BROYDEN FINE-TUNING")
println("="^100)
using Printf

# Print header
println(@sprintf("%-15s | %-20s | %-20s | %-15s | %-15s", 
        "Body", "Newtonian Miss (AU)", "Relativistic Miss (AU)", "N/R Ratio", "Improvement"))
println("-"^100)

# Print results for each body
for i in 1:length(body_names_completed)
    ratio = miss_distances_N[i] / miss_distances_R[i]
    
    # Determine improvement level
    improvement = if miss_distances_N[i] < 1e-10 && miss_distances_R[i] < 1e-10
        "Excellent"
    elseif miss_distances_N[i] < 1e-8 && miss_distances_R[i] < 1e-8
        "Very Good"
    elseif miss_distances_N[i] < 1e-6 && miss_distances_R[i] < 1e-6
        "Good"
    else
        "Moderate"
    end
    
    println(@sprintf("%-15s | %-20.6e | %-20.6e | %-15.6f | %-15s", 
            body_names_completed[i], miss_distances_N[i], miss_distances_R[i], ratio, improvement))
end

println("-"^100)

# Print summary statistics
println("\nSUMMARY STATISTICS:")
println("==================")
println(@sprintf("Total bodies analyzed: %d", length(body_names_completed)))
println(@sprintf("Average Newtonian miss: %.6e AU", mean(miss_distances_N)))
println(@sprintf("Average Relativistic miss: %.6e AU", mean(miss_distances_R)))
println(@sprintf("Average N/R ratio: %.6f", mean(miss_distances_N ./ miss_distances_R)))

# Find best and worst cases
best_newt_idx = argmin(miss_distances_N)
best_rel_idx = argmin(miss_distances_R)
worst_newt_idx = argmax(miss_distances_N)
worst_rel_idx = argmax(miss_distances_R)

println(@sprintf("\nBest Newtonian result: %s (%.6e AU)", 
        body_names_completed[best_newt_idx], miss_distances_N[best_newt_idx]))
println(@sprintf("Best Relativistic result: %s (%.6e AU)", 
        body_names_completed[best_rel_idx], miss_distances_R[best_rel_idx]))
println(@sprintf("Worst Newtonian result: %s (%.6e AU)", 
        body_names_completed[worst_newt_idx], miss_distances_N[worst_newt_idx]))
println(@sprintf("Worst Relativistic result: %s (%.6e AU)", 
        body_names_completed[worst_rel_idx], miss_distances_R[worst_rel_idx]))

println("\n" * "="^100)
println("FINE-TUNING COMPLETE - All trajectories optimized using Broyden method")
println("Target: Proxima Centauri at distance ≈ 4.24 light-years")
println("Integration time: ", tcl, " time units")
println("Integration tolerance: ", tols)
println("="^100)

#-----------------------------------------------------------------------
#   GENERATE TEXT FILE WITH INITIAL DATA FOR EXTERNAL USE
#-----------------------------------------------------------------------

println("\nGenerating comprehensive initial data file...")

open("spacecraft_initial_data_finetuned.txt", "w") do file
    println(file, "="^100)
    println(file, "COMPREHENSIVE SPACECRAFT INITIAL DATA - FINE-TUNED FOR ALL CELESTIAL BODIES")
    println(file, "Generated: ", Dates.now())
    println(file, "Precision: Double64 ($(precision(tpfl(1.0))) bits)")
    println(file, "="^100)
    
    println(file, "\nSPACECRAFT PARAMETERS:")
    println(file, "=====================")
    println(file, "Mass (solar masses): ", mchip)
    println(file, "Initial position (geometric units): ", qchip)
    println(file, "Original velocity: ", Vst)
    println(file, "Original speed: ", norm(Vst), " c")
    println(file, "Target: Proxima Centauri")
    println(file, "Target position: ", flat_space_spacecraft_final)
    
    println(file, "\nCELESTIAL BODY MASSES AND POSITIONS:")
    println(file, "====================================")
    println(file, "Sun:")
    println(file, "  Mass: ", msol, " solar masses")
    println(file, "  Position: ", qsol)
    println(file, "  Momentum: ", psol)
    
    println(file, "\nProxima Centauri:")
    println(file, "  Mass: ", mProx, " solar masses")
    println(file, "  Position: ", qProx)
    println(file, "  Relativistic momentum: ", pProx)
    println(file, "  Newtonian momentum: ", pProxN)
    
    println(file, "\nAlpha Centauri A+B:")
    println(file, "  Mass: ", malpha, " solar masses")
    println(file, "  Position: ", qalpha)
    println(file, "  Velocity: ", valpha)
    println(file, "  Lorentz factor: ", γalpha)
    println(file, "  Relativistic momentum: ", palpha)
    println(file, "  Newtonian momentum: ", palphaN)
    
    println(file, "\nEarth:")
    println(file, "  Mass: ", mEarth, " solar masses")
    println(file, "  Position: ", qEarth)
    println(file, "  Velocity: ", vEarth)
    println(file, "  Lorentz factor: ", γEarth)
    println(file, "  Relativistic momentum: ", pEarth)
    println(file, "  Newtonian momentum: ", pEarthN)
    
    println(file, "\nJupiter:")
    println(file, "  Mass: ", mjup, " solar masses")
    println(file, "  Position: ", qjup)
    println(file, "  Velocity: ", vjup)
    println(file, "  Lorentz factor: ", γjup)
    println(file, "  Relativistic momentum: ", pjup)
    println(file, "  Newtonian momentum: ", pjupN)
    
    println(file, "\nMoon:")
    println(file, "  Mass: ", mMoon, " solar masses")
    println(file, "  Position: ", qMoon)
    println(file, "  Velocity: ", vMoon)
    println(file, "  Lorentz factor: ", γMoon)
    println(file, "  Relativistic momentum: ", pMoon)
    println(file, "  Newtonian momentum: ", pMoonN)
    
    println(file, "\nMars:")
    println(file, "  Mass: ", mMars, " solar masses")
    println(file, "  Position: ", qMars)
    println(file, "  Velocity: ", vMars)
    println(file, "  Lorentz factor: ", γMars)
    println(file, "  Relativistic momentum: ", pMars)
    println(file, "  Newtonian momentum: ", pMarsN)
    
    println(file, "\nFINE-TUNED NEWTONIAN VELOCITIES:")
    println(file, "================================")
    println(file, "Sun:           ", Vchip_sun_newt_FT)
    println(file, "Alpha Centauri: ", Vchip_alpha_newt_FT)
    println(file, "Jupiter:       ", Vchip_jup_newt_FT)
    println(file, "Earth:         ", Vchip_earth_newt_FT)
    println(file, "Proxima:       ", Vchip_prox_newt_FT)
    println(file, "Moon:          ", Vchip_moon_newt_FT)
    println(file, "Mars:          ", Vchip_mars_newt_FT)
    
    println(file, "\nFINE-TUNED RELATIVISTIC VELOCITIES:")
    println(file, "====================================")
    println(file, "Sun:           ", Vchip_sun_rel_FT)
    println(file, "Alpha Centauri: ", Vchip_alpha_rel_FT)
    println(file, "Jupiter:       ", Vchip_jup_rel_FT)
    println(file, "Earth:         ", Vchip_earth_rel_FT)
    println(file, "Proxima:       ", Vchip_prox_rel_FT)
    println(file, "Moon:          ", Vchip_moon_rel_FT)
    println(file, "Mars:          ", Vchip_mars_rel_FT)
    
    println(file, "\nFINAL MISS DISTANCES (AU):")
    println(file, "===========================")
    for i in 1:length(body_names_completed)
        println(file, @sprintf("%-15s: Newtonian = %.6e AU, Relativistic = %.6e AU", 
                body_names_completed[i], miss_distances_N[i], miss_distances_R[i]))
    end
    
    println(file, "\nINTEGRATION PARAMETERS:")
    println(file, "======================")
    println(file, "Time span: ", params.tspan)
    println(file, "Integration time: ", tcl, " time units")
    println(file, "Absolute tolerance: ", params.atol)
    println(file, "Relative tolerance: ", params.rtol)
    println(file, "Integrator: ", params.integrator)
    
    println(file, "\nCONSTANTS AND UNITS:")
    println(file, "====================")
    println(file, "Speed of light (c): ", c)
    println(file, "Astronomical Unit (AU): ", AU)
    println(file, "Speed of light in MKS (cMKS): ", cMKS)
    
    println(file, "\n" * "="^100)
    println(file, "END OF INITIAL DATA FILE")
    println(file, "All values are in geometric units unless otherwise specified")
    println(file, "Precision type: ", typeof(mchip))
    println(file, "="^100)
end

println("Initial data saved to: spacecraft_initial_data_finetuned.txt")
println("File contains all celestial body parameters and fine-tuned velocities")
println("Ready for use in external applications with full Double64 precision")
