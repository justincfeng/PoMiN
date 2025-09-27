#-----------------------------------------------------------------------
#
#   NEWTONIAN VS RELATIVISTIC COMPARISON WITH FINE-TUNING
#   Spacecraft trajectory to Proxima Centauri
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
#   RATIO OF CLOSEST APPROACH TO IMPACT PARAMETER TERMS
#-----------------------------------------------------------------------
"""
    rcb_ratio(b,v,M,tpfl=eltype(mass),G=one(tpfl))

Calculate the ratio of closest approach to impact parameter for a given 
trajectory.
"""
function rcb_ratio(b,v,M,tpfl=eltype(M),G=one(tpfl))
    l=one(tpfl)
    TWO=tpfl(2)
    FOUR=tpfl(4)
    return ( l , (G*M)/(b*v^2) , (G^2*M^2*(l-FOUR*v^2))/(TWO*b^2*v^4) )
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CLOSEST APPROACH
#-----------------------------------------------------------------------
"""
    closest_approach(xa,xb,va,vb)

Calculate the closest approach between two trajectories.
"""
function closest_approach(xa,xb,va,vb)
    Δx = xb-xa
    Δv = vb-va
    return sqrt( dot(Δx,Δx) - dot(Δx,Δv)^2/dot(Δv,Δv) )
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CLOSEST APPROACH VECTOR
#-----------------------------------------------------------------------
"""
    closest_approach_vec(xa,xb,va,vb)

Calculate the closest approach between two trajectories.
"""
function closest_approach_vec(xa,xb,va,vb)
    Δx = xb-xa
    Δv = vb-va
    tcl = - dot(Δx,Δv)/dot(Δv,Δv)
    return Δx + tcl*Δv
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MISS DISTANCE DECOMPOSITION
#-----------------------------------------------------------------------
"""
    miss_decomposition(miss_vector, velocity_vector)

Decompose miss distance into longitudinal (along velocity) and transverse 
(perpendicular to velocity) components.

Returns: (longitudinal_miss, transverse_miss, total_miss)
"""
function miss_decomposition(miss_vec, vel_vec)
    # Normalize velocity vector
    vel_unit = vel_vec / norm(vel_vec)
    
    # Longitudinal component (projection onto velocity direction)
    longitudinal = dot(miss_vec, vel_unit)
    
    # Transverse component (perpendicular to velocity)
    longitudinal_vec = longitudinal * vel_unit
    transverse_vec = miss_vec - longitudinal_vec
    transverse = norm(transverse_vec)
    
    # Total miss distance
    total = norm(miss_vec)
    
    return (longitudinal, transverse, total)
end #-------------------------------------------------------------------

"""
    transverse_miss_distance(position_final, target_position, velocity_initial)

Calculate the transverse miss distance for a trajectory that may stop near target.
"""
function transverse_miss_distance(pos_final, pos_target, vel_initial)
    miss_vec = pos_final - pos_target
    _, transverse, _ = miss_decomposition(miss_vec, vel_initial)
    return transverse
end #-------------------------------------------------------------------

"""
    longitudinal_miss_distance(position_final, target_position, velocity_initial)

Calculate the longitudinal miss distance for a trajectory that may stop near target.
"""
function longitudinal_miss_distance(pos_final, pos_target, vel_initial)
    miss_vec = pos_final - pos_target
    longitudinal, _, _ = miss_decomposition(miss_vec, vel_initial)
    return longitudinal
end #-------------------------------------------------------------------

function scattering_correction(v,b,M,tpfl=eltype(M),G=one(tpfl))
    ν1 = one(tpfl)
    ν2 = tpfl(2)
    ν3 = tpfl(3)
    ν4 = tpfl(4)
    δϕ1 = (ν2 * M * (ν1 + ν1*v^2))/(b * v^2)
    δϕ2 = (ν3 * tpfl(π) * M^2 * (ν4 + v^2))/(ν4 * b^2 * v^4)
    return (δϕ1,δϕ2)
end #-------------------------------------------------------------------

function miss_estimates(ΔX,v,b,M,tpfl=eltype(M),G=one(tpfl))
    δϕ1,δϕ2 = scattering_correction(v,b,M,tpfl,G)
    return (ΔX*δϕ1,ΔX*δϕ2)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0
vst0 = norm(Vst)

γst0 = one(tpfl)/sqrt(one(tpfl)-(vst0/c)^2)
pchip_rel0 = (mchip*γst0) .* Vst
pchip_newt0 = mchip .* Vst

Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt0 = pomin.setup_single_particle(mchip, qchip, pchip_newt0, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR PROXIMA CENTAURI
#-----------------------------------------------------------------------

# Mass, position, momentum
mProx   = tpfl(0.1221)
qProx   = Xpx0
pProx   = mProx .* γVpx
pProxN  = mProx .* Vpx

bProx   = closest_approach(qProx,qchip,Vpx,Vst)
bProx_vec = closest_approach_vec(qProx,qchip,Vpx,Vst)
ΔVstpx   = Vst - Vpx

rcbProx = rcb_ratio(bProx,norm(ΔVstpx),mProx)
missTHProx   = miss_estimates(ΔXst,norm(ΔVstpx),bProx,mProx)

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
vSol = tpfl.([0.0, 0.0, 0.0])

bSol        = closest_approach(qsol,qchip,vSol,Vst)
ΔXSoli      = norm(qsol-qchip)
missTHSol   = miss_estimates(ΔXst,vst,bSol,msol)

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

bAlpha      = closest_approach(qalpha,qchip,valpha,Vst)
ΔXalphai    = norm(qalpha-qchip)
missTHAlpha = miss_estimates(ΔXst,vst,bAlpha,malpha)

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

bEarth      = closest_approach(qEarth,qchip,vEarth,Vst)
ΔXEarti     = norm(qEarth-qchip)
missTHEarth = miss_estimates(ΔXst,vst,bEarth,mEarth)

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

bJup      = closest_approach(qjup,qchip,vjup,Vst)
ΔXJupi    = norm(qjup-qchip)
missTHJup = miss_estimates(ΔXst,vst,bJup,mjup)

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

bMoon      = closest_approach(qMoon,qchip,vMoon,Vst)
ΔXMooni    = norm(qMoon-qchip)
missTHMoon = miss_estimates(ΔXst,vst,bMoon,mMoon)

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

bMars      = closest_approach(qMars,qchip,vMars,Vst)
ΔXMarsi    = norm(qMars-qchip)
missTHMars = miss_estimates(ΔXst,vst,bMars,mMars)

PMars = pomin.setup_single_particle(mMars, qMars, pMars, tpfl)
PMarsN = pomin.setup_single_particle(mMars, qMars, pMarsN, tpfl)

#-----------------------------------------------------------------------
#
#   RUNS
#
#-----------------------------------------------------------------------

# Integration parameters
tspan = (tpfl(0), tcl*tpfl(1.1))
tols  = tpfl(1e-17)

params = pomin.ParametersJulia(tspan, integrator="Vern7", 
                                     atol=tols, rtol=tols)

# Get Proxima's final position (target)
XstF,XpxF,XbF,tcl_from_target = pfs
qtar = XbF(tcl)
flat_space_spacecraft_final = XstF(tcl)  # Flat space spacecraft position for comparison

# Open file for miss distance output
output_file = "miss_distance_results.txt"
file = open(output_file, "w")

# Write header information
println(file, "="^100)
println(file, "MISS DISTANCE ANALYSIS - FINE-TUNED INITIAL DATA (FTID)")
println(file, "="^100)
println(file, "Generated: ", Dates.now())
println(file, "Target: Proxima Centauri (displaced by ", norm(bv)/AU, " AU)")
println(file, "Integration time: ", tcl, " time units")
println(file, "Precision: Double64")
println(file, "="^100)
println(file)

#-----------------------------------------------------------------------
#   SUN
#-----------------------------------------------------------------------

solN    = pomin.solve( Psol, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( Psol, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHSol[1]
dmissPMHO   = missTHSol[2]

# Store results for summary table
bratio_sun = bSol/ΔXSoli
dmisstar_sun = dmisstar
dmisstarN_sun = dmisstarN
fmisstar_sun = fmisstar
fmisstarN_sun = fmisstarN
dmissPMcomp_sun = dmissPMcomp
dmissPMHO_sun = dmissPMHO

println("SUN RESULTS:")
println("Impact parameter: ", bSol, " (", bSol/AU, " AU)")
println("Initial distance: ", ΔXSoli, " (", ΔXSoli/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "SUN RESULTS:")
println(file, "-----------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bSol, " (", bSol/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   ALPHA CENTAURI
#-----------------------------------------------------------------------

solN    = pomin.solve( Palpha, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( Palpha, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHAlpha[1]
dmissPMHO   = missTHAlpha[2]

bAlphaf     = closest_approach(zendN[1:3],zendN[7:9],vend,γV2v(zend[4:6] ./ malpha))

# Store results for summary table
bratio_alpha = bAlphaf/norm(zendN[7:9]-zendN[1:3])
dmisstar_alpha = dmisstar
dmisstarN_alpha = dmisstarN
fmisstar_alpha = fmisstar
fmisstarN_alpha = fmisstarN
dmissPMcomp_alpha = dmissPMcomp
dmissPMHO_alpha = dmissPMHO

println("ALPHA CENTAURI RESULTS:")
println("Impact parameter: ", bAlpha, " (", bAlpha/AU, " AU)")
println("Initial distance: ", ΔXalphai, " (", ΔXalphai/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "ALPHA CENTAURI RESULTS:")
println(file, "----------------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bAlpha, " (", bAlpha/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   JUPITER
#-----------------------------------------------------------------------

solN    = pomin.solve( Pjup, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( Pjup, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHJup[1]
dmissPMHO   = missTHJup[2]

# Store results for summary table
bratio_jup = bJup/ΔXJupi
dmisstar_jup = dmisstar
dmisstarN_jup = dmisstarN
fmisstar_jup = fmisstar
fmisstarN_jup = fmisstarN
dmissPMcomp_jup = dmissPMcomp
dmissPMHO_jup = dmissPMHO

println("JUPITER RESULTS:")
println("Impact parameter: ", bJup, " (", bJup/AU, " AU)")
println("Initial distance: ", ΔXJupi, " (", ΔXJupi/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "JUPITER RESULTS:")
println(file, "---------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bJup, " (", bJup/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   EARTH
#-----------------------------------------------------------------------

solN    = pomin.solve( PEarth, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( PEarth, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHEarth[1]
dmissPMHO   = missTHEarth[2]

# Store results for summary table
bratio_earth = bEarth/ΔXEarti
dmisstar_earth = dmisstar
dmisstarN_earth = dmisstarN
fmisstar_earth = fmisstar
fmisstarN_earth = fmisstarN
dmissPMcomp_earth = dmissPMcomp
dmissPMHO_earth = dmissPMHO

println("EARTH RESULTS:")
println("Impact parameter: ", bEarth, " (", bEarth/AU, " AU)")
println("Initial distance: ", ΔXEarti, " (", ΔXEarti/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "EARTH RESULTS:")
println(file, "-------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bEarth, " (", bEarth/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   PROXIMA
#-----------------------------------------------------------------------

solN    = pomin.solve( PProx, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( PProx, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = rcbProx[2]*bProx
dmissPMHO   = rcbProx[3]*bProx

# Store results for summary table
bratio_prox = bProx/norm(zendN[7:9]-zendN[1:3])
dmisstar_prox = dmisstar
dmisstarN_prox = dmisstarN
fmisstar_prox = fmisstar
fmisstarN_prox = fmisstarN
dmissPMcomp_prox = dmissPMcomp
dmissPMHO_prox = dmissPMHO

println("PROXIMA RESULTS:")
println("Impact parameter: ", bProx, " (", bProx/AU, " AU)")
println("Initial distance: ", ΔXst, " (", ΔXst/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "PROXIMA RESULTS:")
println(file, "-------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bProx, " (", bProx/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   MOON
#-----------------------------------------------------------------------

solN    = pomin.solve( PMoon, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( PMoon, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHMoon[1]
dmissPMHO   = missTHMoon[2]

# Store results for summary table
bratio_moon = bMoon/ΔXMooni
dmisstar_moon = dmisstar
dmisstarN_moon = dmisstarN
fmisstar_moon = fmisstar
fmisstarN_moon = fmisstarN
dmissPMcomp_moon = dmissPMcomp
dmissPMHO_moon = dmissPMHO

println("MOON RESULTS:")
println("Impact parameter: ", bMoon, " (", bMoon/AU, " AU)")
println("Initial distance: ", ΔXMooni, " (", ΔXMooni/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")

# Write to file
println(file, "MOON RESULTS:")
println(file, "------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bMoon, " (", bMoon/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   MARS
#-----------------------------------------------------------------------

solN    = pomin.solve( PMars, params; testparticles=Pchip_newt0, 
                       Newtonian=true )
sol     = pomin.solve( PMars, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)

qmisstarN   = zendN[7:9] - qtar
vendN       = zendN[10:12] ./ mchip

dmisstarN   = closest_approach(zendN[7:9], qtar, vendN, Vpx)
fmisstarN   = norm(qmisstarN)

qmisstar    = zend[7:9] - qtar
vend        = γV2v(zend[10:12] ./ mchip)

dmisstar = closest_approach(zend[7:9], qtar, vend, Vpx)
fmisstar = norm(qmisstar)

dmissPMcomp = missTHMars[1]
dmissPMHO   = missTHMars[2]

# Store results for summary table
bratio_mars = bMars/ΔXMarsi
dmisstar_mars = dmisstar
dmisstarN_mars = dmisstarN
fmisstar_mars = fmisstar
fmisstarN_mars = fmisstarN
dmissPMcomp_mars = dmissPMcomp
dmissPMHO_mars = dmissPMHO

println("MARS RESULTS:")
println("Impact parameter: ", bMars, " (", bMars/AU, " AU)")
println("Initial distance: ", ΔXMarsi, " (", ΔXMarsi/AU, " AU)")
println("Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println("PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println("PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println("Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")  

# Write to file
println(file, "MARS RESULTS:")
println(file, "------------")
println(file, "Newtonian closest approach: ", dmisstarN, " (", dmisstarN/AU, " AU)")
println(file, "PoMiN closest approach: ", dmisstar, " (", dmisstar/AU, " AU)")
println(file, "PM Miss estimate: ", dmissPMcomp, " (", dmissPMcomp/AU, " AU)")
println(file, "Higher order miss distance: ", dmissPMHO, " (", dmissPMHO/AU, " AU)")
println(file, "Impact parameter: ", bMars, " (", bMars/AU, " AU)")
println(file)


#-----------------------------------------------------------------------
#   COLLECT MISS DISTANCES AND CREATE TABLE
#-----------------------------------------------------------------------

# Collect all results for summary table
results = [
    ("Sun", bratio_sun, fmisstar_sun/AU, fmisstarN_sun/AU, dmisstar_sun/AU, dmissPMcomp_sun/AU, (dmissPMHO_sun/AU)),
    ("Alpha Centauri", bratio_alpha, fmisstar_alpha/AU, fmisstarN_alpha/AU, dmisstar_alpha/AU, dmissPMcomp_alpha/AU, (dmissPMHO_alpha/AU)),
    ("Jupiter", bratio_jup, fmisstar_jup/AU, fmisstarN_jup/AU, dmisstar_jup/AU, dmissPMcomp_jup/AU, dmissPMHO_jup/AU),
    ("Earth", bratio_earth, fmisstar_earth/AU, fmisstarN_earth/AU, dmisstar_earth/AU, dmissPMcomp_earth/AU, (dmissPMHO_earth/AU)),
    ("Proxima", bratio_prox, fmisstar_prox/AU, fmisstarN_prox/AU, dmisstar_prox/AU, dmissPMcomp_prox/AU, dmissPMHO_prox/AU),
    ("Moon", bratio_moon, fmisstar_moon/AU, fmisstarN_moon/AU, dmisstar_moon/AU, dmissPMcomp_moon/AU, dmissPMHO_moon/AU),
    ("Mars", bratio_mars, fmisstar_mars/AU, fmisstarN_mars/AU, dmisstar_mars/AU, dmissPMcomp_mars/AU, dmissPMHO_mars/AU)
]

# Write summary table to file
println(file, "="^100)
println(file, "COMPREHENSIVE SUMMARY TABLE - COMMA SEPARATED VALUES")
println(file, "="^100)
println(file)

# CSV Header
println(file, "Body\t\t\tB/D_ratio\t\tEndpt Dist PoMiN (AU)\tEndpt Dist Newtonian (AU)\tClosest App PoMiN (AU)\tPM Miss Estimate (AU)\tHigher Order Miss (AU)")

# CSV Data rows
for (body, bratio, fmiss_rel, fmiss_newt, dmiss_rel, pm_est, ho_est) in results
    println(file, @sprintf("%s\t\t\t%.6e\t\t\t%.6e\t\t\t%.6e\t\t\t%.6e\t\t\t%.6e\t\t\t%.6e", body, bratio, fmiss_rel, fmiss_newt, dmiss_rel, pm_est, ho_est))
end

println(file)
println(file, "Column Definitions:")
println(file, "- Body: Celestial body name")
println(file, "- B/D_ratio: Ratio of impact parameter to initial distance between spacecraft and celestial body")
println(file, "- Endpoint Dist PoMiN: Distance between spacecraft and target at final time (relativistic)")
println(file, "- Endpoint Dist Newtonian: Distance between spacecraft and target at final time (Newtonian)")
println(file, "- Closest Approach PoMiN: Closest approach distance using final velocities (relativistic)")
println(file, "- PM Miss Estimate: Post-Minkowskian theoretical estimate")
println(file, "- Higher Order Miss: Higher-order relativistic corrections")
println(file)
println(file, "="^100)
println(file, "ANALYSIS COMPLETE")
println(file, "="^100)

# Close the output file
close(file)
println("Miss distance results written to: ", output_file)

# Console summary
println("\n" * "="^80)
println("MISS DISTANCE ANALYSIS COMPLETE")
println("="^80)
println("Results saved to: ", output_file)
println("All miss distances calculated using fine-tuned velocities")
println("Target: Proxima Centauri with ", norm(bv)/AU, " AU displacement")
println("="^80)
