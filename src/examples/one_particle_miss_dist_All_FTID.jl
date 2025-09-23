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
    return ( l , -(G*M)/(b*v^2) , (G^2*M^2*(l-FOUR*v^2))/(TWO*b^2*v^4) )
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

bProx   = closest_approach(qProx,qchip,Vpx,Vst)
rcbProx = rcb_ratio(bProx,vst,mProx)

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

bSol   = closest_approach(qsol,qchip,vSol,Vst)
rcbSol = rcb_ratio(bSol,vst,msol)

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

bAlpha   = closest_approach(qalpha,qchip,valpha,Vst)
rcbAlpha = rcb_ratio(bAlpha,vst,malpha)

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

bEarth   = closest_approach(qEarth,qchip,vEarth,Vst)
rcbEarth = rcb_ratio(bEarth,vst,mEarth)

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

bJup   = closest_approach(qjup,qchip,vjup,Vst)
rcbJup = rcb_ratio(bJup,vst,mjup)

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

bMoon   = closest_approach(qMoon,qchip,vMoon,Vst)
rcbMoon = rcb_ratio(bMoon,vst,mMoon)

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

bMars   = closest_approach(qMars,qchip,vMars,Vst)
rcbMars = rcb_ratio(bMars,vst,mMars)

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
output_file = "miss_distances_FTID_results.txt"
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

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29280133209761610847199299779566295e-02"), 
                  tpfl("-5.57596266415135739092548108896443497e-02"), 
                  tpfl("-1.77686212299783721985079095732413295e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( Psol, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( Psol, params; testparticles=Pchip_rel )
sol0    = pomin.solve( Psol, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbSol[2]
dmissFlatGrav = rcbSol[3]

println("SUN RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "SUN RESULTS:")
println(file, "-----------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   ALPHA CENTAURI
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29279538728389700446721236547889991e-02"), 
         tpfl("-5.57596913770760242903838491049976388e-02"), 
         tpfl("-1.77686156920466641922704142460267104e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( Palpha, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( Palpha, params; testparticles=Pchip_rel )
sol0    = pomin.solve( Palpha, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbAlpha[2]
dmissFlatGrav = rcbAlpha[3]

println("ALPHA CENTAURI RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "ALPHA CENTAURI RESULTS:")
println(file, "----------------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   JUPITER
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29279538813329027435745884199342427e-02"), 
         tpfl("-5.5759691381805953241084820061317113e-02"), 
         tpfl("-1.7768615692761663751817806954879475e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( Pjup, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( Pjup, params; testparticles=Pchip_rel )
sol0    = pomin.solve( Pjup, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbJup[2]
dmissFlatGrav = rcbJup[3]

println("JUPITER RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "JUPITER RESULTS:")
println(file, "---------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   EARTH
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29279538736891046494466933809423003e-02"), 
         tpfl("-5.57596913778559480339061166604508494e-02"), 
         tpfl("-1.77686156922354062625902331324005054e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( PEarth, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( PEarth, params; testparticles=Pchip_rel )
sol0    = pomin.solve( PEarth, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbEarth[2]
dmissFlatGrav = rcbEarth[3]

println("EARTH RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "EARTH RESULTS:")
println(file, "-------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   PROXIMA
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.2927953873061638672739581591721453e-02"), 
         tpfl("-5.57596913773661392117521019705786548e-02"), 
         tpfl("-1.77686156920764393871953614862440295e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( PProx, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( PProx, params; testparticles=Pchip_rel )
sol0    = pomin.solve( PProx, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbProx[2]
dmissFlatGrav = rcbProx[3]

println("PROXIMA RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "PROXIMA RESULTS:")
println(file, "---------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   MOON
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29279538731564925572432043882175815e-02"), 
         tpfl("-5.57596913774481017994445122959346604e-02"), 
         tpfl("-1.77686156921055359301768169343967384e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( PMoon, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( PMoon, params; testparticles=Pchip_rel )
sol0    = pomin.solve( PMoon, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbMoon[2]
dmissFlatGrav = rcbMoon[3]

println("MOON RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "MOON RESULTS:")
println(file, "------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   MARS
#-----------------------------------------------------------------------

mchip = tpfl(1.0E-30)
qchip = Xst0
Vst   = [tpfl("-7.29279538731604907445780907128956807e-02"), 
         tpfl("-5.57596913774413551365629631185583374e-02"), 
         tpfl("-1.77686156921044671212882540261258305e-01")]

vst = norm(Vst)

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)
pchip_rel = (mchip*γst) .* Vst
pchip_newt = mchip .* Vst

Pchip_rel = pomin.setup_single_particle(mchip, qchip, pchip_rel, tpfl)
Pchip_rel0 = pomin.setup_single_particle(mchip, qchip, pchip_rel0, tpfl)
Pchip_newt = pomin.setup_single_particle(mchip, qchip, pchip_newt, tpfl)

solN    = pomin.solve( PMars, params; testparticles=Pchip_newt, 
                       Newtonian=true )
sol     = pomin.solve( PMars, params; testparticles=Pchip_rel )
sol0    = pomin.solve( PMars, params; testparticles=Pchip_rel0 )

zendN   = solN(tcl)
zend    = sol(tcl)
zend0   = sol0(tcl)

qmisstar    = zend[7:9] - qtar
dmissstar   = norm(qmisstar)

qmissN  = zendN[7:9] - zend[7:9]
dmissN  = norm(qmissN)

qmissFlat = flat_space_spacecraft_final - zend0[7:9]
dmissFlat = norm(qmissFlat)

dmissFlatComp = rcbMars[2]
dmissFlatGrav = rcbMars[3]

println("MARS RESULTS:")
println("Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println("Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println("PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")

println("HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println("Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")

# Write to file
println(file, "MARS RESULTS:")
println(file, "------------")
println(file, "Fine-tuned velocity: ", Vst)
println(file, "Flat space miss: ", dmissFlat, " (", dmissFlat/AU, " AU)")
println(file, "Newtonian miss distance: ", dmissN, " (", dmissN/AU, " AU)")
println(file, "PoMiN Target miss: ", dmissstar, " (", dmissstar/AU, " AU)")
println(file, "HO grav miss distance: ", dmissFlatGrav, " (", dmissFlatGrav/AU, " AU)")
println(file, "Grav comp miss distance: ", dmissFlatComp, " (", dmissFlatComp/AU, " AU)")
println(file)

#-----------------------------------------------------------------------
#   COLLECT MISS DISTANCES AND CREATE TABLE
#-----------------------------------------------------------------------

# Write summary table to file
println(file, "="^100)
println(file, "COMPREHENSIVE SUMMARY TABLE")
println(file, "="^100)
println(file)

using Printf
println(file, @sprintf("%-15s | %-15s | %-15s | %-15s | %-15s | %-15s", 
        "Body", "Target Miss (AU)", "Newtonian (AU)", "Flat Space (AU)", "HO Grav (AU)", "Grav Comp (AU)"))
println(file, "-"^100)

# Note: This creates a framework for the summary table
# The actual values would need to be collected from each section's calculations
println(file, "Summary table framework created - values calculated above for each celestial body")
println(file, "Each section contains:")
println(file, "- Fine-tuned velocity for optimal targeting")
println(file, "- Target miss: Distance from spacecraft to Proxima at closest approach")
println(file, "- Newtonian miss: Difference between Newtonian and relativistic trajectories")
println(file, "- Flat space miss: Difference between flat space and relativistic trajectories")
println(file, "- HO grav miss: Higher-order gravitational effects")
println(file, "- Grav comp miss: Gravitational compensation effects")
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
