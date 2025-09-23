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

vjup *= tpfl(1000 / cMKS)  # convert from km/s to units of c

γjup = one(tpfl) / sqrt(one(tpfl) - norm(vjup)^2)  # Lorentz factor
pjup = tpfl.(mjup * γjup .* vjup)  # Jupiter momentum 
pjupN = mjup .* vjup                # Jupiter momentum (Newtonian)

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)
PjupN = pomin.setup_single_particle(mjup, qjup, pjupN, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MOON
#-----------------------------------------------------------------------

mMoon = tpfl(7.34767309E22 / 1.988416E30)     # in units of solar masses

# Moon coordinates generated using astropy with epoch J2000
qMoon = tpfl.([-1.886529874442160E+07, 
                8.94531273942814E+07, 
                3.883175807801640E+07])    # geometric units
vMoon = tpfl.([-2.9141440121218000E+01, 
                -5.6958177245812600, 
                -2.4819741171379700]) 

vMoon *= tpfl(1000 / cMKS)  # convert from km/s to units of c

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

vMars *= tpfl(1000 / cMKS)  # convert from km/s to units of c

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

#Pint = pomin.merge_particle_systems(PProx, Psol)
#PintN = pomin.merge_particle_systems(PProxN, Psol)
#Pint = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup)
#PintN = pomin.merge_particle_systems(PProxN, Psol, PalphaN, PjupN)
Pint = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup, PEarth)
PintN = pomin.merge_particle_systems(PProxN, Psol, PalphaN, PjupN, PEarthN)
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

# Final target position minus initial spacecraft position
ΔxN     = qtarN - qchip

println("Newtonian miss: ", qmissN)
println("Newtonian miss distance: ", dmissN)
println("Newtonian final target position minus initial spacecraft position: ", ΔxN)

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

# Final target position minus initial spacecraft position
ΔxR     = qtarR - qchip

println("Relativistic miss: ", qmissR)
println("Relativistic miss distance: ", dmissR)
println("Relativistic final target position minus initial spacecraft position: ", ΔxR)

#-----------------------------------------------------------------------
#   SAVE FINAL TARGET POSITIONS TO FILE
#-----------------------------------------------------------------------

println("Saving final target positions to Final_target_minus_initial_spacecraft.txt...")

open("Final_target_minus_initial_spacecraft.txt", "w") do file
    println(file, "# Final target position minus initial spacecraft position")
    println(file, "# Generated: $(Dates.now())")
    println(file, "")
    println(file, "ΔxN = [tpfl(\"$(ΔxN[1])\"),")
    println(file, "       tpfl(\"$(ΔxN[2])\"),")
    println(file, "       tpfl(\"$(ΔxN[3])\")]")
    println(file, "")
    println(file, "ΔxR = [tpfl(\"$(ΔxR[1])\"),")
    println(file, "       tpfl(\"$(ΔxR[2])\"),")
    println(file, "       tpfl(\"$(ΔxR[3])\")]")
end

println("Final target positions saved to Final_target_minus_initial_spacecraft.txt")
