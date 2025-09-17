include("../pomin.jl")
using .pomin
using DelimitedFiles

include("target_proxima.jl")

q_starchip_Func = pfs[1]
q_target_Func = pfs[3]
time_closest_approach = pfs[4]


#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT 
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0

vst = norm(Vst)
γst = one(tpfl) / sqrt(one(tpfl) - (vst / c)^2)

# Momentum scaled by same factor as mass to preserve velocity
pchip = (mchip * γst) .* (Vst)

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
#   INITIAL DATA SETUP FOR MARS
#-----------------------------------------------------------------------

mMars = tpfl(3.2271058587E-07)  # Mars mass in solar masses
qMars = tpfl.([1.38356874 * AU, -0.00120916 * AU, -0.03786079 * AU])  # Mars position in geometric units (from astropy, J2000)

vMars = tpfl.([1.17347755650600, 2.390740193498500E+01, 1.0934201867474400E+01])  # in km/s, from astropy, J2000

vMars *= 1000 / cMKS  # convert from km/s to units of c

γMars = one(tpfl) / sqrt(one(tpfl) - norm(vMars)^2)  # Lorentz factor
pMars = tpfl.(mMars * γMars * vMars)  # Mars momentum 

PMars = pomin.setup_single_particle(mMars, qMars, pMars, tpfl)

Mars_Spacecraft_angle = rad2deg(acos(dot(qMars,qchip)/(norm(qMars)*norm(qchip))))
println("Angle between Sun-Mars and Sun-spacecraft = ",Mars_Spacecraft_angle," degs")

#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(PMars)

#-----------------------------------------------------------------------
#   SETUP TEST PARTICLE SYSTEM
#-----------------------------------------------------------------------

test_particle_system = pomin.merge_particle_systems(Pchip, PProx)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tcl * tpfl(1.1)) # 22 years
δ = tpfl(7.31E+8) # about one hour
tols = tpfl(1e-16)

params = pomin.ParametersJulia(tspan, integrator="Vern9",
    atol=tols, rtol=tols)

#-----------------------------------------------------------------------
#   RUN SOLVER
#-----------------------------------------------------------------------

println("Running solver...")
sol = pomin.solveT(main_particle_system, params, test_particle_system)

println("Time of closest approach = ",time_closest_approach)

Zend = sol(time_closest_approach)

# structure of Zend: [ q_int[1], q_int[2], q_int[3], p_int[1], p_int[2], p_int[3], q_chip[1], q_chip[2], q_chip[3], q_Prox[1], q_Prox[2], q_Prox[3],
#                      p_chip[1], p_chip[2], p_chip[3], p_Prox[1], p_Prox[2], p_Prox[3]  ]
q_starchip = Zend[7:9]
q_proxima = Zend[10:12]
q_target = q_proxima + bv
miss_dist = norm(q_starchip - q_target) / AU

println("Miss distance: ",miss_dist," AU")
