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
#   INITIAL DATA SETUP FOR ALPHA CENTAURI A + B
#-----------------------------------------------------------------------

malpha = tpfl(2.0429)

# Alpha A+B barycenter position and velocity given by Kervella et al 2017
qalpha = tpfl.([-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13])         # in geometric units
palpha = tpfl.([-6.363558578130740E-05, 1.508126516235040E-04, 1.475021586009610E-04])          # in geometric units

Palpha = pomin.setup_single_particle(malpha, qalpha, palpha, tpfl)


#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(Palpha)

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
