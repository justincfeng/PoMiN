include("../pomin.jl")
using .pomin
using DelimitedFiles

include("../core/initial_data/idgen.jl")

include("target_proxima.jl")

q_starchip_Func = pfs[1]
q_target_Func = pfs[3]
time_closest_approach = pfs[4]


#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT (EXCEPT MOMENTUM)
#-----------------------------------------------------------------------

# Mass and position
mchip = tpfl(1.0E-30)
qchip = Xst0
γst = one(tpfl) / sqrt(one(tpfl) - tpfl(0.2)^2)


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
#   INITIAL DATA SETUP FOR SUN
#-----------------------------------------------------------------------

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR JUPITER (astropy)
#-----------------------------------------------------------------------

mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([3.99442362 * AU, 2.73345644 * AU, 1.07451704 * AU])  # Jupiter position in geometric units (from astropy, J2000)

vjup = tpfl.([-7.8875477321829800, 1.0175858402135900E+01, 4.5538873116399600])  # in km/s, from astropy, J2000

vjup *= 1000 / cMKS  # convert from km/s to units of c

γjup = one(tpfl) / sqrt(one(tpfl) - norm(vjup)^2)  # Lorentz factor
pjup = tpfl.(mjup * γjup * vjup)  # Jupiter momentum 

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)


Jupiter_Spacecraft_angle = rad2deg(acos(dot(qjup, qchip) / (norm(qjup) * norm(qchip))))
println("Angle between Sun-Jup and Sun-spacecraft = ", Jupiter_Spacecraft_angle, " degs")


#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR ALPHA CENTAURI A + B
#-----------------------------------------------------------------------

malpha = 2.0429

# Alpha A+B barycenter position and velocity given by Kervella et al 2017
qalpha = tpfl.([-1.045245216607860E+13, -8.74487747435090E+12, -2.442178248634550E+13])         # in geometric units
palpha = tpfl.([-6.363558578130740E-05, 1.508126516235040E-04, 1.475021586009610E-04])          # in geometric units

Palpha = pomin.setup_single_particle(malpha, qalpha, palpha, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR EARTH
#-----------------------------------------------------------------------

mEarth = 3.0033693739E-06

# Earth coordinates generated using astropy with epoch J2000
qEarth = tpfl.([-1.8667826140E+07, 8.9633743993E+07, 3.8883291714E+07])         # geometric units
pEarth = tpfl.([-2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11])       # geometric units

PEarth = pomin.setup_single_particle(mEarth, qEarth, pEarth, tpfl)

#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(Psol, PProx, Palpha, PEarth, Pjup)


#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tcl * tpfl(1.1)) # 22 years
δ = tpfl(7.31E+8) # about one hour
tols = tpfl(1e-16)

params = pomin.ParametersJulia(tspan, integrator="Vern9",
    atol=tols, rtol=tols)



#-----------------------------------------------------------------------
#   EXPERIMENT INPUT PARAMETERS
#-----------------------------------------------------------------------
theta_tol = 0.00001
N = 2

#-----------------------------------------------------------------------
#   RUN EXPERIMENT
#-----------------------------------------------------------------------

for i = 1:N
    println("\n\nTrial #",i)

    (theta_deg, init_unit_vec) = generateUnitVectorWithinToleranceAngle(theta_tol, q_starchip_Func(0), q_target_Func(time_closest_approach), tpfl)

    println("Initial angle from base vector = ",theta_deg, " degs")

    # set length of init velocity vector to be 0.2 c
    init_vel = 0.2 * init_unit_vec

    # Set up starchip
    pchip = (mchip * γst) .* (init_vel)
    Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)
    test_particle_system = Pchip

    # run solver
    println("Running solver...")
    sol = pomin.solveT(main_particle_system, params, test_particle_system)

    Zend = sol(time_closest_approach)
    q_starchip = Zend[31:33]
    q_proxima = Zend[4:6]
    q_target = q_proxima + bv
    miss_dist = norm(q_starchip - q_target) / AU

    println("Miss distance: ",miss_dist," AU")

    expected_geometric_miss_dist = deg2rad(theta_deg) * 268210   # expected miss distance in AU

    println("Expected miss distance based on geometry: ",expected_geometric_miss_dist)

    if i == 1
        open("miss_distances.csv", "w") do io
            writedlm(io, [theta_deg miss_dist expected_geometric_miss_dist], ',')
        end
    else
        open("miss_distances.csv", "a") do io
            writedlm(io, [ theta_deg miss_dist expected_geometric_miss_dist], ',')
        end
    end

end