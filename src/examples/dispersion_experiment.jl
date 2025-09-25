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
# qchip = Double64[-2.23588844590237e7, 8.68066486466342e7, 2.9882578786605e7]

Vchip_rel_FT  = [tpfl("-7.29280133630346695175238301027084351e-02"), 
                 tpfl("-5.57596267171519947482140178590898884e-02"), 
                 tpfl("-1.77686212323486749509362995970259857e-01")]
Vchip_newt_FT = [tpfl("-7.2928011848008478283629726326358021e-02"), 
                 tpfl("-5.57596297969013443215818900874537475e-02"), 
                 tpfl("-1.77686212073603716257005595684594549e-01")]

# Final target position minus initial spacecraft position

ΔxN = [tpfl("-9.90791377535429553805304070804749628e+12"),
       tpfl("-7.57545201007502618274337340834229983e+12"),
       tpfl("-2.41402511341614011030784420306057395e+13")]

ΔxR = [tpfl("-9.90791377535429379640585529915263168e+12"),
       tpfl("-7.57545201007501210952698262636613967e+12"),
       tpfl("-2.41402511341613904686121086679115484e+13")]


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
# println("Angle between Sun-Jup and Sun-spacecraft = ", Jupiter_Spacecraft_angle, " degs")

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
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(PProx, Psol, Palpha, Pjup, PEarth)

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
# target_distance = norm(ΔxR)   # for relativistic
target_distance = norm(ΔxN)     # for Newtonian
target_distance_in_m = target_distance * Msol2m

# size_of_target_disk_in_m = tpfl(5E8)
# theta_tol = rad2deg(size_of_target_disk_in_m / target_distance_in_m)

theta_tol = tpfl(1E-7)
size_of_target_disk_in_m = deg2rad(theta_tol) * target_distance_in_m

# baseVector = Vchip_rel_FT     # for relativistic
baseVector = Vchip_newt_FT      # for Newtonian

println("Using theta_tol = ",theta_tol)
println("Size of target disk in m = ",size_of_target_disk_in_m)

N = 20



#-----------------------------------------------------------------------
#   RUN EXPERIMENT
#-----------------------------------------------------------------------

for i = 1:N
    println("\n\nTrial #",i)

    (theta_deg, init_unit_vec) = generateUnitVectorWithinToleranceAngle(theta_tol, baseVector, target_distance, tpfl)

    println("Initial angle from base vector = ",theta_deg, " degs")

    # set length of init velocity vector to be length of base velocity vector
    init_vel = init_unit_vec * norm(baseVector)

    # Set up starchip
    γst = one(tpfl) / sqrt(one(tpfl) - (norm(init_vel))^2)
    # pchip = mchip * γst .* init_vel     # for relativistic
    pchip = mchip * init_vel     # for Newtonian
    Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)
    test_particle_system = Pchip

    # run solver
    println("Running solver...")
    # sol = pomin.solve(main_particle_system, params; testparticles=test_particle_system)     # for relativistic
    sol = pomin.solve(main_particle_system, params; testparticles=test_particle_system, Newtonian=true)     # for Newtonian

    Zend = sol(time_closest_approach)
    
    nInt = length(main_particle_system.m)
    ZendInt = Zend[1:6*nInt]
    ZendTest = Zend[6*nInt+1:end]

    q_starchip = ZendTest[1:3]
    q_proxima = ZendInt[1:3]
    q_target = q_proxima + bv
    miss_dist = norm(q_starchip - q_target) / AU

    println("Miss distance: ",miss_dist," AU")

    expected_geometric_miss_dist = deg2rad(theta_deg) * target_distance / AU   # expected miss distance in AU

    println("Expected miss distance based on geometry: ",expected_geometric_miss_dist," AU")

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