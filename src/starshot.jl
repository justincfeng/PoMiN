using DoubleFloats
using Random
using Distributions

include("broyden.jl")

function closestApproach(v_init::RealVec)
    println(stderr,"closestApproach() called")
    # Mass of Sun, Earth, Spacechip, Proxima, Alpha A, Alpha B in M
    masses = Double64[1, 3.0033693739E-06, 1.0057832537E-33, 0.1221, 1.0788, 0.9092] 

    # Make sure spacecraft velocity has length of 0.2
    v_init = v_init / norm(v_init) * 0.2
    # Find spacecraft momentum p = gamma m v
    p_init = masses[3]*v_init/sqrt(1-0.2^2)

    # Initial Z:
    Z_init = [  0, 0, 0,                                                    # Sun pos
                -1.8667826140E+07, 8.9633743993E+07, 3.8883291714E+07,      # Earth pos
                -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,      # Spacechip pos (Moon dist (semimajor axis) from Earth)
                -9.9157694997E+12, -7.5891180191E+12, -2.4171267258E+13,    # Proxima pos
                -1.0514396564E+13, -8.7920419314E+12, -2.4558076740E+13,    # Alpha A pos
                -1.0514148951E+13, -8.7899759755E+12, -2.4558922285E+13,    # Alpha B pos
                0, 0, 0,                                                    # Sun momentum
                -2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11,    # Earth momentum
                p_init[1], p_init[2], p_init[3],                            # Spacechip momentum
                -3.8172314441E-06, 9.0494364315E-06, 8.9932012189E-06,      # Proxima momentum
                -3.2681576581E-05, 8.2784394610E-05, 7.2543810801E-05,      # Alpha A momentum
                -2.9658611095E-05, 6.6362559735E-05, 6.7386459079E-05       # Alpha B momentum
                ]

    tspan = (Double64(0), Double64(1.6018E+14)) # 25 years
    # tspan = (Double64(0), Double64(1.922E+11)) # 10 days
    δ = Double64(1.922E+10) # about one day

    no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
    adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
    maxit = 100000

    sol = hrkintegrator(3, 6, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)
    # sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan)

    println(stderr,"integrator done")
    #printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
    minvec = zeros(3)
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # get q vector for particle 3 (spacecraft) at timestep tn
        q1 = Z2q(sol.N, sol.d, 3, sol.z[tn])
        # get q vector for particle 4 (Proxima Centauri) at timestep tn
        q2 = Z2q(sol.N, sol.d, 4, sol.z[tn])
        # get distance and save it if it's the minimum distance
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end

function closestApproach_2body_Proxima_test_particle()
    println(stderr,"closestApproach_2body_Proxima_test_particle() called")

    # Mass of Spacechip, Proxima test particle
    masses = Double64[1.0057832537E-33, 1E-40] 

    # Initial Z:
    Z_init = [  -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,                  # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
                -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos (Kervella 2017)
                -7.486220592724260E-35, -5.723861062118610E-35, -1.823989847481890E-34, # Spacechip momentum (misses Proxima by 0.0009 AU in absence of gravity)
                -3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45       # Proxima momentum (Kervella 2017)
                ]

    tspan = (Double64(0), Double64(1.410E+14)) # 22 years
    δ = Double64(1E7) # less than a minute
    start_dist = 2.718112140945220E+13  # Proxima-Spacechip distance at time 0
    min_dist = 1E10 # about 100 AU -- once within this distance it will use min_dt
    max_dt = Double64(7.3E8) # about an hour
    min_dt = Double64(1E7) # less than a minute
    adapt = (δ, z, zdot) -> dt_lin_int(δ, z, zdot, start_dist, min_dist, max_dt, min_dt)
    no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
    maxit = 15000000

    sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
    # sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan)

    println(stderr,"integrator done")
    printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
    minvec = zeros(3)
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # get q vector for particle 1 (spacecraft) at timestep tn
        q1 = Z2q(sol.N, sol.d, 1, sol.z[tn])
        # get q vector for particle 2 (Proxima Centauri) at timestep tn
        q2 = Z2q(sol.N, sol.d, 2, sol.z[tn])
        # get distance and save it if it's the minimum distance
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end

function closestApproach_2body_w_MW()
    println(stderr,"closestApproach_2body_w_MW() called")

    # Mass of Spacechip, Proxima test particle
    masses = Double64[1.0057832537E-33, 1E-40] 

    # Initial Z:
    Z_init = [  -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,                  # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
                -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos (Kervella 2017)
                -7.486220592724260E-35, -5.723861062118610E-35, -1.823989847481890E-34, # Spacechip momentum (misses Proxima by 0.0009 AU in absence of gravity)
                -3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45       # Proxima momentum (Kervella 2017)
                ]

    tspan = (Double64(0), Double64(1.410E+14)) # 22 years
    δ = Double64(7.31E+8) # about one hour
    start_dist = 2.718112140945220E+13  # Proxima-Spacechip distance at time 0
    min_dist = 1E10 # about 100 AU -- once within this distance it will use min_dt
    max_dt = Double64(7.3E8) # about an hour
    min_dt = Double64(1E7) # less than a minute
    adapt = (δ, z, zdot) -> dt_lin_int(δ, z, zdot, start_dist, min_dist, max_dt, min_dt)
    no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
    maxit = 10000000

    sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)

    println(stderr,"integrator done")
    printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
    minvec = zeros(3)
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # get q vector for particle 1 (spacecraft) at timestep tn
        q1 = Z2q(sol.N, sol.d, 1, sol.z[tn])
        # get q vector for particle 2 (Proxima Centauri) at timestep tn
        q2 = Z2q(sol.N, sol.d, 2, sol.z[tn])
        # get distance and save it if it's the minimum distance
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end

function closestApproach_2body_massive_Proxima()
    println(stderr,"closestApproach_2body_massive_Proxima() called")

    # Mass of Spacechip, Proxima
    masses = Double64[1.0057832537E-33, ] 

    # Initial Z:
    Z_init = [  -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,                  # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
                -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos (Kervella 2017)
                -7.486220592724260E-35, -5.723861062118610E-35, -1.823989847481890E-34, # Spacechip momentum (misses Proxima by 0.0009 AU in absence of gravity)
                       # Proxima momentum (Kervella 2017)
                ]

    tspan = (Double64(0), Double64(1.410E+14)) # 22 years
    δ = Double64(1E7) # less than a minute
    start_dist = 2.718112140945220E+13  # Proxima-Spacechip distance at time 0
    min_dist = 1E10 # about 100 AU -- once within this distance it will use min_dt
    max_dt = Double64(7.3E8) # about an hour
    min_dt = Double64(1E7) # less than a minute
    adapt = (δ, z, zdot) -> dt_lin_int(δ, z, zdot, start_dist, min_dist, max_dt, min_dt)
    maxit = 10000000

    sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, adapt, tspan, maxit)
    # sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan)

    println(stderr,"integrator done")
    printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
    minvec = zeros(3)
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # get q vector for particle 1 (spacecraft) at timestep tn
        q1 = Z2q(sol.N, sol.d, 1, sol.z[tn])
        # get q vector for particle 2 (Proxima Centauri) at timestep tn
        q2 = Z2q(sol.N, sol.d, 2, sol.z[tn])
        # get distance and save it if it's the minimum distance
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end


function closestApproach_2plus1(mass::Double64, q_init::RealVec, p_init::RealVec)
    println(stderr,"closestApproach_2plus1() called")

    # Mass of Spacechip, Proxima test particle, and 3rd body
    masses = Double64[1.0057832537E-33, 1E-40, mass] 

    # Initial Z:
    Z_init = [  -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,                  # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
                -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos (Kervella 2017)
                q_init[1], q_init[2], q_init[3],                                        # 3rd body pos
                -7.486220592724260E-35, -5.723861062118610E-35, -1.823989847481890E-34, # Spacechip momentum (misses Proxima by 0.0009 AU in absence of gravity)
                -3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45,      # Proxima momentum (Kervella 2017)
                p_init[1], p_init[2], p_init[3]                                         # 3rd body momentum
                ]

    tspan = (Double64(0), Double64(1.410E+14)) # 22 years
    δ = Double64(7.31E+8) # about one hour

    no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
    adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
    maxit = 1000000

    sol = hrkintegrator(3, 3, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
    # sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan)

    println(stderr,"integrator done")
    printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
    minvec = zeros(3)
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # get q vector for particle 1 (spacecraft) at timestep tn
        q1 = Z2q(sol.N, sol.d, 1, sol.z[tn])
        # get q vector for particle 2 (Proxima Centauri) at timestep tn
        q2 = Z2q(sol.N, sol.d, 2, sol.z[tn])
        # get distance and save it if it's the minimum distance
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end

function generateInitVelocityDiskProxima(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using random points in a disk centered on Proxima")

    ProximaDistance = norm(Double64[-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]    # Proxima pos from Kervella 2017
                           -
                           Double64[-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]                # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
    )
    diskRadius = ProximaDistance * Double64(tan(deg2rad(theta_tol)))
    println(stderr, "proxima distance = ", ProximaDistance * 1.47669196951425 * 6.6845871226706E-09, " AU")
    println(stderr, "disk Radius = ", diskRadius * 1.47669196951425 * 6.6845871226706E-09, " AU")

    U = baseVector / norm(baseVector)

    v = nothing
    while true
        # generate vector that is linearly independent of U
        Y = [rand(Double64), rand(Double64), rand(Double64)]
        # ensure Y is not collinear with U
        while abs(dot(U,Y)/(norm(U)*norm(Y))) == 1
            Y = [rand(Double64), rand(Double64), rand(Double64)]
        end

        # make a vector v that is orthogonal to U by subtracting the part of Y that is parallel to U
        a = dot(Y,U)/norm(U)^2
        v = Y - a * U
        println(stderr,"dot(v,U) = ",dot(v,U))

        # ensure v and Y are orthogonal, otherwise repeat
        if dot(v,U) == 0
            break
        end
    end

    # find a vector w that's orthogonal to U and v
    w = cross(U,v)

    # scale v and w to match the disk radius
    V = v / norm(v) * diskRadius
    W = w / norm(w) * diskRadius

    # are U, V, and W all mutually orthogonal?
    println(stderr, "U dot V = ", dot(U, V))
    println(stderr, "U dot W = ", dot(U, W))
    println(stderr, "V dot W = ", dot(V, W))

    # generate random point (x,y) in the unit disk
    x = nothing
    y = nothing
    while true
        x = rand(Double64)
        y = rand(Double64)
        if norm([x,y]) <= diskRadius
            break
        end
    end

    # construct vector Z in the disk centered on Proxima
    Z = x * V + y * W

    # construct init velocity vector via vector addition of U that reaches Proxima plus Z
    init_vel = U * ProximaDistance + Z

    # make init velocity a unit vector
    init_vel = init_vel / norm(init_vel)

    # find angle between baseVector and init_vel
    cosine_theta = dot(init_vel, baseVector) / norm(baseVector)     # init_vel already unit length
    theta_rad = acos(cosine_theta)
    theta_deg = rad2deg(theta_rad)

    println(stderr, "For initial velocity ", init_vel, " angle with base vector is ", theta_deg)

    return theta_deg, init_vel

end

function generateInitVelocityMonteCarloDiskOrigin(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method with disk centered on origin and then rotated to proper orientation and translated to Proxima")

    # compute 1/2 of side length of square using theta_tol
    ProximaDistance = norm(Double64[-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]    # Proxima pos from Kervella 2017
                           -
                           Double64[-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]                # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
    )
    sideLength = ProximaDistance * tan(deg2rad(theta_tol))
    println(stderr, "proxima distance = ", ProximaDistance * 1.47669196951425 * 6.6845871226706E-09, " AU")
    println(stderr, "sidelength = ", sideLength * 1.47669196951425 * 6.6845871226706E-09, " AU")

    # find random point in x-y plane that is uniformly sampled from a disk with radius sideLength
    x = nothing
    y = nothing
    while true
        x = sideLength * (rand() * 2 - 1)
        y = sideLength * (rand() * 2 - 1)
        if norm([x, y]) <= sideLength
            break
        end
    end
    z = 0
    println(stderr, "x, y = ",x,", ",y," where sideLength = ",sideLength," M")

    # now we need to rotate the disk to align its normal vector with baseVector

    # negative of the normalized baseVector will be new normal vector
    new_normal = -baseVector / norm(baseVector)

    # current normal vector points along z-axis
    curr_normal = [0, 0, 1]

    # cos of rotation angle comes from dot product
    cos_theta = dot(new_normal, curr_normal)        # both vectors are unit length
    println(stderr,"cos_theta of rotation = ",cos_theta)

    # axis of rotation is perpendicular to both normals
    rotation_axis = cross(new_normal, curr_normal)
    rotation_axis = rotation_axis / norm(rotation_axis)

    # construct Rodrigues rotation matrix
    sin_theta = sqrt(1 - cos_theta * cos_theta)
    r = rotation_axis   # shorthand
    c = cos_theta   # shorthand
    s = sin_theta   # shorthand
    R = [   r[1] * r[1] * (1-c) + c         r[1] * r[2] * (1-c) - r[3]*s     r[1] * r[3] * (1-c) + r[2]*s    ;
            r[2] * r[1] * (1-c) + r[3]*s    r[2] * r[2] * (1-c) + c          r[2] * r[3] * (1-c) - r[1]*s    ;
            r[3] * r[1] * (1-c) - r[2]*s    r[3] * r[2] * (1-c) + r[1]*s     r[3] * r[3] * (1-c) + c    
    ]

    # rotate x, y, z
    X_old = Double64[x; y; z]       # put random values in a column vector
    println(stderr,"X_old = ", X_old)
    X_rot = R * X_old       # multiply by rotation matrix
    println(stderr, "X_rot = ", X_rot)
    println(stderr, "Rotated through cos_theta = ", dot(X_old,X_rot)/(norm(X_old)*norm(X_rot)))
    # compare new vector and baseVector -- they should be parallel
    println(stderr,"Norm of X_rot cross new_normal (should be 0): ",norm(cross(X_rot,new_normal)))

    # disk will be translated so that it's centered on (x_0, y_0, z_0) which is location of Proxima
    x_0 = -9.90338347925777E+12
    y_0 = -7.58520275447461E+12
    z_0 = -2.41494965557852E+13

    # translate x, y, z and transpose result back to a row vector
    X_new = [x_0 + X_rot[1], y_0 + X_rot[2], z_0 + X_rot[3]]

    # construct init velocity vector via vector addition of baseVector that reaches Proxima plus new vector that lies in the rotated & translated disk
    init_vel = baseVector / norm(baseVector) * ProximaDistance + X_new

    # make init velocity a unit vector
    init_vel = init_vel / norm(init_vel)

    # find angle between baseVector and init_vel
    cosine_theta = dot(init_vel, baseVector) / norm(baseVector)     # init_vel already unit length
    theta_rad = acos(cosine_theta)
    theta_deg = rad2deg(theta_rad)

    println(stderr, "For initial velocity ", init_vel, " angle with base vector is ", theta_deg)

    return theta_deg, init_vel

end

function generateInitVelocityMonteCarloDiskProxima(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method with disk centered on Proxima")

    # compute side length of square using theta_tol
    ProximaDistance = norm([-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]    # Proxima pos from Kervella 2017
                           - [-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]                # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
    )
    sideLength = ProximaDistance * tan(deg2rad(theta_tol))
    println("sidelength = ", sideLength * 1.47669196951425 * 6.6845871226706E-09, " AU")

    # [A,B,C] will be the normal to the plane
    A = baseVector[1] / norm(baseVector)
    B = baseVector[2] / norm(baseVector)
    C = baseVector[3] / norm(baseVector)

    # square will be centered on (x_0, y_0, z_0) which is location of Proxima
    x_0 = -9.90338347925777E+12
    y_0 = -7.58520275447461E+12
    z_0 = -2.41494965557852E+13

    # select random point inside square and discard any that don't lie within the inscribed disk
    x = 1
    y = 1
    while true
        x = rand(Uniform(-sideLength, sideLength))
        y = rand(Uniform(-sideLength, sideLength))
        if norm([x, y]) <= sideLength
            break
        end
    end
    # translate random point to be inside square centered on Proxima
    x += x_0
    y += y_0

    # find z such that (x,y,z) lies on the plane centered on Proxima and normal to baseVector
    z = z_0 + (A * (x_0 - x) + B * (y_0 - y)) / C

    # construct init velocity vector
    init_vel = baseVector / norm(baseVector) * ProximaDistance + [x, y, z]

    # make init velocity a unit vector
    init_vel = init_vel / norm(init_vel)

    # find angle between baseVector and init_vel
    cosine_theta = dot(init_vel, baseVector) / norm(baseVector)     # init_vel already unit length
    theta_rad = acos(cosine_theta)
    theta_deg = rad2deg(theta_rad)

    println(stderr, "For initial velocity ", init_vel, " angle with base vector is ", theta_deg)

    return theta_deg, init_vel

end

function generateInitVelocityMonteCarloProxima(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method with sphere around Proxima")

    # adjust baseVector length to equal distance to Proxima in units of M (solar mass)
    ProximaDistance = norm(   [-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]    # Proxima pos from Kervella 2017
                            - [-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]                # Spacechip pos (Moon dist (semimajor axis) from Earth in +y dir)
                        ) 
    baseVector = baseVector / norm(baseVector) * ProximaDistance

    # compute side length of cube using theta_tol
    sideLength = ProximaDistance * tan(deg2rad(theta_tol))
    println("Sidelength = ", sideLength * 1.47669196951425 * 6.6845871226706E-09, " AU")

    # select random point inside cube and discard any that don't lie within the inscribed sphere
    x = 1
    y = 1
    z = 1
    while true
        x = rand(Uniform(-sideLength, sideLength))
        y = rand(Uniform(-sideLength, sideLength))
        z = rand(Uniform(-sideLength, sideLength))
        if norm([x,y,z]) <= sideLength
            break
        end
    end

    # init velocity vector is baseVector plus random vector inside sphere
    init_vel = baseVector + [x,y,z]

    # make init velocity a unit vector
    init_vel = init_vel / norm(init_vel)
    
    # find angle between baseVector and init_vel
    cosine_theta = dot(init_vel, baseVector) / norm(baseVector)     # init_vel already unit length
    theta_rad = acos(cosine_theta)
    theta_deg = rad2deg(theta_rad)

    println(stderr, "For initial velocity ", init_vel, " angle with base vector is ", theta_deg)

    return theta_deg, init_vel
end

function generateInitVelocityMonteCarloLaunchPoint(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method with sphere around launch point")

    # normalize baseVector
    baseVector = baseVector / norm(baseVector)

    i = 0
    while true
        # generate random point on unit sphere using method of Marsaglia
        x = 1
        y = 1
        while true
            # pick x and y randomly from (-1,1) until x^2 + y^2 < 1
            x = rand() * 2 - 1
            y = rand() * 2 - 1
            if x^2 + y^2 < 1 
                break
            end
        end
        Xsquare = x^2 + y^2
        a = 2 * x * sqrt(1-Xsquare)
        b = 2 * y * sqrt(1-Xsquare)
        c = 1 - 2 * Xsquare
        rand_v = [a, b, c]
        i += 1

        # normalize random vector
        rand_v = rand_v / norm(rand_v)

        # find angle between baseVector and random vector
        cosine_theta = dot(rand_v,baseVector) # both unit vectors
        theta_rad = acos(cosine_theta)
        theta_deg = rad2deg(theta_rad)

        # compare to tolerance angle
        if theta_deg < theta_tol
            println(stderr,"For initial velocity ", rand_v, " angle with base vector is ", theta_deg, " (after ", i, " tries)")
            return theta_deg, rand_v
        end
    end
end

function generateInitVelocityBullseye(theta_tol, baseVector::RealVec)

    # use N rings
    N = 10
    num_per_ring = floor(1000/N)
    theta_1 = theta_tol / sqrt(N)

    # convert baseVector to spherical coordinates
    x_0 = baseVector[1]
    y_0 = baseVector[2]
    z_0 = baseVector[3]
    theta_0 = acos(z_0)
    phi_0 = sign(y_0) * acos(x_0 / sqrt(x_0^2 + y_0^2))
    
    # store sign of phi_0 for later
    sign_phi = sign(phi_0)

    # populate "bullseye" (central circle)
    for j in 1:num_per_ring
        
        # unfinished

    end

    # populate each ring
    for i in 2:N
        prev_theta = theta_1 * sqrt(i-1)
        next_theta = theta_1 * sqrt(i)
        rand_theta = prev_theta + rand() * (next_theta - prev_theta)

        # unfinished
    end

end


function generateInitVelocityHatBox(theta_tol, baseVector::RealVec)

    # generates initial velocity vector within theta_tol of baseVector
    # using method based on Archimedes' Hat-Box Theorem
    # see https://mathworld.wolfram.com/SpherePointPicking.html equations (1) and (2)
    # and https://mathworld.wolfram.com/ArchimedesHat-BoxTheorem.html 
    # and https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php 
    # and https://math.stackexchange.com/questions/711594/uniform-sampling-from-part-of-sphere-surface 

    # normalize baseVector
    baseVector = baseVector / norm(baseVector)

    # convert baseVector to spherical coordinates
    x_0 = baseVector[1]
    y_0 = baseVector[2]
    z_0 = baseVector[3]
    theta_0 = acos(z_0)
    phi_0 = sign(y_0) * acos(x_0 / sqrt(x_0^2 + y_0^2))
    
    # store sign of phi_0 for later
    sign_phi = sign(phi_0)

    # println("theta_0 = ", theta_0, " phi_0 = ", phi_0)

    # find theta min/max and phi min/max
    delta = deg2rad(theta_tol)
    theta_min = theta_0 - delta
    theta_max = theta_0 + delta
    phi_min = phi_0 - delta
    phi_max = phi_0 + delta

    # println("phi_min = ", phi_min, " phi_max = ", phi_max)

    # find u min/max and v min/max
    u_min = theta_min / (2*pi)
    u_max = theta_max / (2*pi)
    v_min = (cos(phi_min) + 1) / 2
    v_max = (cos(phi_max) + 1) / 2

    # println("u_min = ", u_min, " u_max = ", u_max)
    # println("v_min = ", v_min, " v_max = ", v_max)

    while true
            # putting this in a loop that repeats until angle is less than tolerance angle
            # because I'm sampling from unit square rather than unit disk, 
            # so there's a chance new vector could end up pointing to one of the corners inside unit square but outside unit disk

        # generate random u and v
        u = u_min + rand() * (u_max - u_min)
        v = v_min + rand() * (v_max - v_min)

        # println("u = ", u, " v = ", v)

        # find theta, phi
        theta = 2 * pi * u
        phi = acos(2*v - 1)

        # ensure phi is still in same quadrant as phi_0
        # NOTE: this won't work if phi ends up in a different quadrant than phi_0 after varying it by random angle!  But it works as long as random
        # angle is small enough (which it typically is in our case) that phi stays in the same quadrant as phi_0
        phi = sign_phi * abs(phi)

        # convert theta, phi to unit vector in cartesian coords
        x = sin(theta) * cos(phi)
        y = sin(theta) * sin(phi)
        z = cos(theta)
        rand_v = [x, y, z]

        # find angle between baseVector and random vector
        cosine_theta = dot(rand_v,baseVector) # both unit vectors
        theta_rad = acos(cosine_theta)
        theta_deg = rad2deg(theta_rad)

        # println("theta_deg = ", theta_deg)

        if theta_deg < theta_tol
            println(stderr,"For initial velocity ", rand_v, " angle with base vector is ", theta_deg)
            return theta_deg, rand_v
        end

    end
end

function testInitVector(theta_tol)

    
    baseVector = Double64[-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13]    # Proxima pos from Kervella 2017
                -Double64[-1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07]                  # Spacechip init pos (Moon dist (semimajor axis) from Earth in +y dir)

    # normalize baseVector
    baseVector = baseVector / norm(baseVector)

    # # test conversion to spherical and back

    # println("Before = ", baseVector)
    # # convert baseVector to spherical coordinates
    # x_0 = baseVector[1]
    # y_0 = baseVector[2]
    # z_0 = baseVector[3]
    # theta_0 = acos(z_0)
    # phi_0 = sign(y_0) * acos(x_0 / sqrt(x_0^2 + y_0^2))
    # # convert theta_0, phi_0 to unit vector in cartesian coords
    # x = sin(theta_0) * cos(phi_0)
    # y = sin(theta_0) * sin(phi_0)
    # z = cos(theta_0)
    # rand_v = [x, y, z]
    # println("After =  ", rand_v)

    # # test all octants

    # for i in 1:20
        # generateInitVelocity(theta_tol, Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01])
        # test all octants
        # generateInitVelocity(theta_tol, [1,1,1])
        # generateInitVelocity(theta_tol, [1,1,-1])
        # generateInitVelocity(theta_tol, [1,-1,1])
        # generateInitVelocity(theta_tol, [1,-1,-1])
        # generateInitVelocity(theta_tol, [-1,1,1])
        # generateInitVelocity(theta_tol, [-1,1,-1])
        # generateInitVelocity(theta_tol, [-1,-1,1])
        # generateInitVelocity(theta_tol, [-1,-1,-1])
    # end

    # test for uniformity

    N = 6
    theta_1 = theta_tol / sqrt(N)
    num_in_ring = zeros(N)
    for i in 1:1000
        angle = generateInitVelocityDiskProxima(theta_tol, baseVector)[1]
        if angle <= theta_1
            num_in_ring[1] += 1
        else
            for j in 2:N
                theta_next = theta_1 * sqrt(j)
                if angle <= theta_next
                    num_in_ring[j] += 1
                    break
                end
            end
        end
    end
    println(num_in_ring)

end

function runExperiment(iterations,tolerance)

    for i in 1:iterations
        println(stderr,"Trial ",i)

        (theta_deg, v_init) = generateInitVelocityDiskProxima(tolerance, Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01])
        
        minvec = closestApproach(v_init)

        distance_in_M = norm(minvec)
        distance_in_AU = distance_in_M * 1.47669196951425 * 6.6845871226706E-09

        println(stderr,"For initial velocity ", v_init, ", closest approach was ", distance_in_AU, " AU")

        println(v_init[1],",",v_init[2],",",v_init[3],",",distance_in_AU)
    end
end

function secantMethod(f, x0, x1, tol, cutoff)
    f0 = f(x0)
    println("f0 = ", f0 * 1.47669196951425 * 6.6845871226706E-09)
    for iterations in 1:cutoff
        f1 = f(x1)
        println("x1 = ",x1)
        println("f1 = ", f1 * 1.47669196951425 * 6.6845871226706E-09)
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        println("x2 = ",x2)
        if abs(f1) <= tol
            return x2
        end
        x0 = x1
        x1 = x2
        f0 = f1
    end
    return x2
end

function BroydenMethod(v_0)
    J = ForwardDiff.jacobian(x -> closestApproach(x), v_0)
    # J = Double64[1.2417723007378808e14 -1.4560336524476594e13 -4.641925143472228e13; -1.456033137057556e13 1.3211341961647948e14 -3.546165106914633e13; -4.64192486563975e13 -3.5461661518072055e13 3.018271672896033e13]
    f_0 = closestApproach(v_0)
    println(stderr,"Init norm F = ", norm(f_0))
    return bsolve(x -> closestApproach(x), J, f_0, v_0)
end