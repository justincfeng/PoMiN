using DoubleFloats

include("broyden.jl")

function closestApproach(v_init::RealVec)
    println(stderr,"closestApproach() called")
    # Mass of Sun, Earth, Spacechip, Proxima, Alpha A, Alpha B in M
    masses = Double64[1, 3.0033693739E-06, 1.0057832537E-33, 0.1221, 1.0788, 0.9092] 

    # Make sure spacecraft velocity has length of 0.2
    v_init = v_init / norm(v_init) * 0.2
    # Find spacecraft momentum p = gamma m v
    p_init = masses[3]*v_init/sqrt(1-0.2^2)

    # Initial Z: Positions first, momenta second
    Z_init = [0, 0, 0, -1.8667826140E+07, 8.9633743993E+07, 3.8883291714E+07, -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07, -9.9157694997E+12, -7.5891180191E+12, -2.4171267258E+13, -1.0514396564E+13, -8.7920419314E+12, -2.4558076740E+13, -1.0514148951E+13, -8.7899759755E+12, -2.4558922285E+13, 0, 0, 0, -2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11, p_init[1], p_init[2], p_init[3], -3.8172314441E-06, 9.0494364315E-06, 8.9932012189E-06, -3.2681576581E-05, 8.2784394610E-05, 7.2543810801E-05, -2.9658611095E-05, 6.6362559735E-05, 6.7386459079E-05]

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
        # get distance and save it to the list of minima
        dist = norm(q1-q2)
        if dist < mindist
            mindist = dist
            minvec = q1-q2
        end
    end
    return minvec
end

function generateInitVelocityMonteCarlo(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method")

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

function generateInitVelocity(theta_tol, baseVector::RealVec)

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
    theta_min = theta_0 - delta/2
    theta_max = theta_0 + delta/2
    phi_min = phi_0 - delta/2
    phi_max = phi_0 + delta/2

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

    # # test conversion to spherical and back
    # baseVector = Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01]
    # # normalize baseVector
    # baseVector = baseVector / norm(baseVector)
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

    for i in 1:20
        generateInitVelocity(theta_tol, Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01])
        # test all octants
        # generateInitVelocity(theta_tol, [1,1,1])
        # generateInitVelocity(theta_tol, [1,1,-1])
        # generateInitVelocity(theta_tol, [1,-1,1])
        # generateInitVelocity(theta_tol, [1,-1,-1])
        # generateInitVelocity(theta_tol, [-1,1,1])
        # generateInitVelocity(theta_tol, [-1,1,-1])
        # generateInitVelocity(theta_tol, [-1,-1,1])
        # generateInitVelocity(theta_tol, [-1,-1,-1])
    end

end

function runExperiment(iterations,tolerance)

    for i in 1:iterations
       (theta_deg, v_init) = generateInitVelocity(tolerance, Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01])
        
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