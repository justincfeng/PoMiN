using DoubleFloats

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

    tspan = (Double64(0), Double64(1.922E+14)) # 30 years
    # tspan = (Double64(0), Double64(1.922E+11)) # 10 days
    δ = Double64(1.922E+10) # about one day

    no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
    adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
    maxit = 100000

    sol = hrkintegrator(3, 6, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)
    
    println(stderr,"integrator done")
    printcsv(sol)

    # find distance of closest approach
    mindist = Double64(1E100)
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
        end
    end
    return mindist
end

function generateInitVelocity(theta_tol, baseVector::RealVec)
    println(stderr, "Finding initial velocity vector using Monte Carlo method")

    # normalize baseVector
    baseVector = baseVector / norm(baseVector)

    while true
        # generate random vector in [-1,1]^3
        x = rand() * 2 - 1
        y = rand() * 2 - 1
        z = rand() * 2 - 1
        rand_v = [x, y, z]

        # normalize random vector
        rand_v = rand_v / norm(rand_v)

        # find angle between baseVector and random vector
        cosine_theta = dot(rand_v,baseVector) # both unit vectors
        theta_rad = acos(cosine_theta)
        theta_deg = rad2deg(theta_rad)

        # compare to tolerance angle
        if theta_deg < theta_tol
            println(stderr,"For initial velocity ", rand_v, " angle with base vector is ", theta_deg)
            return theta_deg, rand_v
        end
    end
end

function runExperiment(iterations)

    for i in 1:iterations
       (theta_deg, v_init) = generateInitVelocity(0.1, Double64[-7.00365851051320E-02, -5.35052039083642E-02, -1.70575701379575E-01])

        distance_in_M = closestApproach(v_init) 
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