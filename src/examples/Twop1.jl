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