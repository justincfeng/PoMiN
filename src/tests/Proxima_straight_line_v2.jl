using DoubleFloats

# use only spacechip and proxima, with proxima having very small mass
masses = Double64[1.0057832537E-33, 1e-40] # Spacechip, Proxima in solar masses

# allow Proxima to travel in straight line (by removing all other gravitational influences)
Z_init = Double64[  -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07,                  # spacechip pos
                    -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos (Kervella 2017)
                    -7.488708327115080E-35, -5.721560589493490E-35, -1.823959902534090E-34, # spacechip mom 
                    -3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45       # Proxima mom (Kervella 2017)
                    ]

tspan = (Double64(0), Double64(1.922E+14)) # 30 years
# tspan = (Double64(0), Double64(1.6017785E+14)) # 25 years
δ = Double64(1E7) # less than a minute
# δ = Double64(1.75E+10) # about one day

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
start_dist = 2.720605E+13  # 4.2465 light years
min_dist = 1E9 # about 10 AU -- once within this distance it will use min_dt
max_dt = Double64(1.75E+10) # about one day
min_dt = Double64(1.2E7) # about a minute
# adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
adapt = (δ, z, zdot) -> dt_lin_int(δ, z, zdot, start_dist, min_dist, max_dt, min_dt)
maxit = 10000000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)