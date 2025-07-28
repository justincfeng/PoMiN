using DoubleFloats

masses = Double64[2.0429, 0.1221]        # the Alpha Centauri binary star and Proxima Centauri with mass in units of solar masses

# Using data from Kervella et al 
Z_init = Double64[  -1.04524521660786E+13, -8.74487747435090E+12, -2.44217824863455E+13,    # Alpha AB pos
                    -9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13,    # Proxima pos
                    -6.36355857813074E-05, 1.50812651623504E-04, 1.47502158600961E-04,      # Alpha AB mom
                    -4.08766299402035E-06, 8.84245935987489E-06, 8.33398125371119E-06 ]     # Proxima mom

tspan = (Double64(0), Double64(3.5e18))  # one orbit ~ 550,000 years ~ 3.524E+18 M
δ = Double64(1.17e13)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
maxit = 1000000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)