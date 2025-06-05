using DoubleFloats

masses = Double64[1.0788, 0.9092, 1]        # Alpha Centauri A & B in units of solar masses

Z_init = Double64[-1.0514396567E+13, -8.7920419041E+12, -2.4558076843E+13, -1.0514148949E+13, -8.7899759136E+12, -2.4558922297E+13, 0, 0, 0, -3.2681576591E-05, 8.2784394612E-05, 7.2543810790E-05, -3.5191057704E-05, 7.8741673392E-05, 7.9956568487E-05, 0, 0, 0]
tspan = (Double64(0), Double64(1.92213424E+15))
δ = Double64(1.92213424E+11)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
maxit = 10000

sol = hrkintegrator(3, 3, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)