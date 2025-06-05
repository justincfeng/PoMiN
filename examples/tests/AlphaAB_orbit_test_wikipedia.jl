using DoubleFloats

masses = Double64[1.0788, 0.9092]        # Alpha Centauri A & B in units of solar masses

Z_init = Double64[-1.008633862E+13, -9.2053185745E+12, -2.4039340207E+13, -1.1010729524E+13, -9.2051242414E+12, -2.4039340207E+13, -5.6748967592E-05, 1.5454926272E-04, 1.3207642456E-04, -7.1413412414E-05, 1.3872133524E-04, 1.2324033123E-04]
tspan = (Double64(0), Double64(1.92213424E+15))
δ = Double64(1.92213424E+11)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
maxit = 10000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)