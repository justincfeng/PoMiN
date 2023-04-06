using DoubleFloats

masses = Double64[1.988, 0.1221]        # the Alpha Centauri binary star and Proxima Centauri with mass in units of solar masses

Z_init = Double64[-1.0422461557E+13, -8.7018418912E+12, -2.4325683942E+13, -9.9157696012E+12, -7.5891180437E+12, -2.4171267294E+13, -6.19207190E-05, 1.47027953E-04, 1.43891876E-04, -3.8172314445E-06, 9.0494364315E-06, 8.9932012183E-06]
tspan = (Double64(0), Double64(3.53E+18))
δ = Double64(3.53E+14)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
maxit = 10000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)