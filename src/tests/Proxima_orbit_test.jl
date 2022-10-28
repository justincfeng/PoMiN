using DoubleFloats

masses = Double64[1.988, 0.1221]        # the Alpha Centauri binary star and Proxima Centauri with mass in units of solar masses

Z_init = Double64[-1.0422461557E+13, -8.7018418912E+12, -2.4325683942E+13, -9.9157696012E+12, -7.5891180437E+12, -2.4171267294E+13, -6.19207190E-05, 1.47027953E-04, 1.43891876E-04, -6.21511557E-05, 1.47340537E-04, 1.46424931E-04]
tspan = (Double64(0), Double64(3.53E+18))
δ = Double64(3.53E+14)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.1)   # last parameter is Courant number
maxit = 10000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

println(soln2csv(sol))