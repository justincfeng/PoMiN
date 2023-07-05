using DoubleFloats

masses = Double64[1, 0.00000300336937]        # Sun and Earth masses in units of solar masses

Z_init = Double64[-1.708859462494220E+17, 0, 4.346342845091530E+14, -1.708859463641260E+17, 6.958323268751290E+07, 4.346343398098600E+14, 4.302976828056160E-05, 8.192334178066620E-04, 2.595128660641620E-05, -9.75E-11, 2.268729561576910E-09, 3.148365038646750E-10]
tspan = (Double64(0), Double64(1.922E+12))
δ = Double64(1.922E+8)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)