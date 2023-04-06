using DoubleFloats

masses = Double64[1, 0.00000300336937]        # Sun and Earth masses in units of solar masses

Z_init = Double64[-1.708859462494220E+17, 0, 4.346342845091530E+14, 4.302976828056160E-05, 8.192334178066610E-04, 2.595128660641620E-05, -1.708859463454590E+17, -2.005051128548710E+07, 4.346343009265680E+14, 2.008478809088120E-10, 2.319118451106750E-09, 3.366825528988520E-10]
tspan = (Double64(0), Double64(7.5E+21))
δ = Double64(1.5E+17)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)