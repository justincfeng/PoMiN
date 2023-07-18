using DoubleFloats

masses = Double64[1]        # one object, the Sun, with mass in units of solar masses

Z_init = Double64[-1.708859462494220E+17 0 4.346342845091530E+14 4.302976828056160E-05 8.192334178066610E-04 2.595128660641620E-05]  # Sun's position and momentum in galactocentric frame using 2019 data
tspan = (Double64(0), Double64(7.5E+21)) # about 1 billion years
δ = Double64(1.5E+17) # about 23,000 years (50,000 time steps)

noadapt = (δ, z, zdot) -> δ  # effectively turns off adaptive time stepping
maxit = 100000

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_MW(3, masses, x, 0, 0, 0), δ, noadapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)