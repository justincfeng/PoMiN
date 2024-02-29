# Tests the Milky Way potential by placing the Sun in orbit around the galactic center

using DoubleFloats

masses = Double64[1]        # one object, the Sun, with mass in units of solar masses

tspan = (Double64(0), Double64(7.5E+21)) # about 1 billion years
δ = Double64(1.5E+17) # about 23,000 years (50,000 time steps)
noadapt = (δ, z, zdot) -> δ  # turns off adaptive time stepping
maxit = 100000


# galactocentric frame

# Z_init = Double64[-1.708859462494220E+17 0 4.346342845091530E+14 4.302976828056160E-05 8.192334178066610E-04 2.595128660641620E-05]  # orbital radius = 8178 pc (2019 data) and speed from astropy
# Z_init = Double64[1.654953934741780E+17 0 0 0 7.563565858618100E-04 0]  # using 2020 VERA data (orbital radius = 7920 pc) 
Z_init = Double64[1.755254114968200E+17 0 0 0 8.072251103795280E-04 0]  # using Irrgang model I (orbital radius = 8400 pc)
# Z_init = [-1.755248460092340E+17, 0, 4.346342845091530E+14, 4.302976828056160E-05, 8.056595739976890E-04, 2.595128660641620E-05]    # astropy data adjusted to match Irrgang orbital radius (8400 pc) and velocity (242 km/s)

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_MW(3, masses, x, 0, 0, 0), δ, noadapt, tspan, maxit)

# heliocentric coordinate system

# Z_init = Double64[0 0 0 4.302976828056160E-05 8.192334178066610E-04 2.595128660641620E-05]  # orbital radius = 8178 pc (2019 data) and speed from astropy
# sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, noadapt, tspan, maxit)

printcsv(sol)