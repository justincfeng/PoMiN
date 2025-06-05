using DoubleFloats

masses = Double64[1, 1E-07]

# uncomment for eccentricity = 0
# Z_init = Double64[0, 0, 0, 4e10, 0, 0, 0, -5e-13, 0, 0, 5e-13, 0]
# tspan = (Double64(0), Double64(16E+16))
# δ = Double64(16E+12)

# uncomment for eccentricity = 0.5
# a = 4e10
# rp = a*0.5      # periapsis = a*(1-e)
# vp = sqrt(3*(masses[1]+masses[2])/a)  # velocity at periapsis = sqrt( (1+e)/(1-e) * (M+m)/a )
# p2 = masses[2] * vp
# Z_init = Double64[0, 0, 0, rp, 0, 0, 0, -1*p2, 0, 0, p2, 0]
# tspan = (Double64(0), Double64(16E+16))
# δ = Double64(16E+12)

# uncomment for eccentricity = 0.9
# a = 4e10
# rp = a*0.1      # periapsis = a*(1-e)
# vp = sqrt(19*(masses[1]+masses[2])/a)  # velocity at periapsis = sqrt( (1+e)/(1-e) * (M+m)/a )
# p2 = masses[2] * vp
# Z_init = Double64[0, 0, 0, rp, 0, 0, 0, -1*p2, 0, 0, p2, 0]
# tspan = (Double64(0), Double64(16E+16))
# δ = Double64(16E+11)

# uncomment for eccentricity = 0.95
a = 4e10
rp = a*0.05      # periapsis = a*(1-e)
vp = sqrt(39*masses[1]/a)  # velocity at periapsis = sqrt( (1+e)/(1-e) * M/a )
p2 = masses[2] * vp
Z_init = Double64[0, 0, 0, rp, 0, 0, 0, -1*p2, 0, 0, p2, 0]
tspan = (Double64(0), Double64(16E+16))
δ = Double64(4E+11)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.001)   # last parameter is Courant number
maxit = 1500000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)