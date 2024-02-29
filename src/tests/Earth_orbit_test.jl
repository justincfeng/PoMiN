using DoubleFloats

masses = Double64[1, 0.00000300336937]        # the Sun and Earth with mass in units of solar masses

#Z_init = Double64[0, 0, 0, -1.866782584728260E+07, 8.96337439732320E+07, 3.888329146643650E+07, 0, 0, 0, -2.983904203718850E-10, -5.038888951024780E-11, -2.184604906902130E-11]  # ICRS coords J2000, init conditions from astropy
Z_init = Double64[0, 0, 0, -1.866782584728260E+07, 8.96337439732320E+07, 3.888329146643650E+07, 0, 8.072251103795280E-04, 0, -2.983904203718850E-10, 2.374006281698500E-09, -2.184604906902130E-11]  # using Irrgang model I for init vel/mom
tspan = (Double64(0), Double64(1.60175E+13)) # about 2.5 years
δ = Double64(3.2035E+8)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 1000000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH_plus_MW(3, masses, x, Double64(1.755254114968200E+17), Double64(0), Double64(0)), δ, no_adapt, tspan, maxit)  # using Irrgang model I
printcsv(sol)

# sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan, 1e-10, 1e-4)
# printODEsoln(sol)

#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)
# sol = hointegrator(Z_init, x -> pomin.dH(3, masses, x), tspan, 1e-6, 1e-3, AutoTsit5(Rodas5()))
# printODEsoln(sol)