using DoubleFloats

masses = Double64[1]        # one object, the Sun, with mass in units of solar masses

Z_init = Double64[-1.6971700079E+17 0 4.3463742829E+14 4.3029768281E-05 8.1923341781E-04 2.5951286606E-05]  # Sun's position and momentum in galactocentric frame
tspan = (Double64(0), Double64(1.922E+13)) # 3 years
δ = Double64(1.922E+9) # about 2 hrs

noadapt = (δ, z, zdot) -> δ  # effectively turns off adaptive time stepping
maxit = 100000

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, noadapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)