using DoubleFloats

numhalves = 8  # how many times to halve the timestep; smallest timestep will be δ divided by 2^(numhalves-1)

tadap = t -> 1.0  # effectively turns off adaptive time stepping since δ will be multiplied by 1 each timestep
maxit = 1000000

# Sun's mass and Mercury's mass
masses = Double64[1 0.000000166]

# Uncomment this section for orbit eccentricity = 0.5
Z_init = Double64[0 0 0 0 -3.27E-12 0 3860700000 0 0 0 3.27E-12 0] # Mercury orbiting Sun w eccentricity=0.5
δ = Double64(6.15E+13)
tspan = (Double64(0), Double64(1.23E+16))

# # Uncomment this section for orbit eccentricity = 0.95
# Z_init = Double64[0 0 0 0 -3.73E-12 0 3860700000 0 0 0 3.73E-12 0] # Mercury orbiting Sun w eccentricity=0.95
# δ = Double64(2.70E+13)
# tspan = (Double64(0), Double64(2.70E+16))

p_results = zeros(Double64,numhalves)

for i in 0:numhalves-1
    divisor = 2^i
    println("t/" * string(divisor))

    timestep = δ / Double64(divisor)
    println("timestep = " * string(timestep))
    sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), timestep, tadap, tspan, maxit)

    # get the last row
    last_row = length(sol.t)
    # verify end time matches what was specified in tspan
    if (sol.t[last_row] == tspan[2])
        # get p^2 for smaller object and save it
        p = Z2p(2, 3, 2, sol.z[last_row])   # n=2 particles, d=3 dims, particle number i=2
        println("final p for particle 2 = " * string(p))
        p_results[i+1] = p[1]*p[1] + p[2]*p[2] + p[3]*p[3]
    else
        error("Error: final time of " * string(sol.t[last_row]) * " doesn't equal expected value " * string(tspan[2]))
        if (sol.t[last_row] < tspan[2])
            error("Is max iterations large enough?")
        end
    end
end
println(p_results)
# calculate Q values
for i in 2:numhalves-1
    Q = (p_results[i-1] - p_results[i]) / (p_results[i] - p_results[i+1])
    println(string(Q))
end