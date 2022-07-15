using DoubleFloats

numhalves = 8  # how many times to halve the timestep; smallest timestep will be δ divided by 2^(numhalves-1)

tadap = t -> 1.0  # effectively turns off adaptive time stepping since δ will be multiplied by 1 each timestep
maxit = 1000000

# Sun's mass and Mercury's mass
masses = Double64[1 0.000000166]

# Uncomment this section for orbit eccentricity = 0.5
Z_init = Double64[0 0 0 3860700000 0 0 0 -3.27E-12 0 0 3.27E-12 0] # Mercury orbiting Sun w eccentricity=0.5
δ = Double64(1)
tspan = (Double64(0), Double64(8))

# # Uncomment this section for orbit eccentricity = 0.95
# Z_init = Double64[0 0 0 3860700000 0 0 0 -3.73E-12 0 0 3.73E-12 0] # Mercury orbiting Sun w eccentricity=0.5
# δ = Double64(1.23E+13)
# tspan = (Double64(0), Double64(2.46E+14))

results = zeros(Double64,numhalves)

for i in 0:numhalves-1
    divisor = 2^i
    println("t/" * string(divisor))

    timestep = δ / Double64(divisor)
    println("timestep = " * string(timestep))
    #sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), timestep, tadap, tspan, maxit)
    sol = hsintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), timestep, Double64(1), tspan, maxit)

    # get the last row
    last_row = length(sol.t)
    # verify end time matches what was specified in tspan
    if (sol.t[last_row] == tspan[2])
        # get q^2 for smaller object and save it
        q = Z2q(2, 3, 2, sol.z[last_row])   # n=2 particles, d=3 dims, particle number i=2
        println("final q for particle 2 = " * string(q))
        #results[i+1] = q[1]*q[1] + q[2]*q[2] + q[3]*q[3]
        results[i+1] = norm(q)
        println("length of sol.t = " * string(length(sol.t)))
        println("final time = " * string(sol.t[last_row]))
    else
        error("Error: final time of " * string(sol.t[last_row]) * " doesn't equal expected value " * string(tspan[2]))
        if (sol.t[last_row] < tspan[2])
            error("Is max iterations large enough?")
        end
    end
end
println(results)
# calculate Q values
for i in 2:numhalves-1
    Q = (results[i-1] - results[i]) / (results[i] - results[i+1])
    println(string(Q))
end