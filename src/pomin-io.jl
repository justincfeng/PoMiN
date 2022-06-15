
# Convenience function that converts 1,2,3 into x,y,z for output readability
function num2xyz(i::Int)
    if i == 1
        return "x"
    elseif i == 2
        return "y"
    elseif i == 3
        return "z"
    end
end

# Converts the soln data structure into CSV data
# CSV data will be organized as follows:
#   q's for particle 1
#   p's for particle 1
#   q's for particle 2
#   p's for particle 2
#   ...
function soln2csv(sol::soln)
    
    #N = Int(length(sol.z[1]) / 6)  # number of particles

    # Write CSV header
    csv_str = "iteration,time,"
    for i in 1:sol.N
        # write q header
        for j in 1:sol.d
            csv_str = string(csv_str, "q", num2xyz(j), i, ",")
        end
        for j in 1:sol.d
            csv_str = string(csv_str, "p", num2xyz(j), i)
            if j < sol.d  # add a comma after each except the last one
                csv_str *= ","
            end
        end
        if i < sol.N  # add a comma after each particle except the last one
            csv_str *= ","
        end
    end
    csv_str *= "\n"
    
    # iterate through particle data at each timestep
    numsteps = length(sol.t)
    for tn in 1:numsteps  # tn is timestep number
        # write timestep number and time
        csv_str = string(csv_str, tn, ",", sol.t[tn], ",")
        for i in 1:sol.N
            # get q values for particle i at timestep t
            q = Z2q(sol.N, sol.d, i, sol.z[tn])
            # write q values for particle i at timestep t
            for j in 1:sol.d
                csv_str = string(csv_str, q[j], ",")
            end

            # get p values for particle i at timestep t
            p = Z2p(sol.N, sol.d, i, sol.z[tn])
            # write p values for particle i at timestep t
            for j in 1:sol.d
                csv_str = string(csv_str, p[j])
                if j < sol.d || i < sol.N # add a comma after each except the last value of the last particle
                    csv_str *= ","
                end
            end
        end
        if tn < numsteps   # add newline after each timestep except the last
            csv_str *= "\n"
        end
    end
    return csv_str
end