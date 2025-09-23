#!/usr/bin/env julia

"""
Script to extract and tabulate miss distances from the output of one_particle_miss_dist_All.jl
Run this after running the main script to get a formatted table.
"""

using Printf

function parse_miss_distances(output_text::String)
    bodies = ["Sun", "Alpha Centauri", "Jupiter", "Earth", "Proxima", "Moon", "Mars"]
    newtonian_distances = Float64[]
    relativistic_distances = Float64[]
    
    lines = split(output_text, '\n')
    
    for line in lines
        if contains(line, "Newtonian miss distance:")
            # Extract the numerical value
            match_result = match(r"Newtonian miss distance:\s*([\d\.e\-\+]+)", line)
            if match_result !== nothing
                push!(newtonian_distances, parse(Float64, match_result.captures[1]))
            end
        elseif contains(line, "Relativistic miss distance:")
            # Extract the numerical value
            match_result = match(r"Relativistic miss distance:\s*([\d\.e\-\+]+)", line)
            if match_result !== nothing
                push!(relativistic_distances, parse(Float64, match_result.captures[1]))
            end
        end
    end
    
    return newtonian_distances, relativistic_distances
end

function create_miss_distance_table(newtonian_distances, relativistic_distances, bodies)
    println("="^80)
    println("MISS DISTANCE SUMMARY TABLE")
    println("="^80)
    println(@sprintf("%-15s | %-18s | %-19s | %-12s", "Body", "Newtonian Miss", "Relativistic Miss", "Ratio (N/R)"))
    println("-"^80)
    
    for i in 1:min(length(bodies), length(newtonian_distances), length(relativistic_distances))
        ratio = newtonian_distances[i] / relativistic_distances[i]
        println(@sprintf("%-15s | %-18.6e | %-19.6e | %-12.6f", 
                bodies[i], newtonian_distances[i], relativistic_distances[i], ratio))
    end
    println("="^80)
end

# Example usage:
println("Miss Distance Table Generator")
println("To use this script:")
println("1. Run: julia one_particle_miss_dist_All.jl > output.txt")
println("2. Then use this script to parse the output")
println()

# If you have output text, you can parse it like this:
# output_text = read("output.txt", String)
# newtonian_distances, relativistic_distances = parse_miss_distances(output_text)
# bodies = ["Sun", "Alpha Centauri", "Jupiter", "Earth", "Proxima", "Moon", "Mars"]
# create_miss_distance_table(newtonian_distances, relativistic_distances, bodies)

# For demonstration, here's what the table would look like:
bodies = ["Sun", "Alpha Centauri", "Jupiter", "Earth", "Proxima", "Moon", "Mars"]
println("Example table format:")
create_miss_distance_table([1e-5, 2e-6, 3e-7, 4e-8, 5e-9, 6e-10, 7e-11], 
                          [1.1e-5, 2.2e-6, 3.3e-7, 4.4e-8, 5.5e-9, 6.6e-10, 7.7e-11], 
                          bodies)
