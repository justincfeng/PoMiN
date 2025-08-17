#-----------------------------------------------------------------------
#
#   3-BODY CLOSEST APPROACH CALCULATION
#   Modernized version compatible with current PoMiN framework
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, DoubleFloats

"""
    closest_approach_3body(third_body_mass::Real, q_init::Vector, p_init::Vector; 
                          tpfl::Type=Double64, simulation_years::Real=22.0)

Calculate closest approach in a 3-body system: spacecraft, Proxima Centauri, and a third body.
This is a modernized version of the legacy Twop1.jl script.

Arguments:
- third_body_mass: Mass of the third body
- q_init: Initial position vector [x, y, z] of the third body
- p_init: Initial momentum vector [px, py, pz] of the third body
- tpfl: Floating point type for precision (default: Double64)
- simulation_years: Duration of simulation in years (default: 22.0)

Returns:
- Named tuple with closest approach distance, vector, and time
"""
function closest_approach_3body(third_body_mass::Real, q_init::Vector, p_init::Vector; 
                               tpfl::Type=Double64, simulation_years::Real=22.0)
    
    println("Starting 3-body closest approach calculation...")
    
    # Convert inputs to specified precision
    third_body_mass = tpfl(third_body_mass)
    q_init = convert(Vector{tpfl}, q_init)
    p_init = convert(Vector{tpfl}, p_init)
    
    #-----------------------------------------------------------------------
    #   PHYSICAL PARAMETERS (from original script)
    #-----------------------------------------------------------------------
    
    # Masses: Spacecraft, Proxima Centauri, Third body
    m_spacecraft = tpfl(1.0057832537E-33)  # Spacecraft mass
    m_proxima = tpfl(1E-40)                # Proxima test particle mass
    m_third = third_body_mass              # Third body mass
    
    # Initial positions
    q_spacecraft = tpfl.([−1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07])
    q_proxima = tpfl.([−9.90338347925777E+12, −7.58520275447461E+12, −2.41494965557852E+13])
    q_third = q_init
    
    # Initial momenta
    p_spacecraft = tpfl.([−7.486220592724260E-35, −5.723861062118610E-35, −1.823989847481890E-34])
    p_proxima = tpfl.([−3.34779933990201E-45, 7.24198145771899E-45, 6.82553747232694E-45])
    p_third = p_init
    
    #-----------------------------------------------------------------------
    #   SET UP 3-BODY SYSTEM
    #-----------------------------------------------------------------------
    
    # Create 3-body particle system directly
    masses = [m_spacecraft, m_proxima, m_third]
    positions = [q_spacecraft, q_proxima, q_third]
    momenta = [p_spacecraft, p_proxima, p_third]
    
    particles = Particles(masses, positions, momenta)
    
    #-----------------------------------------------------------------------
    #   INTEGRATION PARAMETERS
    #-----------------------------------------------------------------------
    
    # Time span (22 years in the original script)
    seconds_per_year = tpfl(3.154e7)  # Approximate seconds per year
    t_end = tpfl(simulation_years) * seconds_per_year * tpfl(1e7)  # Scale factor from original
    t_start = tpfl(0.0)
    
    # Set up integration parameters with high precision
    params = ParametersJulia((t_start, t_end), 
                           integrator="Vern9",
                           atol=1e-20, 
                           rtol=1e-20,
                           Nrec=1000)  # Save many points for closest approach analysis
    
    #-----------------------------------------------------------------------
    #   SOLVE THE SYSTEM
    #-----------------------------------------------------------------------
    
    println("Integrating 3-body system...")
    @time solution = solve(particles, params)
    println("Integration completed with $(length(solution.t)) time steps")
    
    #-----------------------------------------------------------------------
    #   FIND CLOSEST APPROACH
    #-----------------------------------------------------------------------
    
    println("Analyzing closest approach...")
    
    # Initialize closest approach tracking
    min_distance = tpfl(1e100)
    min_vector = zeros(tpfl, 3)
    min_time = tpfl(0.0)
    min_index = 1
    
    # Extract positions for all time steps
    n_steps = length(solution.t)
    
    for i in 1:n_steps
        # Get positions of spacecraft (particle 1) and Proxima (particle 2)
        q_spacecraft_t = solution.u[i][1:3]      # First 3 components: spacecraft position
        q_proxima_t = solution.u[i][4:6]         # Next 3 components: Proxima position
        
        # Calculate separation vector and distance
        separation_vector = q_spacecraft_t - q_proxima_t
        distance = norm(separation_vector)
        
        # Update minimum if this is closer
        if distance < min_distance
            min_distance = distance
            min_vector = separation_vector
            min_time = solution.t[i]
            min_index = i
        end
    end
    
    #-----------------------------------------------------------------------
    #   RESULTS
    #-----------------------------------------------------------------------
    
    println("Closest approach analysis complete:")
    println("  Minimum distance: $(min_distance)")
    println("  Time of closest approach: $(min_time)")
    println("  Step index: $(min_index) / $(n_steps)")
    
    # Return comprehensive results
    return (
        distance = min_distance,
        vector = min_vector,
        time = min_time,
        step_index = min_index,
        total_steps = n_steps,
        solution = solution
    )
end

#-----------------------------------------------------------------------
#   EXAMPLE USAGE
#-----------------------------------------------------------------------

# Example: Third body with moderate mass at distant position
if abspath(PROGRAM_FILE) == @__FILE__
    # Third body parameters
    third_mass = 1e-30  # Moderate mass
    third_pos = [1e13, 1e13, 1e13]  # Distant position
    third_mom = [1e-40, -1e-40, 0]  # Small momentum
    
    # Run calculation
    result = closest_approach_3body(third_mass, third_pos, third_mom, 
                                  tpfl=Double64, simulation_years=22.0)
    
    println("\nFinal Results:")
    println("Closest approach distance: $(result.distance)")
    println("Separation vector: $(result.vector)")
    println("Time of closest approach: $(result.time)")
end
