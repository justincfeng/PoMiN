#include "pomin.h"
#include "RK4.h"

int main(int argc, char** argv)
{
    string input_data;

    // Parse command line arguments

    // Get input
    iostream my_input = new iostream();
    string description, delim;
    double start_time, end_time, timestep, grav, courant = 1.0;
    int interations;
    int N;  // number of particles
    
    my_input >> description >> delim >> start_time >> delim >> end_time >> delim >> timestep >> delim >> iterations >> delim >> courant >> "," >> grav >> "," >> N;
    
    // Construct pomin
    Pomin * pomin = new Pomin(N);

    // Populate pomin object with input data
    double m, qx, qy, qz, px, py, pz;
    for(int i=0; i<N; i++)
    {
        my_input >> m >> "," >> qx >> "," >> ;
        pomin->setParticle(i, m, qx, qy, qz, px, py, pz);
    }
    
    // Construct integrator
    RK4 rk4(pomin, courant);

    // Create output file

    // Solve equations
    //rk4.solve(ofstream);
}
