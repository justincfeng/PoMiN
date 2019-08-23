#include "integrator.h"

class RK4 : Integrator
{
    public:
    
    RK4(Pomin* pomin, double courantNum);
    double getTimeStep();
    void solve(ofstream output);

};
