#include"pomin.h"
#include<fstream>

class Integrator
{
    public:
        Integrator(Pomin* pomin, double courantNum = 1.0) : m_pomin(pomin), courant(courantNum)
        { }
    
    virtual double getTimeStep() = 0;
    virtual void solve(ofstream output) = 0;
    
    protected:
    
    Pomin* m_pomin;
    double courant;
    double timestep;
};
