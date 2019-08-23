#include "RK4.h"

RK4::RK4(Pomin* pomin, double courantNum) : Integrator(pomin, courantNum)
{
}

double RK4::getTimeStep()
{
    /************************************************************
     * CALCULATE TIMESTEP
     ************************************************************
    double ats = m_pomin->getSmallesta();
    double bts = m_pomin->getSmallestb();
    if (courant!=0) {
        double qdotsmallest = sqrt((part[bts].qdot[0]-part[ats].qdot[0]) * (part[bts].qdot[0]-part[ats].qdot[0]) \
                                   + (part[bts].qdot[1]-part[ats].qdot[1]) * (part[bts].qdot[1]-part[ats].qdot[1]) \
                                   + (part[bts].qdot[2]-part[ats].qdot[2]) * (part[bts].qdot[2]-part[ats].qdot[2]));
        
        double adapted_step = courant * (m_pomin->getRsmallest()/qdotsmallest);
        
        // Ensure the timestep is non0 and less than the maximum size.
        if (adapted_step > 0){
            if(adapted_step < maxtimestep){
                timestep = adapted_step;
            }
            else{
                timestep=maxtimestep;
            }
        }
    }
     */
    return 0;
}

void RK4::solve(ofstream output)
{
    
}
