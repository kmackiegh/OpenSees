/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for inSimulationation on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        

// Written: MHS 
//
// Description: This file contains the SimulationStorage class implementation

#include <SimulationStorage.h>
#include <Information.h>
#include <string.h>

SimulationStorage::SimulationStorage()
{
    beta = 0;
    pf = 0;
    cov = 0;
    qbar = 0;
    rstdv = 0;
}


SimulationStorage::~SimulationStorage()
{
}


const char *
SimulationStorage::getClassType(void) const
{
    return "Simulation";
}


int 
SimulationStorage::setVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"betaSimulation") == 0) {

    }
    else if (strcmp(variable,"pfSimulation") == 0) {
        
    }
    else if (strcmp(variable,"covSimulation") == 0) {
        
    }
    else if (strcmp(variable,"qbarSimulation") == 0) {
        
    }
    else if (strcmp(variable,"rstdvSimulation") == 0) {
        
    }
    
    else {
        opserr << "SimulationStorage:: unknown variable " << variable <<
            " in setVariable()" << endln;
    }
        
    return 0;
}


int 
SimulationStorage::getVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"betaSimulation") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(beta);
    }
    else if (strcmp(variable,"pfSimulation") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(pf);
    }
    else if (strcmp(variable,"covSimulation") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(cov);
    }
    else if (strcmp(variable,"qbarSimulation") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(qbar);
    }
    else if (strcmp(variable,"rstdvSimulation") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(rstdv);
    }
    
    else {
        opserr << "SimulationStorage:: unknown variable " << variable <<
            " in getVariable()" << endln;
    }
    
    return 0;
}
