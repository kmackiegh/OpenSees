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
** file 'COPYRIGHT'  in main directory for information on usage and   **
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
// Description: This file contains the FOSMStorage class implementation

#include <FOSMStorage.h>
#include <Information.h>
#include <string.h>

FOSMStorage::FOSMStorage()
{
    alpha = 0;
    gradientX = 0;
    beta = 0;
}


FOSMStorage::~FOSMStorage()
{
    if (alpha != 0)
        delete alpha;
    if (gradientX != 0)
        delete gradientX;
}


const char *
FOSMStorage::getClassType(void) const
{
    return "FOSM";
}


int 
FOSMStorage::setVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"alphaFOSM") == 0) {
        alpha = new Vector(*(theInfo.theVector));
    }
    
    else if (strcmp(variable,"gradientXFOSM") == 0) {
        gradientX = new Vector(*(theInfo.theVector));
    }
    
    else if (strcmp(variable,"betaFOSM") == 0) {

    }
    
    else {
        opserr << "FOSMStorage:: unknown variable " << variable <<
            " in setVariable()" << endln;
    }
        
    return 0;
}


int 
FOSMStorage::getVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"alphaFOSM") == 0) {
        if (alpha != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*alpha);
        }
        else
            return -1;
    }
    
    else if (strcmp(variable,"gradientXFOSM") == 0) {
        if (gradientX != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*gradientX);
        }
        else
            return -1;
    }
    
    else if (strcmp(variable,"betaFOSM") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(beta);
    }
    
    else {
        opserr << "FOSMStorage:: unknown variable " << variable <<
            " in getVariable()" << endln;
    }
    
    return 0;
}
