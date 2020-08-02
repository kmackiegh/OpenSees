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

// Written: KRM
// Created: Aug 2019
// Revision: A
//
// Description: This file contains the class implementation for BilinearTS.
//
// What: "@(#) BilinearTS.C, revA"

#include <string.h>

#include <BilinearTS.h>
#include <BilinearTS2D.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


void *
OPS_NewBilinearTS(void)
{
    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 4) {
        opserr << "Want: nDMaterial BilinearTS tag? delc? sigc? lamcr? <cmult?>" << endln;
        return 0;
    }

    int iData[1];
    double dData[3];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial BilinearTS \n";
        return 0;
    }

    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial BilinearTS: " << iData[0] <<"\n";
        return 0;
    }
    
    double cmult = 1.0;
    numData = 1;
    if (numArgs == 5) {
        if (OPS_GetDouble(&numData, &cmult) != 0) {
            opserr << "WARNING invalid cmult: nDMaterial BilinearTS: " << cmult <<"\n";
            return 0;
        }
    }

    theMaterial = new BilinearTS(iData[0], dData[0], dData[1], dData[2], cmult);

    return theMaterial;
}



BilinearTS::BilinearTS
(int tag, int classTag, double dc, double sc, double lc, double cm)
  :NDMaterial(tag, classTag), 
delc(dc), sigc(sc), lamcr(lc), cmult(cm)
{
    // do some input checks
    if (delc < 0)
        delc = fabs(delc);
    if (sigc < 0)
        sigc = fabs(sigc);
    if (lamcr < 0 || lamcr > 1)
        lamcr = 0.5;
    if (cmult <= 0)
        cmult = 1.0;
    
    // derived properties
}

BilinearTS::BilinearTS
(int tag, double dc, double sc, double lc, double cm)
  :NDMaterial(tag, ND_TAG_BilinearTS),
delc(dc), sigc(sc), lamcr(lc), cmult(cm)
{
    // derived properties
}

BilinearTS::~BilinearTS()
{
	
}

double
BilinearTS::getRho()
{ 
    return 0;
}

NDMaterial*
BilinearTS::getCopy (const char *type)
{
    if (strcmp(type,"2D") == 0 || strcmp(type,"2d") == 0) {
        BilinearTS2D *theModel;
        theModel = new BilinearTS2D (this->getTag(), delc, sigc, lamcr, cmult);

        return theModel;
    }

    // Handle other cases
    else
        return NDMaterial::getCopy(type);
}

int
BilinearTS::setTrialStrain (const Vector &v)
{
    opserr << "BilinearTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
BilinearTS::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "BilinearTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
BilinearTS::setTrialStrainIncr (const Vector &v)
{
    opserr << "BilinearTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
BilinearTS::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "BilinearTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
BilinearTS::getTangent (void)
{
    opserr << "BilinearTS::getTangent -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Matrix *ret = new Matrix();
    return *ret;
}

const Matrix&
BilinearTS::getInitialTangent (void)
{
    return this->getTangent();
}

const Vector&
BilinearTS::getStress (void)
{
    opserr << "BilinearTS::getStress -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
BilinearTS::getStrain (void)
{
    opserr << "BilinearTS::getStrain -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
BilinearTS::getState (void)
{
    opserr << "BilinearTS::getState -- subclass responsibility\n";
    exit(-1);
    
    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

int
BilinearTS::commitState (void)
{
    opserr << "BilinearTS::commitState -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
BilinearTS::revertToLastCommit (void)
{
    opserr << "BilinearTS::revertToLastCommit -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
BilinearTS::revertToStart (void)
{
    opserr << "BilinearTS::revertToStart -- subclass responsibility\n";
    exit(-1);
    return -1;
}

NDMaterial*
BilinearTS::getCopy (void)
{
    opserr << "BilinearTS::getCopy -- subclass responsibility\n";
    exit(-1);
    return 0;
}

const char*
BilinearTS::getType (void) const
{
    opserr << "BilinearTS::getType -- subclass responsibility\n";
    exit(-1);

    return 0;
}

int
BilinearTS::getOrder (void) const
{
    opserr << "BilinearTS::getOrder -- subclass responsibility\n";
    exit(-1);
    return -1;
}

Response*
BilinearTS::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	const char *matType = this->getType();
    
	output.tag("NdMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());
    
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"slip") == 0 || strcmp(argv[0],"deformation") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0],"state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int BilinearTS::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStrain();
			return 0;
		case 3:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getState();
			return 0;
		default:
			return -1;
	}
}

int
BilinearTS::updateState (const Information &matInfo)
{
    Vector *delp = matInfo.theVector;
    
    // Felipe, need some logic on how to update sigc and delc for Bilinear
    opserr << "BilinearTS::updateState recevied delp = " << delp << endln;
    
    return 0;
}

int
BilinearTS::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(5);

    data(0) = this->getTag();
    data(1) = delc;
    data(2) = sigc;
    data(3) = lamcr;
    data(4) = cmult;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "BilinearTS::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int
BilinearTS::recvSelf (int commitTag, Channel &theChannel,
				    FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static Vector data(5);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "BilinearTS::recvSelf -- could not recv Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    delc = data(1);
    sigc = data(2);
    lamcr = data(3);
    cmult = data(4);
    
    return res;
}

void
BilinearTS::Print (OPS_Stream &s, int flag)
{
	s << "Bilinear Traction Slip Material Model" << endln;
	s << "\tdelc: " << delc << ", sigc: " << sigc << endln;
    s << "\tlamcr: " << lamcr << endln;

	return;
}

void
BilinearTS::Shear_Envlp (double Delt, double Deln,
                            double &Tt, double &ETt, double &ETn)
{
    // shear monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delc;
    double Delnbar = Deln/delc;
    double lameff = sqrt(Deltbar*Deltbar+Delnbar*Delnbar);
    
    if (lameff <= lamcr) {
        Tt = sigc/lamcr*Deltbar;
        ETn = 0;
        ETt = sigc/lamcr/delc;
    } else {
        Tt = sigc/lameff*Deltbar*(1-lameff)/(1-lamcr);
        ETn = -delc*sigc/(1-lamcr)*Deltbar/delc/lameff*Delnbar/delc/lameff - (1-lameff)*delc*sigc/(1-lamcr)*Deltbar*Delnbar/delc/delc/pow(lameff,3);
        ETt = -delc*sigc/(1-lamcr)*pow(Deltbar/delc/lameff,2) + (1-lameff)*delc*sigc/(1-lamcr)*(1/lameff/delc/delc - Deltbar*Deltbar/delc/delc/pow(lameff,3));
    }
    return;
}

void
BilinearTS::Normal_Envlp (double Delt, double Deln,
                             double &Tn, double &ENt, double &ENn)
{
    // shear monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delc;
    double Delnbar = Deln/delc;
    double lameff = sqrt(Deltbar*Deltbar+Delnbar*Delnbar);
    
    if (lameff <= lamcr) {
        Tn = sigc/lamcr*Delnbar;
        ENn = 0;
        ENt = sigc/lamcr/delc;
    } else {
        Tn = sigc/lameff*Delnbar*(1-lameff)/(1-lamcr);
        ENn = -delc*sigc/(1-lamcr)*pow(Delnbar/delc/lameff,2) + (1-lameff)*delc*sigc/(1-lamcr)*(1/lameff/delc/delc - Delnbar*Delnbar/delc/delc/pow(lameff,3));
        ENt = -delc*sigc/(1-lamcr)*Deltbar/delc/lameff*Delnbar/delc/lameff - (1-lameff)*delc*sigc/(1-lamcr)*Deltbar*Delnbar/delc/delc/pow(lameff,3);
    }
    
    // compression normal multiplier
    // NYI for bilinear
    double cloc = 1.0;
    if (Deln < 0)
        cloc = cmult;
    
    return;
}
