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
// Description: This file contains the class implementation for ElasticTS.
//
// What: "@(#) ElasticTS.C, revA"

#include <string.h>

#include <ElasticTS.h>
#include <ElasticTS2D.h>

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
OPS_NewElasticTS(void)
{
    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 5) {
        opserr << "Want: nDMaterial ElasticTS $tag $delt $deln $tau_max $sig_max" << endln;
        return 0;
    }

    int iData[1];
    double dData[4];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial ElasticTS \n";
        return 0;
    }

    numData = numArgs - 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial ElasticTS: " << iData[0] <<"\n";
        return 0;
    }

    theMaterial = new ElasticTS(iData[0], dData[0], dData[1], dData[2], dData[3]);

    return theMaterial;
}



ElasticTS::ElasticTS
(int tag, int classTag, double d1, double d2, double s1, double s2)
  :NDMaterial(tag, classTag), 
delt(d1), deln(d2), tau_max(s1), sig_max(s2)
{
    // do some input checks
    if (delt < 0)
        delt = fabs(delt);
    if (deln < 0)
        deln = fabs(deln);
    if (tau_max < 0)
        tau_max = fabs(tau_max);
    if (sig_max < 0)
        sig_max = fabs(sig_max);
    
    // derived properties
    kit = tau_max/delt;
    kin = sig_max/deln;
}

ElasticTS::ElasticTS
(int tag, double d1, double d2, double s1, double s2)
  :NDMaterial(tag, ND_TAG_ElasticTS),
delt(d1), deln(d2), tau_max(s1), sig_max(s2)
{
    // derived properties
    kit = tau_max/delt;
    kin = sig_max/deln;
}

ElasticTS::~ElasticTS()
{
	
}

double
ElasticTS::getRho()
{ 
    return 0;
}

NDMaterial*
ElasticTS::getCopy (const char *type)
{
    if (strcmp(type,"2D") == 0 || strcmp(type,"2d") == 0) {
        ElasticTS2D *theModel;
        theModel = new ElasticTS2D (this->getTag(), delt, deln, tau_max, sig_max);

        return theModel;
    }

    // Handle other cases
    else
        return NDMaterial::getCopy(type);
}

int
ElasticTS::setTrialStrain (const Vector &v)
{
    opserr << "ElasticTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticTS::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "ElasticTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticTS::setTrialStrainIncr (const Vector &v)
{
    opserr << "ElasticTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticTS::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "ElasticTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
ElasticTS::getTangent (void)
{
    opserr << "ElasticTS::getTangent -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Matrix *ret = new Matrix();
    return *ret;
}

const Matrix&
ElasticTS::getInitialTangent (void)
{
    return this->getTangent();
}

const Vector&
ElasticTS::getStress (void)
{
    opserr << "ElasticTS::getStress -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
ElasticTS::getStrain (void)
{
    opserr << "ElasticTS::getStrain -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
ElasticTS::getState (void)
{
    opserr << "ElasticTS::getState -- subclass responsibility\n";
    exit(-1);
    
    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

int
ElasticTS::commitState (void)
{
    opserr << "ElasticTS::commitState -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticTS::revertToLastCommit (void)
{
    opserr << "ElasticTS::revertToLastCommit -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
ElasticTS::revertToStart (void)
{
    opserr << "ElasticTS::revertToStart -- subclass responsibility\n";
    exit(-1);
    return -1;
}

NDMaterial*
ElasticTS::getCopy (void)
{
    opserr << "ElasticTS::getCopy -- subclass responsibility\n";
    exit(-1);
    return 0;
}

const char*
ElasticTS::getType (void) const
{
    opserr << "ElasticTS::getType -- subclass responsibility\n";
    exit(-1);

    return 0;
}

int
ElasticTS::getOrder (void) const
{
    opserr << "ElasticTS::getOrder -- subclass responsibility\n";
    exit(-1);
    return -1;
}

Response*
ElasticTS::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	const char *matType = this->getType();
    
	output.tag("NdMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());
    
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int ElasticTS::getResponse (int responseID, Information &matInfo)
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
ElasticTS::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(5);

    data(0) = this->getTag();
    data(1) = delt;
    data(2) = deln;
    data(3) = tau_max;
    data(4) = sig_max;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTS::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int
ElasticTS::recvSelf (int commitTag, Channel &theChannel,
				    FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static Vector data(5);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTS::recvSelf -- could not recv Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    delt = data(1);
    deln = data(2);
    tau_max = data(3);
    sig_max = data(4);
    
    return res;
}

void
ElasticTS::Print (OPS_Stream &s, int flag)
{
	s << "Elastic Traction Slip Material Model" << endln;
	s << "\tdelt: " << delt << ", deln: " << deln << endln;
    s << "\ttau_max: " << tau_max << ", sig_max: " << sig_max << endln;
    s << "\tkit: " << kit << ", kin: " << kin << endln;

	return;
}

void
ElasticTS::Shear_Envlp (double Delt, double Deln,
                        double &Tt, double &ETt, double &ETn)
{
    // shear monotonic envelope function
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    Tt = kit*Delt;
    ETn = 0;
    ETt = kit;
    
    return;
}

void
ElasticTS::Normal_Envlp (double Delt, double Deln,
                         double &Tn, double &ENt, double &ENn)
{
    // shear monotonic envelope function
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    Tn = kin*Deln;
    ENn = kin;
    ENt = 0;
    
    return;
}

