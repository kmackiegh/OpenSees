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
// Description: This file contains the class implementation for ExponentialTS.
//
// What: "@(#) ExponentialTS.C, revA"

#include <string.h>

#include <ExponentialTS.h>
#include <ExponentialTS2D.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits>

void *
OPS_NewExponentialTS(void)
{
    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 9) {
        opserr << "Want: nDMaterial ExponentialTS tag? delt? deln? tau_max? sig_max? lambda? alpha? beta? cmult?" << endln;
        return 0;
    }

    int iData[1];
    double dData[8];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial ExponentialTS \n";
        return 0;
    }

    numData = numArgs - 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial ExponentialTS: " << iData[0] <<"\n";
        return 0;
    }

    theMaterial = new ExponentialTS(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7]);

    return theMaterial;
}



ExponentialTS::ExponentialTS
(int tag, int classTag, double d1, double d2, double s1, double s2, double l, double a, double b, double k)
  :NDMaterial(tag, classTag), 
delt(d1), deln(d2), tau_max(s1), sig_max(s2), lambda(l), alpha(a), beta(b), cmult(k)
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
    if (lambda < 0)
        lambda = 1.0;
    if (alpha < 0)
        alpha = 1.0;
    if (beta <= 0)
        beta = 1.0;
    if (cmult <= 0)
        cmult = 1;
    
    // derived properties
    phit = delt/lambda*tau_max * exp(1);
    phin = deln/alpha*sig_max * exp(1);


}

ExponentialTS::ExponentialTS
(int tag, double d1, double d2, double s1, double s2, double l, double a, double b, double k)
  :NDMaterial(tag, ND_TAG_ExponentialTS),
delt(d1), deln(d2), tau_max(s1), sig_max(s2), lambda(l), alpha(a), beta(b), cmult(k)
{
    // derived properties
    phit = delt/lambda*tau_max * exp(1);
    phin = deln/alpha*sig_max * exp(1);
}

ExponentialTS::~ExponentialTS()
{
	
}

double
ExponentialTS::getRho()
{ 
    return 0;
}

NDMaterial*
ExponentialTS::getCopy (const char *type)
{
    if (strcmp(type,"2D") == 0 || strcmp(type,"2d") == 0) {
        ExponentialTS2D *theModel;
        theModel = new ExponentialTS2D (this->getTag(), delt, deln, tau_max, sig_max, lambda, alpha, beta, cmult);

        return theModel;
    }

    // Handle other cases
    else
        return NDMaterial::getCopy(type);
}

int
ExponentialTS::setTrialStrain (const Vector &v)
{
    opserr << "ExponentialTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ExponentialTS::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "ExponentialTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ExponentialTS::setTrialStrainIncr (const Vector &v)
{
    opserr << "ExponentialTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ExponentialTS::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "ExponentialTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
ExponentialTS::getTangent (void)
{
    opserr << "ExponentialTS::getTangent -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Matrix *ret = new Matrix();
    return *ret;
}

const Matrix&
ExponentialTS::getInitialTangent (void)
{
    return this->getTangent();
}

const Vector&
ExponentialTS::getStress (void)
{
    opserr << "ExponentialTS::getStress -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
ExponentialTS::getStrain (void)
{
    opserr << "ExponentialTS::getStrain -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
ExponentialTS::getState (void)
{
    opserr << "ExponentialTS::getState -- subclass responsibility\n";
    exit(-1);
    
    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

int
ExponentialTS::commitState (void)
{
    opserr << "ExponentialTS::commitState -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ExponentialTS::revertToLastCommit (void)
{
    opserr << "ExponentialTS::revertToLastCommit -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
ExponentialTS::revertToStart (void)
{
    opserr << "ExponentialTS::revertToStart -- subclass responsibility\n";
    exit(-1);
    return -1;
}

NDMaterial*
ExponentialTS::getCopy (void)
{
    opserr << "ExponentialTS::getCopy -- subclass responsibility\n";
    exit(-1);
    return 0;
}

const char*
ExponentialTS::getType (void) const
{
    opserr << "ExponentialTS::getType -- subclass responsibility\n";
    exit(-1);

    return 0;
}

int
ExponentialTS::getOrder (void) const
{
    opserr << "ExponentialTS::getOrder -- subclass responsibility\n";
    exit(-1);
    return -1;
}

Response*
ExponentialTS::setResponse (const char **argv, int argc, OPS_Stream &output)
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

int ExponentialTS::getResponse (int responseID, Information &matInfo)
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
ExponentialTS::updateState (const Information &matInfo)
{
    Vector *delp = matInfo.theVector;

	// store in a local variable
	//elmf = *matInfo.theVector;
    
    // Felipe, need some logic on how to update sigc and delc for Elastic (or not at all)
    //opserr << "ElasticTS::updateState recevied force = " << elmf << endln;
    
    return 0;
}

int
ExponentialTS::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(9);

    data(0) = this->getTag();
    data(1) = delt;
    data(2) = deln;
    data(3) = tau_max;
    data(4) = sig_max;
    data(5) = lambda;
    data(6) = alpha;
    data(7) = beta;
    data(8) = cmult;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ExponentialTS::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int
ExponentialTS::recvSelf (int commitTag, Channel &theChannel,
				    FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static Vector data(9);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ExponentialTS::recvSelf -- could not recv Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    delt = data(1);
    deln = data(2);
    tau_max = data(3);
    sig_max = data(4);
    lambda = data(5);
    alpha = data(6);
    beta = data(7);
    cmult = data(8);
    
    return res;
}

void
ExponentialTS::Print (OPS_Stream &s, int flag)
{
	s << "Exponential Traction Slip Material Model" << endln;
	s << "\tdelt: " << delt << ", deln: " << deln << endln;
    s << "\ttau_max: " << tau_max << ", sig_max: " << sig_max << endln;
    s << "\tlambda: " << lambda << ", alpha: " << alpha << endln;
    s << "\tphit: " << phit << ", phin: " << phin << endln;
    s << "\tbeta: " << beta << ", cmult: " << cmult << endln;

	return;
}

void
ExponentialTS::Shear_Envlp (double Delt, double Deln,
                            double &Tt, double &ETt, double &ETn)
{
    // shear monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delt;
    double Delnbar = Deln/deln;
	
	//Delnbar = 0;
	if (Delnbar < 0)
		Delnbar = 0;

    double e1 = exp(-alpha*Delnbar);
    double e2 = exp(-lambda*fabs(Deltbar));
    
	double sgn = 1;	
    if (Delt < 0)
	   sgn = -1;

    Tt = pow(lambda,2)*phit/delt * (1+alpha*Delnbar) * Deltbar * e1*e2;
    ETn = -pow(alpha,2)*pow(lambda,2)*phit/(deln*delt) * e1*e2 * Delnbar*Deltbar;
    ETt = pow(lambda,2)*phit/pow(delt,2) * e1*e2 * (1+alpha*Delnbar) * (1-sgn*lambda*Deltbar);

    return;
}

void
ExponentialTS::Normal_Envlp (double Delt, double Deln,
                             double &Tn, double &ENt, double &ENn)
{
    // normal monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delt;
    double Delnbar = Deln/deln;
	
	//Deltbar = 0;
	if (Delnbar < 0)
		Deltbar = 0;
    
    double e1 = exp(-alpha*Delnbar);
    double e2 = exp(-lambda*fabs(Deltbar));
	
    Tn = pow(alpha,2)*phin/deln * Delnbar * (1+lambda*fabs(Deltbar)) * e1*e2;
    ENn = pow(alpha,2)*phin/pow(deln,2) * e1*e2 * (1+lambda*fabs(Deltbar))*(1-alpha*Delnbar);
	ENt = -pow(alpha,2)*pow(lambda,2)*phin/(deln*delt) * e1*e2 * Delnbar*Deltbar;

    // compression normal multiplier
    if (Deln < 0) {
		//Tn = cmult * pow(alpha,2)*phin/deln * Delnbar;
        // isn't the line below incorrect, should be deln squared?
		ENn = cmult * pow(alpha,2)*phin/pow(deln,2);
        Tn = ENn * Deln;
        // Felipe is the line below necessary? You already said Deltbar is 0 if Delnbar < 0?
		//ENt = cmult * 0;
    }
    
    // smoothe the tangent in the compression transition
    double crit_dwidth = 0.2;
    if (Delnbar > -crit_dwidth && Delnbar < crit_dwidth) {
        double temp = ENn;
        ENn = (1 - erf(2*Delnbar/(2*crit_dwidth)))/2 * temp*(cmult-1) + temp;
    }

    return;
}

