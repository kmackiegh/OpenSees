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
#include <algorithm>

void *
OPS_NewElasticTS(void)
{
    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 7) {
        opserr << "Want: nDMaterial ElasticTS tag? delt? deln? tau_max? sig_max? phi_ang? cmult?" << endln;
        return 0;
    }

    int iData[1];
    double dData[6];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial ElasticTS \n";
        return 0;
    }

    numData = numArgs - 1;;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial ElasticTS: " << dData[0] <<"\n";
        return 0;
    }

    theMaterial = new ElasticTS(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);

    return theMaterial;
}



ElasticTS::ElasticTS
(int tag, int classTag, double d1, double d2, double s1, double s2, double f, double cm)
  :NDMaterial(tag, classTag), 
delt(d1), deln(d2), tau_max(s1), sig_max(s2), phi_ang(f), cmult(cm)
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
    if (phi_ang < 0)
        phi_ang = fabs(phi_ang);
    if (cmult <= 0)
        cmult = 1.0;
    
    // derived properties
    //kit = tau_max/delt;
    kin = sig_max/deln;
}

ElasticTS::ElasticTS
(int tag, double d1, double d2, double s1, double s2, double f, double cm)
  :NDMaterial(tag, ND_TAG_ElasticTS),
delt(d1), deln(d2), tau_max(s1), sig_max(s2), phi_ang(f), cmult(cm)
{
    // derived properties
    //kit = tau_ult/delt;
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
        theModel = new ElasticTS2D (this->getTag(), delt, deln, tau_max, sig_max, phi_ang, cmult);

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
	else if (strcmp(argv[0],"slip") == 0 || strcmp(argv[0],"deformation") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0],"state") == 0)
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
ElasticTS::updateState (const Information &matInfo)
{
    //Vector *delp = matInfo.theVector;

	// store in a local variable
	elmf = *matInfo.theVector;
    
    // Felipe, need some logic on how to update sigc and delc for Elastic (or not at all)
   // opserr << "ElasticTS::updateState recevied force = " << elmf << endln;
    
    return 0;
}

int
ElasticTS::retrieveState (const Vector &elmf)
{   
    return 0;
}

int
ElasticTS::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(7);

    data(0) = this->getTag();
    data(1) = delt;
    data(2) = deln;
    data(3) = tau_max;
    data(4) = sig_max;
	data(5) = phi_ang;
	data(6) = cmult;

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

    static Vector data(7);

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
	phi_ang = data(5);
	cmult = data(6);
    
    return res;
}

void
ElasticTS::Print (OPS_Stream &s, int flag)
{
	s << "Elastic Traction Slip Material Model" << endln;
	s << "\tdelt: " << delt << ", deln: " << deln << endln;
    s << "\ttau_max: " << tau_max << ", sig_max: " << sig_max << endln;
    s << "\tphi_ang: " << phi_ang << ", kin: " << kin << endln;

	return;
}

void
ElasticTS::Shear_Envlp (double Delt, double Deln,
                        double &Tt, double &ETt, double &ETn)
{
    // shear monotonic envelope function
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
	double kit = 0; 
	update_ShearEnergy(kit);

    Tt = kit*Delt;
    ETn = 0;
    ETt = kit;

	//opserr << "ElasticTS::Shear_Envlp computed kit = " << kit << endln;
    
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

	// compression:
	if (Deln < 0) {
		Tn = cmult*kin*Deln;
		ENn = cmult*kin*Deln;
		ENt = 0;
	}
    
    return;
}

void
ElasticTS::update_ShearEnergy (double &kit_new)
{
	// shear monotonic update function
	// takes current kit and returns it updated
	// Mohr-Coulomb envelope tau = c + sigma*tan(phi_ang);
	
	//opserr << "ShearEnergy elmf = " << elmf << endln;	

	double pi = 3.14159265359;
	double phi = pi*phi_ang/180;

	double sig_max = 0;
	double sig_min = 0;
	double sigma = 0;
	if (elmf.Size() != 0) {
		//int i = intp(1);
		double area = elmf(2);
		
		sig_max = std::max(elmf(0),elmf(1))/area;
		sig_min = std::min(elmf(0),elmf(1))/area;
		sigma = elmf(1)/area;
	}

	// tension cut-off:
	if (sig_min > 0)
		sig_min = 0;
	if (sig_max > 0)
		sig_min = 0;
	if (sigma > 0)
		sigma = 0;
	
	// failure parameters
	double fang = pi/2 + phi;									
	//double fsig = (sig_max+sig_min)/2 + (sig_max-sig_min)/2 * cos(2*fang);
	double fsig = sigma;

	// Compression cap:
	// none for elasticTS

	double tau_ult = tau_max;					
	if (fsig <= 0) {										// if compression
        tau_ult = tau_max - fsig*tan(phi);					// tau @failure
        //opserr << "fsigma = " << fsig << " ftau = " << tau_ult/tau_max << endln;
	}

    kit_new = tau_ult/delt;

	return;
}
