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
// Description: This file contains the class implementation for PDRExponentialTS.
//
// What: "@(#) PDRExponentialTS.C, revA"

#include <string.h>

#include <PDRExponentialTS.h>
#include <PDRExponentialTS2D.h>

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
OPS_NewPDRExponentialTS(void)
{
    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 9) {
        opserr << "Want: nDMaterial PDRExponentialTS tag? delt? deln? tau_max? sig_max? phi_ang? lambda? alpha? sig_cap? beta? cmult?" << endln;
        return 0;
    }

    int iData[1];
    double dData[10];

    int numData = 1;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: nDMaterial PDRExponentialTS \n";
        return 0;
    }

    numData = numArgs - 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid data: nDMaterial PDRExponentialTS: " << iData[0] <<"\n";
        return 0;
    }

    theMaterial = new PDRExponentialTS(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9]);

    return theMaterial;
}



PDRExponentialTS::PDRExponentialTS
(int tag, int classTag, double d1, double d2, double s1, double s2, double fp, double l, double a, double sc, double b, double k)
  :NDMaterial(tag, classTag), 
delt(d1), deln(d2), tau_max(s1), sig_max(s2), phi_ang(fp), lambda(l), alpha(a), sig_cap(sc), beta(b), cmult(k)
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
    if (phi_ang < 0)
        phi_ang = fabs(phi_ang);
    if (sig_cap == 0)
        sig_cap = -1e6*sig_max;
    if (sig_cap > 0)
        sig_cap = -sig_cap;
    if (beta <= 0)
        beta = 1.0;
    if (cmult < 1.0)
        cmult = 1.0;
    
    // derived properties
    phit = delt/lambda*tau_max * exp(1);
    phin = deln/alpha*sig_max * exp(1);

    sigma_cmp = 0;
    
}

PDRExponentialTS::PDRExponentialTS
(int tag, double d1, double d2, double s1, double s2, double fp, double l, double a, double sc, double b, double k)
  :NDMaterial(tag, ND_TAG_PDRExponentialTS),
delt(d1), deln(d2), tau_max(s1), sig_max(s2), phi_ang(fp), lambda(l), alpha(a), sig_cap(sc), beta(b), cmult(k)
{
    // derived properties
    phit = delt/lambda*tau_max * exp(1);
    phin = deln/alpha*sig_max * exp(1);

    sigma_cmp = 0;
}

PDRExponentialTS::~PDRExponentialTS()
{
    
}

double
PDRExponentialTS::getRho()
{ 
    return 0;
}

NDMaterial*
PDRExponentialTS::getCopy (const char *type)
{
    if (strcmp(type,"2D") == 0 || strcmp(type,"2d") == 0) {
        PDRExponentialTS2D *theModel;
        theModel = new PDRExponentialTS2D (this->getTag(), delt, deln, tau_max, sig_max, phi_ang, lambda, alpha, sig_cap, beta, cmult);

        return theModel;
    }

    // Handle other cases
    else
        return NDMaterial::getCopy(type);
}

int
PDRExponentialTS::setTrialStrain (const Vector &v)
{
    opserr << "PDRExponentialTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
PDRExponentialTS::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "PDRExponentialTS::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
PDRExponentialTS::setTrialStrainIncr (const Vector &v)
{
    opserr << "PDRExponentialTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
PDRExponentialTS::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "PDRExponentialTS::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
PDRExponentialTS::getTangent (void)
{
    opserr << "PDRExponentialTS::getTangent -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Matrix *ret = new Matrix();
    return *ret;
}

const Matrix&
PDRExponentialTS::getInitialTangent (void)
{
    return this->getTangent();
}

const Vector&
PDRExponentialTS::getStress (void)
{
    opserr << "PDRExponentialTS::getStress -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
PDRExponentialTS::getStrain (void)
{
    opserr << "PDRExponentialTS::getStrain -- subclass responsibility\n";
    exit(-1);

    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

const Vector&
PDRExponentialTS::getState (void)
{
    opserr << "PDRExponentialTS::getState -- subclass responsibility\n";
    exit(-1);
    
    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

int
PDRExponentialTS::commitState (void)
{
    opserr << "PDRExponentialTS::commitState -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
PDRExponentialTS::revertToLastCommit (void)
{
    opserr << "PDRExponentialTS::revertToLastCommit -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
PDRExponentialTS::revertToStart (void)
{
    opserr << "PDRExponentialTS::revertToStart -- subclass responsibility\n";
    exit(-1);
    return -1;
}

NDMaterial*
PDRExponentialTS::getCopy (void)
{
    opserr << "PDRExponentialTS::getCopy -- subclass responsibility\n";
    exit(-1);
    return 0;
}

const char*
PDRExponentialTS::getType (void) const
{
    opserr << "PDRExponentialTS::getType -- subclass responsibility\n";
    exit(-1);

    return 0;
}

int
PDRExponentialTS::getOrder (void) const
{
    opserr << "PDRExponentialTS::getOrder -- subclass responsibility\n";
    exit(-1);
    return -1;
}

Response*
PDRExponentialTS::setResponse (const char **argv, int argc, OPS_Stream &output)
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

int PDRExponentialTS::getResponse (int responseID, Information &matInfo)
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
PDRExponentialTS::updateState (const Information &matInfo)
{
    elmf = *matInfo.theVector;
    
    double area = 0; 
    if (elmf.Size() != 0) {
        area = elmf(2);
        sigma_cmp = elmf(0)/area;
    }

    // tension cut-off:
    if (sigma_cmp > 0)
        sigma_cmp = 0;

    // Compression cap:
    if (sigma_cmp < sig_cap)
        sigma_cmp = sig_cap;

    //opserr << "\tsigma = " << elmf(0) << " area = " << area << " = " << sigma_cmp << endln;
        
    return 0;
}

int
PDRExponentialTS::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(11);

    data(0) = this->getTag();
    data(1) = delt;
    data(2) = deln;
    data(3) = tau_max;
    data(4) = sig_max;
    data(5) = phi_ang;
    data(6) = lambda;
    data(7) = alpha;
    data(8) = sig_cap;
    data(9) = beta;
    data(10) = cmult;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "PDRExponentialTS::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int
PDRExponentialTS::recvSelf (int commitTag, Channel &theChannel,
                    FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static Vector data(11);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "PDRExponentialTS::recvSelf -- could not recv Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    delt = data(1);
    deln = data(2);
    tau_max = data(3);
    sig_max = data(4);
    phi_ang = data(5);
    lambda = data(6);
    alpha = data(7);
    sig_cap = data(8);
    beta = data(9);
    cmult = data(10);
    
    return res;
}

void
PDRExponentialTS::Print (OPS_Stream &s, int flag)
{
    s << "Exponential Traction Slip Material Model" << endln;
    s << "\tdelt: " << delt << ", deln: " << deln << endln;
    s << "\ttau_max: " << tau_max << ", sig_max: " << sig_max << ", sig_cap: " << sig_cap << endln;
    s << "\tlambda: " << lambda << ", alpha: " << alpha << endln;
    s << "\tphi_ang: " << phi_ang << ", phin: " << phin << endln;
    s << "\tbeta: " << beta << endln;

    return;
}

void
PDRExponentialTS::Frict_Envlp (double Delt, double &Tr, double &ETr)
{
    // frcition monotonic function
    // takes current Delt and returns Tr and damage
    // Coulomb's Law envelope tau = sigma*tan(phi_ang);
    
    double Deltbar = Delt/delt;
    double e2 = exp(-lambda*fabs(Deltbar));
    
    double sgn = 1; 
    if (Delt < 0) {
       sgn = -1;
    }

    double pi = acos(-1.0);
    double phi = pi*phi_ang/180;

    Tr = 0;
    ETr = 0;                    
    if (sigma_cmp < 0) {                                     // if compression
        double d = sgn * (1 - e2);
        Tr = fabs(sigma_cmp)*tan(phi) * d;
        ETr = lambda * fabs(sigma_cmp)*tan(phi);
        //opserr << "\tsigma = " << sigma_cmp << " d = " << d <<  endln;
        //opserr << "\tTr = " << Tr << " ETr = " << ETr <<  endln;
    }

    return;
}

void
PDRExponentialTS::Shear_Envlp (double Delt, double Deln,
                            double &Tt, double &ETt, double &ETn)
{
    // shear monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delt;
    double Delnbar = Deln/deln;

    double sgn = 1; 
    if (Delt < 0)
       sgn = -1;
    
    if (Delnbar < 0)
        Delnbar = 0;

    double e1 = exp(-alpha*Delnbar);
    double e2 = exp(-lambda*fabs(Deltbar));

    double Tr = 0;
    double ETr = 0;
    Frict_Envlp(Delt, Tr, ETr);

    Tt = pow(lambda,2)*phit/delt * (1+alpha*Delnbar) * Deltbar * e1*e2 + Tr;
    ETn = -pow(alpha,2)*pow(lambda,2)*phit/(deln*delt) * e1*e2 * Delnbar*Deltbar;
    ETt = pow(lambda,2)*phit/pow(delt,2) * e1*e2 * (1+alpha*Delnbar) * (1-sgn*lambda*Deltbar) + ETr;

    return;
}

void
PDRExponentialTS::Normal_Envlp (double Delt, double Deln,
                             double &Tn, double &ENt, double &ENn)
{
    // normal monotonic envelope function
    // unloading and reloading reflected elsewhere
    // takes current Deln and Delt and returns Tt and dTt/Deln, dTt/Delt
    
    double Deltbar = Delt/delt;
    double Delnbar = Deln/deln;
    
    double sgn = 1; 
    if (Delt < 0)
       sgn = -1;
    if (Delnbar < 0)
        Deltbar = 0;
    
    double e1 = exp(-alpha*Delnbar);
    double e2 = exp(-lambda*fabs(Deltbar));
    
    Tn = pow(alpha,2)*phin/deln * Delnbar * (1+lambda*fabs(Deltbar)) * e1*e2;
    ENn = pow(alpha,2)*phin/pow(deln,2) * e1*e2 * (1+lambda*fabs(Deltbar))*(1-alpha*Delnbar);
    ENt = -pow(alpha,2)*pow(lambda,2)*phin/(deln*delt) * e1*e2 * Delnbar*Deltbar;

    // compression normal multiplier
    if (Deln < 0) {
        ENn = cmult * pow(alpha,2)*phin/pow(deln,2);
        Tn = ENn * Deln;
    }
    
    // smoothe the tangent in the compression transition
    double crit_dwidth = (1/alpha)*0.2;
    if (Delnbar > -crit_dwidth && Delnbar < crit_dwidth) {
        double temp = ENn;
        ENn = (1 - erf(2*Delnbar/crit_dwidth))/2 * temp*(cmult-1) + temp;
    }

    return;
}

