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
                                                                        
#include <PDRExponentialTS2D.h>
#include <Channel.h>
#include <Matrix.h>

#include <math.h>
#include <float.h>
#include <limits>

Vector PDRExponentialTS2D::stress(2);
Matrix PDRExponentialTS2D::tangent(2,2);
Vector PDRExponentialTS2D::state(1);

PDRExponentialTS2D::PDRExponentialTS2D
(int tag, double d1, double d2, double s1, double s2, double fp, double l, double a, double fr, double sc, double b, double k) :
 PDRExponentialTS (tag, ND_TAG_PDRExponentialTS2D,
                d1, d2, s1, s2, fp, l, a, fr, sc, b, k),
 sigma(2), Tstress(2), D(2,2), epsilon(2),
 Cepsilon(2), Cstress(2)
{
    this->initialize();
}

PDRExponentialTS2D::PDRExponentialTS2D():
 PDRExponentialTS (0, ND_TAG_PDRExponentialTS2D,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
sigma(2), Tstress(2), D(2,2), epsilon(2),
Cepsilon(2), Cstress(2)
{
    this->initialize();
}

PDRExponentialTS2D::~PDRExponentialTS2D ()
{

}

int
PDRExponentialTS2D::initialize(void)
{
    // initialize local storage
    sigma.Zero();
    Tstress.Zero();
    D.Zero();
    epsilon.Zero();
    
    Cepsilon.Zero();
    Cstress.Zero();
    
    // initialize stress
    sigt = 0;
    sign = 0;
    ETt = 0;
    ETn = 0;
    ENt = 0;
    ENn = 0;
    Shear_Envlp(0,0,sigt,ETt,ETn);
    Normal_Envlp(0,0,sign,ENt,ENn);
    
    // populate matrices
    sigma(0) = sigt;
    sigma(1) = sign;
    D(0,0) = ETt;
    D(0,1) = ETn;
    D(1,0) = ENt;
    D(1,1) = ENn;
    
    // history init
    delmax = 0;
    
    return 0;
}

int
PDRExponentialTS2D::setTrialStrain (const Vector &strain)
{
    epsilon = strain;
    Vector deps = epsilon - Cepsilon;
    
    // trial stress
    Tstress = Cstress + D*deps;
    
    // effective displacement check
    double deleff = sqrt(beta*beta*strain(0)*strain(0)+strain(1)*strain(1));
    //opserr << "LOADING (?) = " << deleff << "rat = " << deleff/delmax << endln;

    if (deleff >= delmax) {
        // loading condition
        Shear_Envlp(strain(0),strain(1),sigt,ETt,ETn);
        Normal_Envlp(strain(0),strain(1),sign,ENt,ENn);
        delmax = deleff;

    } else {
        // unloading or reloading
        double ct = delmax/deleff;
        double tsigt, tsign;
        double tETt, tETn, tENt, tENn;
        Shear_Envlp(strain(0)*ct,strain(1)*ct,tsigt,tETt,tETn);
        Normal_Envlp(strain(0)*ct,strain(1)*ct,tsign,tENt,tENn);

        // factors added to tangent, using old way
        double deffdn = strain(1)/deleff;
        double deffdt = strain(0)*beta*beta/deleff;
        ETt = 1/delmax*(deffdt*tsigt + deleff*tETt);
        ETn = 1/delmax*(deffdn*tsigt + deleff*tETn);
        ENt = 1/delmax*(deffdt*tsign + deleff*tENt);
        ENn = 1/delmax*(deffdn*tsign + deleff*tENn);

       
        // correct above stress values
        sigt = tsigt/ct;
        sign = tsign/ct;

	//opserr << "UNLOADING -RELOADING = " << deleff << "rat = " << deleff/delmax << endln;
    }
    
    // store in vector and matrix form from above
    sigma(0) = sigt;
    sigma(1) = sign;
    D(0,0) = ETt;
    D(0,1) = ETn;
    D(1,0) = ENt;
    D(1,1) = ENn;
    
    return 0;
}

int
PDRExponentialTS2D::setTrialStrain (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrain(strain) ;
}

int
PDRExponentialTS2D::setTrialStrainIncr (const Vector &strain)
{
    static Vector newStrain(2);
    newStrain = epsilon + strain;
    
    return this->setTrialStrain(newStrain);
}

int
PDRExponentialTS2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrainIncr(strain);
}

const Matrix&
PDRExponentialTS2D::getTangent (void)
{
    tangent = D;
    return tangent;
}

const Matrix&
PDRExponentialTS2D::getInitialTangent (void)
{
    sigt = 0;
    sign = 0;
    ETt = 0;
    ETn = 0;
    ENt = 0;
    ENn = 0;
    Shear_Envlp(0,0,sigt,ETt,ETn);
    Normal_Envlp(0,0,sign,ENt,ENn);
    
    // populate matrices
    D(0,0) = ETt;
    D(0,1) = ETn;
    D(1,0) = ENt;
    D(1,1) = ENn;
    
    tangent = D;
    return tangent;
}

const Vector&
PDRExponentialTS2D::getStress (void)
{
    stress = sigma;
    return stress;
}

const Vector&
PDRExponentialTS2D::getStrain (void)
{
    return epsilon;
}

const Vector&
PDRExponentialTS2D::getState (void)
{
    // store quantities in output vector
    state(0) = delmax;
    
    return state;
}

int
PDRExponentialTS2D::commitState (void)
{
    Cepsilon = epsilon;
    Cstress = sigma;
    
    return 0;
}

int
PDRExponentialTS2D::revertToLastCommit (void)
{
    epsilon = Cepsilon;
    sigma = Cstress;
    
    return 0;
}

int
PDRExponentialTS2D::revertToStart (void)
{
    this->initialize();
    
    return 0;
}

NDMaterial*
PDRExponentialTS2D::getCopy (void)
{
    PDRExponentialTS2D *theCopy =
        new PDRExponentialTS2D (this->getTag(), delt,deln,tau_max,sig_max,phi_ang,lambda,alpha,phi_res,sig_cap,beta,kcmp);
  
    theCopy->sigma = sigma;
    theCopy->D = D;
    theCopy->epsilon = epsilon;
    
    theCopy->Cepsilon = Cepsilon;
    theCopy->Cstress = Cstress;
    
    return theCopy;
}

const char*
PDRExponentialTS2D::getType (void) const
{
    return "2D";
}

int
PDRExponentialTS2D::getOrder (void) const
{
  return 2;
}

int 
PDRExponentialTS2D::sendSelf(int commitTag, Channel &theChannel)
{
    opserr << "PDRExponentialTS2D::sendSelf()" << endln;
    static Vector data(6);
  
    // note this is incomplete, should send other vectors as well?
    data(0) = this->getTag();
    data(1) = Cepsilon(0);
    data(2) = Cepsilon(1);
    data(3) = Cstress(0);
    data(4) = Cstress(1);
    data(5) = delmax;
  
    int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "PDRExponentialTS2D::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int 
PDRExponentialTS2D::recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker)
{
    opserr << "PDRExponentialTS2D::recvSelf()" << endln;
    static Vector data(6);
  
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "PDRExponentialTS2D::sendSelf -- could not send Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    Cepsilon(0) = data(1);
    Cepsilon(1) = data(2);
    Cstress(2) = data(3);
    Cstress(3) = data(4);
    delmax = data(5);

    epsilon = Cepsilon;
    stress = Cstress;

    return res;
}
