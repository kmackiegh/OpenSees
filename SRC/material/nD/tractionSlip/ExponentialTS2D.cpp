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
                                                                        
#include <ExponentialTS2D.h>
#include <Channel.h>
#include <Matrix.h>

#include <math.h>
#include <float.h>

Vector ExponentialTS2D::stress(2);
Matrix ExponentialTS2D::tangent(2,2);
Vector ExponentialTS2D::state(4);

ExponentialTS2D::ExponentialTS2D
(int tag, double d1, double d2, double s1, double s2, double b) :
 ExponentialTS (tag, ND_TAG_ExponentialTS2D,
                d1, d2, s1, s2, b),
 sigma(2), Tstress(2), D(2,2), epsilon(2),
 Cepsilon(2), Cstress(2)
{
    this->initialize();
}

ExponentialTS2D::ExponentialTS2D():
 ExponentialTS (0, ND_TAG_ExponentialTS2D,
                0.0, 0.0, 0.0, 0.0, 0.0),
sigma(2), Tstress(2), D(2,2), epsilon(2),
Cepsilon(2), Cstress(2)
{
    this->initialize();
}

ExponentialTS2D::~ExponentialTS2D ()
{

}

int
ExponentialTS2D::initialize(void)
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
ExponentialTS2D::setTrialStrain (const Vector &strain)
{
    epsilon = strain;
    Vector deps = epsilon - Cepsilon;
    
    // trial stress
    Tstress = Cstress + D*deps;
    
    // effective displacement check
    double deleff = sqrt(beta*beta*strain(0)*strain(0)+strain(1)*strain(1));
    
    if (deleff >= delmax) {
        // loading condition
        Shear_Envlp(strain(0),strain(1),sigt,ETt,ETn);
        Normal_Envlp(strain(0),strain(1),sign,ENt,ENn);
        delmax = deleff;
    
    } else {
        // unloading or reloading
        double ct = delmax/deleff;
        double tETt, tETn, tENt, tENn;
        Shear_Envlp(strain(0)*ct,strain(1)*ct,sigt,tETt,tETn);
        Normal_Envlp(strain(0)*ct,strain(1)*ct,sign,tENt,tENn);
        
        // factors added to tangent
        ETt = beta*beta*strain(0)/delmax/deleff*sigt + strain(1)*strain(1)/deleff/deleff*tETt - beta*beta*strain(0)*strain(1)/deleff/deleff*tETn;
        ETn = strain(1)/delmax/deleff*sigt + beta*beta*strain(0)*strain(0)/deleff/deleff*tETn - strain(0)*strain(1)/deleff/deleff*tETt;
        ENt = beta*beta*strain(0)/delmax/deleff*sign + strain(1)*strain(1)/deleff/deleff*tENt - beta*beta*strain(0)*strain(1)/deleff/deleff*tENn;
        ENn = strain(1)/delmax/deleff*sign + beta*beta*strain(0)*strain(0)/deleff/deleff*tENn - strain(0)*strain(1)/deleff/deleff*tENt;
        
        // correct above stress values
        sigt = sigt/ct;
        sign = sign/ct;
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
ExponentialTS2D::setTrialStrain (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrain(strain) ;
}

int
ExponentialTS2D::setTrialStrainIncr (const Vector &strain)
{
    static Vector newStrain(2);
    newStrain = epsilon + strain;
    
    return this->setTrialStrain(newStrain);
}

int
ExponentialTS2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrainIncr(strain);
}

const Matrix&
ExponentialTS2D::getTangent (void)
{
    tangent = D;
    return tangent;
}

const Matrix&
ExponentialTS2D::getInitialTangent (void)
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
ExponentialTS2D::getStress (void)
{
    stress = sigma;
    return stress;
}

const Vector&
ExponentialTS2D::getStrain (void)
{
    return epsilon;
}

const Vector&
ExponentialTS2D::getState (void)
{
    double pi = acos(-1.0);
    
    // NYI
    
    return state;
}

int
ExponentialTS2D::commitState (void)
{
    Cepsilon = epsilon;
    Cstress = sigma;
    
    return 0;
}

int
ExponentialTS2D::revertToLastCommit (void)
{
    epsilon = Cepsilon;
    sigma = Cstress;
    
    return 0;
}

int
ExponentialTS2D::revertToStart (void)
{
    this->initialize();
    
    return 0;
}

NDMaterial*
ExponentialTS2D::getCopy (void)
{
    ExponentialTS2D *theCopy =
        new ExponentialTS2D (this->getTag(), delt,deln,tau_max,sig_max,beta);
  
    theCopy->sigma = sigma;
    theCopy->D = D;
    theCopy->epsilon = epsilon;
    
    theCopy->Cepsilon = Cepsilon;
    theCopy->Cstress = Cstress;
    
    return theCopy;
}

const char*
ExponentialTS2D::getType (void) const
{
    return "2D";
}

int
ExponentialTS2D::getOrder (void) const
{
  return 2;
}

int 
ExponentialTS2D::sendSelf(int commitTag, Channel &theChannel)
{
    opserr << "ExponentialTS2D::sendSelf()" << endln;
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
        opserr << "ExponentialTS2D::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int 
ExponentialTS2D::recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker)
{
    opserr << "ExponentialTS2D::recvSelf()" << endln;
    static Vector data(6);
  
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ExponentialTS2D::sendSelf -- could not send Vector\n";
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
