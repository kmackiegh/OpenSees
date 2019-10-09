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
                                                                        
#include <ElasticTS2D.h>
#include <Channel.h>
#include <Matrix.h>

#include <math.h>
#include <float.h>

Vector ElasticTS2D::stress(2);
Matrix ElasticTS2D::tangent(2,2);
Vector ElasticTS2D::state(1);

ElasticTS2D::ElasticTS2D
(int tag, double d1, double d2, double s1, double s2) :
 ElasticTS (tag, ND_TAG_ElasticTS2D,
                d1, d2, s1, s2),
 sigma(2), Tstress(2), D(2,2), epsilon(2),
 Cepsilon(2), Cstress(2)
{
    this->initialize();
}

ElasticTS2D::ElasticTS2D():
 ElasticTS (0, ND_TAG_ElasticTS2D,
                0.0, 0.0, 0.0, 0.0),
sigma(2), Tstress(2), D(2,2), epsilon(2),
Cepsilon(2), Cstress(2)
{
    this->initialize();
}

ElasticTS2D::~ElasticTS2D ()
{

}

int
ElasticTS2D::initialize(void)
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
    
    return 0;
}

int
ElasticTS2D::setTrialStrain (const Vector &strain)
{
    epsilon = strain;
    Vector deps = epsilon - Cepsilon;
    
    // trial stress
    Tstress = Cstress + D*deps;
    
    // loading condition
    Shear_Envlp(strain(0),strain(1),sigt,ETt,ETn);
    Normal_Envlp(strain(0),strain(1),sign,ENt,ENn);
    
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
ElasticTS2D::setTrialStrain (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrain(strain) ;
}

int
ElasticTS2D::setTrialStrainIncr (const Vector &strain)
{
    static Vector newStrain(2);
    newStrain = epsilon + strain;
    
    return this->setTrialStrain(newStrain);
}

int
ElasticTS2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrainIncr(strain);
}

const Matrix&
ElasticTS2D::getTangent (void)
{
    tangent = D;
    return tangent;
}

const Matrix&
ElasticTS2D::getInitialTangent (void)
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
ElasticTS2D::getStress (void)
{
    stress = sigma;
    return stress;
}

const Vector&
ElasticTS2D::getStrain (void)
{
    return epsilon;
}

const Vector&
ElasticTS2D::getState (void)
{
    // store quantities in output vector
    state(0) = 0;
    
    return state;
}

int
ElasticTS2D::commitState (void)
{
    Cepsilon = epsilon;
    Cstress = sigma;
    
    return 0;
}

int
ElasticTS2D::revertToLastCommit (void)
{
    epsilon = Cepsilon;
    sigma = Cstress;
    
    return 0;
}

int
ElasticTS2D::revertToStart (void)
{
    this->initialize();
    
    return 0;
}

NDMaterial*
ElasticTS2D::getCopy (void)
{
    ElasticTS2D *theCopy =
        new ElasticTS2D (this->getTag(), delt,deln,tau_max,sig_max);
  
    theCopy->sigma = sigma;
    theCopy->D = D;
    theCopy->epsilon = epsilon;
    
    theCopy->Cepsilon = Cepsilon;
    theCopy->Cstress = Cstress;
    
    return theCopy;
}

const char*
ElasticTS2D::getType (void) const
{
    return "2D";
}

int
ElasticTS2D::getOrder (void) const
{
  return 2;
}

int 
ElasticTS2D::sendSelf(int commitTag, Channel &theChannel)
{
    opserr << "ElasticTS2D::sendSelf()" << endln;
    static Vector data(5);
  
    // note this is incomplete, should send other vectors as well?
    data(0) = this->getTag();
    data(1) = Cepsilon(0);
    data(2) = Cepsilon(1);
    data(3) = Cstress(0);
    data(4) = Cstress(1);
  
    int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTS2D::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int 
ElasticTS2D::recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker)
{
    opserr << "ElasticTS2D::recvSelf()" << endln;
    static Vector data(5);
  
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTS2D::sendSelf -- could not send Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    Cepsilon(0) = data(1);
    Cepsilon(1) = data(2);
    Cstress(2) = data(3);
    Cstress(3) = data(4);

    epsilon = Cepsilon;
    stress = Cstress;

    return res;
}
