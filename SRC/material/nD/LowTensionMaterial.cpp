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
                                                                        
// $Revision: 1.25 $                                                              
// $Date: 2009-01-29 00:42:03 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/LowTensionMaterial.cpp,v $                                                                
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for LowTensionMaterial.
//
// What: "@(#) LowTensionMaterial.C, revA"

#include <string.h>

#include <LowTensionMaterial.h>
#include <LowTensionPlaneStress.h>
//#include <LowTensionPlaneStrain2D.h>
//#include <LowTensionAxiSymm.h>
//#include <LowTensionThreeDimensional.h>
//#include <LowTensionPlateFiber.h>
//#include <LowTensionBeamFiber.h>
//#include <LowTensionBeamFiber2d.h>

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
OPS_NewLowTensionMaterial(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 13) {
    opserr << "Want: nDMaterial LowTension $tag $E $Eh $Es $nu $fc $ftc $ftm $shr $fres $fcu $epscu $rat <$rho> <$angle_limit>" << endln;
    return 0;	
  }
  
  int iData[1];
  double dData[14];
  dData[12] = 0.0;
  dData[13] = 100.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial LowTension \n";
    return 0;
  }
  
  numData = numArgs - 1;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial LowTension: " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new LowTensionMaterial(iData[0], 
	dData[0], dData[1], dData[2],
	dData[3], dData[4], dData[5],
	dData[6], dData[7], dData[8],
    dData[9], dData[10], dData[11],
    dData[12], dData[13]);
  
  return theMaterial;
}



LowTensionMaterial::LowTensionMaterial
(int tag, int classTag, double e, double eh, double es,
 double nuxy, double fcmax, double cr, double fmax,
 double sr, double fr, double fcult, double ecu, double lam,
 double r, double ang_lim)
  :NDMaterial(tag, classTag), 
E(e), Eh(eh), Es(es),
nu(nuxy), fc(fcmax), ftc(cr), ftm(fmax),
shr(sr), fres(fr), fcu(fcult), epscu(ecu), rat(lam),
rho(r), tlim(ang_lim)
{

}

LowTensionMaterial::LowTensionMaterial
(int tag, double e, double eh, double es,
 double nuxy, double fcmax, double cr, double fmax,
 double sr, double fr, double fcult, double ecu, double lam,
 double r, double ang_lim)
  :NDMaterial(tag, ND_TAG_LowTension), 
E(e), Eh(eh), Es(es),
nu(nuxy), fc(fcmax), ftc(cr), ftm(fmax),
shr(sr), fres(fr), fcu(fcult), epscu(ecu), rat(lam),
rho(r), tlim(ang_lim)
{

}

LowTensionMaterial::~LowTensionMaterial()
{
	
}

double
LowTensionMaterial::getRho() 
{ 
  return rho;
}

NDMaterial*
LowTensionMaterial::getCopy (const char *type)
{
  if (strcmp(type,"PlaneStress") == 0 || strcmp(type,"PlaneStress2D") == 0) {
    LowTensionPlaneStress *theModel;
      theModel = new LowTensionPlaneStress (this->getTag(), E, Eh, Es,
                                            nu, fc, ftc, ftm, shr, fres,
                                            fcu, epscu, rat, rho, tlim);

    return theModel;
  }

  // Handle other cases
  else
    return NDMaterial::getCopy(type);
}

int
LowTensionMaterial::setTrialStrain (const Vector &v)
{
    opserr << "LowTensionMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
LowTensionMaterial::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "LowTensionMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
LowTensionMaterial::setTrialStrainIncr (const Vector &v)
{
    opserr << "LowTensionMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
LowTensionMaterial::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "LowTensionMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
LowTensionMaterial::getTangent (void)
{
  opserr << "LowTensionMaterial::getTangent -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Matrix *ret = new Matrix();
  return *ret;
}

const Matrix&
LowTensionMaterial::getInitialTangent (void)
{
  return this->getTangent();
}

const Vector&
LowTensionMaterial::getStress (void)
{
  opserr << "LowTensionMaterial::getStress -- subclass responsibility\n";
  exit(-1);
    
  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

const Vector&
LowTensionMaterial::getStrain (void)
{
  opserr << "LowTensionMaterial::getStrain -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

const Vector&
LowTensionMaterial::getState (void)
{
    opserr << "LowTensionMaterial::getState -- subclass responsibility\n";
    exit(-1);
    
    // Just to make it compile
    Vector *ret = new Vector();
    return *ret;
}

int
LowTensionMaterial::commitState (void)
{
  opserr << "LowTensionMaterial::commitState -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
LowTensionMaterial::revertToLastCommit (void)
{
  opserr << "LowTensionMaterial::revertToLastCommit -- subclass responsibility\n";
  exit(-1);
    
  return -1;
}

int
LowTensionMaterial::revertToStart (void)
{
  opserr << "LowTensionMaterial::revertToStart -- subclass responsibility\n";
  exit(-1);
  return -1;
}

NDMaterial*
LowTensionMaterial::getCopy (void)
{
  opserr << "LowTensionMaterial::getCopy -- subclass responsibility\n";
  exit(-1);
  return 0;
}

const char*
LowTensionMaterial::getType (void) const
{
  opserr << "LowTensionMaterial::getType -- subclass responsibility\n";
  exit(-1);	

  return 0;
}

int
LowTensionMaterial::getOrder (void) const
{
  opserr << "LowTensionMaterial::getOrder -- subclass responsibility\n";
  exit(-1);
  return -1;
}

Response*
LowTensionMaterial::setResponse (const char **argv, int argc, OPS_Stream &output)
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

int LowTensionMaterial::getResponse (int responseID, Information &matInfo)
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
LowTensionMaterial::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(15);

    data(0) = this->getTag();
    data(1) = E;
    data(2) = Eh;
    data(3) = Es;
    data(4) = nu;
    data(5) = fc;
    data(6) = ftc;
    data(7) = ftm;
    data(8) = shr;
    data(9) = fres;
    data(10) = fcu;
    data(11) = epscu;
    data(12) = rat;
    data(13) = rho;
    data(14) = tlim;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "LowTensionMaterial::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int
LowTensionMaterial::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static Vector data(15);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "LowTensionMaterial::recvSelf -- could not recv Vector\n";
        return res;
    }

    this->setTag((int)data(0));
    E = data(1);
    Eh = data(2);
    Es = data(3);
    nu = data(4);
    fc = data(5);
    ftc = data(6);
    ftm = data(7);
    shr = data(8);
    fres = data(9);
    fcu = data(10);
    epscu = data(11);
    rat = data(12);
    rho = data(13);
    tlim = data(14);

    return res;
}

void
LowTensionMaterial::Print (OPS_Stream &s, int flag)
{
	s << "Low Tension Material Model" << endln;
	s << "\tE: " << E << ", Eh: " << Eh << ", Es: " << Es << endln;
	s << "\tnu: " << nu << endln;
    s << "\tfc: " << fc << ", fcu: " << fcu << endln;
    s << "\tepscu: " << epscu << ", ratio: " << rat << endln;
	s << "\tftc: " << ftc << ", ftm: " << ftm << endln;
	s << "\tshr: " << shr << endln;
    s << "\tfres: " << fres << endln;
	s << "\trho: " << rho << endln;
    if ( tlim <= 90 && tlim >= -90 )
        s << "\tlimit on the cracking angle:  " << tlim << endln;

	return;
}

void
LowTensionMaterial::Tens_Envlp (double epsc, double &sigc, double &Ect)
{
    /*--------------------------------------------------------------------------
     ! monotonic envelope of concrete in tension from Concrete02
     ! returned variables
     !   sigc  = stress corresponding to eps
     !   Ect  = tangent concrete modulus
     !-------------------------------------------------------------------------*/
    
    double eps0 = ftc/E;
    double eps1 = eps0 + (ftm-ftc)/Eh;
    double epsu = eps1 + ftm/fabs(Es);
    
    if (epsc <= eps0) {
        sigc = epsc*E;
        Ect  = E;
    } else if (epsc <= eps1) {
        sigc = ftc + (epsc-eps0)*Eh;
        Ect = Eh;
    } else {
        if (epsc <= epsu) {
            Ect  = Es;
            sigc = ftm + Es*(epsc-eps1);
        } else {
            Ect  = 1.0e-10;
            sigc = 0.0;
        }
        
        // implement a residual stress if requested
        if (sigc < fres) {
            sigc = fres;
            Ect = 1.0e-10;
        }
    }
    return;
}

void
LowTensionMaterial::Comp_Envlp (double epsc, double &sigc, double &Ect)
{
    /*--------------------------------------------------------------------------
     ! monotonic envelope of concrete in compression from Concrete02
     ! returned variables
     !   sigc  = stress corresponding to eps
     !   Ect  = tangent concrete modulus
     !-------------------------------------------------------------------------*/
    
    // this would be for bilinear behavior
    //double eps0 = -fc/E;
    
    //if (epsc >= eps0) {
    //    sigc = epsc*E;
    //    Ect  = E;
    //} else {
    //    Ect  = 1.0e-10;
    //    sigc = -fc;
    //}
    
    // this is from Concrete02
    double epsc0  = 2.0*fc/E;
    double ratLocal = epsc/epsc0;
    // modify to get bilinear backbone, ratLocal not needed anymore
    epsc0 = fc/E;
    
    if (epsc >= epsc0) {
        // from Concrete02
        sigc = fc*ratLocal*(2.0-ratLocal);
        Ect  = E*(1.0-ratLocal);
        
        // modified
        sigc = E*epsc;
        Ect = E;
        
    } else {
        //   linear descending branch between epsc0 and epscu
        if (epsc > epscu) {
            sigc = (fcu-fc)*(epsc-epsc0)/(epscu-epsc0)+fc;
            Ect  = (fcu-fc)/(epscu-epsc0);
        } else {
            // flat friction branch for strains larger than epscu
            sigc = fcu;
            Ect  = 1.0e-10;
        }
    }
    
    return;
}

