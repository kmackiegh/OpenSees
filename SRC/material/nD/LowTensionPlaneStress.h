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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticOrthotropicPlaneStress2D.h,v $

#ifndef LowTensionPlaneStress_h
#define LowTensionPlaneStress_h

// Written: fmk
// Created: 10/11
//
// Description: 
//
// What: "@(#) LowTensionPlaneStress.h, revA"

#include <LowTensionMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class LowTensionPlaneStress : public LowTensionMaterial
{
  public:
    LowTensionPlaneStress(int tag, double E, double Eh, double Es,
                          double nu, double fc, double ftc, double ftm,
                          double shr, double fres, double fcu, double epscu, double rat,
                          double rho, double tlim);
    LowTensionPlaneStress();
    ~LowTensionPlaneStress();

    const char *getClassType(void) const {return "LowTensionPlaneStress";};

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);
    
    const Vector &getStress (void);
    const Vector &getStrain (void);
    const Vector &getState (void);
    
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    

  protected:

  private:
    int initialize(void);
    int rotate_principal(void);
    int update_principal(void);
    int update_tangent(void);
    int setUniaxialStrain(int indx, double &sig, double &e);
    
    static Vector stress;	// Stress vector ... class-wide for returns
    static Matrix tangent;	// Elastic constants ... class-wide for returns
    static Vector state;    // vector for recorders ... class-wide for returns
    
    // local plane stress storage
    Vector sigma;           // stress vector local
    Matrix D;               // stiffness matrix local
    Vector epsilon;         // Trial strains
    Vector sigprin;         // principal stress
    Vector epsprin;         // principal strain
    
    // committed storage
    Vector Cepsilon;	    // Committed strain
    Vector Cstress;         // Committed stress
    Vector Csigprin;        // Committed principal stress
    Vector Cepsprin;        // Committed principal strain
    Vector Ctangent;        // Committed principal tangent stiffness
    
    // storage for recording purpses
    Vector crack_point;     // strain state at which cracking occurred
    Vector crack_strain;    // principal strain after cracking
    
    double E1;
    double E2;
    double G;
    double G12;
    double nu12;
    double nu21;
    double rnu;
    double rG;
    
    double sig1;
    double sig2;
    double tau;
    
    double th;
    double th_cr;
    int cracker;
    
    // hstvP : Concrete HISTORY VARIABLES last committed step from Concrete02
    Vector ecminP;  //  hstP(1)
    Vector deptP;   //  hstP(2)
    
    // hstv : Concrete HISTORY VARIABLES  current step from Concrete02
    Vector ecmin;
    Vector dept;
    
};

#endif
