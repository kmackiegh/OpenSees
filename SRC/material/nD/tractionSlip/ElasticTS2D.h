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

#ifndef ElasticTS2D_h
#define ElasticTS2D_h

// Written: KRM
// Created: 09/2019
//
// Description: 
//
// What: "@(#) ElasticTS2D.h, revA"

#include <ElasticTS.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticTS2D : public ElasticTS
{
  public:
    ElasticTS2D(int tag, double d1, double d2, double s1, double s2);
    ElasticTS2D();
    ~ElasticTS2D();

    const char *getClassType(void) const {return "ElasticTS2D";};

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
    
    static Vector stress;	// Stress vector ... class-wide for returns
    static Matrix tangent;	// Elastic constants ... class-wide for returns
    static Vector state;    // vector for recorders ... class-wide for returns
    
    // local plane stress storage
    Vector sigma;           // stress vector local
    Vector Tstress;         // trial stress vector
    Matrix D;               // stiffness matrix local
    Vector epsilon;         // Trial strains
    
    // committed storage
    Vector Cepsilon;	    // Committed strain
    Vector Cstress;         // Committed stress
    
    // some local storage for passing into envelopes
    double sigt;
    double sign;
    double ETt;
    double ETn;
    double ENt;
    double ENn;
    
};

#endif