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
                                                                        
#ifndef BilinearTS2D_h
#define BilinearTS2D_h

// Written: krm
// Created: 1/2020
//
// Description: 
//
// What: "@(#) BilinearTS2D.h, revA"

#include <BilinearTS.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class BilinearTS2D : public BilinearTS
{
  public:
    BilinearTS2D(int tag, double dc, double sc, double lc, double cm);
    BilinearTS2D();
    ~BilinearTS2D();

    const char *getClassType(void) const {return "BilinearTS2D";};

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
    
    // history variables
    double delmax;
    
};

#endif
