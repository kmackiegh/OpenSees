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
                                                                        
// $Revision: 1.13 $
// $Date: 2006-09-05 21:21:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/LowTensionMaterial.h,v $
                                                                        
                                                                        
#ifndef LowTensionMaterial_h
#define LowTensionMaterial_h

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class LowTensionMaterial : public NDMaterial
{
  public:
    // Only called by subclasses to pass their tags to NDMaterialModel
    LowTensionMaterial (int tag, int classTag, 
                        double E, double Eh, double Es,
                        double nu, double fc, double ftc, double ftm,
                        double shr, double fres, double fcu, double epscu, double rat,
                        double rho, double tlim);

    // Called by clients
    LowTensionMaterial (int tag, 
                        double E, double Eh, double Es,
                        double nu, double fc, double ftc, double ftm,
                        double shr, double fres, double fcu, double epscu, double rat,
                        double rho, double tlim);

    // For parallel processing
    LowTensionMaterial (void);

    virtual ~LowTensionMaterial (void);

    virtual const char *getClassType(void) const {return "LowTensionMaterial";};

    virtual double getRho( ) ;

    virtual int setTrialStrain (const Vector &v);
    virtual int setTrialStrain (const Vector &v, const Vector &r);
    virtual int setTrialStrainIncr (const Vector &v);
    virtual int setTrialStrainIncr (const Vector &v, const Vector &r);
    virtual const Matrix &getTangent (void);
    virtual const Matrix &getInitialTangent (void);
    virtual const Vector &getStress (void);
    virtual const Vector &getStrain (void);
    virtual const Vector &getState (void);

    virtual int commitState (void);
    virtual int revertToLastCommit (void);
    virtual int revertToStart (void);
    
    // Create a copy of material parameters AND state variables
    // Called by GenericSectionXD
    virtual NDMaterial *getCopy (void);

    // Create a copy of just the material parameters
    // Called by the continuum elements
    virtual NDMaterial *getCopy (const char *type);

    // Return a string indicating the type of material model
    virtual const char *getType (void) const;

    virtual int getOrder (void) const;
    
    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag = 0);
    void Tens_Envlp(double epsc, double &sigc, double &Ect);
    void Comp_Envlp(double epsc, double &sigc, double &Ect);
    
  protected:
    double E;
    double Eh;
    double Es;
    double nu;
    double fc;
    double ftc;
    double ftm;
    double shr;
    double fres;
    
    // properties from Concrete02
    double fcu;   // stress at ultimate (crushing) strain
    double epscu; // ultimate (crushing) strain
    double rat;   // ratio between unloading slope at epscu and original slope
    
    // optional arguments
    double rho; //mass per unit 3D volume
    double tlim; //limits cracking angle to be within small range of this angle (in degrees between -90 and 90)
    
  private:

};


#endif
