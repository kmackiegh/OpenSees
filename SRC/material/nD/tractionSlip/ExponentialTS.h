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
                                                                
#ifndef ExponentialTS_h
#define ExponentialTS_h

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ExponentialTS : public NDMaterial
{
  public:
    // Only called by subclasses to pass their tags to NDMaterialModel
    ExponentialTS (int tag, int classTag,
                   double d1, double d2, double s1, double s2, double l, double a, double b, double cmult);

    // Called by clients
    ExponentialTS (int tag,
                   double d1, double d2, double s1, double s2, double l, double a, double b, double k);

    // For parallel processing
    ExponentialTS (void);

    virtual ~ExponentialTS (void);

    virtual const char *getClassType(void) const {return "ExponentialTS";};

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
    int updateState (const Information &matInfo);

    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag = 0);
    void Shear_Envlp (double Delt, double Deln, double &Tt, double &ETt, double &ETn);
    void Normal_Envlp (double Delt, double Deln, double &Tn, double &ENt, double &ENn);
    
  protected:
    // passed as arguments
    double delt;
    double deln;
    double tau_max;
    double sig_max;
    double lambda;
    double alpha;
    double beta;
    double cmult;

    // derived
    double phit;
    double phin;
    
  private:

};


#endif
