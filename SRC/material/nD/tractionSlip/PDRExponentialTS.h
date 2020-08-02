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
                                                                
#ifndef PDRExponentialTS_h
#define PDRExponentialTS_h

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class PDRExponentialTS : public NDMaterial
{
  public:
    // Only called by subclasses to pass their tags to NDMaterialModel
    PDRExponentialTS (int tag, int classTag,
                   double d1, double d2, double s1, double s2, double fp, double l, double a, double fr, double sc, double b, double k);

    // Called by clients
    PDRExponentialTS (int tag,
                   double d1, double d2, double s1, double s2, double fp, double l, double a, double fr, double sc, double b, double k);

    // For parallel processing
    PDRExponentialTS (void);

    virtual ~PDRExponentialTS (void);

    virtual const char *getClassType(void) const {return "PDRExponentialTS";};

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
    int updateState (const Information &matInformation);

    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag = 0);
    void Shear_Envlp (double Delt, double Deln, double &Tt, double &ETt, double &ETn);
    void Normal_Envlp (double Delt, double Deln, double &Tn, double &ENt, double &ENn);
/**/
	void update_ShearEnergy (double &tau_ult, double &tau_res, double &phit_new);
	void residual_slip (double Delnbar, double tau_ult, double tau_res, double &del_res);
/**/
    
  protected:
	// passed as element info
	Vector elmf;

    // passed as arguments
    double delt;
    double deln;
    double tau_max;
    double sig_max;
    double phi_ang;
    double lambda;
    double alpha;
    double phi_res;
    double sig_cap;
    double beta;
    double kcmp;

    // derived
    double del_res;
	//double phit;
    double phin;
    
  private:

};


#endif
