#ifndef Inerter_h
#define Inerter_h

// this probably started as Reza's inerter uniaxialMaterial
// since modified by Apostolakis and Mackie

#include <UniaxialMaterial.h>

class Inerter : public UniaxialMaterial
{
public:
    Inerter(int tag, int Inerter_type, double C);
    Inerter();
    ~Inerter();
    
    const char *getClassType(void) const {return "Inerter";};
    
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStrainRate(void);
    double getStress(void);
    
    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);
    double getCurrentTime(void);
    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    
protected:
    
private:
    // Fixed Input Material Variables
    double C;
    int Inerter_type;
    
    // Trial State Variables
    double Tstrain; // Trial Strain
    double TVel;   // Trial Velocity
    double TAcc;   // Trial Acceleration
    double Tstress; // Trial Stress
//    double CurrentTime;
    
    // Committeed State Variables
    double Cstrain;  // Committed Strain
    double CVel;    // Committed Velocity
    double CAcc;    // Committed Acceleration
    double Cstress;  // Committed Stress
//    double CommitTime;
};

#endif
