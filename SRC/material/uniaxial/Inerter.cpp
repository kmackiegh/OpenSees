#include <Inerter.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <string.h>

#include <Domain.h>
#include <OPS_Globals.h>
#include <elementAPI.h>


void *
OPS_Inerter()
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    //
    // parse the input line for the material parameters
    //
    int iData[2];
    double dData[1];

    int numData;
    numData = 2;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial Inerter tag" << endln;
        return 0;
    }

    numData = 1;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid C uniaxialMaterial Inerter\n";
        return 0;
    }

    //
    // create a new material
    //
    theMaterial = new Inerter(iData[0], iData[1], dData[0]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type Inerter\n";
        return 0;
    }

    // return the material
    return theMaterial;
}

Inerter::Inerter(int tag, int Inerter_type, double c)
:UniaxialMaterial(tag, 0),
C(c), Inerter_type(Inerter_type), Tstrain(0.0), TVel(0.0), TAcc(0.0), Tstress(0.0), Cstrain(0.0), CVel(0.0), CAcc(0.0), Cstress(0.0)
{
//    this->revertToStart();
}

Inerter::Inerter()
:UniaxialMaterial(0, 0),
C(0.0), Inerter_type(0), Tstrain(0.0), TVel(0.0), TAcc(0.0), Tstress(0.0), Cstrain(0.0), CVel(0.0), CAcc(0.0), Cstress(0.0)
{
//    this->revertToStart();
}

Inerter::~Inerter()
{
    // does nothing
}

int
Inerter::setTrialStrain(double strain, double strainRate)
{
//    this->revertToLastCommit();

    Tstrain = strain;
    TVel = strainRate;
    
    double nextStrain = 0;
    nextStrain = Tstrain + (TVel*ops_Dt);
    
//    double alfa =1;
//    double alfaM =1;
//    double alfaF =0.9;
//
//    double gamma = 1.5-alpha;
//    double beta = 0.25*(2-alpha)*(2-alpha);
    
//    double gamma = 0.5+alfaM-alfaF;
//    double beta = 0.25*(1+alfaM-alfaF)*(1+alfaM-alfaF);
    
    if (Inerter_type == 1) {
        if ( fabs(TVel) < 1.0e-9 ) {
            TAcc = 0.0;
        } else {
//            TAcc = (TVel-CVel)/ops_TheActiveDomain;
            TAcc = ((1/(0.25*ops_Dt*ops_Dt))*(Tstrain-Cstrain))-((1/(0.25*ops_Dt))*CVel)-CAcc;
//            TAcc = ((1/(0.5*ops_Dt))*(TVel-CVel))-CAcc;
//            TAcc = (nextStrain-(2*Tstrain)+Cstrain)/(ops_Dt*ops_Dt);
//            TAcc = ((1/(beta*ops_Dt*ops_Dt))*(Tstrain-Cstrain))-((1/(beta*ops_Dt))*CVel)-((0.5-beta)/beta)*CAcc;

        }

		Tstress = C*TAcc;

	} else {
            TAcc = (TVel-CVel)/ops_Dt;
//            TAcc = ((1/(0.25*ops_Dt*ops_Dt))*(Tstrain-Cstrain))-((1/(0.25*ops_Dt))*CVel)-CAcc;
//            TAcc = ((1/(0.5*ops_Dt))*(TVel-CVel))-CAcc;
//            TAcc = (nextStrain-(2*Tstrain)+Cstrain)/(ops_Dt*ops_Dt);
//            TAcc = ((1/(beta*ops_Dt*ops_Dt))*(Tstrain-Cstrain))-((1/(beta*ops_Dt))*CVel)-((0.5-beta)/beta)*CAcc;

		if (TAcc/CVel > 0) {
			Tstress = C*TAcc;
		} else {
            TAcc = 0.0;
			Tstress = 0.0;
		}
	}

    return 0;
}

double
Inerter::getStress(void)
{
    if (Inerter_type == 1) {
        return Tstress;
    } else {
        if (TAcc/CVel > 0) {
            return Tstress;
        } else {
            return 0.0;
        }
    }
}

double
Inerter::getTangent(void)
{
//    return Tstress/TAcc;
    return 0.0;
    
}

double
Inerter::getInitialTangent(void)
{
//    return C;
    return 0.0;
}

double
Inerter::getDampTangent(void)
{
    
//    return 0.0;
//    return C;
    if (Inerter_type == 1) {
       return (C/ops_Dt);
    } else {
        if (TAcc/CVel > 0) {
            return C/ops_Dt;
        } else {
            return 0.0;
        }
    }
}

double
Inerter::getStrain(void)
{
    if (Inerter_type == 1) {
        return Tstrain;
    } else {
        if (TAcc/CVel > 0) {
            return Tstrain;
        } else {
            return 0.0;
        }
    }
}

double
Inerter::getStrainRate(void)
{
    if (Inerter_type == 1) {
        return TVel;
    } else {
        if (TAcc/CVel > 0) {
            return TVel;
        } else {
            return 0.0;
        }
    }
}

int
Inerter::commitState(void)
{
    Cstrain = Tstrain;
    CVel = TVel;
    Cstress = Tstress;
    CAcc = TAcc;

    return 0;
}

int
Inerter::revertToLastCommit(void)
{
    Tstrain = Cstrain;
    TVel = CVel;
    Tstress = Cstress;
    TAcc = CAcc;

    return 0;
}

int
Inerter::revertToStart(void)
{
    Tstrain = 0.0;
    TVel = 0.0;
    Tstress = 0.0;
    TAcc = 0.0;
    
    Cstrain = 0.0;
    CVel = 0.0;
    Cstress = 0.0;
    CAcc = 0.0;

    return 0;
}

UniaxialMaterial *
Inerter::getCopy(void)
{
    Inerter *theCopy = new Inerter(this->getTag(), Inerter_type, C);

    theCopy->Tstrain = Tstrain;
    theCopy->TVel = TVel;
    theCopy->Tstress = Tstress;
    theCopy->TAcc = TAcc;
    theCopy->Cstrain = Cstrain;
    theCopy->CVel = CVel;
    theCopy->Cstress = Cstress;
    theCopy->CAcc = CAcc;

    return theCopy;
}


int
Inerter::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(7);
    data(0) = this->getTag();
    data(1) = C;
    data(2) = Inerter_type;
    data(3) = Tstrain;
    data(4) = TVel;
    data(5) = Tstress;
    data(6) = TAcc;
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "Inerter::sendSelf() - failed to send data\n";

    return res;
}

int
Inerter::recvSelf(int cTag, Channel &theChannel,
                  FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(7);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "Inerter::recvSelf() - failed to receive data\n";
        C = 0;
        this->setTag(0);
    }
    else {
        this->setTag((int)data(0));
        C = data(1);
        Inerter_type = (int)data(2);
        Tstrain = data(3);
        TVel = data(4);
        Tstress = data(5);
        TAcc = data(6);

        Cstrain = Tstrain;
        CVel = TVel;
        Cstress = Tstress;
        CAcc = TAcc;
    }

    return res;
}


void
Inerter::Print(OPS_Stream &s, int flag)
{
    s << "Inerter tag: " << this->getTag() << endln;
    s << "  C: " << C << endln;
    s << "  Inerter_type: " << Inerter_type << endln;
}

int
Inerter::setParameter(const char **argv, int argc, Parameter &param)
{

    if (strcmp(argv[0],"C") == 0) {
        param.setValue(C);
        return param.addObject(1, this);
    }
    if (strcmp(argv[0],"Inerter_type") == 0) {
        param.setValue(Inerter_type);
        return param.addObject(2, this);
    }
    return -1;
}


int
Inerter::updateParameter(int parameterID, Information &info)
{
    switch(parameterID) {
        case 1:
            C = info.theDouble;
            return 0;
        case 2:
            Inerter_type = (int)info.theDouble;
            return 0;
        default:
            return -1;
    }
}

