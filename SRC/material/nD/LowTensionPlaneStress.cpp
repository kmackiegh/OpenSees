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
                                                                        
#include <LowTensionPlaneStress.h>           
#include <Channel.h>
#include <Matrix.h>

#include <math.h>
#include <float.h>

Vector LowTensionPlaneStress::stress(3);
Matrix LowTensionPlaneStress::tangent(3,3);
Vector LowTensionPlaneStress::state(6);

LowTensionPlaneStress::LowTensionPlaneStress
(int tag, double E, double Eh, double Es,
 double nu, double fc, double ftc, double ftm,
 double shr, double fres, double fcu, double epscu, double rat,
 double rho, double tlim) :
 LowTensionMaterial (tag, ND_TAG_LowTensionPlaneStress,
                     E, Eh, Es, nu, fc, ftc, ftm, shr, fres, fcu, epscu, rat, rho, tlim),
 sigma(3), D(3,3), epsilon(3), sigprin(3), epsprin(3),
 Cepsilon(3), Cstress(3), Csigprin(3), Cepsprin(3), Ctangent(3),
 crack_point(3), crack_strain(2),
 ecminP(2), deptP(2), ecmin(2), dept(2)
{
    this->initialize();
}

LowTensionPlaneStress::LowTensionPlaneStress():
 LowTensionMaterial (0, ND_TAG_LowTensionPlaneStress,
                     0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
sigma(3), D(3,3), epsilon(3), sigprin(2), epsprin(2),
Cepsilon(3), Cstress(3), Csigprin(3), Cepsprin(3), Ctangent(3),
crack_point(3), crack_strain(2),
ecminP(2), deptP(2), ecmin(2), dept(2)
{
    this->initialize();
}

LowTensionPlaneStress::~LowTensionPlaneStress ()
{

}

int
LowTensionPlaneStress::initialize(void)
{
    // initialize local storage
    sigma.Zero();
    D.Zero();
    epsilon.Zero();
    sigprin.Zero();
    epsprin.Zero();
    
    Cepsilon.Zero();
    Cstress.Zero();
    Csigprin.Zero();
    Cepsprin.Zero();
    Ctangent.Zero();
    
    crack_point.Zero();
    crack_strain.Zero();
    
    double d = E/(1.0-nu*nu); 
    D(0,0) = d;
    D(0,1) = nu*d;
    D(1,0) = D(0,1);
    D(1,1) = d;
    D(2,2) = 0.5*(1.0-nu)*d;
    
    // set the new orthotropic constants
    G = E/(2.0*(1.0+nu));
    E1 = E;
    E2 = E;
    G12 = G;
    nu12 = nu;
    nu21 = nu;
    rnu = 1;
    rG = shr;
    
    // initialize stress
    sig1 = 0;
    sig2 = 0;
    tau = 0;
    
    // Concrete02 storage
    ecminP.Zero();
    deptP.Zero();
    ecmin.Zero();
    dept.Zero();
    
    // initialize angles, th_cr needs to be outside of -pi
    th = 0.0;
    th_cr = -1.5*3.1416;
    cracker = 0;
    
    // check stress values for consistency
    if (ftc > ftm) {
        opserr << "LowTensionPlaneStress::initialize()" << endln;
        opserr << "cannot have ftc > ftm, setting them to be equal" << endln;
        ftc = ftm;
    }
    
    if (fres > ftm) {
        opserr << "LowTensionPlaneStress::initialize()" << endln;
        opserr << "cannot have fres > ftm, setting them to be equal" << endln;
        fres = ftm;
    }
    
    // just for book keeping, make sure compression values are negative
    if (fc > 0)
        fc = -fc;
    if (fcu > 0)
        fcu = -fcu;
    if (epscu > 0)
        epscu = -epscu;
    
    // Es must be negative
    if (Es > 0)
        Es = -Es;
    
    return 0;
}

int
LowTensionPlaneStress::setTrialStrain (const Vector &strain)
{
    epsilon = strain;
    Vector deps = epsilon - Cepsilon;
    
    if (cracker == 0) {
        // uncracked
        sigma = D * epsilon;
        
        // compute trial stress using elastic tangent
        //double eps0 = deps(0);
        //double eps1 = deps(1);
        //double eps2 = deps(2);
        
        //sigma(0) = Cstress(0) + D(0,0)*eps0 + D(0,1)*eps1 + D(0,2)*eps2;
        //sigma(1) = Cstress(1) + D(1,0)*eps0 + D(1,1)*eps1 + D(1,2)*eps2;
        //sigma(2) = Cstress(2) + D(2,0)*eps0 + D(2,1)*eps1 + D(2,2)*eps2;
        
        // update principal stress and strain
        this->rotate_principal();
        
        // send principal directions to the constitutive models
        //setUniaxialStrain(0,sig1,E1);
        //setUniaxialStrain(1,sig2,E2);
        
        // update crack status
        this->update_principal();
        
    } else {
        // cracked
        
        // update principal stress and strain
        this->rotate_principal();
        
        // send principal directions to the constitutive models
        setUniaxialStrain(0,sig1,E1);
        setUniaxialStrain(1,sig2,E2);
    
        // first principal direction
        nu12 = E1/E*nu;
        if (nu12 < 0)
            nu12 = 0;
        G12 = G*rG;
        
        // second principal direction
        nu21 = E2/E*nu;
        if (nu21 < 0)
            nu21 = 0;
        G12 = G*rG;
        
        // add back shear stress
        tau = G12*epsprin(2);
        
        // update tangent and stress
        this->update_tangent();
    }
    
    // save tangent terms for later
    Ctangent(0) = E1;
    Ctangent(1) = E2;
    Ctangent(2) = G12;
    
    return 0;
}

int
LowTensionPlaneStress::setTrialStrain (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrain(strain) ;
}

int
LowTensionPlaneStress::setTrialStrainIncr (const Vector &strain)
{
    static Vector newStrain(3);
    newStrain = epsilon + strain;
    
    return this->setTrialStrain(newStrain);
}

int
LowTensionPlaneStress::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
    return this->setTrialStrainIncr(strain);
}

const Matrix&
LowTensionPlaneStress::getTangent (void)
{
    tangent = D;
    return tangent;
}

const Matrix&
LowTensionPlaneStress::getInitialTangent (void)
{
    return this->getTangent();
}

const Vector&
LowTensionPlaneStress::getStress (void)
{
    stress = sigma;
    return stress;
}

const Vector&
LowTensionPlaneStress::getStrain (void)
{
    return epsilon;
}

const Vector&
LowTensionPlaneStress::getState (void)
{
    double pi = acos(-1.0);
    
    // compute principal cracking strain if cracked (otherwise returns zero)
    if (th_cr < -pi) {
        crack_strain.Zero();
    } else {
        // crack_strain is in the principal axes but technically should also have the shear term
        Vector deps = epsilon - crack_point;
        double st = sin(th_cr);
        double ct = cos(th_cr);
        
        crack_strain(0) = deps(0)*ct*ct + deps(1)*st*st + deps(2)*st*ct;
        crack_strain(1) = deps(0)*st*st + deps(1)*ct*ct - deps(2)*st*ct;
        //crack_strain(2) = 2.0*(deps(1)-deps(0))*st*ct + deps(2)*(ct*ct-st*st);
    }
    
    // store quantities in output vector
    state(0) = sigprin(0);
    state(1) = sigprin(1);
    state(2) = th_cr;
    state(2) = th;
    state(3) = crack_strain(0);
    state(4) = crack_strain(1);
    state(3) = epsprin(0);
    state(4) = epsprin(1);
    
    return state;
}

int
LowTensionPlaneStress::commitState (void)
{
    Cepsilon = epsilon;
    Cstress = sigma;
    Csigprin = sigprin;
    Cepsprin = epsprin;
    
    ecminP = ecmin;
    deptP = dept;
    
    return 0;
}

int
LowTensionPlaneStress::revertToLastCommit (void)
{
    epsilon = Cepsilon;
    sigma = Cstress;
    sigprin = Csigprin;
    epsprin = Cepsprin;
    
    ecmin = ecminP;
    dept = deptP;
    
    // note that the constant E1, E2, etc may have been modified, don't have a way to
    // return them here directly other than recomputing the stress terms
    
    return 0;
}

int
LowTensionPlaneStress::revertToStart (void)
{
    this->initialize();
    
    return 0;
}

NDMaterial*
LowTensionPlaneStress::getCopy (void)
{
    LowTensionPlaneStress *theCopy =
        new LowTensionPlaneStress (this->getTag(), E, Eh, Es,
                               nu, fc, ftc, ftm, shr, fres, fcu, epscu, rat,
                               rho, tlim);
  
    theCopy->sigma = sigma;
    theCopy->D = D;
    theCopy->epsilon = epsilon;
    theCopy->sigprin = sigprin;
    theCopy->epsprin = epsprin;
    
    theCopy->Cepsilon = Cepsilon;
    theCopy->Cstress = Cstress;
    theCopy->Csigprin = Csigprin;
    theCopy->Cepsprin = Cepsprin;
    theCopy->Ctangent = Ctangent;
    
    theCopy->crack_point = crack_point;
    theCopy->crack_strain = crack_strain;
    
    theCopy->ecminP = ecminP;
    theCopy->deptP = deptP;
    theCopy->ecmin = ecmin;
    theCopy->dept = dept;
    
    // see about copying other values
    
    return theCopy;
}

const char*
LowTensionPlaneStress::getType (void) const
{
    return "PlaneStress";
}

int
LowTensionPlaneStress::getOrder (void) const
{
  return 3;
}

int 
LowTensionPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
    opserr << "LowTensionPlaneStress::sendSelf()" << endln;
    static Vector data(21);
  
    // note this is incomplete, should send other vectors as well?
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
    data(15) = Cepsilon(0);
    data(16) = Cepsilon(1);
    data(17) = Cepsilon(2);
    data(18) = Cepsilon(3);
    data(19) = Cepsilon(4);
    data(20) = Cepsilon(5);
  
    int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "LowTensionPlaneStress::sendSelf -- could not send Vector\n";
        return res;
    }

    return res;
}

int 
LowTensionPlaneStress::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
    opserr << "LowTensionPlaneStress::recvSelf()" << endln;
    static Vector data(21);
  
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "LowTensionPlaneStress::sendSelf -- could not send Vector\n";
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
    epsilon(0)=data(15);
    epsilon(1)=data(16);
    epsilon(2)=data(17);
    epsilon(3)=data(18);
    epsilon(4)=data(19);
    epsilon(5)=data(20);

    Cepsilon = epsilon;

    return res;
}

int
LowTensionPlaneStress::rotate_principal(void)
{
    // current stress theta on interval [-pi,+pi] radians
    th = 0.5*atan2(sigma(2),(sigma(0)-sigma(1))/2.0);
    
    if (cracker == 0) {
        // uncracked
        double detr = sqrt( 0.25*pow(sigma(0)-sigma(1),2) + sigma(2)*sigma(2) );
        sigprin(0) = (sigma(0)+sigma(1))/2.0 + detr;
        sigprin(1) = (sigma(0)+sigma(1))/2.0 - detr;
        sigprin(2) = 0;
        
        // through isotropic behavior
        epsprin(0) = sigprin(0)/E;
        epsprin(1) = sigprin(1)/E;
        epsprin(2) = 0;
        //opserr << "epsprin from sigprin = (" << epsprin(0) << ", " << epsprin(1) << ") and from actual computation (";
        
        // through actual plane strain
        //detr = sqrt( 0.25*pow(epsilon(0)-epsilon(1),2) + 0.25*epsilon(2)*epsilon(2) );
        //epsprin(0) = (epsilon(0)+epsilon(1))/2.0 + detr;
        //epsprin(1) = (epsilon(0)+epsilon(1))/2.0 - detr;
        //opserr << epsprin(0) << ", " << epsprin(1) << ")" << endln;
        
    } else {
        // cracked
        double st = sin(th_cr);
        double ct = cos(th_cr);
        
        sigprin(0) = sigma(0)*ct*ct + sigma(1)*st*st + 2.0*sigma(2)*st*ct;
        sigprin(1) = sigma(0)*st*st + sigma(1)*ct*ct - 2.0*sigma(2)*st*ct;
        sigprin(2) = -(sigma(0)-sigma(1))*st*ct + sigma(2)*(ct*ct-st*st);
        
        // these are just rotated stresses, not actually principal (there is a shear too)
        //sigprin(0) = (sigma(0)+sigma(1))/2.0 + (sigma(0)-sigma(1))/2.0*cos(2*th_cr) + sigma(2)*sin(2*th_cr);
        //sigprin(1) = (sigma(0)+sigma(1))/2.0 - (sigma(0)-sigma(1))/2.0*cos(2*th_cr) - sigma(2)*sin(2*th_cr);
        
        epsprin(0) = epsilon(0)*ct*ct + epsilon(1)*st*st + epsilon(2)*st*ct;
        epsprin(1) = epsilon(0)*st*st + epsilon(1)*ct*ct - epsilon(2)*st*ct;
        epsprin(2) = 2.0*(epsilon(1)-epsilon(0))*st*ct + epsilon(2)*(ct*ct-st*st);
        //opserr << "after cracking sigprin = (" << sigprin(0) << ", " << sigprin(1) << ") and epsprin = (" << epsprin(0) << ", " << epsprin(1) << ")" << endln;
        //epsprin(0) = (epsilon(0)+epsilon(1))/2.0 + (epsilon(0)-epsilon(1))/2.0*cos(2*th_cr) + epsilon(2)*sin(2*th_cr);
        //epsprin(1) = (epsilon(0)+epsilon(1))/2.0 - (epsilon(0)-epsilon(1))/2.0*cos(2*th_cr) - epsilon(2)*sin(2*th_cr);
    }
    
    return 0;
}

int
LowTensionPlaneStress::update_principal(void)
{
    double pi = acos(-1.0);
    
    // compute cracking angle to check if there was a limit specified
    // any tlim is considered out of range if > 90 or < -90. Tolerance
    // is 4 degrees on either side of tlim. Note only tlim is in degrees.
    int allow_change = 1;
    if (fabs(tlim) <= 90.0) {
        if ( fabs(tlim - th*180.0/pi) < 4.0 )
            allow_change = 1;
        else
            allow_change = 0;
    }
    
    // this could be cracker equals 0
    if (th_cr < -pi) {
        // uncracked
        if (sigprin(0) > ftc && allow_change == 1) {
            // store information at cracking
            th_cr = th;
            cracker = 1;
            
            opserr << "cracked at theta = " << th_cr << " and principal stress (" << sigprin(0)
                   << ", " << sigprin(1) << ")" << endln;
            crack_point(0) = epsilon(0);
            crack_point(1) = epsilon(1);
            crack_point(2) = epsilon(2);
            
            // update the uniaxial constitutive models with this point
            Cepsprin.Zero();
            Csigprin.Zero();
            ecminP = ecmin;
            deptP = dept;
            
            //opserr << "updating uniaxial with previous sigprin = " << sigprin << " and epsprin = " << epsprin;
            double sigc;
            double Ect;
            setUniaxialStrain(0,sigc,Ect);
            //opserr << " and after has sig1 = " << sigc;
            setUniaxialStrain(1,sigc,Ect);
            //opserr << " and sig2 = " << sigc << endln;
        }
    }
    
    return 0;
}

int
LowTensionPlaneStress::update_tangent(void)
{
    // update cracked terms
    D.Zero();
    D(0,0) = 1.0/E1;
    D(0,1) = -nu21/E2;
    D(1,0) = -nu12/E1;
    D(1,1) = 1.0/E2;
    D(2,2) = 1.0/G12;

    D.Invert(D);
    
    // transformation matrix
    Matrix TT(3,3);
    double cc=cos(th_cr);
    double ss=sin(th_cr);
    
    if (cracker == 0) {
        cc=cos(th);
        ss=sin(th);
    }
    
    TT(0,0) = cc*cc;
    TT(0,1) = ss*ss;
    TT(0,2) = 2.0*ss*cc;
    TT(1,0) = ss*ss;
    TT(1,1) = cc*cc;
    TT(1,2) = -2.0*ss*cc;
    TT(2,0) = -ss*cc;
    TT(2,1) = ss*cc;
    TT(2,2) = cc*cc-ss*ss;

    Matrix temp = D*TT;

    // change back to the global system
    TT.Invert(TT);
	D = TT * temp;
    
    // stress terms
    Vector temps(3);
    temps(0) = sig1;
    temps(1) = sig2;
    temps(2) = tau;
    sigma = TT * temps;
    
    return 0;
}

int
LowTensionPlaneStress::setUniaxialStrain(int indx, double &sig, double &e)
{
    // indx is the index into the vectors of stored quantities in each principal direction
    double ec0 = E;
    
    // retrieve concrete history variables
    ecmin(indx) = ecminP(indx);
    dept(indx) = deptP(indx);
    
    // calculate current strain
    double eps = epsprin(indx);
    double deps = eps - Cepsprin(indx);
    
    // default values
    sig = Csigprin(indx) + ec0 * deps;
    e = ec0;
    
    //opserr << endln << "setTrial on " << indx << " with eps = " << eps << ", and ecmin = " << ecmin(indx) << " and dept = " << dept(indx) << endln;
    
    //if (fabs(deps) < DBL_EPSILON) {
    //    sig = Csigprin(indx);
    //    e = ec0;
    //    return 0;
    //}
    
    // if the current strain is less than the smallest previous strain
    // call the monotonic envelope in compression and reset minimum strain
    if (eps < ecmin(indx)) {
        //opserr << "   compression envelope with ecmin = " << ecmin(indx) << endln;
        this->Comp_Envlp(eps, sig, e);
        ecmin(indx) = eps;
    
    } else {
        // else, if the current strain is between the minimum strain and ept
        // (which corresponds to zero stress) the material is in the unloading-
        // reloading branch and the stress remains between sigmin and sigmax
        double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
        double sigmr = ec0 * epsr;
        
        // calculate the previous minimum stress sigmm from the minimum
        // previous strain ecmin and the monotonic envelope in compression
        double sigmm;
        double dumy;
        this->Comp_Envlp(ecmin(indx), sigmm, dumy);
        
        // calculate current reloading slope Er (Eq. 2.35 in EERC Report)
        // calculate the intersection of the current reloading slope Er
        // with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report)
        double er = (sigmm - sigmr) / (ecmin(indx) - epsr);
        double ept = ecmin(indx) - sigmm / er;
        
        if (eps <= ept) {
            //opserr << "   eps <= ept (=" << ept << ") branch, which is compression side unload/reload" << endln;
            double sigmin = sigmm + er * (eps - ecmin(indx));
            double sigmax = er * 0.5 * (eps - ept);
            sig = Csigprin(indx) + ec0 * deps;
            e = ec0;
            //opserr << "   epsr = " << epsr << ", sigmr = " << sigmr << ", sigmm = " << sigmm << ", er = " << er
            //    << ", sigmin = " << sigmin << ", sigmax = " << sigmax << endln;
            
            if (sig <= sigmin) {
                //opserr << "      inside sigmin with sig = " << sig << endln;
                sig = sigmin;
                e = er;
            }
            if (sig >= sigmax) {
                //opserr << "      inside sigmax with sig = " << sig << endln;
                sig = sigmax;
                e = 0.5 * er;
            }
            
        } else {
            // else, if the current strain is between ept and epn
            // (which corresponds to maximum remaining tensile strength)
            // the response corresponds to the reloading branch in tension
            // Since it is not saved, calculate the maximum remaining tensile
            // strength sicn (Eq. 2.43 in EERC Report)
            double epn = ept + dept(indx);
            double sicn;
            //opserr << "   else branch with epn = " << epn << " and ept = " << ept << endln;
            
            if (eps <= epn) {
                //opserr << "      else branch plus eps <= epn" << endln;
                this->Tens_Envlp(dept(indx), sicn, e);
                //e = ec0;
                //sig = sicn - e * (epn - eps);
                //if (sig < 0) {
                //    sig = 0;
                //    e = 1.0e-10;
                //}
                
                if (dept(indx) > ftc/E) {
                    e = sicn / dept(indx);
                } else {
                    e = ec0;
                }
                sig = e * (eps - ept);
                
            } else {
                //opserr << "      else branch plus else" << endln;
                // else, if the current strain is larger than epn the response
                // corresponds to the tensile envelope curve shifted by ept
                double epstmp = eps - ept;
                this->Tens_Envlp(epstmp, sig, e);
                dept(indx) = eps - ept;
            }
        }
    }
    
    //opserr << endln << "setTrial yields sig = " << sig << " and tangent = " << e << endln;
    
    return 0;
}
