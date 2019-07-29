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
                                                                        
// $Revision: 1.27 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/InerterElement/InerterElement.cpp,v $

// Written: KRM
// Created: 7/2019
// Revision: A
//
// Description: This file contains the implementation for the InerterElement class.
//
// What: "@(#) InerterElement.C, revA"

#include "InerterElement.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <vector>

// initialise the class wide variables
Matrix InerterElement::InerterElementM6(6,6);
Matrix InerterElement::InerterElementM12(12,12);
Vector InerterElement::InerterElementV6(6);
Vector InerterElement::InerterElementV12(12);

void* OPS_InerterElement()
{
    int ndm = OPS_GetNDM();
    // first scan the command line to obtain eleID, iNode, jNode, and the
    // orientation of ele xPrime and yPrime not along the global x and y axis
    
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 7) {
        opserr << "WARNING too few arguments " <<
            "want - element InerterElement eleTag? iNode? jNode? inerterType? C? " <<
            "<-orient x1? x2? x3? y1? y2? y3?>\n";

        return 0;
    }

    // eleTag, iNode, jNode, inerterType
    int idata [4];
    numdata = 4;
    if (OPS_GetIntInput(&numdata,idata) < 0) {
        opserr << "WARNING: failed to get integer data\n";
        return 0;
    }
    
    // C
    double ddata [1];
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata,ddata) < 0) {
        opserr << "WARNING: failed to get double data\n";
        return 0;
    }

    // create the vectors for the element orientation
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

    // finally check the command line to see if user specified orientation
    int doRayleighDamping = 0;
    const char* type;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        if (strcmp(type,"-doRayleigh") == 0) {
            doRayleighDamping = 1;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                numdata = 1;
                if (OPS_GetIntInput(&numdata,&doRayleighDamping) < 0) {
                    opserr<<"WARNING: invalid integer\n";
                    return 0;
                }
            }
        } else if (strcmp(type,"-orient") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 6) {
                opserr<<"WARNING: insufficient orient values\n";
                return 0;
            }
            numdata = 3;
            if (OPS_GetDoubleInput(&numdata,&x(0)) < 0) {
                opserr<<"WARNING: invalid double input\n";
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata,&y(0)) < 0) {
                opserr<<"WARNING: invalid double input\n";
                return 0;
            }
        }
    }

    Element *theEle = 0;
    theEle = new InerterElement(idata[0], ndm, idata[1], idata[2], idata[3],
                                ddata[0], x, y, doRayleighDamping);
    
    return theEle;
}


//  Constructor
InerterElement::InerterElement(int tag, int dim, int Nd1, int Nd2,
            int iType, double Cin,
            const Vector &x, const Vector &yp, int doRayleigh)
 :Element(tag,0),
  connectedExternalNodes(2),
  dimension(dim), numDOF(0), transformation(3,3), useRayleighDamping(doRayleigh),
  theMatrix(0), theVector(0)
{
    // check type
    inerterType = iType;
    if (inerterType != 1 || inerterType != 2) {
        opserr << "InerterElement::InerterElement invalid inerterType input = " << inerterType << endln;
        exit(-1);
    }
    
    // check constant
    C = Cin;
    if (C < 0)
        C = fabs(C);

    // establish the connected nodes and set up the transformation matrix for orientation
    this->setUp( Nd1, Nd2, x, yp);
    
}


//   Constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
InerterElement::InerterElement(void)
  :Element(0,0),
  connectedExternalNodes(2),
  dimension(0), numDOF(0), transformation(3,3), useRayleighDamping(0),
  theMatrix(0), theVector(0)
{
    
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)
        opserr << "FATAL ZeroLength::ZeroLength - failed to create an ID of correct size\n";
    
}


//  Destructor:
//  delete must be invoked on any objects created by the object
InerterElement::~InerterElement()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to

}


int
InerterElement::getNumExternalNodes(void) const
{
    return 2;
}


const ID &
InerterElement::getExternalNodes(void)
{
    return connectedExternalNodes;
}


Node **
InerterElement::getNodePtrs(void)
{
  return theNodes;
}


int
InerterElement::getNumDOF(void)
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the InerterElement element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
InerterElement::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    // set default values for error conditions
    numDOF = 3;
    theMatrix = &InerterElementM6;
    theVector = &InerterElementV6;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ( theNodes[0] == 0 || theNodes[1] == 0 ) {
      if (theNodes[0] == 0) 
        opserr << "WARNING InerterElement::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING InerterElement::setDomain() - Nd2: " << Nd2 << " does not exist in ";

      opserr << "model for InerterElement ele: " << this->getTag() << endln;

      return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING InerterElement::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
                "have differing dof at ends for InerterElement " << this->getTag() << endln;
      return;
    }	

    // Check that length is zero within tolerance
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;
    
    vm = (v1<v2) ? v2 : v1;

    if (L > LENTOL*vm)
      opserr << "WARNING InerterElement::setDomain(): Element " << this->getTag() << " has L= " << L <<
                ", which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set the number of dof for element and set matrix and vector pointer
    if (dimension == 2 && dofNd1 == 3) {
        numDOF = 6;
        theMatrix = &InerterElementM6;
        theVector = &InerterElementV6;
        elemType  = D2N6;
    }
    else if (dimension == 3 && dofNd1 == 6) {
        numDOF = 12;
        theMatrix = &InerterElementM12;
        theVector = &InerterElementV12;
        elemType  = D3N12;
    }
    else {
        opserr << "WARNING InerterElement::setDomain cannot handle " << dimension <<
            "dofs at nodes in " << dofNd1 << " d problem\n";
        return;
    }
    
    // create the basic deformation-displacement transformation matrix for the element
    this->setTran1d( elemType );
   
}   	 


int
InerterElement::commitState()
{
    int code=0;

    // call element commitState to do any base class stuff
    if ((code = this->Element::commitState()) != 0) {
      opserr << "InerterElement::commitState () - failed in base class";
    }    

    return code;
}

int
InerterElement::revertToLastCommit()
{
    int code=0;
    return code;
}


int
InerterElement::revertToStart()
{   
    int code=0;
    return code;
}


int
InerterElement::update(void)
{
    double curd, curv, cura;

    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector diffD  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getTrialVel();
    const Vector& vel2  = theNodes[1]->getTrialVel();
    Vector diffV = vel2-vel1;
    const Vector& acc1  = theNodes[0]->getTrialAccel();
    const Vector& acc2  = theNodes[1]->getTrialAccel();
    Vector diffA = acc2-acc1;
    
    // loop over dofs
    for (int mat=0; mat<numDOF/2; mat++) {
        // compute strain and rate; set as current trial for material
        curd = this->computeCurrentStrain1d(mat,diffD);
        curv = this->computeCurrentStrain1d(mat,diffV);
        cura = this->computeCurrentStrain1d(mat,diffA);
        
        // set trial curd, curv, cura on element
        Tstress = 0;
        
        if (mat == 0) {
            if (inerterType == 1) {
                if (fabs(curv) > 1.0e-9)
                    Tstress = C*cura;
                
            } else if (inerterType == 2) {
                if (cura/curv > 0)
                    Tstress = C*cura;
                
            }
        }
    }

    return 0;
}

const Matrix &
InerterElement::getTangentStiff(void)
{
    double E = 0;

    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    
    // loop over dofs
    Matrix& tran = *t1d;
    for (int mat=0; mat<numDOF/2; mat++) {
        // get tangent for material
        //E = theMaterial1d[mat]->getTangent();

        // compute contribution of material to tangent matrix
        for (int i=0; i<numDOF; i++)
            for(int j=0; j<i+1; j++)
                stiff(i,j) +=  tran(mat,i) * E * tran(mat,j);

    }
    
    // complete symmetric stiffness matrix
    for (int i=0; i<numDOF; i++)
        for(int j=0; j<i; j++)
            stiff(j,i) = stiff(i,j);

    return stiff;
}


const Matrix &
InerterElement::getInitialStiff(void)
{
    double E = 0;

    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    
    // loop over dofs
    Matrix& tran = *t1d;
    for (int mat=0; mat<numDOF/2; mat++) {
        // get tangent for material
        //E = theMaterial1d[mat]->getInitialTangent();

        // compute contribution of material to tangent matrix
        for (int i=0; i<numDOF; i++)
            for(int j=0; j<i+1; j++)
                stiff(i,j) +=  tran(mat,i) * E * tran(mat,j);
        
    }

    // complete symmetric stiffness matrix
    for (int i=0; i<numDOF; i++)
        for(int j=0; j<i; j++)
            stiff(j,i) = stiff(i,j);

    return stiff;
}
    

const Matrix &
InerterElement::getDamp(void)
{
    // damp is a reference to the matrix holding the damping matrix
    Matrix& damp = *theMatrix;

    // zero damping matrix
    damp.Zero();

    // get Rayleigh damping matrix
    if (useRayleighDamping == 1) {
        damp = this->Element::getDamp();

    } else {
        // loop over dofs and add their damping tangents
        double eta = 0;
        Matrix& tran = *t1d;
        for (int mat=0; mat<numDOF/2; mat++) {
            // get tangent for material
            //eta = theMaterial1d[mat]->getDampTangent();
            
            // reza had this depend on acc/vel > 0 for type 2
            if (mat == 0)
                eta = C/ops_Dt;

            // compute contribution of material to tangent matrix
            for (int i=0; i<numDOF; i++)
                for(int j=0; j<i+1; j++)
                    damp(i,j) +=  tran(mat,i) * eta * tran(mat,j);

        }
    }

    // complete symmetric damping matrix
    for (int i=0; i<numDOF; i++)
        for(int j=0; j<i; j++)
            damp(j,i) = damp(i,j);

    return damp;
}


const Matrix &
InerterElement::getMass(void)
{
    // no mass
    theMatrix->Zero();
    return *theMatrix;
}


void 
InerterElement::zeroLoad(void)
{
    // does nothing now
}


int 
InerterElement::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr << "InerterElement::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
    return -1;
}


int 
InerterElement::addInertiaLoadToUnbalance(const Vector &accel)
{
    // does nothing as element has no mass yet!
    return 0;
}


const Vector &
InerterElement::getResistingForce()
{
    double force;

    // zero the residual
    theVector->Zero();

    // loop over dofs
    for (int mat=0; mat<numDOF/2; mat++) {
        // get resisting force for material
        //force = theMaterial1d[mat]->getStress();
        force = Tstress;

        // compute residual due to resisting force
        for (int i=0; i<numDOF; i++)
            (*theVector)(i)  += (*t1d)(mat,i) * force;

    }

    return *theVector;
}


const Vector &
InerterElement::getResistingForceIncInertia()
{	
    // this already includes damping forces from materials
    this->getResistingForce();

    // add the damping forces from rayleigh damping
    if (useRayleighDamping == 1) {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
            *theVector += this->getRayleighDampingForces();
        }
    }

    return *theVector;
}


int
InerterElement::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// InerterElement packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments

	// Make one size bigger so not a multiple of 3, otherwise will conflict
	// with classTags ID
	static ID idData(8);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = connectedExternalNodes(0);
	idData(4) = connectedExternalNodes(1);
    idData(5) = inerterType;
    idData(6) = C;
	idData(7) = useRayleighDamping;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "InerterElement::sendSelf -- failed to send ID data\n";
	  return res;
	}

	// Send the 3x3 direction cosine matrix, have to send it since it is only set
	// in the constructor and not setDomain()
	res += theChannel.sendMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
        opserr <<  "InerterElement::sendSelf -- failed to send transformation Matrix\n";
        return res;
	}

	return res;
}


int
InerterElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // InerterElement creates an ID, receives the ID and then sets the
  // internal data with the data in the ID

  static ID idData(8);

  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "InerterElement::recvSelf -- failed to receive ID data\n";
			    
    return res;
  }

  res += theChannel.recvMatrix(dataTag, commitTag, transformation);
  if (res < 0) {
    opserr << "InerterElement::recvSelf -- failed to receive transformation Matrix\n";
			    
    return res;
  }

  this->setTag(idData(0));
  dimension = idData(1);
  numDOF = idData(2);
  connectedExternalNodes(0) = idData(3);
  connectedExternalNodes(1) = idData(4);
    inerterType = idData(5);
    C = idData(6);
  useRayleighDamping = idData(7);
  
  return res;
}


int
InerterElement::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0 )
       return 0;

    static Vector v1(3);
    static Vector v2(3);

    float d1 = 1.0;

    if (displayMode == 1 || displayMode == 2) {

        theNodes[0]->getDisplayCrds(v1, fact);
        theNodes[1]->getDisplayCrds(v2, fact);
        d1 = Tstress;
        
    } else {
        theNodes[0]->getDisplayCrds(v1, 0.);
        theNodes[1]->getDisplayCrds(v2, 0.);

    }
    if (v1 != v2)
      return theViewer.drawLine(v1, v2, d1, d1);	
    else
      return theViewer.drawPoint(v1, d1, 10);	
}


void
InerterElement::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain= 0.0;
    double force = 0.0;
     
    for (int i=0; i<numDOF; i++)
	(*theVector)(i) = (*t1d)(0,i)*force;
    
    if (flag == OPS_PRINT_CURRENTSTATE) { // print everything
        s << "Element: " << this->getTag();
        s << " type: InerterElement  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1) << endln;
        s << " C: " << C << " for inerter type: " << inerterType << endln;
    }
     
    else if (flag == 1) {
        s << this->getTag() << "  " << strain << "  ";
    }

}

Response*
InerterElement::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","InerterElement");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    char outputData[10];

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            char outputData[20];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(numDOF));

    } else if ((strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) ||
	       (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)) {

        for (int i=0; i<numDOF/2; i++) {
            sprintf(outputData,"P%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, Vector(numDOF/2));

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	       strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"basicDeformation") == 0) {

            for (int i=0; i<numDOF/2; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 3, Vector(numDOF/2));

    } else if (strcmp(argv[0],"basicStiffness") == 0) {

            for (int i=0; i<numDOF/2; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 13, Matrix(numDOF/2,numDOF/2));

    } else if ((strcmp(argv[0],"defoANDforce") == 0) ||
        (strcmp(argv[0],"deformationANDforces") == 0) ||
        (strcmp(argv[0],"deformationsANDforces") == 0)) {
      
      int i;
      for (i=0; i<numDOF/2; i++) {
	sprintf(outputData,"e%d",i+1);
	output.tag("ResponseType",outputData);
      }
      for (i=0; i<numDOF/2; i++) {
	sprintf(outputData,"P%d",i+1);
	output.tag("ResponseType",outputData);
	    }
      theResponse = new ElementResponse(this, 4, Vector(numDOF));
      
    } else if ((strcmp(argv[0],"dampingForces") == 0) || (strcmp(argv[0],"rayleighForces") == 0)) {
            theResponse = new ElementResponse(this, 15, Vector(numDOF));
    }

    output.endTag();

    return theResponse;
}

int 
InerterElement::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    case 15:
      theVector->Zero();
      if (useRayleighDamping == 1) {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            *theVector += this->getRayleighDampingForces();      
      }
      return eleInformation.setVector(*theVector);

    case 2:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < numDOF/2; i++) {
                //(*(eleInformation.theVector))(i) = theMaterial1d[i]->getStress();
            }
        }
        return 0;

    case 3:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < numDOF/2; i++) {
                //(*(eleInformation.theVector))(i) = theMaterial1d[i]->getStrain();
            }
        }
        return 0;

    case 13:
        if (eleInformation.theMatrix != 0) {
            for (int i = 0; i < numDOF/2; i++) {
                //(*(eleInformation.theMatrix))(i,i) = theMaterial1d[i]->getTangent();
            }
        }
        return 0;

    case 4:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < numDOF/2; i++) {
                //(*(eleInformation.theVector))(i) = theMaterial1d[i]->getStrain();
                //(*(eleInformation.theVector))(i+numDOF/2) = theMaterial1d[i]->getStress();
            }
        }
        return 0;      

    default:
        return -1;
    }
}

int
InerterElement::setParameter(const char **argv, int argc, Parameter &param)
{
  int result = -1;  

  if (argc < 1)
    return -1;

  // need to setparameter on C and other things in the class
    
  return result;
}


// Establish the external nodes and set up the transformation matrix
// for orientation
void
InerterElement::setUp( int Nd1, int Nd2,
		   const Vector &x,
		   const Vector &yp )
{ 
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL InerterElement::setUp - failed to create an ID of correct size\n";
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    // check that vectors for orientation are correct size
    if ( x.Size() != 3 || yp.Size() != 3 )
	opserr << "FATAL InerterElement::setUp - incorrect dimension of orientation vectors\n";

    // establish orientation of element for the transformation matrix
    // z = x cross yp
    Vector z(3);
    z(0) = x(1)*yp(2) - x(2)*yp(1);
    z(1) = x(2)*yp(0) - x(0)*yp(2);
    z(2) = x(0)*yp(1) - x(1)*yp(0);

    // y = z cross x
    Vector y(3);
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);

    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();

    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0) {
      opserr << "FATAL InerterElement::setUp - invalid vectors to constructor\n";
    }
    
    // create transformation matrix of direction cosines
    for (int i=0; i<3; i++ ) {
        transformation(0,i) = x(i)/xn;
        transformation(1,i) = y(i)/yn;
        transformation(2,i) = z(i)/zn;
    }

}


// Set basic deformation-displacement transformation matrix
void
InerterElement::setTran1d(Etype elemType)
{
    enum Dtype { TRANS, ROTATE };
    
    int   indx, dir;
    Dtype dirType;
    
    // Create 1d transformation matrix
    t1d = new Matrix(numDOF/2,numDOF);
    
    if (t1d == 0)
	opserr << "FATAL InerterElement::setTran1d - can't allocate 1d transformation matrix\n";
    
    // Use reference for convenience and zero matrix.
    Matrix& tran = *t1d;
    tran.Zero();
    
    // loop over materials, setting row in tran for each material depending on dimensionality of element
    for ( int i=0; i<numDOF/2; i++ ) {
	
	dir  = i;	        // direction 0 to 5;
	indx = dir % 3;		// direction 0, 1, 2 for axis of translation or rotation
	
	// set direction type to translation or rotation
	dirType = (dir<3) ? TRANS : ROTATE;
	
	// now switch on dimensionality of element
	switch (elemType) {
	  	    
	  case D2N6:
	    if (dirType == TRANS) {
            tran(i,3) = transformation(indx,0);
	        tran(i,4) = transformation(indx,1);
	        tran(i,5) = 0.0;
	    } else if (dirType == ROTATE) {
            tran(i,3) = 0.0;
            tran(i,4) = 0.0;
            tran(i,5) = transformation(indx,2);
	    }
	    break;
		    
	  case D3N12:
	    if (dirType == TRANS) {
            tran(i,6)  = transformation(indx,0);
	        tran(i,7)  = transformation(indx,1);
	        tran(i,8)  = transformation(indx,2);
            tran(i,9)  = 0.0;
            tran(i,10) = 0.0;
            tran(i,11) = 0.0;
	    } else if (dirType == ROTATE) {
            tran(i,6)  = 0.0;
	        tran(i,7)  = 0.0;
	        tran(i,8)  = 0.0;
            tran(i,9)  = transformation(indx,0);
            tran(i,10) = transformation(indx,1);
            tran(i,11) = transformation(indx,2);
	    }
	    break;
		 
	} // end switch
	
	// fill in first half of transformation matrix with negative sign
	
	for (int j=0; j < numDOF/2; j++ )
	    tran(i,j) = -tran(i,j+numDOF/2);
	
    } // end loop over dofs
}
		     

// Compute current strain for 1d material mat
// dispDiff are the displacements of node 2 minus those of node 1
double
InerterElement::computeCurrentStrain1d( int mat,
				    const Vector& dispDiff ) const
{
    double strain = 0.0;

    for (int i=0; i<numDOF/2; i++){
        strain += -dispDiff(i) * (*t1d)(mat,i);
    }

    return strain;
}


void
InerterElement::updateDir(const Vector& x, const Vector& y)
{
	this->setUp(connectedExternalNodes(0), connectedExternalNodes(1), x, y);
	this->setTran1d(elemType);
}
