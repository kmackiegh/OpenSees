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
                                                                        
// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/CohesiveZoneQuad/CohesiveZoneQuad.cpp,v $

// Written: KRM
// Created: Aug 2019
//
// Description: This file contains the class definition for CohesiveZoneQuad.

#include <CohesiveZoneQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>

void* OPS_CohesiveZoneQuad()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 2) {
        opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with CohesiveZoneQuad element\n";
        return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element CohesiveZoneQuad eleTag? iNode? jNode? kNode? lNode? thk? matTag? intType? <-vecn n1? n2? n3?>\n";
        return 0;
    }

    // CohesiveZoneQuadId, iNode, jNode, kNode, lNode
    int idata[5];
    int num = 5;
    if (OPS_GetIntInput(&num,idata) < 0) {
        opserr << "WARNING: invalid integer inputs\n";
        return 0;
    }

    double thk = 1.0;
    num = 1;
    if (OPS_GetDoubleInput(&num,&thk) < 0) {
        opserr << "WARNING: invalid double inputs\n";
        return 0;
    }

    // import material tag and check pointer
    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
        opserr << "WARNING: invalid matTag\n";
        return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matTag;
        opserr << "\nCohesiveZoneQuad element: " << idata[0] << endln;
        return 0;
    }
    
    // get integration type and check range
    int integType;
    num = 1;
    if (OPS_GetIntInput(&num,&integType) < 0) {
        opserr << "WARNING: invalid intType\n";
        return 0;
    }
    
    if (integType < 1 || integType > 3) {
        opserr << "WARNING intType must be between 1 and 3";
        opserr << " for CohesiveZoneQuad element: " << idata[0] << endln;
        return 0;
    }
    
    // optional outward normal vector
    Vector vecn(3);
    vecn(0) = 0; vecn(1) = 1; vecn(2) = 0;
    const char* type;
    
    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        if (strcmp(type, "-vecn") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING: insufficient arguments after -vecn\n";
                return 0;
            }
            num = 3;
            if (OPS_GetDoubleInput( &num,&vecn(0) ) < 0) {
                opserr << "WARNING invalid -vecn value for ele  " << idata[0] << endln;
                return 0;
            }
        }
    }
    
    // easier to make sure it's a unit vector here
    vecn = vecn / vecn.Norm();

    // create element
    return new CohesiveZoneQuad(idata[0],idata[1],idata[2],idata[3],idata[4],
			                *mat,thk,integType,vecn);
}


Matrix CohesiveZoneQuad::K(8,8);
Vector CohesiveZoneQuad::P(8);
double CohesiveZoneQuad::shp[3][2];
double CohesiveZoneQuad::pts[2];
double CohesiveZoneQuad::wts[2];


CohesiveZoneQuad::CohesiveZoneQuad(int tag, int nd1, int nd2, int nd3, int nd4,
			   NDMaterial &m, double t, int it, const Vector _vec)
:Element (tag, ELE_TAG_CohesiveZoneQuad),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), delp(4), thickness(t), vecn(_vec), ae(4,8), Ki(0)
{
    // default node indexing
    int i;
    for (i=0; i < 4; i++)
        indx[i] = i;
    
    // integration scheme
    if (it == 1) {
        // Gauss
        pts[0] = -0.5773502691896258;
        pts[1] =  0.5773502691896258;
        wts[0] = 1.0;
        wts[1] = 1.0;
    } else if (it == 2) {
        // Lobatto or Newton-Cotes
        pts[0] = -1.0;
        pts[1] =  1.0;
        wts[0] = 1.0;
        wts[1] = 1.0;
    } else if (it == 3) {
        // Chebyshev NYI ----------------
        pts[0] = -0.5773502691896258;
        pts[1] =  0.5773502691896258;
        wts[0] = 1.0;
        wts[1] = 1.0;
    }

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[2];
    
    if (theMaterial == 0) {
        opserr << "CohesiveZoneQuad::CohesiveZoneQuad - failed allocate material model pointer\n";
        exit(-1);
    }

    for (i = 0; i < 2; i++) {
        // Get copies of the material model for each integration point
        theMaterial[i] = m.getCopy("2D");
            
        // Check allocation
        if (theMaterial[i] == 0) {
            opserr << "CohesiveZoneQuad::CohesiveZoneQuad -- failed to get a copy of material model\n";
            exit(-1);
        }
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    
    for (i=0; i<4; i++)
        theNodes[i] = 0;
}

CohesiveZoneQuad::CohesiveZoneQuad()
:Element (0,ELE_TAG_CohesiveZoneQuad),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), delp(4), thickness(0.0), vecn(0), ae(4,8), Ki(0)
{
    // integration scheme defaults to Gauss, this may be able to initialize to zero
    pts[0] = -0.5773502691896258;
    pts[1] =  0.5773502691896258;
    wts[0] = 1.0;
    wts[1] = 1.0;
  
    for (int i=0; i<4; i++) {
        theNodes[i] = 0;
        indx[i] = i;
    }
}

CohesiveZoneQuad::~CohesiveZoneQuad()
{    
  for (int i = 0; i < 2; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;

  if (Ki != 0)
    delete Ki;
}

int
CohesiveZoneQuad::getNumExternalNodes() const
{
    return 4;
}

const ID&
CohesiveZoneQuad::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
CohesiveZoneQuad::getNodePtrs(void)
{
    return theNodes;
}

int
CohesiveZoneQuad::getNumDOF()
{
    return 8;
}

void
CohesiveZoneQuad::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        theNodes[2] = 0;
        theNodes[3] = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
        //opserr << "FATAL ERROR CohesiveZoneQuad (tag: %d), node not found in domain",
        //	this->getTag());
        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
        //opserr << "FATAL ERROR CohesiveZoneQuad (tag: %d), has differing number of DOFs at its nodes",
        //	this->getTag());
        return;
    }
    this->DomainComponent::setDomain(theDomain);

    // check the outward normal vector to find the face it's perpendicular to
    Vector l12 = theNodes[1]->getCrds() - theNodes[0]->getCrds();
    Vector l23 = theNodes[2]->getCrds() - theNodes[1]->getCrds();
    Vector l34 = theNodes[3]->getCrds() - theNodes[2]->getCrds();
    Vector l41 = theNodes[0]->getCrds() - theNodes[3]->getCrds();
    
    Vector crossk(4);
    crossk(0) = vecn(0)*l12(1) - vecn(1)*l12(0);
    crossk(1) = vecn(0)*l23(1) - vecn(1)*l23(0);
    crossk(2) = vecn(0)*l34(1) - vecn(1)*l34(0);
    crossk(3) = vecn(0)*l41(1) - vecn(1)*l41(0);
    
	// vface = 1 swapped for vface = 3 to account for rotation of IP
    int vface = 0;
    double kcrit = crossk(0);
    if ( crossk(1) > kcrit) {
        vface = 3;
        kcrit = crossk(1);
    }
    if ( crossk(2) > kcrit) {
        vface = 2;
        kcrit = crossk(2);
    }
    if ( crossk(3) > kcrit) {
        vface = 1;
        kcrit = crossk(3);
    }
    
    // store ordering of displacement dofs
    Vector dsp_order(4);
    for (int kl = 0; kl < 4; kl++) {
        if (vface == 0 || vface == 2)
            dsp_order(kl) = kl;
        if (vface == 1 || vface == 3)
            dsp_order(kl) = 1-kl;
        
        if (dsp_order(kl) < 0)
            dsp_order(kl) = dsp_order(kl)+4;
    }

	// rotation matrix to get the normal and shear dofs lined up
    Matrix rot(4,4);
    for (int kl = 0; kl < 4; kl++) {
		if (vface == 0 || vface == 2)
	        rot(kl,dsp_order(kl)) = vface-1;
		
		if (vface == 1)
			rot(kl,dsp_order(kl)) = dsp_order(kl)-kl;

		if (vface == 3)
			rot(kl,dsp_order(kl)) = kl-dsp_order(kl);
	}
    
    // renumber node indices so outward normal always points in positive eta direction
    vface = vface - 2;
    for (int kl = 0; kl < 4; kl++) {
        if (vface < 0)
            vface = vface+4;
        if (vface > 3)
            vface = vface-4;
        
        indx[vface] = kl;
        vface++;
    }
    
    //opserr << "final node index order in natural system is " << indx[0] << " " << indx[1] << " " << indx[2] << " " << indx[3] << endln;
    
    // transformation matrix
    ae.Zero();
    for (int kl = 0; kl < 4; kl++) {
        ae(kl,kl) = -1;
        int rowno = kl+2;
        if (rowno > 3)
            rowno = rowno - 4;
        ae(rowno,kl+4) = 1;
    }
	
    
    

	//opserr << "dsp_order: " << dsp_order;
	//opserr << "ae: " << ae;
	//opserr << "rot: " << rot << endln;
	
    
    // and the final form for transformation after rotation is multiplied
    Matrix temp = ae;
    temp.addMatrixProduct(0.0,rot,ae,1.0);
    ae = temp;
	
	//opserr << "ae: rot*ae = " << ae;
	
}

int
CohesiveZoneQuad::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "CohesiveZoneQuad::commitState () - failed in base class";
    }
    
    //const Vector &eleVec = delp;
    Information eleData;
    //eleData.setVector(delp);
    Vector telp(3);
    double dvol;
    
    // Loop over the integration points and commit the material states
    for (int i = 0; i < 2; i++) {
        // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i]);
        dvol *= (wts[i]);

        telp(0) = delp(2*i);
        telp(1) = delp(2*i+1);
        telp(2) = dvol;
        eleData.setVector(telp);
        
        // first call sends local element forces to integration point
        retVal += theMaterial[i]->updateState(eleData);
        retVal += theMaterial[i]->commitState();
    }

    return retVal;
}

int
CohesiveZoneQuad::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 2; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
CohesiveZoneQuad::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 2; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
CohesiveZoneQuad::update()
{
    // use reordered indices from indx (old way)
	//const Vector &disp1 = theNodes[indx[0]]->getTrialDisp();
	//const Vector &disp2 = theNodes[indx[1]]->getTrialDisp();
	//const Vector &disp3 = theNodes[indx[2]]->getTrialDisp();
	//const Vector &disp4 = theNodes[indx[3]]->getTrialDisp();
    
    // combine all nodal displacements for transformation (new way)
    Vector disp(8);
    for (int beta = 0; beta < 4; beta++) {
        const Vector &dispi = theNodes[indx[beta]]->getTrialDisp();
        disp(beta*2) = dispi(0);
        disp(beta*2+1) = dispi(1);
    }
	
    // initialize relative displacements
	Vector delu(4);
    delu.Zero();

	// relative displacements, this is the transformation A * ue (old way)
	//delu(0) = disp4(0)-disp1(0);	delu(1) = disp4(1)-disp1(1);
	//delu(2) = disp3(0)-disp2(0);	delu(3) = disp3(1)-disp2(1);
    
    // actually use transformation (new way)
    delu.addMatrixVector(1.0,ae,disp,1.0);
   	
	//opserr << "   matrix disp: " << disp;
	//opserr << "   matrix ae: " << ae;
	//opserr << "   matrix delu: " << delu;

	int ret = 0;
    Vector dsp(2);

	// Loop over the integration points
	for (int i = 0; i < 2; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i]);

        // Interpolate displacement
		dsp.Zero();
		for (int beta = 0; beta < 2; beta++) {
			dsp(0) += shp[2][beta]*delu(beta*2);
			dsp(1) += shp[2][beta]*delu(beta*2+1);
		}

		// Set the material displacements
		ret += theMaterial[i]->setTrialStrain(dsp);
        //opserr << "   CZQuad IP" << i+1 << " has prescribed slips of " << dsp << endln;
	}

	return ret;
}


const Matrix&
CohesiveZoneQuad::getTangentStiff()
{

	double dvol;
    Matrix delk(4,4);
    delk.Zero();

	// Loop over the integration points
	for (int i = 0; i < 2; i++) {

        // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i]);
        dvol *= (thickness*wts[i]);

        // Get the material tangent
        const Matrix &D = theMaterial[i]->getTangent();

        // Perform numerical integration
        //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
        //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

        double D00 = D(0,0); double D01 = D(0,1);
        double D10 = D(1,0); double D11 = D(1,1);
        //opserr << "inside getTangent, IP returning stiffness " << D;

        for (int alpha = 0; alpha < 2; alpha++) {
            int ia = alpha*2;
            
            for (int beta = 0; beta < 2; beta++) {
                int ib = beta*2;
                
                delk(ia,ib) += dvol * shp[2][alpha] * D00 * shp[2][beta];
                delk(ia,ib+1) += dvol * shp[2][alpha] * D01 * shp[2][beta];
                delk(ia+1,ib) += dvol * shp[2][alpha] * D10 * shp[2][beta];
                delk(ia+1,ib+1) += dvol * shp[2][alpha] * D11 * shp[2][beta];
            }
        }
	}
    //opserr << delk;
    
    // need A^T * k * A and then indx mapping to nodal dofs
    K.Zero();
    Matrix Ktemp = K;
    Ktemp.addMatrixTripleProduct(1.0,ae,delk,1.0);
    
    for (int alpha = 0; alpha < 4; alpha++) {
        int ia = indx[alpha]*2;
        
        for (int beta = 0; beta < 4; beta++) {
            int ib = indx[beta]*2;
            
            K(ia,ib) = Ktemp(alpha*2,beta*2);
            K(ia,ib+1) = Ktemp(alpha*2,beta*2+1);
            K(ia+1,ib) = Ktemp(alpha*2+1,beta*2);
            K(ia+1,ib+1) = Ktemp(alpha*2+1,beta*2+1);
        }
    }
    
	//opserr << "tangent stiffness = " << K;
	
	return K;
}


const Matrix&
CohesiveZoneQuad::getInitialStiff()
{
    if (Ki != 0)
        return *Ki;

    double dvol;
    Matrix delk(4,4);
    delk.Zero();

    // Loop over the integration points
    for (int i = 0; i < 2; i++) {

        // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i]);
        dvol *= (thickness*wts[i]);

        // Get the material tangent
        const Matrix &D = theMaterial[i]->getInitialTangent();

        double D00 = D(0,0); double D01 = D(0,1);
        double D10 = D(1,0); double D11 = D(1,1);

        for (int alpha = 0; alpha < 2; alpha++) {
            int ia = alpha*2;
            
            for (int beta = 0; beta < 2; beta++) {
                int ib = beta*2;
                
                delk(ia,ib) += dvol * shp[2][alpha] * D00 * shp[2][beta];
                delk(ia,ib+1) += dvol * shp[2][alpha] * D01 * shp[2][beta];
                delk(ia+1,ib) += dvol * shp[2][alpha] * D10 * shp[2][beta];
                delk(ia+1,ib+1) += dvol * shp[2][alpha] * D11 * shp[2][beta];
            }
        }
    }
    
    // need A^T * k * A and then indx mapping to nodal dofs
    K.Zero();
    Matrix Ktemp = K;
    Ktemp.addMatrixTripleProduct(1.0,ae,delk,1.0);
    
    for (int alpha = 0; alpha < 4; alpha++) {
        int ia = indx[alpha]*2;
        
        for (int beta = 0; beta < 4; beta++) {
            int ib = indx[beta]*2;
            
            K(ia,ib) = Ktemp(alpha*2,beta*2);
            K(ia,ib+1) = Ktemp(alpha*2,beta*2+1);
            K(ia+1,ib) = Ktemp(alpha*2+1,beta*2);
            K(ia+1,ib+1) = Ktemp(alpha*2+1,beta*2+1);
        }
    }

    Ki = new Matrix(K);
    return K;
}

const Matrix&
CohesiveZoneQuad::getMass()
{
	K.Zero();
	return K;
}

void
CohesiveZoneQuad::zeroLoad(void)
{
	Q.Zero();
  	return;
}

int 
CohesiveZoneQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "CohesiveZoneQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int 
CohesiveZoneQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
    // does nothing as element has no mass yet!
    return 0;
}

const Vector&
CohesiveZoneQuad::getResistingForce()
{
	double dvol;
    delp.Zero();

	// Loop over the integration points
	for (int i = 0; i < 2; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i]);
		dvol *= (thickness*wts[i]);

		// Get material stress response
		const Vector &sigma = theMaterial[i]->getStress();
        //opserr << "   CZQuad IP" << i+1 << " returned stresses of " << sigma;

		// Perform numerical integration on internal force
		//P = P + (N^ sigma) * intWt * detJ;
        for (int alpha = 0; alpha < 2; alpha++) {
            int ia = alpha*2;
            delp(ia) += dvol*shp[2][alpha]*sigma(0);
			delp(ia+1) += dvol*shp[2][alpha]*sigma(1);
		}
	}
	
    // need A^T * q and then indx mapping to nodal dofs
    P.Zero();
    Vector Ptemp = P;
    Ptemp.addMatrixTransposeVector(1.0,ae,delp,1.0);

    for (int alpha = 0; alpha < 4; alpha++) {
        int ia = indx[alpha]*2;
        P(ia) = Ptemp(alpha*2);
        P(ia+1) = Ptemp(alpha*2+1);
    }
    
	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);
	//opserr << "Resisting force = " << P;

	return P;
}

const Vector&
CohesiveZoneQuad::getResistingForceIncInertia()
{
	// Compute the current resisting force
	this->getResistingForce();

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P += this->getRayleighDampingForces();

	return P;
}

int
CohesiveZoneQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(2);
  data(0) = this->getTag();
  data(1) = thickness;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING CohesiveZoneQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  
  // Now quad sends the ids of its materials
  int matDbTag;
  static ID idData(8);
  
  int i;
  for (i = 0; i < 2; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+2) = matDbTag;
  }
  
  idData(4) = connectedExternalNodes(0);
  idData(5) = connectedExternalNodes(1);
  idData(6) = connectedExternalNodes(2);
  idData(7) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING CohesiveZoneQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 2; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING CohesiveZoneQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
CohesiveZoneQuad::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(2);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING CohesiveZoneQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);

  static ID idData(8);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING CohesiveZoneQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(4);
  connectedExternalNodes(1) = idData(5);
  connectedExternalNodes(2) = idData(6);
  connectedExternalNodes(3) = idData(7);
  
  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[2];
    if (theMaterial == 0) {
      opserr << "CohesiveZoneQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 2; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+2);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	    opserr << "CohesiveZoneQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	    return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "CohesiveZoneQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 2; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+2);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	    delete theMaterial[i];
	    theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	    if (theMaterial[i] == 0) {
          opserr << "CohesiveZoneQuad::recvSelf() - material " << i << "failed to create\n";
	      return -1;
	    }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "CohesiveZoneQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }
  
  return res;
}

void
CohesiveZoneQuad::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    s << "#CohesiveZoneQuad\n";
    
    int i;
    const int numNodes = 4;
    const int nstress = 1;
    
    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = 2;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;
  }
  
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nCohesiveZoneQuad, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (Tt Tn)" << endln;
	for (int i = 0; i < 2; i++)
		s << "\t\tIntegration point " << i+1 << ": " << theMaterial[i]->getStress();
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"CohesiveZoneQuad\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
      s << connectedExternalNodes(1) << ", ";
      s << connectedExternalNodes(2) << ", ";
      s << connectedExternalNodes(3) << "], ";
      s << "\"thickness\": " << thickness << ", ";
      s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int
CohesiveZoneQuad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

    // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 2 we will plot material stresses otherwise 0.0
    static Vector values(2);

    for (int j=0; j<2; j++)
	   values(j) = 0.0;

    if (displayMode < 3 && displayMode > 0) {
        for (int i=0; i<2; i++) {
            const Vector &stress = theMaterial[i]->getStress();
            values(i) = stress(displayMode-1);
        }
    }

    // now  determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end3Crd = theNodes[2]->getCrds();	
    const Vector &end4Crd = theNodes[3]->getCrds();	

    static Matrix coords(4,3);

    if (displayMode >= 0) {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        const Vector &end3Disp = theNodes[2]->getDisp();
        const Vector &end4Disp = theNodes[3]->getDisp();

        for (int i = 0; i < 2; i++) {
            coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
            coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
            coords(2,i) = end3Crd(i) + end3Disp(i)*fact;
            coords(3,i) = end4Crd(i) + end4Disp(i)*fact;
        }
        
    } else {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        const Matrix &eigen3 = theNodes[2]->getEigenvectors();
        const Matrix &eigen4 = theNodes[3]->getEigenvectors();
        
        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 2; i++) {
                coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
                coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
                coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
            }
        } else {
            for (int i = 0; i < 2; i++) {
                coords(0,i) = end1Crd(i);
                coords(1,i) = end2Crd(i);
                coords(2,i) = end3Crd(i);
                coords(3,i) = end4Crd(i);
            }
        }
    }

    int error = 0;

    // finally we  the element using drawPolygon
    error += theViewer.drawPolygon (coords, values, this->getTag());

    return error;
}

Response*
CohesiveZoneQuad::setResponse(const char **argv, int argc,
			  OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","CohesiveZoneQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, P);
  }   

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 2) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("xi",pts[pointNum-1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    }
  }
    
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<2; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("xi",pts[i]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","Tt");
      output.tag("ResponseType","Tn");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(4));
  }

  else if ((strcmp(argv[0],"strain") ==0) || (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<2; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("xi",pts[i]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","Dt");
      output.tag("ResponseType","Dn");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(4));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int 
CohesiveZoneQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {
    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {
    // Loop over the integration points
    static Vector stresses(4);
    int cnt = 0;
    for (int i = 0; i < 2; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      cnt += 2;
    }
    
    return eleInfo.setVector(stresses);
      
  } else if (responseID == 4) {
    // Loop over the integration points
    static Vector stresses(4);
    int cnt = 0;
    for (int i = 0; i < 2; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      cnt += 2;
    }

    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}

int
CohesiveZoneQuad::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    int res = -1;

    // quad pressure loading
    if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {
        if (argc < 3)
            return -1;

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= 2)
            return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
        else
            return -1;
    }

    // otherwise it could be just a forall material parameter
    else {
        int matRes = res;
        for (int i=0; i<2; i++) {
            matRes =  theMaterial[i]->setParameter(argv, argc, param);

            if (matRes != -1)
                res = matRes;
        }
    }

    return res;
}
    
int
CohesiveZoneQuad::updateParameter(int parameterID, Information &info)
{
    int res = -1;
    int matRes = res;
    switch (parameterID) {
        case -1:
            return -1;

        case 1:
            for (int i = 0; i<2; i++) {
                matRes = theMaterial[i]->updateParameter(parameterID, info);
            }
            if (matRes != -1) {
                res = matRes;
            }
            return res;

        default:
            return -1;
    }
}

double CohesiveZoneQuad::shapeFunction(double xi)
{
	const Vector &nd1Crds = theNodes[indx[0]]->getCrds();
	const Vector &nd2Crds = theNodes[indx[1]]->getCrds();
	const Vector &nd3Crds = theNodes[indx[2]]->getCrds();
	const Vector &nd4Crds = theNodes[indx[3]]->getCrds();

    // retaining only xi for interface
    double eta = 0;
	double oneMinuseta = 1.0-eta;
    double onePluseta = 1.0+eta;
    double oneMinusxi = 1.0-xi;
    double onePlusxi = 1.0+xi;

	shp[2][0] = 0.5*oneMinusxi;	    // N_1
	shp[2][1] = 0.5*onePlusxi;		// N_2

    // but to get the Jacobian we retain the original shape functions to get element dimensions
    // note that by sampling with eta = 0, we cannot get correct Jacobian for distorted element
    double J[2][2];
    
    J[0][0] = 0.25 * (-nd1Crds(0)*oneMinuseta + nd2Crds(0)*oneMinuseta +
                nd3Crds(0)*(onePluseta) - nd4Crds(0)*(onePluseta));

    J[0][1] = 0.25 * (-nd1Crds(0)*oneMinusxi - nd2Crds(0)*onePlusxi +
                nd3Crds(0)*onePlusxi + nd4Crds(0)*oneMinusxi);

    J[1][0] = 0.25 * (-nd1Crds(1)*oneMinuseta + nd2Crds(1)*oneMinuseta +
                nd3Crds(1)*onePluseta - nd4Crds(1)*onePluseta);

    J[1][1] = 0.25 * (-nd1Crds(1)*oneMinusxi - nd2Crds(1)*onePlusxi +
                nd3Crds(1)*onePlusxi + nd4Crds(1)*oneMinusxi);

    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    //opserr << "detJ at xi = " << xi << " is " << detJ << endln;
    
    // to make loops easier, multiply by 2/h0 here
	//double detJ = J[0]/2;
    
	// shape function derivatives not yet completed properly
    shp[0][0] = 0;	        // N_1,1
    shp[0][1] = 0;		    // N_2,1
	
    shp[1][0] = 0;	        // N_1,2
    shp[1][1] = 0;		    // N_2,2

    return detJ;
}

