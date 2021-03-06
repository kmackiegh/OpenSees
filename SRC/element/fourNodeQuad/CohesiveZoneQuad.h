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
                                                                        
// $Revision: 1.15 $
// $Date: 2009-08-07 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/CohesiveZoneQuad/CohesiveZoneQuad.h,v $
                                                                        
// Written: KRM
// Created: Aug 2019
//
// Description: This file contains the class definition for CohesiveZoneQuad.

#ifndef CohesiveZoneQuad_h
#define CohesiveZoneQuad_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class CohesiveZoneQuad : public Element
{
  public:
    CohesiveZoneQuad(int tag, int nd1, int nd2, int nd3, int nd4,
		 NDMaterial &m, double t, int it, const Vector vin = 0);
    CohesiveZoneQuad();
    ~CohesiveZoneQuad();

    const char *getClassType(void) const {return "CohesiveZoneQuad";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getMass(void);    

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);

    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:
    
  private:
    NDMaterial **theMaterial; // pointer to the ND material objects
    ID connectedExternalNodes; // Tags of quad nodes
    Node *theNodes[4];

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads
    Vector delp;            // local forces after element integration
    Vector slp;         	// local slips after element integration

    double thickness;	        // Element thickness
    Vector vecn;                // outward normal vector
    Matrix ae;                  // transformation matrix
    int indx[4];                // node numbering index
    static double shp[3][2];	// Stores shape functions and derivatives (overwritten)
    static double pts[2];   	// Stores quadrature points
    static double wts[2];		// Stores quadrature weights

    double shapeFunction(double xi);

    Matrix *Ki;
};

#endif

