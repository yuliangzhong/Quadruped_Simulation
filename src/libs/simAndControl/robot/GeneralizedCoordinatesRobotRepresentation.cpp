#include <robot/GeneralizedCoordinatesRobotRepresentation.h>
#include <robot/RBJoint.h>
#include <utils/utils.h>


namespace crl {

GeneralizedCoordinatesRobotRepresentation::
    GeneralizedCoordinatesRobotRepresentation(Robot *a) {
    robot = a;
    resize(q, getDOFCount());
    syncGeneralizedCoordinatesWithRobotState();
}

//returns the total number of degrees of freedom in the reduced representation
int GeneralizedCoordinatesRobotRepresentation::getDOFCount() const {
    //we will have 3 for root position, 3 for root orientation, and then 1 for each joint,
    //assuming they're all hinge joints, which is what this data structure is designed for
    return 6 + (int)robot->jointList.size();
}

//returns the index of the generalized coordinate corresponding to this joint
int GeneralizedCoordinatesRobotRepresentation::getQIdxForJoint(
    RBJoint *joint) const {
    //the root of the robot has a nullptr joint, and as such this should correspond to the first qIdx of the root DOFs
    if (joint == NULL) return 5;
    return 6 + joint->jIndex;
}

//returns the index of the generalized coordinate corresponding to the index of this joint
int GeneralizedCoordinatesRobotRepresentation::getQIdxForJointIdx(
    int jIdx) const {
    return 6 + jIdx;
}

//returns a pointer to the joint corresponding to this generalized coordinate index.
//If the index corresponds to a root DOF, this method will return NULL.
RBJoint *GeneralizedCoordinatesRobotRepresentation::getJointForQIdx(
    int qIdx) const {
    if (qIdx < 6) return NULL;
    return robot->jointList[qIdx - 6];
}

/**
* In the tree-like hierarchy of joints/dofs, this method returns the parent index of the
* dof corresponding to qIdx
*/
int GeneralizedCoordinatesRobotRepresentation::getParentQIdxOf(int qIdx) {
    if (qIdx < 6) return qIdx - 1;
    return getQIdxForJoint(getJointForQIdx(qIdx)->parent->pJoint);
}

// updates q and qDot given current state of robot
void GeneralizedCoordinatesRobotRepresentation::
    syncGeneralizedCoordinatesWithRobotState() {
    // write out the position of the root...
    RobotState state(robot);
    P3D pos = state.getPosition();
    Quaternion orientation = state.getOrientation();

    q[0] = pos.x;
    q[1] = pos.y;
    q[2] = pos.z;

    // Root
    computeEulerAnglesFromQuaternion(orientation, getQAxis(5), getQAxis(4),
                                     getQAxis(3), q[5], q[4], q[3]);

    // Now go through each joint, and decompose it as appropriate...
    for (uint i = 0; i < robot->jointList.size(); i++) {
        int qIdx = getQIdxForJointIdx(i);
        // Only 1-dof hinge joints
        computeRotationAngleFromQuaternion(state.getJointRelativeOrientation(i),
                                           getQAxis(qIdx), q[qIdx]);
    }
}

// returns the axis corresponding to the indexed generalized coordinate,
// expressed in local coordinates
V3D GeneralizedCoordinatesRobotRepresentation::getQAxis(int qIndex) const {
    if (qIndex >= 0 || qIndex < 6) {
        // the first three are the translational dofs of the body
        if (qIndex == 0) return V3D(1, 0, 0);
        if (qIndex == 1) return V3D(0, 1, 0);
        if (qIndex == 2) return V3D(0, 0, 1);
        if (qIndex == 3) return RBGlobals::worldUp;  // y - yaw
        if (qIndex == 4)
            return RBGlobals::worldUp.cross(robot->forward);  // x - pitch
        if (qIndex == 5) return robot->forward;               // z - roll
    }

    return getJointForQIdx(qIndex)->rotationAxis;
}

void GeneralizedCoordinatesRobotRepresentation::
    syncRobotStateWithGeneralizedCoordinates() {
    RobotState rs(robot);
    getReducedRobotState(rs);
    robot->setState(&rs);
}

// given the current state of the generalized representation, output the reduced
// state of the robot
void GeneralizedCoordinatesRobotRepresentation::getReducedRobotState(
    RobotState &state) {
    // set the position, velocity, rotation and angular velocity for the root
    state.setPosition(P3D(0, 0, 0) + getQAxis(0) * q[0] + getQAxis(1) * q[1] +
                      getQAxis(2) * q[2]);
    state.setOrientation(getOrientationFor(robot->root));

    for (uint i = 0; i < robot->jointList.size(); i++) {
        int qIdx = getQIdxForJointIdx(i);
        Quaternion jointOrientation =
            getRotationQuaternion(q[qIdx], getQAxis(qIdx));

        state.setJointRelativeOrientation(jointOrientation, i);
    }
    // and done...
}

// sets the current q values
void GeneralizedCoordinatesRobotRepresentation::setQ(const dVector &qNew) {
    assert(q.size() == qNew.size());
    // NOTE: we don't update the angular velocities. The assumption is that the
    // correct behavior is that the joint relative angular velocities don't
    // change, although the world relative values of the rotations do
    q = qNew;
}

// gets the current q values
void GeneralizedCoordinatesRobotRepresentation::getQ(dVector &q_copy) {
    q_copy = q;
}

void GeneralizedCoordinatesRobotRepresentation::getQFromReducedState(
    const RobotState &rs, dVector &q_copy) {
    dVector q_old = q;

    RobotState oldState(robot);
    robot->setState((RobotState *)&rs);
    syncGeneralizedCoordinatesWithRobotState();
    getQ(q_copy);
    robot->setState(&oldState);

    setQ(q_old);
}

/**
    pLocal is expressed in the coordinate frame of the link that pivots about DOF qIdx.
    This method returns the point in the coordinate frame of the parent of qIdx after
    the DOF rotation has been applied.
*/
P3D GeneralizedCoordinatesRobotRepresentation::
    getCoordsInParentQIdxFrameAfterRotation(int qIndex, const P3D &pLocal) {
    // if qIndex <= 2, this q is a component of position of the base. 
    if (qIndex <= 2) return pLocal;
    
    // My code start------------------------------------------
    RBJoint *Joint = getJointForQIdx(qIndex);
    V3D rotate = rotateVec(V3D(pLocal)-V3D(Joint->cJPos), q[qIndex], getQAxis(qIndex));
    return Joint->pJPos + getP3D(rotate);
    // My code end------------------------------------------

    // TODO: Ex.1 Forward Kinematics
    // this is a subfunction for getWorldCoordinates() and compute_dpdq()
    // return the point in the coordinate frame of the parent of qIdx after
    // the DOF rotation has been applied.
    // 
    // Hint:
    // - use rotateVec(const V3D &v, double alpha, const V3D &axis) to get a vector 
    // rotated around axis by angle alpha.
}

// returns the world coordinates for point p, which is specified in the local
// coordinates of rb (relative to its COM): p(q)
P3D GeneralizedCoordinatesRobotRepresentation::getWorldCoordinates(const P3D &p,
                                                                   RB *rb) {
    // TODO: Ex.1 Forward Kinematics
    // implement subfunction getCoordsInParentQIdxFrameAfterRotation() first.
    //
    // Hint: you may want to use the following functions
    // - getQIdxForJoint()
    // - getParentQIdxOf()
    // - getCoordsInParentQIdxFrameAfterRotation() 

    
    //My code start------------------------------------------
    if(rb->pJoint == NULL) // if rb is the base, then no such pJoint!
    {
        V3D trueP = rotateVec(V3D(p), q[3],V3D(0, 1, 0));
        V3D pp = trueP + V3D(q[0],q[1],q[2]);
        return getP3D(pp);
    }
    RBJoint *J = rb->pJoint; // if rb is not the base, then pJoint exists!
    int Jindex = getQIdxForJoint(J);
    P3D temp_p = p;
    while(true)
    {
        if(getParentQIdxOf(Jindex)<=5) break;
        temp_p = getCoordsInParentQIdxFrameAfterRotation(Jindex, temp_p);
        Jindex = getParentQIdxOf(Jindex);
    }
    temp_p = getCoordsInParentQIdxFrameAfterRotation(Jindex,temp_p);
    V3D trueP = rotateVec(V3D(temp_p), q[3],V3D(0, 1, 0));
    V3D pp = trueP + V3D(q[0],q[1],q[2]);
    return getP3D(pp);
    //My code end------------------------------------------
    
}   


// returns the global orientation associated with a specific dof q...
Quaternion GeneralizedCoordinatesRobotRepresentation::getWorldRotationForQ(
    int qIndex) {
    Quaternion qRes = Quaternion::Identity();
    // 2 here is the index of the first translational DOF of the root -- these
    // dofs do not contribute to the orientation of the rigid bodies...
    while (qIndex > 2) {
        qRes = getRelOrientationForQ(qIndex) * qRes;
        qIndex = getParentQIdxOf(qIndex);
    }
    return qRes;
}

Quaternion GeneralizedCoordinatesRobotRepresentation::getRelOrientationForQ(
    int qIndex) {
    if (qIndex < 3) return Quaternion::Identity();
    return getRotationQuaternion(q[qIndex], getQAxis(qIndex));
}

// this is a somewhat slow function to use if we must iterate through multiple
// rigid bodies...
V3D GeneralizedCoordinatesRobotRepresentation::getWorldCoordsAxisForQ(
    int qIndex) {
    if (qIndex < 3) return getQAxis(qIndex);
    return getWorldRotationForQ(qIndex) * getQAxis(qIndex);
}

// returns the world-relative orientation for rb
Quaternion GeneralizedCoordinatesRobotRepresentation::getOrientationFor(
    RB *rb) {
    int qIndex = getQIdxForJoint(rb->pJoint);
    return getWorldRotationForQ(qIndex);
}

// computes the jacobian dp/dq that tells you how the world coordinates of p
// change with q. p is expressed in the local coordinates of rb
void GeneralizedCoordinatesRobotRepresentation::compute_dpdq(const P3D &p, RB *rb, Matrix &dpdq) 
{
    resize(dpdq, 3, (int)q.size());

    Vector3d I_r_BP = V3D(getWorldCoordinates(p,rb)-V3D(q[0],q[1],q[2]));
    Matrix C_IB(3,3);
    C_IB(0,0) = cos(q[3]); C_IB(0,1) = 0; C_IB(0,2) = sin(q[3]);
    C_IB(1,0) = 0; C_IB(1,1) = 1; C_IB(1,2) = 0;
    C_IB(2,0) = -sin(q[3]); C_IB(2,1) = 0; C_IB(2,2) = cos(q[3]);
    Matrix C_BI(3,3);
    C_BI(0,0) = cos(q[3]); C_BI(0,1) = 0; C_BI(0,2) = -sin(q[3]);
    C_BI(1,0) = 0; C_BI(1,1) = 1; C_BI(1,2) = 0;
    C_BI(2,0) = sin(q[3]); C_BI(2,1) = 0; C_BI(2,2) = cos(q[3]);
    Vector3d p_base = C_BI * I_r_BP;

    V3D v = V3D(-sin(q[3])*p_base[0]+cos(q[3])*p_base[2],0,-cos(q[3])*p_base[0]-sin(q[3])*p_base[2]);

    dpdq.block(0,0,3,3).setIdentity();
    dpdq.block(0,3,3,1) = v;
    dpdq.block(0,4,3,2).setZero();

//------------------Jeo in Robot Dynamics---------------------------
    RBJoint* J = rb->pJoint;
    int JIndex = getQIdxForJoint(J);
    P3D pp = p;
    P3D I_p_IE = getWorldCoordinates(pp,rb);

    while(J!=NULL){
        V3D I_n_i = getWorldCoordsAxisForQ(JIndex);
        P3D I_p_Ii = getWorldCoordinates(J->cJPos,J->child);
        dpdq.block(0,JIndex,3,1) = I_n_i.cross(V3D(I_p_IE-I_p_Ii));
        J = J->parent->pJoint;
        JIndex = getQIdxForJoint(J);
    }
}

// estimates the linear jacobian dp/dq using finite differences
void GeneralizedCoordinatesRobotRepresentation::estimate_linear_jacobian(
    const P3D &p, RB *rb, Matrix &dpdq) 
{
    resize(dpdq, 3, (int)q.size());

    for (int i = 0; i < q.size(); i++)
    {
        double val = q[i];
        double h = 0.0001;

        // DONE: Ex. 2-1 Inverse Kinematics - Jacobian by Finite Difference
        // compute Jacobian matrix dpdq_i by FD and fill dpdq
        q[i] = val + h;
        P3D p_p = getWorldCoordinates(p,rb);  // DONE: fixed this: p(qi + h);

        q[i] = val - h;
        P3D p_m = getWorldCoordinates(p,rb);  // DONE: fixed this: p(qi - h)

        V3D dpdq_i = V3D(p_p-p_m);  // DONE: fixed this: compute derivative dp(q)/dqi

        // set Jacobian matrix components
        dpdq(0, i) = dpdq_i[0]/2/h;
        dpdq(1, i) = dpdq_i[1]/2/h;
        dpdq(2, i) = dpdq_i[2]/2/h;

        // finally, we don't want to change q[i] value. back to original value.
        q[i] = val;
    }
}

}  // namespace crl
