#pragma once

#include <robot/GeneralizedCoordinatesRobotRepresentation.h>
#include <robot/Robot.h>

namespace crl {

struct IK_EndEffectorTargets {
    RB *rb = nullptr;
    P3D p;       // local coordinates of end effector in rb's frame
    P3D target;  // target position in world frame
};

class IK_Solver {
public:
    IK_Solver(Robot *robot) : robot(robot) {}

    ~IK_Solver(void) {}

    /**
     * add IK end effector target to solver. Specify the end effector point p, which 
     * is specified in the local coordinates of rb and its target expressed in world frame.
     */
    void addEndEffectorTarget(RB *rb, P3D p, P3D target) {
        endEffectorTargets.push_back(IK_EndEffectorTargets());
        endEffectorTargets.back().rb = rb;
        endEffectorTargets.back().p = p;
        endEffectorTargets.back().target = target;
    }

    void solve(int nSteps = 10) {
        GeneralizedCoordinatesRobotRepresentation gcrr(robot);
        
        for (uint i = 0; i < nSteps; i++) {
            dVector q;
            gcrr.getQ(q);

            // get current generalized coordinates of the robots

            // TODO: Ex.2-2 Inverse Kinematics
            //
            // update generalized coordinates of the robot by solving IK.

            // remember, we don't update base pose since we assume it's already at
            // the target position and orientation
            
            // TODO: here, compute deltaq using GD, Newton, or Gauss-Newton.
            // end effector targets are stored in endEffectorTargets vector.
            //
            // Hint:
            // - use gcrr.estimate_linear_jacobian(p, rb, dpdq) function for Jacobian matrix.
            // - if you already implemented analytic Jacobian, you can use gcrr.compute_dpdq(const P3D &p, RB *rb, Matrix &dpdq)
            // - don't forget we use only last q.size() - 6 columns (use block(0,6,3,q.size() - 6) function)
            // - when you compute inverse of the matrix, use ldlt().solve() instead of inverse() function. this is numerically more stable.
            //   see https://eigen.tuxfamily.org/dox-devel/group__LeastSquares.html

            // TODO: your implementation should be here.
            dVector deltaq(q.size() - 6); // fixed base problem
            deltaq.setZero(); // for each iteration, set deltaq = 0
            //My code start--------------------------------------------
            double tau = 1; //step rate; to keep linearize valid
            Matrix J(3,q.size());
            for(int it = 0; it < 4; it++)
            {
                P3D p = endEffectorTargets[it].p;
                RB *rb = endEffectorTargets[it].rb;
                gcrr.estimate_linear_jacobian(p, rb, J); // Compute Jacobian numerically
                //gcrr.compute_dpdq(p, rb, J); // Compute Jacobian analytically
                Vector3d v = V3D(gcrr.getWorldCoordinates(p,rb) - endEffectorTargets[it].target);
                Matrix Jacobian = J.block(0,6,3,q.size() - 6);
                Matrix JT = Jacobian.transpose();
                deltaq = deltaq -tau*(JT*Jacobian).ldlt().solve(JT*v);
                
                //-----------Analytical Jacobian Checker-----------------//
                Matrix AJ(3,q.size());
                gcrr.compute_dpdq(p, rb, AJ);
                float a = (AJ-J).norm();
                printf("Delta norm is %f\n",a);
                //-----------Analytical Jacobian Checker-----------------//
                
            }
            //My code end----------------------------------------------

            q.tail(q.size() - 6) += deltaq;

            // now update gcrr with q
            gcrr.setQ(q);
        }

        gcrr.syncRobotStateWithGeneralizedCoordinates();

        // clear end effector targets
        // we will add targets in the next step again.
        endEffectorTargets.clear();
    }

private:
    Robot *robot;
    std::vector<IK_EndEffectorTargets> endEffectorTargets;
};

}  // namespace crl