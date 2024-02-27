#include "elasticStretchingForceShell.h"
#include "time_steppers/baseTimeStepper.h"

elasticStretchingForceShell::elasticStretchingForceShell(const shared_ptr<softRobots>& m_soft_robots) :
                                               baseForce(m_soft_robots)
{
    f.setZero(3);
    Jss.setZero(6,6);
    Id3 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;
}

elasticStretchingForceShell::~elasticStretchingForceShell()
{
    ;
}

void elasticStretchingForceShell::computeForce(double dt)
{
    int shell_limb_idx = 0;
    for (const auto& shell_limb : soft_robots->shell_limbs) {

        for (int i = 0; i < shell_limb->ne; i++)
        { 
            // if (shell_limb->isEdgeJoint[i]) continue;
            epsX = shell_limb->edge_len(i) / shell_limb->ref_len(i) - 1.0;
            f = shell_limb->EA * (shell_limb->tangent).row(i) * epsX;

            // cout<<"stretching force being added to the force vector: "<<-f;

            int n1 = shell_limb->EdgeIsBet[i][0] ; // node 1 number
            int n2 = shell_limb->EdgeIsBet[i][1] ; // node 2 number
            for (int k = 0; k < 3; k++)
            {
                ind = 3*n1 + k;
                stepper->addForceShell(ind, - f[k], shell_limb_idx); // subtracting elastic force


                ind = 3*n2 + k;
                stepper->addForceShell(ind, f[k], shell_limb_idx); // adding elastic force
            } 
        }
        shell_limb_idx++;
    }
}

void elasticStretchingForceShell::computeForceAndJacobian(double dt)
{
    computeForce(dt);

    int shell_limb_idx = 0;
    for (const auto& shell_limb : soft_robots->shell_limbs) {
        for (int i = 0; i < shell_limb->ne; i++)
        {
            int node0ind = shell_limb->EdgeIsBet[i][0] ; // node 1 number
            int node1ind = shell_limb->EdgeIsBet[i][1] ; // node 2 number

            int ind[6] = { 3 * node0ind, 3 * node0ind + 1, 3 * node0ind + 2, \
            3 * node1ind, 3 * node1ind + 1, 3 * node1ind + 2 };

            // if (shell_limb->isEdgeJoint[i]) continue;
            len = shell_limb->edge_len[i];
            refLength = shell_limb->ref_len[i];

            dxx(0) = shell_limb->x(3*node1ind+0) - shell_limb->x(3*node0ind+0);
            dxx(1) = shell_limb->x(3*node1ind+1) - shell_limb->x(3*node0ind+1);
            dxx(2) = shell_limb->x(3*node1ind+2) - shell_limb->x(3*node0ind+2);

            u = dxx;
            v = u.transpose();
            M0= shell_limb->EA * ((1/refLength - 1/len) * Id3 + (1/len) * (u*v) / (u.norm() * u.norm()));

            Jss.block(0,0,3,3) = M0;
            Jss.block(3,3,3,3) = M0;
            Jss.block(3,0,3,3) = -M0;
            Jss.block(0,3,3,3) = -M0;

            // cout<<"stretching jacobian being added to the jaconian matrix: "<<-Jss;

            // confirm if this is correctly done
            int p_counter = 0, q_counter = 0;
            for (int j : ind) {
                q_counter = 0;
                for (int k : ind) {
                    stepper->addJacobianShell(j, k, Jss(p_counter, q_counter), shell_limb_idx);
                    q_counter++;
                }
                p_counter++;
            }
        }
        shell_limb_idx++;
    }
}
