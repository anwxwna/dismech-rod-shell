#include "gravityForce.h"
#include "time_steppers/baseTimeStepper.h"

gravityForce::gravityForce(const shared_ptr<softRobots>& soft_robots, Vector3d g_vector) :
                                           baseForce(soft_robots), g_vector(g_vector)
{
    setGravity();
}

gravityForce::~gravityForce() = default;


void gravityForce::computeForce(double dt)
{
    int limb_idx = 0;
    for (const auto& limb : soft_robots->limbs) {
        mass_gravity = mass_gravities[limb_idx];
        for (int i = 0; i < limb->ndof; i++)
        {
            if (limb->isDOFJoint[i]) continue;
            stepper->addForce(i, -mass_gravity[i], limb_idx); // subtracting gravity force
        }
        limb_idx++;
    }

    // TODO: store these values like above
    double force;
    for (const auto & joint : soft_robots->joints) {
        for (int i = 0; i < 3; i++) {
            force = g_vector[i] * joint->mass;
            stepper->addForce(4*joint->joint_node+i, -force, joint->joint_limb);
        }
    }

    // shell
    int shell_limb_idx = 0;
    for (const auto& shell_limb : soft_robots->shell_limbs) {
        mass_gravity = mass_gravities[shell_limb_idx];
        // cout<<"gravity force being added to the force vector: "<<mass_gravity<<endl;
        for (int i = 0; i < shell_limb->ndof; i++)
        {
            // if (limb->isDOFJoint[i]) continue;
            stepper->addForceShell(i, -mass_gravity[i], shell_limb_idx); // subtracting gravity force
        }
        shell_limb_idx++;
    }
}

void gravityForce::computeForceAndJacobian(double dt)
{
    computeForce(dt);
}

void gravityForce::setGravity()
{
    for (const auto& limb : soft_robots->limbs) {
        mass_gravity = VectorXd::Zero(limb->ndof);
        for (int i = 0; i < limb->nv; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                int ind = 4*i + k;
                mass_gravity[ind] = g_vector[k] * limb->mass_array[ind];
            }
        }
        mass_gravities.push_back(mass_gravity);
    }

    // shell
    for (const auto& shell_limb : soft_robots->shell_limbs) {
        mass_gravity = VectorXd::Zero(shell_limb->ndof);
        for (int i = 0; i < shell_limb->nv; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                int ind = 3*i + k;
                mass_gravity[ind] = g_vector[k] * shell_limb->mass_array[ind];
            }
        }
        mass_gravities.push_back(mass_gravity);
        cout<<"total gravity force = "<<mass_gravity.sum()<<endl;
    }
    
}
