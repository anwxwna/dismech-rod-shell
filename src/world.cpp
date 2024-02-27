#include "world.h"


world::world(const shared_ptr<softRobots>& soft_robots,
             const shared_ptr<forceContainer>& forces,
             const simParams& sim_params) :
             soft_robots(soft_robots), forces(forces),
             time_step(0), curr_time(0.0), total_time(sim_params.sim_time) {

    
    // Declare inner elastic forces. 
    // add shell elastic forces
    if (sim_params.structure_shell){
        forces->addForce(make_shared<elasticStretchingForceShell>(soft_robots));
        forces->addForce(make_shared<elasticBendingForceShell>(soft_robots));
    }
    else { // add rod elastic forces
        forces->addForce(make_shared<elasticStretchingForce>(soft_robots));
        forces->addForce(make_shared<elasticBendingForce>(soft_robots));
        forces->addForce(make_shared<elasticTwistingForce>(soft_robots));
    }
    cout<< "no. of forces after adding elastic forces: " << forces->forces.size()<<endl;

    // Declare inertial force. Should be avoided for explicit methods
    if (sim_params.nis != FORWARD_EULER && sim_params.nis != VERLET_POSITION) {
        forces->addForce(make_shared<inertialForce>(soft_robots));
    }

    cout<< "no. of forces after adding inertial force: " << forces->forces.size()<<endl;

    // Set up the time stepper
    switch(sim_params.nis) {
        case FORWARD_EULER:
            stepper = make_shared<forwardEuler>(soft_robots, forces, sim_params);
            break;
        case VERLET_POSITION:
            stepper = make_shared<verletPosition>(soft_robots, forces, sim_params);
            break;
        case BACKWARD_EULER:
            cout<<"selected backwardEuler"<<endl;
            stepper = make_shared<backwardEuler>(soft_robots, forces, sim_params, PARDISO_SOLVER);
            cout<<"constructed backwardEuler"<<endl;
            break;
        case IMPLICIT_MIDPOINT:
            stepper = make_shared<implicitMidpoint>(soft_robots, forces, sim_params, PARDISO_SOLVER);
            break;
    }
    cout<<"after stepper declaration"<<endl;

    stepper->initStepper();

    cout<<"after the initStepper() step"<<endl;

    
    // for (const auto& limb : soft_robots->limbs) cout<<"rod-limb no. of vertices: "<<limb->nv; // this kind of thing doesn't give errors

    // cout<< soft_robots->limbs[0]->nv; //issue

    // cout<<"size of offsets for rod limbs set for the time-stepper linked with force0:"<< forces->forces[0]->stepper->offsets.size()<<endl; gives segmentation fault

    if (sim_params.enable_2d_sim) {
        for (const auto& limb : soft_robots->limbs) limb->enable2DSim();
    }

    // Update boundary conditions
    updateCons();

    cout<<"updateCons() done"<<endl;

    // Allocate every thing to prepare for the first iteration
    stepper->updateSystemForNextTimeStep();

    cout<<"stepper ->updateSystemForNextTimeStep() done"<<endl;
}


world::~world() = default;


void world::updateCons()
{
    for (const auto &limb : soft_robots->limbs)
        limb->updateMap();
    for (const auto &shell_limb : soft_robots->shell_limbs)
        shell_limb->updateMap();

    cout<<"we are here"<<endl;
    stepper->update();
}


void world::updateTimeStep() {
    curr_time += stepper->stepForwardInTime();
    time_step++;
}


void world::printSimData()
{
    auto cf = forces->cf;
    auto ff = forces->ff;
    if (cf && ff) {
        if (cf->getNumCollisions() > 0) {
            printf("time: %.4f | iters: %i | con: %i | min_dist: %.6f | floor_con: %i | f_min_dist: %.6f\n",
                   curr_time, stepper->iter, cf->getNumCollisions(), cf->getMinDist(), ff->num_contacts, ff->min_dist);
        }
        else {
            printf("time: %.4f | iters: %i | con: %i | min_dist: %s | floor_con: %i | f_min_dist: %.6f\n",
                   curr_time, stepper->iter, 0, "N/A", ff->num_contacts, ff->min_dist);
        }
    }
    else if (cf) {
        if (cf->getNumCollisions() > 0) {
            printf("time: %.4f | iters: %i | con: %i | min_dist: %.6f\n",
                   curr_time, stepper->iter, cf->getNumCollisions(), cf->getMinDist());
        }
        else {
            printf("time: %.4f | iters: %i | con: %i | min_dist: %s\n",
                   curr_time, stepper->iter, 0, "N/A");
        }
    }
    else if (ff) {
        printf("time: %.4f | iters: %i | floor_con: %i | f_min_dist: %.6f\n",
               curr_time, stepper->iter, ff->num_contacts, ff->min_dist);
    }
    else {
        printf("time: %.4f | iters: %i\n",
               curr_time, stepper->iter);
    }
}


bool world::simulationRunning() const {
    if (curr_time < total_time)
        return true;
    else
    {
        if (verbosity) cout << "Completed simulation." << endl;
        return false;
    }
}


double world::getCoordinate(int i, int limb_idx)
{
    return soft_robots->limbs[limb_idx]->x[i];
}

double world::getShellCoordinate(int i, int shell_limb_idx)
{
    return soft_robots->shell_limbs[shell_limb_idx]->x[i];
}


VectorXd world::getM1(int i, int limb_idx)
{
    return soft_robots->limbs[limb_idx]->m1.row(i);
}


VectorXd world::getM2(int i, int limb_idx)
{
    return soft_robots->limbs[limb_idx]->m2.row(i);
}


int world::getTimeStep() const
{
    return time_step;
}


double world::getCurrentTime() const
{
    return curr_time;
}
