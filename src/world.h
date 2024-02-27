#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"
#include "robotDescription.h"

#include "mechanics/softRobots.h"
#include "mechanics/forceContainer.h"

// include inner force classes
#include "mechanics/inner_forces/inertialForce.h"
#include "mechanics/inner_forces/elasticStretchingForce.h"
#include "mechanics/inner_forces/elasticBendingForce.h"
#include "mechanics/inner_forces/elasticTwistingForce.h"

// include shell elastic forces
#include "mechanics/inner_forces/elasticStretchingForceShell.h"
#include "mechanics/inner_forces/elasticBendingForceShell.h"

// include time stepper
#include "time_steppers/forwardEuler.h"
#include "time_steppers/verletPosition.h"
#include "time_steppers/backwardEuler.h"
#include "time_steppers/implicitMidpoint.h"


class world
{
public:
    world(const shared_ptr<softRobots>& soft_robots,
          const shared_ptr<forceContainer>& forces,
          const simParams& sim_params);
    ~world();
    void updateTimeStep();
    double getCoordinate(int i, int limb_idx);
    double getShellCoordinate(int i, int shell_limb_idx);
    VectorXd getM1(int i, int limb_idx);
    VectorXd getM2(int i, int limb_idx);

    double getCurrentTime() const;
    bool simulationRunning() const;
    int getTimeStep() const;
    void printSimData();

    shared_ptr<softRobots> soft_robots;

private:
    shared_ptr<forceContainer> forces;
    shared_ptr<baseTimeStepper> stepper;
    

    int time_step;
    double curr_time;
    double total_time;

    void updateCons();
};

#endif
