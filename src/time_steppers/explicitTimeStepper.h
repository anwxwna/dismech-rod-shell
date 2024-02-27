#ifndef EXPLICITTIMESTEPPER_H
#define EXPLICITTIMESTEPPER_H

#include "baseTimeStepper.h"

class explicitTimeStepper : public baseTimeStepper
{
public:
    explicitTimeStepper(const shared_ptr<softRobots>& soft_robots,
                        const shared_ptr<forceContainer>& forces,
                        const simParams& sim_params);
    ~explicitTimeStepper() override;

    void prepSystemForIteration() override;
    void integrator() override;
    void addJacobian(int ind1, int ind2, double p, int limb_indx) override;
    void addJacobianShell(int ind1, int ind2, double p, int shell_limb_indx) override;
    void addJacobianShell(int ind1, int ind2, double p, int shell_limb_indx1, int shell_limb_indx2) override;
    void addJacobian(int ind1, int ind2, double p, int limb_indx1, int limb_idx2) override;

protected:
    vector<VectorXd> inverse_masses;

private:
    void constructInverseMassVector();

};


#endif
