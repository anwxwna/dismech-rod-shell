#ifndef ELASTICSTRETCHINGFORCESHELL_H
#define ELASTICSTRETCHINGFORCESHELL_H

#include "mechanics/baseForce.h"

class baseTimeStepper;

class elasticStretchingForceShell : public baseForce
{
public:
    elasticStretchingForceShell(const shared_ptr<softRobots>& m_soft_robots);
    ~elasticStretchingForceShell() override;
    void computeForce(double dt) override;
    void computeForceAndJacobian(double dt) override;

private:
    double len, refLength;
    double epsX;
    Vector3d u;
    Vector3d dxx;
    Vector3d f;
    Matrix3d Id3;
    Matrix3d M0;
    Matrix<double,1,3> v;
    Matrix<double,6,6> Jss;
    
    double EA;
    int ind, ind1, ind2;	
};

#endif
