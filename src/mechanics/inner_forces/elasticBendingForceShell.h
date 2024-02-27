#ifndef ELASTICBENDINGFORCESHELL_H
#define ELASTICBENDINGFORCESHELL_H

#include "mechanics/baseForce.h"

class baseTimeStepper;

class elasticBendingForceShell : public baseForce
{
public:
    elasticBendingForceShell(const shared_ptr<softRobots>& m_soft_robots);
    ~elasticBendingForceShell() override;
    void computeForce(double dt) override;
    void computeForceAndJacobian(double dt) override;

private:
    void initValues(Vector3d m_e0, Vector3d m_e1 , Vector3d m_e2, Vector3d m_e3, Vector3d m_e4);
    void computeGradTheta(VectorXd& gradTheta);
    void computeHessTheta(MatrixXd& hessTheta);
    double signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n);
    Matrix3d s(Matrix3d &mat);

    vector<double> EIMatrices;
    vector<int> hinge;
    double theta;
    MatrixXd hessTheta;
    VectorXd gradTheta;
    MatrixXd Jbb;
    VectorXd f;

    Vector3d m_e0, m_e1, m_e2, m_e3, m_e4;

    Vector3d m_nn1;
    Vector3d m_nn2;
    Vector3d m_m1;
    Vector3d m_m2;
    Vector3d m_m3;
    Vector3d m_m4;
    Vector3d m_m01;
    Vector3d m_m02;


    double m_cosA1;
    double m_cosA2;
    double m_cosA3;
    double m_cosA4;

    double m_sinA1;
    double m_sinA2;
    double m_sinA3;
    double m_sinA4;

    double m_psi;
    double m_zeta;
    double m_xi;

    double m_h1;
    double m_h2;
    double m_h3;
    double m_h4;
    double m_h01;
    double m_h02;

    Matrix3d Id3;
};

#endif