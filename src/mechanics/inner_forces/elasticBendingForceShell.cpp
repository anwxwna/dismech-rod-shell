#include "elasticBendingForceShell.h"
#include "time_steppers/baseTimeStepper.h"

elasticBendingForceShell::elasticBendingForceShell(const shared_ptr<softRobots>& m_soft_robots) :
                                         baseForce(m_soft_robots)
{
    Id3 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    for (const auto& shell_limb : soft_robots->shell_limbs) {
        double EI = shell_limb->EI;
        EIMatrices.push_back(EI);

        // int nh = shell_limb->nh;
        // gradthetas.push_back(make_shared<MatrixXd>(MatrixXd::Zero(nh, 12)));
    }

    // for (const auto& joint : soft_robots->joints) {
    //     int nb = joint->num_bending_combos;
    //     gradKappa1s.push_back(make_shared<MatrixXd>(MatrixXd::Zero(nb, 11)));
    //     gradKappa2s.push_back(make_shared<MatrixXd>(MatrixXd::Zero(nb, 11)));
    // }

    Jbb = MatrixXd::Zero(12,12);

    f = VectorXd::Zero(12);
}

elasticBendingForceShell::~elasticBendingForceShell()
{
    ;
}


// void elasticBendingForceShell::initValues() {}
void elasticBendingForceShell::initValues(Vector3d m_e0, Vector3d m_e1 , Vector3d m_e2, Vector3d m_e3, Vector3d m_e4) {
    m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    m_cosA3 = -m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    m_cosA4 = -m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    m_sinA1 = (m_e0.cross(m_e1)).norm() / (m_e0.norm() * m_e1.norm());
    m_sinA2 = (m_e0.cross(m_e2)).norm() / (m_e0.norm() * m_e2.norm());
    m_sinA3 = (-m_e0.cross(m_e3)).norm() / (m_e0.norm() * m_e3.norm());
    m_sinA4 = (-m_e0.cross(m_e4)).norm() / (m_e0.norm() * m_e4.norm());

    m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / (m_nn1.norm());
    m_nn2 = -m_e0.cross(m_e4);
    m_nn2 = m_nn2 / (m_nn2.norm());

    m_m1 = (m_nn1).cross(m_e1 / m_e1.norm());
    m_m2 = -(m_nn2).cross(m_e2 / m_e2.norm());
    m_m3 = -(m_nn1).cross(m_e3 / m_e3.norm());
    m_m4 = (m_nn2).cross(m_e4 / m_e4.norm());
    m_m01 = -(m_nn1).cross(m_e0 / m_e0.norm());
    m_m02 = (m_nn2).cross(m_e0 / m_e0.norm());

    m_h1 = m_e0.norm() * m_sinA1;
    m_h2 = m_e0.norm() * m_sinA2;
    m_h3 = m_e0.norm() * m_sinA3; // no negative (this has negative in plates_shells matlab working code) 
    m_h4 = m_e0.norm() * m_sinA4; // no negative (this has negative in plates_shells matlab working code) 
    m_h01 = m_e1.norm() * m_sinA1;
    m_h02 = m_e2.norm() * m_sinA2;

}

void elasticBendingForceShell::computeGradTheta(VectorXd& gradTheta) {
    gradTheta.segment(0, 3) = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4;
    gradTheta.segment(3, 3) = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2;
    gradTheta.segment(6, 3) = - m_nn1 / m_h01;
    gradTheta.segment(9, 3) = - m_nn2 / m_h02;
}

void elasticBendingForceShell::computeHessTheta(MatrixXd& hessTheta) {
    Matrix3d M331 = m_cosA3 / (m_h3 * m_h3) * m_m3 * m_nn1.transpose();
    Matrix3d M311 = m_cosA3 / (m_h3 * m_h1) * m_m1 * m_nn1.transpose();
    Matrix3d M131 = m_cosA1 / (m_h1 * m_h3) * m_m3 * m_nn1.transpose();
    Matrix3d M3011 = m_cosA3 / (m_h3 * m_h01) * m_m01 * m_nn1.transpose();
    Matrix3d M111 = m_cosA1 / (m_h1 * m_h1) * m_m1 * m_nn1.transpose();
    Matrix3d M1011 = m_cosA1 / (m_h1 * m_h01) * m_m01 * m_nn1.transpose();

    Matrix3d M442 = m_cosA4 / (m_h4 * m_h4) * m_m4 * m_nn2.transpose();
    Matrix3d M422 = m_cosA4 / (m_h4 * m_h2) * m_m2 * m_nn2.transpose();
    Matrix3d M242 = m_cosA2 / (m_h2 * m_h4) * m_m4 * m_nn2.transpose();
    Matrix3d M4022 = m_cosA4 / (m_h4 * m_h02) * m_m02 * m_nn2.transpose();
    Matrix3d M222 = m_cosA2 / (m_h2 * m_h2) * m_m2 * m_nn2.transpose();
    Matrix3d M2022 = m_cosA2 / (m_h2 * m_h02) * m_m02 * m_nn2.transpose();

    Matrix3d B1 = 1 / (pow(m_e0.norm(), 2)) * m_nn1 * m_m01.transpose();
    Matrix3d B2 = 1 / (pow(m_e0.norm(), 2)) * m_nn2 * m_m02.transpose();

    Matrix3d N13 = 1 / (m_h01 * m_h3) * m_nn1 * m_m3.transpose();
    Matrix3d N24 = 1 / (m_h02 * m_h4) * m_nn2 * m_m4.transpose();
    Matrix3d N11 = 1 / (m_h01 * m_h1) * m_nn1 * m_m1.transpose();
    Matrix3d N22 = 1 / (m_h02 * m_h2) * m_nn2 * m_m2.transpose();
    Matrix3d N101 = 1 / (m_h01 * m_h01) * m_nn1 * m_m01.transpose();
    Matrix3d N202 = 1 / (m_h02 * m_h02) * m_nn2 * m_m02.transpose();

    hessTheta.block(0,0, 3,3) = s(M331) - B1 + s(M442) - B2;
    hessTheta.block(0,3, 3,3) = M311 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
    hessTheta.block(0,6, 3,3) = M3011 - N13;
    hessTheta.block(0,9, 3,3) = M4022 - N24;
    hessTheta.block(3,3, 3,3) = s(M111) - B1 + s(M222) - B2;
    hessTheta.block(3,6, 3,3) = M1011 - N11;
    hessTheta.block(3,9, 3,3) = M2022 - N22;
    hessTheta.block(6,6, 3,3) = -s(N101);
    hessTheta.block(9,9, 3,3) = -s(N202);

    // symmetric matrix
    hessTheta.block(3,0, 3,3) = (hessTheta.block(0,3, 3,3)).transpose();
    hessTheta.block(6,0, 3,3) = (hessTheta.block(0,6, 3,3)).transpose();
    hessTheta.block(9,0, 3,3) = (hessTheta.block(0,9, 3,3)).transpose();
    hessTheta.block(6,3, 3,3) = (hessTheta.block(3,6, 3,3)).transpose();
    hessTheta.block(9,3, 3,3) = (hessTheta.block(3,9, 3,3)).transpose();
}

void elasticBendingForceShell::computeForce(double dt)
{
    int shell_limb_idx = 0;
    int ci = 0;
    for (const auto& shell_limb : soft_robots->shell_limbs) {
        
        for (int i = 1; i < shell_limb->nh; i++)
        {
            gradTheta = VectorXd::Zero(12);
            hessTheta = MatrixXd::Zero(12,12);

            hinge = shell_limb->HingeIsBet[i];

            int node0_ind =  hinge[0];
            int node1_ind =  hinge[1];
            int node2_ind =  hinge[2];
            int node3_ind =  hinge[3];

            int ind[12] = { \
            3 * node0_ind, 3 * node0_ind + 1, 3 * node0_ind + 2, \
            3 * node1_ind, 3 * node1_ind + 1, 3 * node1_ind + 2, \
            3 * node2_ind, 3 * node2_ind + 1, 3 * node2_ind + 2, \
            3 * node3_ind, 3 * node3_ind + 1, 3 * node3_ind + 2 };

            Vector3d node0, node1, node2, node3;
            node0 = shell_limb->getVertex(node0_ind);
            node1 = shell_limb->getVertex(node1_ind);
            node2 = shell_limb->getVertex(node2_ind);
            node3 = shell_limb->getVertex(node3_ind);

            m_e0 = node1 - node0;
            m_e1 = node2 - node0;
            m_e2 = node3 - node0;
            m_e3 = node2 - node1;
            m_e4 = node3 - node1;

            initValues(m_e0, m_e1, m_e2, m_e3, m_e4);

            computeGradTheta(gradTheta);
            computeHessTheta(hessTheta);

            f = EIMatrices[shell_limb_idx] * (shell_limb->phi(i) - shell_limb->phi_bar(i)) * gradTheta;

            // check here
            int j = 0;
            for (int k : ind) {
                stepper->addForceShell(k, f[j], shell_limb_idx); // subtracting elastic force (+f or -f??)
                j++;
            }
        }
        shell_limb_idx++;
    }
}


void elasticBendingForceShell::computeForceAndJacobian(double dt)
{
    int shell_limb_idx = 0;
    int ci = 0;
    for (const auto& shell_limb : soft_robots->shell_limbs) {
        
        for (int i = 1; i < shell_limb->nh; i++)
        {
            gradTheta = VectorXd::Zero(12);
            hessTheta = MatrixXd::Zero(12,12);

            hinge = shell_limb->HingeIsBet[i];

            int node0_ind =  hinge[0];
            int node1_ind =  hinge[1];
            int node2_ind =  hinge[2];
            int node3_ind =  hinge[3];

            int ind[12] = { \
            3 * node0_ind, 3 * node0_ind + 1, 3 * node0_ind + 2, \
            3 * node1_ind, 3 * node1_ind + 1, 3 * node1_ind + 2, \
            3 * node2_ind, 3 * node2_ind + 1, 3 * node2_ind + 2, \
            3 * node3_ind, 3 * node3_ind + 1, 3 * node3_ind + 2 };

            Vector3d node0, node1, node2, node3;
            node0 = shell_limb->getVertex(node0_ind);
            node1 = shell_limb->getVertex(node1_ind);
            node2 = shell_limb->getVertex(node2_ind);
            node3 = shell_limb->getVertex(node3_ind);

            m_e0 = node1 - node0;
            m_e1 = node2 - node0;
            m_e2 = node3 - node0;
            m_e3 = node2 - node1;
            m_e4 = node3 - node1;

            initValues(m_e0, m_e1, m_e2, m_e3, m_e4);

            computeGradTheta(gradTheta);
            computeHessTheta(hessTheta);

            // cout<<"gradTheta is: "<<gradTheta<<endl;
            // cout<<"phi is: "<< shell_limb->phi(i) << endl;
            // cout<<"phi bar is: "<< shell_limb->phi_bar(i) << endl;

            // cout<< "Hinge Angles are:"<<endl;

            f = EIMatrices[shell_limb_idx] * (shell_limb->phi(i) - shell_limb->phi_bar(i)) * gradTheta;

            int j = 0;
            for (int k : ind) {
                stepper->addForceShell(k, f[j], shell_limb_idx); // subtracting elastic force (+f or -f??)
                j++;
            }
            
            Jbb = EIMatrices[shell_limb_idx] * (gradTheta*gradTheta.transpose() + (shell_limb->phi(i))*hessTheta);
            // cout<<"bending force being added to the force vector: "<<-f<<endl;
            // cout<<"bending jacobian being added to the jacobian matrix: "<< -Jbb;

            int p_counter = 0, q_counter = 0;
            for (int j : ind) {
                q_counter = 0;
                for (int k : ind) {
                    stepper->addJacobianShell(j, k, Jbb(p_counter, q_counter), shell_limb_idx); // (+Jbb or -Jbb?)
                    q_counter++;
                }
                p_counter++;
            }   
        }
        shell_limb_idx++;
    }
}

double elasticBendingForceShell::signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n)
{
    // Compute the angle between two vectors
    Vector3d w = u.cross(v);
    double angle = atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0)
        return -angle;
    else
        return angle;
}
Matrix3d elasticBendingForceShell::s(Matrix3d &mat) {
    return (mat + mat.transpose());
}