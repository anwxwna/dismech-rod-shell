#include "elasticShell.h"


elasticShell::elasticShell(int limb_idx, const vector<Vector3d>& nodes, const vector<vector<int>>& Face_Nodes, double rho, double thickness,
               double youngs_modulus, double poisson_ratio, double mu) :
                       limb_idx(limb_idx), ndof(nodes.size()*3), nv(nodes.size()), nf(Face_Nodes.size()),
                       all_nodes(nodes), Face_Nodes(Face_Nodes), rho(rho), thickness(thickness), youngM(youngs_modulus), poisson_ratio(poisson_ratio), mu(mu)
{
    // vector<int> vec {1,2,3};
    // EdgeIsBet.push_back(vec);
    // Calculate all geometrical parameters and initialize dof vector, velocity vector, undeformed edge lengths, undeformed hinge angles
    setup();
    // cout<<limb_idx<<endl;
}

void elasticShell::setup()
{
    int edge_index = 0;
    int hinge_index = 0;

    vector<int> third_node;

    for (int c = 0; c < nf; c++) {
        int node1_number = (Face_Nodes)[c][0];
        int node2_number = (Face_Nodes)[c][1];
        int node3_number = (Face_Nodes)[c][2];

        Vector3d node1_position = (all_nodes)[node1_number];
        Vector3d node2_position = (all_nodes)[node2_number];
        Vector3d node3_position = (all_nodes)[node3_number];

        Vector3d edge1 = node3_position - node2_position;
        Vector3d edge2 = node1_position - node3_position;
        Vector3d edge3 = node2_position - node1_position;

        Vector3d face_normal = edge3.cross(edge1);

        Vector3d face_unit_norm = face_normal.normalized();

        normals.push_back(face_unit_norm);

        vector<int> edge1_between{ node2_number, node3_number };
        vector<int> edge1_between_negative{ node3_number, node2_number };

        vector<int> edge2_between{ node3_number, node1_number };
        vector<int> edge2_between_negative{ node1_number, node3_number };

        vector<int> edge3_between{ node1_number, node2_number };
        vector<int> edge3_between_negative{ node2_number, node1_number };

        Face_Edges.push_back(std::vector<int>());
        sign_face_edges.push_back(std::vector<int>());

        //////////////////////////////////////////////////////////////
        /* Edge1 */
        bool bool_edge1_between_absent = true;
        bool bool_edge1_between_negative_absent = true;
        //int edge1_between_present = 0;
        int Locb1_between = -1;

        // check if edge1 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet.size(); i++) {
            if (edge1_between == (EdgeIsBet)[i]) {
                bool_edge1_between_absent = false;
                Locb1_between = i;
            }
            else if (edge1_between_negative == (EdgeIsBet)[i]) {
                bool_edge1_between_negative_absent = false;
                Locb1_between = i;
            }
        }
        // if not present, it is just an edge for now
        if (bool_edge1_between_absent && bool_edge1_between_negative_absent) {
            Edges.push_back(edge1);
            EdgeIsBet.push_back(edge1_between);
            third_node.push_back(node1_number);

            Face_Edges[c].push_back(edge_index);
            edge_avg_normals.push_back(face_unit_norm);
            (sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb1_between];
            int third_node_new = node1_number;

            HingeIsBet.push_back({ node2_number, node3_number, third_node_old, third_node_new });
            (Face_Edges)[c].push_back(Locb1_between);
            edge_avg_normals.push_back(((edge_avg_normals)[Locb1_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge1_between_absent) {
                (sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge1_between_negative_absent) {
                (sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }

        }
        //////////////////////////////////////////////////////////////
        /* Edge2 */
        bool bool_edge2_between_absent = true;
        bool bool_edge2_between_negative_absent = true;
        int Locb2_between = -1;

        // check if edge2 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet.size(); i++) {
            if (edge2_between == (EdgeIsBet)[i]) {
                bool_edge2_between_absent = false;
                Locb2_between = i;
            }
            else if (edge2_between_negative == (EdgeIsBet)[i]) {
                bool_edge2_between_negative_absent = false;
                Locb2_between = i;
            }
        }

        // if not present, it is just an edge for now
        if (bool_edge2_between_absent && bool_edge2_between_negative_absent) {
            Edges.push_back(edge2);
            EdgeIsBet.push_back(edge2_between);
            third_node.push_back(node2_number);

            (Face_Edges)[c].push_back(edge_index);
            edge_avg_normals.push_back(face_unit_norm);
            (sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb2_between];
            int third_node_new = node2_number;

            HingeIsBet.push_back({ node3_number, node1_number, third_node_old, third_node_new });
            (Face_Edges)[c].push_back(Locb2_between);
            edge_avg_normals.push_back(((edge_avg_normals)[Locb2_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge2_between_absent) {
                (sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge2_between_negative_absent) {
                (sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }
        }
        //////////////////////////////////////////////////////////////
        /* Edge3 */
        bool bool_edge3_between_absent = true;
        bool bool_edge3_between_negative_absent = true;
        int Locb3_between = -1;

        // check if edge2 is already present in edge_is_between array
        for (int i = 0; i < EdgeIsBet.size(); i++) {
            if (edge3_between == (EdgeIsBet)[i]) {
                bool_edge3_between_absent = false;
                Locb3_between = i;
            }
            else if (edge3_between_negative == (EdgeIsBet)[i]) {
                bool_edge3_between_negative_absent = false;
                Locb3_between = i;
            }
        }

        // if not present, it is just an edge for now
        if (bool_edge3_between_absent && bool_edge3_between_negative_absent) {
            Edges.push_back(edge3);
            EdgeIsBet.push_back(edge3_between);
            third_node.push_back(node3_number);

            (Face_Edges)[c].push_back(edge_index);
            edge_avg_normals.push_back(face_unit_norm);
            (sign_face_edges)[c].push_back(1);
            edge_index++;
        }
        // if present, its a hinge
        else {
            int third_node_old = third_node[Locb3_between];
            int third_node_new = node3_number;

            HingeIsBet.push_back({ node1_number, node2_number, third_node_old, third_node_new });
            (Face_Edges)[c].push_back(Locb3_between);
            edge_avg_normals.push_back(((edge_avg_normals)[Locb3_between] + face_unit_norm) * 0.5);
            hinge_index++;

            if (!bool_edge3_between_absent) {
                (sign_face_edges)[c].push_back(1);
            }
            else if (!bool_edge3_between_negative_absent) {
                (sign_face_edges)[c].push_back(-1);
            }
            else {
                cout << "error in edge sign finding" << endl;
            }
        }
    }

    ne = edge_index;
    nh = hinge_index;

    // check for debugging: check if size of edges and hinges vecs are correct
    if (EdgeIsBet.size() != ne || HingeIsBet.size() != nh) {
        cout << "error in edge/hinges vectors" << endl;
    }
    // Initialize dof vector
        x = VectorXd::Zero(ndof);
    for (int i = 0; i < nv; i++)
    {
        x(3 * i) = all_nodes[i](0);
        x(3 * i + 1) = all_nodes[i](1);
        x(3 * i + 2) = all_nodes[i](2);
    }
    x0 = x;

    // Initialize velocity vector
    u = VectorXd::Zero(ndof);
    u0 = u;

    // code below is same as elasticRod::setup()
    // We will start off with an unconstrained system
    ncons = 0;
    uncons = ndof;
    unique_dof = ndof;
    isConstrained = new int[ndof];
    isDOFJoint = new int[ndof];
    isNodeJoint = new int[nv];
    isEdgeJoint = new int[ne];
    DOFoffsets = new int[ndof];
    for (int i = 0; i < ndof; i++)
    {
        isConstrained[i] = 0;
        isDOFJoint[i] = 0;
        DOFoffsets[i] = 0;
    }
    for (int i = 0; i < nv; i++)
    {
        isNodeJoint[i] = 0;
        pair<int, int> non_joint{i, limb_idx};
        joint_ids.push_back(non_joint);
    }
    for (int i = 0; i < ne; i++)
    {
        isEdgeJoint[i] = 0;
    }

    // Setup the map from free dofs to all dof
    unconstrainedMap = new int[uncons]; // maps xUncons to x
    fullToUnconsMap = new int[ndof];
    setupMap();

    // compute reference lengths
    setReferenceLength();
    // compute natural hinge-angles
    setPhiBar();
    // set mass array
    setMass();
    // set tangent
    tangent = MatrixXd::Zero(ne, 3);
    computeTangent();

    // compute edge length
    edge_len = VectorXd(ne);
    computeEdgeLen();
    // compute phi's
    phi = VectorXd(ne);
    computeHingeAngles();

    // compute elastic stiffness
    computeElasticStiffness();

    return;
}

// destructor
elasticShell::~elasticShell()
{
    delete[] isConstrained;
    delete[] unconstrainedMap;
    delete[] fullToUnconsMap;
    delete[] isDOFJoint;
    delete[] isNodeJoint;
    delete[] isEdgeJoint;
    delete[] DOFoffsets; 
}

void elasticShell::computeElasticStiffness()
{
    double mean_edge_len =  accumulate(ref_len.begin(), ref_len.end(), 0.0)/ne; //m (avg edge length)
    
    EI = 2 * youngM * pow(thickness, 3) / (sqrt(3) * 12);;
    EA = 0.5 * sqrt(3) * youngM * thickness * pow(mean_edge_len, 2);
}

double elasticShell::getHingeAngle(int k)
{ 
    //         x2
    //         /\
    //        /  \
    //     e1/    \e3
    //      /  t0  \
    //     /        \
    //    /    e0    \
    //  x0------------x1
    //    \          /
    //     \   t1   /
    //      \      /
    //     e2\    /e4
    //        \  /
    //         \/
    //         x3
    //
    // Edge orientation: e0,e1,e2 point away from x0
    //                      e3,e4 point away from x1
    int n0 = HingeIsBet[k][0];
    int n1 = HingeIsBet[k][1];
    int n2 = HingeIsBet[k][2];
    int n3 = HingeIsBet[k][3];

    Vector3d node0, node1, node2, node3;
    node0 << x(3*n0), x(3*n0 + 1), x(3*n0 + 2); 
    node1 << x(3*n1), x(3*n1 + 1), x(3*n1 + 2); 
    node2 << x(3*n2), x(3*n2 + 1), x(3*n2 + 2); 
    node3 << x(3*n3), x(3*n3 + 1), x(3*n3 + 2); 

    Vector3d m_e0, m_e1, m_e2;
    m_e0 = node1 - node0;
    m_e1 = node2 - node0;
    m_e2 = node3 - node0;

    Vector3d normal0 = m_e0.cross(m_e1);
    Vector3d normal1 = m_e2.cross(m_e0);

    // cout<<"normal of face0: "<<normal0<<endl;
    // cout<<"normal of face1: "<<normal0<<endl;

    double theta = signedAngle(normal0, normal1, m_e0);

    return theta;
}

void elasticShell::addJoint(int node_num, bool remove_dof, int joint_node, int joint_limb)
{
    if (remove_dof)
    {
        unique_dof -= 3;
        //        uncons -= 3;
        for (int i = 4 * node_num + 3; i < ndof; i++)
        {
            DOFoffsets[i] -= 3;
        }

        isDOFJoint[4 * node_num] = 1;
        isDOFJoint[4 * node_num + 1] = 1;
        isDOFJoint[4 * node_num + 2] = 1;

        if (node_num == 0)
        {
            isNodeJoint[0] = 1;
            isEdgeJoint[0] = 1;
            joint_ids[0] = pair<int, int>(joint_node, joint_limb);
        }
        else if (node_num == nv - 1)
        {
            isNodeJoint[nv - 1] = 1;
            isEdgeJoint[nv - 2] = 1;
            joint_ids[nv - 1] = pair<int, int>(joint_node, joint_limb);
        }
        else
        {
            throw runtime_error("Tried removing dofs at the mid point of an edge.");
        }
    }
    else
    {
        isDOFJoint[4 * node_num] = 2;
        isDOFJoint[4 * node_num + 1] = 2;
        isDOFJoint[4 * node_num + 2] = 2;

        // NOTE: Might be able to delete this
        if (node_num == 0)
        {
            isNodeJoint[0] = 2;
            isEdgeJoint[0] = 2;
        }
        else if (node_num == nv - 1)
        {
            isNodeJoint[nv - 1] = 2;
            isEdgeJoint[nv - 2] = 2;
        }
        else
        {
            isNodeJoint[node_num] = 2;
            isEdgeJoint[node_num - 1] = 2;
            isEdgeJoint[node_num] = 2;
        }
    }
}

void elasticShell::updateMap()
{
    ncons = 0;
    for (int i = 0; i < ndof; i++)
    {
        if (isConstrained[i] > 0 || isDOFJoint[i] == 1)
        {
            ncons++;
        }
    }
    uncons = ndof - ncons;

    delete[] unconstrainedMap;
    delete[] fullToUnconsMap;
    // Setup the map from free dofs to all dof
    unconstrainedMap = new int[uncons]; // maps xUncons to x
    fullToUnconsMap = new int[ndof];
    setupMap();
}

void elasticShell::freeVertexBoundaryCondition(int k)
{
    isConstrained[3 * k] = 0;
    isConstrained[3 * k + 1] = 0;
    isConstrained[3 * k + 2] = 0;
}

// functions to handle boundary condition
void elasticShell::setVertexBoundaryCondition(Vector3d position, int k)
{
    isConstrained[3 * k] = 1;
    isConstrained[3 * k + 1] = 1;
    isConstrained[3 * k + 2] = 1;
    // Store in the constrained dof vector
    x(3 * k) = position(0);
    x(3 * k + 1) = position(1);
    x(3 * k + 2) = position(2);
    updateMap();
}

int elasticShell::getIfConstrained(int k) const
{
    return isConstrained[k];
}

void elasticShell::setMass () 
{

    mass_array = VectorXd::Zero(ndof);

    for (int i = 0; i<nf; i++){
        int node1ind = (Face_Nodes)[i][0];
        int node2ind = (Face_Nodes)[i][1];
        int node3ind = (Face_Nodes)[i][2];

        Vector3d node1_position = (all_nodes)[node1ind];
        Vector3d node2_position = (all_nodes)[node2ind];
        Vector3d node3_position = (all_nodes)[node3ind];

        Vector3d edge1 = node3_position - node2_position;
        Vector3d edge2 = node1_position - node3_position;
        Vector3d edge3 = node2_position - node1_position;

        Vector3d face_normal = edge3.cross(edge1);

        double face_A = 0.5 * face_normal.norm();
        double Mface = rho * face_A * thickness;
        int face_i_node_dofs[9] = {3*node1ind, 3*node1ind + 1, 3*node1ind + 2, \
        3*node2ind, 3*node2ind + 1, 3*node2ind + 2, 3*node3ind, 3*node3ind + 1, 3*node3ind + 2};
        
        for (int j : face_i_node_dofs){
            mass_array(j) += Mface/3 ;
        }
    }
}

void elasticShell::setReferenceLength()
{
    // This function is only run once at sim initilization
    ref_len = VectorXd(ne);
    for (int i = 0; i < ne; i++){
        ref_len(i) = (Edges[i].norm());
    }
}

void elasticShell::setPhiBar()
{   // This function is only run once at sim initilization
// cout<<"setPhiBar is called! "<<endl;
    phi_bar = VectorXd(nh);
    for (int j = 0; j < nh; j++){
        phi_bar(j) = getHingeAngle(j);
    }
}

Vector3d elasticShell::getVertex(int k)
{
    return x.segment<3>(3*k);
}

Vector3d elasticShell::getPreVertex(int k)
{
    return x0.segment<3>(3*k);
}

Vector3d elasticShell::getVelocity(int k)
{
    return u.segment<3>(3*k);
}

Vector3d elasticShell::getTangent(int k)
{
    return tangent.row(k);
}

void elasticShell::computeEdgeLen()
{
    for (int i = 0; i < ne; i++)
    {
        int node0ind = EdgeIsBet[i][0];
        int node1ind = EdgeIsBet[i][1];
        edge_len(i) = (x.segment(3 * node1ind, 3) - x.segment(3 * node0ind, 3)).norm();
    }
}

void elasticShell::computeHingeAngles()
{ 
    // cout<<"computeHingeAngles is called!"<<endl;
    for (int j = 0; j < nh; j++){
        phi(j) = getHingeAngle(j);
    }
}

void elasticShell::computeTangent()
{
    for (int i = 0; i < ne; i++)
    {
        int node0ind = EdgeIsBet[i][0];
        int node1ind = EdgeIsBet[i][1];
        tangent.row(i) = x.segment(3 * node1ind, 3) - x.segment(3 * node0ind, 3);
        tangent.row(i) = tangent.row(i)/((tangent.row(i)).norm());
    }
}

double elasticShell::signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n)
{
    // Compute the angle between two vectors
    Vector3d w = u.cross(v);
    double angle = atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0)
        return -angle;
    else
        return angle;
}

// void elasticShell::rotateAxisAngle(Vector3d &v, const Vector3d &z, const double &theta)
// {
//     // Compute the vector when it rotates along another vector into certain angle
//     if (theta != 0) // if theta=0, v = v
//     {
//         double cs, ss;
//         cs = cos(theta);
//         ss = sin(theta);
//         v = cs * v + ss * z.cross(v) + z.dot(v) * (1.0 - cs) * z;
//     }
// }

void elasticShell::setupMap()
{
    int c = 0;
    for (int i = 0; i < ndof; i++)
    {
        if (isConstrained[i] == 0 && isDOFJoint[i] != 1)
        {
            unconstrainedMap[c] = i;
            fullToUnconsMap[i] = c;
            c++;
        }
    }
}

void elasticShell::prepareForIteration()
{
    computeEdgeLen();
    computeHingeAngles();
    computeTangent();
}

double elasticShell::updateNewtonX(double *dx, int offset, double alpha)
{   int ind;
    double max_dx = 0;
    double curr_dx = 0;
    for (int c = 0; c < uncons; c++)
    {
        ind = unconstrainedMap[c]; 

        x[ind] -= alpha * dx[offset + c];
        curr_dx = abs(dx[offset + c]);
        if (curr_dx > max_dx) {
            max_dx = curr_dx;
        }
    }
    return max_dx;
}

void elasticShell::updateGuess(double weight, double dt)
{
    int ind;
    for (int c = 0; c < uncons; c++)
    {
        ind = unconstrainedMap[c];
        x[ind] = x0[ind] + weight * u[ind] * dt;
    }
}