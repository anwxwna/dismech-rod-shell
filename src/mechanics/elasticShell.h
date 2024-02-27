#ifndef ELASTICSHELL_H
#define ELASTICSHELL_H

#include "eigenIncludes.h"

class elasticShell
{
    // NOTE: probably could move more stuff to private
    public:
    elasticShell(int limb_idx, const vector<Vector3d>& nodes, const vector<vector<int> >& Face_Nodes, double rho, double thickness,
               double youngs_modulus, double poisson_ratio, double mu);
    ~elasticShell();


    // utility functions
    void setVertexBoundaryCondition(Vector3d position, int k);
    void prepareForIteration();
    double updateNewtonX(double *dx, int offset, double alpha=1.0);
    void updateGuess(double weight, double dt);

    int limb_idx;

    // utility functions
    Vector3d getVertex(int k);
    Vector3d getPreVertex(int k);
    Vector3d getVelocity(int k);
    Vector3d getTangent(int k);
    double getHingeAngle(int k);

    // Elastic stiffness values
    double poisson_ratio;
    double youngM; // Young's and shear modulus
    double EA;             // stretching stiffness
    double EI;             // bending stiffness

    int nv;                    // number of vertices
    int ne;                    // number of edges
    int nh;                    // number of hinges
    int nf;                    // number of faces
    int ndof;                  // number of degrees of freedom = 3*nv 
    int ncons;                 // number of constrained dof
    int uncons;                // number of unconstrained dof
    double rho;                // density
    double thickness;          // cross-sectional thickness of the shell
    vector<double> dm; // mass of each face
    double mu; // friction coefficient
    
    // Edge length
    VectorXd edge_len;
    // reference lengths
    VectorXd ref_len;
    // hinge angles
    VectorXd phi;
    // intial reference hinge angles
    VectorXd phi_bar;
    

    vector<Vector3d> all_nodes;
    vector<vector<int> > Face_Nodes;
    vector<Vector3d> Edges; 
    vector<vector<int> > Face_Edges;
    vector<vector<int> > sign_face_edges; 
    vector<vector<int> > EdgeIsBet; 
    vector<vector<int> > HingeIsBet;
    vector<Vector3d> edge_avg_normals;


    // dof vector before time step
    VectorXd x0;
    // dof vector after time step
    VectorXd x;
    // dof vector for line search /***** what is this?
    VectorXd x_ls;
    // velocity vector
    VectorXd u;
    VectorXd u0;

    // Face Normals
    vector<Vector3d> normals;
    // Tangents
    MatrixXd tangent;
    MatrixXd tangent_old;

    // lumped mass
    VectorXd mass_array;
    // lumped mass at unconstrained dofs
    VectorXd mUncons;

    // boundary conditions
    int *isConstrained;
    int getIfConstrained(int k) const;
    int *unconstrainedMap;
    int *fullToUnconsMap;

    void updateMap();

    void freeVertexBoundaryCondition(int k);

    void addJoint(int node_num, bool remove_dof, int joint_node, int joint_limb);
    int unique_dof;
    int *isDOFJoint;
    int *isNodeJoint;
    int *isEdgeJoint;
    int *DOFoffsets;
    vector<pair<int, int>> joint_ids;

private:
    void setupMap(); //********* not understood this clearly, kept as is for shell for now.

    // NOTE: perhaps move these to util.h later?
    void setup();
    void computeElasticStiffness();
    void setMass();
    void setReferenceLength();
    void setPhiBar();
    void computeTangent();
    void computeEdgeLen();
    void computeHingeAngles();
    static double signedAngle(const Vector3d &u, const Vector3d &v, const Vector3d &n);
};

#endif
