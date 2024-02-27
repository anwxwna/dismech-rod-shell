#include "robotDescription.h"

extern ofstream logging_output_file;  // defined in main.cpp
/*
 * Dynamic Cantilever Example
 *
 * Define your soft robot structure(s), boundary conditions,
 * custom external forces, and loggers in the function below.
 */

void get_robot_description(int argc, char** argv,
                           const shared_ptr<softRobots>& soft_robots,
                           const shared_ptr<forceContainer>& forces,
                           shared_ptr<worldLogger>& logger,
                           simParams& sim_params) {

    sim_params.dt = 1e-4;
    sim_params.sim_time = 2;
    sim_params.dtol = 1e-3;
    sim_params.show_mat_frames = false;
    sim_params.render_scale = 0.2;
    // sim_params.render = false;
    sim_params.enable_2d_sim = false;
    sim_params.structure_shell = true;
    sim_params.line_search = false; // this is having some issues in case of shell, need to check

    double thickness = 0.001; // 1mm
    double young_mod = 1e5;
    double density = 1200;
    double poisson = 0.5;
    double friction_coeff = 0;

    // Read the input txt file into:
    // Nodes
    // Face_Nodes

    vector<Vector3d> Nodes;
    vector<vector<int> > Face_Nodes;

    // string inputFileName = "examples/shell_hemisphere/input hemisphere for c++.txt"; // hemispherical shell
    // string inputFileName = "examples/shell_hemisphere/input hemisphere denser for c++.txt"; // hemispherical shell denser
    string inputFileName = "examples/shell_hemisphere/input hemisphere for c++ less dense.txt"; // hemispherical shell less dense

    inputProcessorFunc(inputFileName, Nodes, Face_Nodes);

    // debugging
    cout << "Node positions:" << endl;
    for (int i = 0; i < Nodes.size(); i++) {
        cout << (Nodes)[i] << endl;
    }
    cout << "Face Nodes:" << endl;
    for (int i = 0; i < Face_Nodes.size(); i++) {
        for (int j = 0; j < (Face_Nodes)[i].size(); j++) {
            cout << (Face_Nodes)[i][j] << ", ";
        }
        cout << ";" << endl;
    }

    // Create a shell beam along the x-y plane
    soft_robots->addShellLimb(Nodes, Face_Nodes, density, thickness, young_mod, poisson, friction_coeff);

    // cout<< "no. of vertices: " << soft_robots->shell_limbs[0]->nv<<endl;
    // cout<< "no. of edges: " << soft_robots->shell_limbs[0]->ne<<endl;
    // cout<< "no. of faces: " << soft_robots->shell_limbs[0]->nf<<endl;
    // cout<< "no. of hinges: " << soft_robots->shell_limbs[0]->nh<<endl;
    // cout<<"thickness is: "<<soft_robots->shell_limbs[0]->thickness<<endl;
    // cout<<"bending stiffness is: "<<soft_robots->shell_limbs[0]->EI<<endl;
    // cout<<"stretching stiffness is: "<<soft_robots->shell_limbs[0]->EA<<endl;

    // Fix center point
    // vector<int> Fixed_node_indices{ 27 };
    vector<int> Fixed_node_indices{ 15 }; // less dense
    // vector<int> Fixed_node_indices{ 35 }; // denser

    for (int i=0; i<Fixed_node_indices.size(); i++){
        soft_robots->lockNodeShell(0,Fixed_node_indices[i]);
    }

    // Add gravity
    Vector3d gravity_vec(0.0, 0.0, -9.8);
    forces->addForce(make_shared<gravityForce>(soft_robots, gravity_vec));


   // Set logger to record nodes
    string logfile_base = "log_files/hemispherical_shell";
    int logging_period = 1;
    logger = make_shared<shellNodeLogger>(logfile_base, convert_float_to_scientific_str(young_mod),
                                        logging_output_file, logging_period);
}

void inputProcessorFunc(const std::string inputFileName, vector<Vector3d>& Nodes, vector<vector<int>>& facenodes) {
    ifstream inputFile(inputFileName);
    string line, typeSpec;

    Nodes.clear();
    facenodes.clear();

    while (getline(inputFile, line)) {
        // Trim leading and trailing spaces
        line.erase(line.find_last_not_of(" \t") + 1);
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty()) {
            continue;
        } else if (line[0] == '*') {
            // Set typeSpec without leading spaces and convert to lowercase for case-insensitive comparison
            typeSpec = line.substr(1);
            transform(typeSpec.begin(), typeSpec.end(), typeSpec.begin(), ::tolower);
        } else {
            istringstream iss(line);
            string token;
            vector<string> dataLine;

            while (getline(iss, token, ',')) {
                dataLine.push_back(token);
            }

            if (dataLine.empty()) {
                continue;
            }
            // change "nodes\r" below to "nodes" if needed
            if (typeSpec == "nodes\r") {
                if (dataLine.size() != 3) {
                    cerr << "Warning. Invalid input for nodes." << endl;
                } else {
                    double x = stod(dataLine[0]);
                    double y = stod(dataLine[1]);
                    double z = stod(dataLine[2]);
                    Nodes.push_back(Vector3d(x, y, z));
                }
            } 
            // change "facenodes\r" below to "facenodes" if needed
            else if (typeSpec == "facenodes\r") {
                if (dataLine.size() != 3) {
                    cerr << "Warning. Invalid input for edges." << endl;
                } else {
                    int node1 = stoi(dataLine[0]);
                    int node2 = stoi(dataLine[1]);
                    int node3 = stoi(dataLine[2]);
                    facenodes.push_back({node1, node2, node3});
                }
            }
        }
    }

    inputFile.close();
}
