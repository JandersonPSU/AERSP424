#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


double ffunc(double t,double &p,double &q, double &r);
void derivative_function(vector<double> x, vector<vector<double>> dcm, vector<double> v_body, double time,vector<double> &xdot,double &Cdot);
double dotProduct(vector<vector<double>> a, vector<vector<double>> b);
  

int main()
{
    const double PI=3.1414926;
    double dt = 0.2; //time step
    vector<double> x0 = {0.0, 0.0, 0.0}; // radians
    vector<double> V_body = {60.0 * 6076.0 / 3600.0, 0.0, 0.0}; // Convert knots to ft/sec
    vector<vector<double>> x; // 2D vector to store state over time
    x.push_back({x0[0], x0[1], x0[2], V_body[0], V_body[1], V_body[2], 0.0, 0.0, 0.0});
    vector<double> xd; //derivatives vector
    vector<double> t = {0.0}; //time vector
    double tmax = 60.0; // Simulation time limit (60 seconds)
    int func_calls = 0; // To track function calls
    double phi = x[0][0];
    double theta = x[0][1];
    double psi = x[0][2];

    // Initialize DCM (Direction Cosine Matrix) as identity matrix (3x3)
    vector<vector<double>> identity_DCM = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    // phi DCM
    vector<vector<double>> DCM_phi = {
        {1.0, 0.0, 0.0},
        {0.0, cos(phi), -sin(phi)},
        {0.0, sin(phi), cos(phi)}
    };

    // theta DCM
    vector<vector<double>> DCM_theta = {
        {cos(theta), 0.0, sin(theta)},
        {0.0, 1, 0.0},
        {-sin(theta), 0, cos(theta)}
    };
    // psi DCM
    vector<vector<double>> DCM_psi = {
        {cos(psi), -sin(psi), 0.0},
        {sin(psi), cos(psi), 1},
        {0.0, 0.0, 1}
    };

    //Multiply the DCMs
    vector<vector<double>> intermediate_DCM;
    vector<vector<double>> DCM;
    intermediate_DCM=multiplyMatrices(identity_DCM,DCM_phi);
    intermediate_DCM=multiplyMatrices(intermediate_DCM,DCM_theta);
    DCM=multiplyMatrices(intermediate_DCM,DCM_psi);

    double current_time;
    
    vector<double> xdot1;
    vector<vector<double>> Cdot1;
    vector<double> xdot2;
    vector<vector<double>> Cdot2;
    vector<double> xdot3;
    vector<vector<double>> Cdot3;
    vector<double> xdot4;
    vector<vector<double>> Cdot4;
    for (int index=0;index<(int)(tmax/dt);index++){
        current_time= index*dt;
        vector<double> xn= x[index];

        derivative_function(xn,DCM,V_body,current_time,xdot1,Cdot1);
        derivative_function(,xdot2,Cdot2);
        derivative_function(,xdot3,Cdot3);
        derivative_function(,xdot4,Cdot4);
        
    }
}





//Function to get pqr at each time step:
double ffunc(double t,double &p,double &q, double &r){
 const double PI=3.1415926;
 p=PI/6;
 q=cos(t*6/PI);
 r=3*sin(t*30/PI);
}


//Function to get derivative:
void derivative_function(vector<double> x, vector<vector<double>> dcm, vector<double> v_body, double time,vector<double> &xdot,vector<vector<double>> &Cdot){
    double phi=x[0];
    double theta=x[1];
    double psi=x[2];

    double p,q,r;
    ffunc(time,p,q,r);
    //Might need to change to a 3x3, with 0s filling empty spaces
    vector<double> pqr_vector={p,q,r};
    
    vector<vector<double>> pqr_matrix = {
        {p, 0.0, 0.0},
        {q, 0.0, 0.0},
        {r, 0.0, 0.0}
    };
    
    //Gimbal Eqn:
    vector<vector<double>> qdot(3, std::vector<double>(3));
    qdot[0][0] = 1;
    qdot[0][1] = tan(theta) * sin(phi);
    qdot[0][2] = tan(theta) * cos(phi);

    qdot[1][0] = 0;
    qdot[1][1] = cos(phi);
    qdot[1][2] = -sin(phi);

    qdot[2][0] = 0;
    qdot[2][1] = sin(phi) / cos(theta);
    qdot[2][2] = cos(phi) / cos(theta);

    vector<vector<double>> qdot_matrix;
    qdot_matrix=multiplyMatrices(qdot,pqr_matrix);

    //Strapdown:
    vector<vector<double>> strapdown={
        {0.0, -r, q},
        {r, 0.0, -p},
        {-q, p, 0.0}
    };
    
    //Cdot solving
    Cdot=multiplyMatrices(dcm,strapdown);

    vector<vector<double>> V_NED;
    //getting v_body into a useable form, MIGHT want to change this later( check if need row or column vector)
    vector<vector<double>> v_body_matrix={
        {v_body[0],0.0, 0.0},
        {v_body[1], 0.0, 0.0},
        {v_body[2], 0.0, 0.0}
    };
    //V_NED solving
    V_NED=multiplyMatrices(dcm,v_body_matrix);
    vector<double> V_NED_Used={V_NED[0][0],V_NED[1][0],V_NED[2][0]};
    //Filling xdot
    xdot.insert(xdot.end(), qdot.begin(), qdot.end());

    // Add zero elements to xdot
    xdot.push_back(0.0);
    xdot.push_back(0.0);
    xdot.push_back(0.0);

    // Add V_NED elements to xdot
    xdot.insert(xdot.end(), V_NED.begin(), V_NED.end());
}

//Function to compute the dot product of two 3x1 vectors 
double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2) {
   
    double result;

    // Calculate the dot product
    for (size_t i = 0; i < 3; ++i) {
        result += vec1[i] * vec2[i];
    }

    return result;
}

//Function to compute 3x3 matrix multiplication:
vector<vector<double>> multiplyMatrices(vector<vector<double>>& mat1, vector<vector<double>>& mat2) {
    // Initialize the result matrix with zeros
    vector<vector<double>> result(3, vector<double>(3, 0.0));

    // Perform matrix multiplication
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}