#ifndef COMMONS_H
#define COMMONS_H

#include <eigen3/Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <vector>
#define orien_thresh 0.3925
#define dist_thresh 500

#define robot_r 95
#define robot_b 320
#define MAX_DIST 8000

inline Eigen::Vector2f  end_point_2_global(Eigen::Vector3f robot_pose,Eigen::Vector2f L_endpoints) {
    double x_robot = robot_pose(0,0);
    double y_robot = robot_pose(1,0);
    double theta_robot = robot_pose(2,0);

    Eigen::Matrix3f Rot_rob_2_gl;
    Rot_rob_2_gl << cos(theta_robot), -1*sin(theta_robot), x_robot, sin(theta_robot), cos(theta_robot), y_robot , 0, 0, 1;
    Eigen::Matrix3f Rot_line_2_rob;
    Rot_line_2_rob << cos(-M_PI/2), -1*sin(-M_PI/2), 0, sin(-M_PI/2), cos(-M_PI/2), 0, 0, 0, 1;

    Eigen::Vector3f L_endpoints_hom;
    L_endpoints_hom << L_endpoints(0), L_endpoints(1), 1;

    Eigen::Vector3f G_endpoints = Rot_rob_2_gl*Rot_line_2_rob*L_endpoints_hom;

    Eigen::Vector2f results;
    results << G_endpoints(0), G_endpoints(1);

    return results;

}

inline Eigen::Vector2f end_point_2_robot(Eigen::Vector2f L_endpoints){
    Eigen::Matrix3f Rot_line_2_rob;
    Rot_line_2_rob << cos(-M_PI/2), -1*sin(-M_PI/2), 0, sin(-M_PI/2), cos(-M_PI/2), 0, 0, 0, 1;
    Eigen::Vector3f L_endpoints_hom;
    L_endpoints_hom << L_endpoints(0), L_endpoints(1), 1;
    Eigen::Vector3f R_endpoints = Rot_line_2_rob*L_endpoints_hom;
    Eigen::Vector2f results;
    results << R_endpoints(0),R_endpoints(1);
    return results;
}

// checked
inline void  p_2_line_match(Eigen::MatrixXf& state_end_1, Eigen::MatrixXf& state_end_2, Eigen::MatrixXf& measure_end_1, Eigen::MatrixXf& measure_end_2,\
Eigen::MatrixXf line, Eigen::MatrixXf endpoint_ref_1, Eigen::MatrixXf endpoint_ref_2, Eigen::MatrixXf end1_test_gl, Eigen::MatrixXf end2_test_gl) {
    int line_num = line.cols();
    for(int i = 0; i<line_num;i++){
        float theta_ref = line(1,i);
        float d_ref = line(0,i);
        //projected points calculation
        state_end_1(0,i) = sin(theta_ref)*(sin(theta_ref)*endpoint_ref_1(0,i)-cos(theta_ref)*endpoint_ref_1(1,i))\
            -cos(theta_ref)*(-1*d_ref);
        state_end_1(1,i) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_ref_1(0,i)+cos(theta_ref)*endpoint_ref_1(1,i))\
            -sin(theta_ref)*(-1*d_ref);

        state_end_2(0,i) = sin(theta_ref)*(sin(theta_ref)*endpoint_ref_2(0,i)-cos(theta_ref)*endpoint_ref_2(1,i))\
            -cos(theta_ref)*(-1*d_ref);
        state_end_2(1,i) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_ref_2(0,i)+cos(theta_ref)*endpoint_ref_2(1,i))\
            -sin(theta_ref)*(-1*d_ref);

        measure_end_1(0,i) = sin(theta_ref)*(sin(theta_ref)*end1_test_gl(0,i)-cos(theta_ref)*end1_test_gl(1,i))\
            -cos(theta_ref)*(-1*d_ref);
        measure_end_1(1,i) = cos(theta_ref)*(-1*sin(theta_ref)*end1_test_gl(0,i)+cos(theta_ref)*end1_test_gl(1,i))\
            -sin(theta_ref)*(-1*d_ref);

        measure_end_2(0,i) = sin(theta_ref)*(sin(theta_ref)*end2_test_gl(0,i)-cos(theta_ref)*end2_test_gl(1,i))\
            -cos(theta_ref)*(-1*d_ref);
        measure_end_2(1,i) = cos(theta_ref)*(-1*sin(theta_ref)*end2_test_gl(0,i)+cos(theta_ref)*end2_test_gl(1,i))\
            -sin(theta_ref)*(-1*d_ref);

    }
}

inline void  p_2_line_state(Eigen::MatrixXf& endpoint_state_new,\
Eigen::VectorXf line, Eigen::MatrixXf endpoint_state) {
    int line_num = endpoint_state.cols();
    for(int i = 0; i<line_num;i++){
        float theta_ref = line(2*i+1,0);
        float d_ref = line(2*i,0);
        //projected points calculation
        endpoint_state_new(0,i) = sin(theta_ref)*(sin(theta_ref)*endpoint_state(0,i)-cos(theta_ref)*endpoint_state(2,i))\
            -cos(theta_ref)*(-1*d_ref);
        endpoint_state_new(2,i) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_state(0,i)+cos(theta_ref)*endpoint_state(2,i))\
            -sin(theta_ref)*(-1*d_ref);
        endpoint_state_new(1,i) = sin(theta_ref)*(sin(theta_ref)*endpoint_state(1,i)-cos(theta_ref)*endpoint_state(3,i))\
            -cos(theta_ref)*(-1*d_ref);
        endpoint_state_new(3,i) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_state(1,i)+cos(theta_ref)*endpoint_state(3,i))\
            -sin(theta_ref)*(-1*d_ref);


    }
}

// checked
inline void line_para_cal(Eigen::MatrixXf d_est,Eigen::MatrixXf theta_est,float& d_H,float& theta_H){
// compute the line parameters by least square 
    Eigen::MatrixXf u_stack(d_est.rows(),2);
    u_stack.col(0) = d_est.array()*theta_est.array().cos();
    u_stack.col(1) = d_est.array()*theta_est.array().sin();
    Eigen::MatrixXf u_vec = u_stack.colwise().mean().transpose(); //2x1 vector
    Eigen::MatrixXf ones = Eigen::MatrixXf::Ones(1,theta_est.rows());
    Eigen::MatrixXf dif_u_p = u_vec*ones-u_stack.transpose();
    Eigen::MatrixXf A = dif_u_p*dif_u_p.transpose();
    Eigen::EigenSolver<Eigen::MatrixXf> es(A);
    Eigen::MatrixXf D = es.eigenvalues().real();
    Eigen::MatrixXf V = es.eigenvectors().real();
    int idx_min_r;
    int idx_min_c;
    float minNorm = D.minCoeff(&idx_min_r,&idx_min_c);
    Eigen::MatrixXf V_min = V.block(0,idx_min_r,V.rows(),1);
    float cos_theta= V_min(0);
    float sin_theta= V_min(1);
    Eigen::MatrixXf uv_product = u_vec.transpose()*V_min;
    float d_theta = uv_product(0,0);
    Eigen::Matrix<float,3,1> coef_est(cos_theta,sin_theta,d_theta);
    if(d_theta<0)
        coef_est = -coef_est;
    theta_H = atan2(coef_est(1),coef_est(0));
    d_H = coef_est(2);
}


inline void AddMatrixToRow(Eigen::MatrixXf& M,Eigen::MatrixXf m){
    if(M.cols() == 0){
        M = m;
    } else{
        Eigen::MatrixXf f(M.rows()+m.rows(),M.cols());
        f.block(0,0,M.rows(),M.cols()) = M;
        f.block(M.rows(),0,m.rows(),m.cols()) = m;
        M = f;
    }
}

inline void AddMatrixToRow(Eigen::MatrixXi& M,Eigen::MatrixXi m){
    if(M.cols() == 0){
        M = m;
    } else {
        Eigen::MatrixXi f(M.rows()+m.rows(),M.cols());
        f.block(0,0,M.rows(),M.cols()) = M;
        f.block(M.rows(),0,m.rows(),m.cols()) = m;
        M = f;
    }

}

inline void AddMatrixToCol(Eigen::MatrixXf& M,Eigen::MatrixXf m){
    if(M.rows() == 0){
        M = m;
    }
    else{
        Eigen::MatrixXf f(M.rows(),M.cols()+m.cols());
        f.block(0,0,M.rows(),M.cols()) = M;
        f.block(0,M.cols(),m.rows(),m.cols()) = m;
        M = f;}
}

inline Eigen::VectorXi sort_vector_int(Eigen::VectorXi v){
    std::vector<int> v1;
    for(int i=0;i<v.rows();i++)
        v1.push_back(v(i));
    std::sort (v1.begin(), v1.begin()+v1.size());
    for(int i=0;i<v.rows();i++)
        v(i) = v1[i];
    return v;
}

inline Eigen::VectorXi unique(Eigen::VectorXi v){
    std::vector<int> v1;
    for(int i=0;i<v.rows();i++)
        v1.push_back(v(i));
    std::sort (v1.begin(), v1.begin()+v1.size());

    int old = -1;
    std::vector<int> v2;

    for(int i=0;i<v.rows();i++) {
        if (v1[i] != old) {
            v2.push_back(v1[i]);
            old = v1[i];
        }
    }

    Eigen::VectorXi result(v2.size());
    for(int i=0;i<v2.size();i++) {
        result(i) = v2[v2.size()-i-1];
    }
    return result;
}

inline void removeColumn(Eigen::MatrixXf& matrix, int colToRemove) {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

inline void removeRow(Eigen::MatrixXf& matrix, int rowToRemove) {
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

#endif //COMMONS_H

