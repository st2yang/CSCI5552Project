#ifndef KALMANFILTER_H
#define KALMANFILTER_H


#include "commons.h"

class KalmanFilter {
public:
    
//     void state_propagation(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_minus,\
//      Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_minus_k_minus, Eigen::MatrixXf R_xytheta, Eigen::VectorXf u, float dt);
  
  void state_propagation(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_minus,\
     Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_minus_k_minus, Eigen::MatrixXf R_vw, Eigen::VectorXf u, float dt);

    void measurement_update(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new,\
    Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_k_minus, Eigen::MatrixXf line_robot, Eigen::MatrixXf R_line, Eigen::VectorXi line_idx);

    void update_mached_line_endpoint(Eigen::MatrixXf& endpoints_in_state_new,\
    Eigen::VectorXf state_all, Eigen::MatrixXf endpoints_in_state, Eigen::VectorXi line_idx, Eigen::MatrixXf endpoints_line_measure);

    void cov_update_new_line_measure(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new,\
    Eigen::MatrixXf Cov_k_k, Eigen::Matrix2f R_new_line, Eigen::VectorXf state_all, Eigen::Vector2f new_line_para);

    void cal_H_mat(Eigen::MatrixXf& H, int& eq_idx,\
     Eigen::Vector3f robot_gl,Eigen::Vector2f line_gl, int line_idx, int num_state);

    void add_line_2_state(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new, Eigen::Vector4f& end_point_gl,\
     Eigen::VectorXf state_all, Eigen::MatrixXf Cov, Eigen::Vector2f x_part, Eigen::Vector2f y_part, Eigen::Matrix2f R_add_single);

};


#endif //KALMANFILTER_H
