#ifndef LINEMATCH_H
#define LINEMATCH_H


#include <eigen3/Eigen/Dense>
#include <iostream>
#include "commons.h"

/*******************************************************************************
 * Class Definitions
 ******************************************************************************/

class LineMatch {
public:

    void performLineMatching(Eigen::MatrixXf& line_robot, Eigen::MatrixXf& R_line, Eigen::VectorXi& line_idx, Eigen::MatrixXf& endpoints_line_measure, Eigen::MatrixXf& line_un_match, Eigen::MatrixXf& endpoints_un_match, Eigen::MatrixXf& R_un_match,\
  Eigen::MatrixXf line_in_state, Eigen::MatrixXf endpoints_in_state, Eigen::MatrixXf line_test, Eigen::MatrixXf endpoints_test ,Eigen::MatrixXf R_test, Eigen::Vector3f robot_pose);

};


#endif //LINEMATCH_H
