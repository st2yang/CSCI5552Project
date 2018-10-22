#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* abs */ /* srand, rand, abs */
#include <math.h>       /* cos */
#include <algorithm>    // std::sort
#include <time.h>       /* time */
#include <map>


class LineExtraction_RANSAC{
/* This code is for line extraction by RANSAC
%
% Input:
% d_est & theta_est: measurement value from 2D laser scanner
% thresh: maxmum distance from points to lines
% min_mum: minimum number of points on a line
% maxIter: maximum number of iterations
% min_rest: minimum number of points left after iterations
%
% Output: 
% line: N*3 matrix, each row represents 3 parameters (cos(theta),sin(theta),d) of a line in
% 2D space (polar coordinate). 
% endpoints: 4*N matrix, each column represents two end points for a
% line(theta from small to large).
% R_mat: 2N*2 matrix, each 2x2 matrix, represent the corresponding line
% parameters' covariance matrix, d_star(the first one), theta_star
*/

public:
  Eigen::Matrix<float,-1,-1> line;
  Eigen::Matrix<float,-1,-1> endpoints_in_line;
  Eigen::Matrix<float,-1,-1> R_mat;
  
  LineExtraction_RANSAC(){;}
  virtual ~LineExtraction_RANSAC() {;}
  
  void Init(Eigen::MatrixXf d_est_int,Eigen::MatrixXf theta_est_int,Eigen::Matrix2f Cov_Laser_int,float thresh_int,int min_num_int,int maxIter_int,int min_rest_int,float theta_thresh_int, float seg_dist_thresh_int, float line_seg_dist_thresh_int);
  void LineExtraction_main();
  void LineExtraction(Eigen::Matrix<float,-1,-1>& endpoints);
  void LineSegment_Diff_method(Eigen::MatrixXf d_inlier, Eigen::MatrixXf theta_inlier,float seg_dist_thresh,int min_num,Eigen::Matrix2f Cov_Laser,Eigen::Matrix<float,-1,-1>& NewLines,Eigen::Matrix<float,-1,-1>& NewEndpoints,Eigen::Matrix<float,-1,-1>& R_total,Eigen::Matrix<int,-1,-1>& EraseIdx);
  void line_para_cal(Eigen::MatrixXf d_est,Eigen::MatrixXf theta_est,float& d_H,float& theta_H);
  Eigen::MatrixXf Ob_Line_Cov_cal(Eigen::MatrixXf d_all, Eigen::MatrixXf theta_all,float d_star, float theta_star,Eigen::MatrixXf Cov_single_point);
  
  
  void removeRow(Eigen::MatrixXf& matrix, int rowToRemove);
  void removeRow(Eigen::MatrixXi& matrix, int rowToRemove);
  void removeColumn(Eigen::MatrixXf& matrix, int colToRemove);
  Eigen::MatrixXi sort_indexes(Eigen::MatrixXf v);
  Eigen::MatrixXi sort_matrix_int(Eigen::MatrixXi v);
  void AddMatrixToRow(Eigen::MatrixXf& M,Eigen::MatrixXf m);
  void AddMatrixToRow(Eigen::MatrixXi& M,Eigen::MatrixXi m);
  void AddMatrixToCol(Eigen::MatrixXf& M,Eigen::MatrixXf m);
  void AddRow(Eigen::MatrixXi& M,int m);
  Eigen::MatrixXi GenerateUniqueInt(int tam, int number);
  bool Contains(std::vector<int> generatedValues,int num);

private:
    Eigen::MatrixXf d_est;
    Eigen::MatrixXf theta_est;
    Eigen::Matrix2f Cov_Laser;
    float thresh;
    int min_num;
    int maxIter;
    int min_rest;
    float theta_thresh;
    float seg_dist_thresh;
    float line_seg_dist_thresh;

};
