#include "LineExtraction_RANSAC.h"


void LineExtraction_RANSAC::Init(Eigen::MatrixXf d_est_int,Eigen::MatrixXf theta_est_int,Eigen::Matrix2f Cov_Laser_int,float thresh_int,int min_num_int,int maxIter_int,int min_rest_int,float theta_thresh_int, float seg_dist_thresh_int, float line_seg_dist_thresh_int){
    d_est = d_est_int;
    theta_est = theta_est_int;
    Cov_Laser = Cov_Laser_int;
    thresh = thresh_int;
    min_num = min_num_int;
    maxIter = maxIter_int;
    min_rest = min_rest_int;
    theta_thresh = theta_thresh_int;
    seg_dist_thresh = seg_dist_thresh_int;
    line_seg_dist_thresh = line_seg_dist_thresh_int;
}


void LineExtraction_RANSAC::LineExtraction_main(){ 
   /*
   % This function is the main call function to do the Line Extraction
   % its outputs are the line parameter, 
   % the end points on the line(not the end points that are usually not on the line)
   % Input: d_est, nx1 vector, the distance measurements for each point
   %        theta_est, nx1 vector, the angle measurements for each point
   % Output: line_para,N*3 matrix, each row represents 3 parameters (cos(theta),sin(theta),d) of a line in
   %                   2D space (polar coordinate).
   %         endpoints_in_line,4xN matrix, each column represents end points
   %                           in the order of (x1,x2,y1,y2)
   %         R_mat, 2Nx2 matrix, each 2x2 matrix represents the covariance of
   %                the estimated line parameters,d_star(the first one), theta_star
   % 
   */
   
   // 1. RANSAC
   Eigen::Matrix<float,-1,-1> endpoints;
   LineExtraction(endpoints);
   endpoints_in_line.resize(endpoints.rows(),endpoints.cols());
   //std::cout << "finish extraction " << std::endl;
   
   Eigen::MatrixXi idx_remove;
   for(int i = 0;i<line.cols();i++){
    // obtain each line parameters
    double cos_theta = cos(line(1,i));
    double sin_theta = sin(line(1,i));
    double d_theta = line(0,i);
    // project the start and end points to the line
    double start_x = endpoints(0,i)*cos(endpoints(1,i));
    double start_y = endpoints(0,i)*sin(endpoints(1,i));
    double end_x = endpoints(2,i)*cos(endpoints(3,i));
    double end_y = endpoints(2,i)*sin(endpoints(3,i));
    
    
    endpoints_in_line(0,i) = sin_theta*(sin_theta*start_x-cos_theta*start_y)-cos_theta*(-1*d_theta);
    endpoints_in_line(2,i) = cos_theta*(-1*sin_theta*start_x+cos_theta*start_y)-sin_theta*1*(-1*d_theta);
    
    endpoints_in_line(1,i) = sin_theta*(sin_theta*end_x-cos_theta*end_y)-cos_theta*(-1*d_theta);
    endpoints_in_line(3,i) = cos_theta*(-1*sin_theta*end_x+cos_theta*end_y)-sin_theta*(-1*d_theta);
    
    if(pow((pow((endpoints_in_line(0,i)-endpoints_in_line(1,i)),2)+pow((endpoints_in_line(2,i)-endpoints_in_line(3,i)),2)),0.5) < line_seg_dist_thresh)
        AddRow(idx_remove,i);
   }

   if(idx_remove.rows()>0){
     Eigen::MatrixXi idx_sort = sort_matrix_int(idx_remove);
     for(int i=idx_sort.rows()-1;i>0;i--){
         int RemoveLineIdx = idx_sort(i,0);
         removeColumn(line,RemoveLineIdx);
         removeColumn(endpoints_in_line,RemoveLineIdx);
         removeRow(R_mat,2*RemoveLineIdx);// remove the highest index first to make sure the index is correct
         removeRow(R_mat,2*RemoveLineIdx+1);
     }    
   }
}

void LineExtraction_RANSAC::LineExtraction(Eigen::Matrix<float,-1,-1>& endpoints){
  srand (time(NULL));
  int iter = 0;
  int outIter = 5;
while ((iter < outIter) and (d_est.rows() > min_rest)){
    int max_inlier_num = 0;
    Eigen::MatrixXi max_inlier_idx;
    for(int inner = 1;inner<=maxIter;inner++){
        
        // 1. randomly select two points, not just two, a small set
        Eigen::MatrixXi idx = GenerateUniqueInt(d_est.rows(),min_num);
	
        // 2.compute the line parameter by Linear Least Square
	float d_H;
	float theta_H;
	Eigen::MatrixXf d_est_select(idx.rows(),1);
	Eigen::MatrixXf theta_est_select(idx.rows(),1);
	for(int select_idx = 0;select_idx<idx.rows();select_idx++){
	   d_est_select(select_idx,0) = d_est(idx(select_idx,0),0);
	   theta_est_select(select_idx,0) = theta_est(idx(select_idx,0),0);}
        line_para_cal(d_est_select,theta_est_select,d_H, theta_H);
        
        // 3. compute distance from other points to this line
        Eigen::Matrix<int,-1,-1> inlier;
        int inlier_num = 0;
        Eigen::MatrixXf d_test = d_est; // use only un-selected points
        Eigen::MatrixXf theta_test = theta_est;
	Eigen::MatrixXi idx_sort = sort_matrix_int(idx);
        for(int row_id=idx_sort.rows()-1;row_id>0;row_id--){
	   int RemoveLineIdx = idx_sort(row_id,0);
           removeRow(d_test, RemoveLineIdx);
           removeRow(theta_test, RemoveLineIdx);
       }
       for(int i = 0;i<d_test.rows();i++){
            float dist = abs(d_test(i,0)*cos(theta_H - theta_test(i,0)) - d_H);
            if(dist < thresh){
                AddRow(inlier,i);
                inlier_num = inlier_num + 1;
            }
        }
        
        // 4. find the best line for the current remaining data points
        
        if(inlier_num >= min_num && inlier_num >=max_inlier_num){
	    max_inlier_idx.resize(inlier.rows()+idx.rows(),inlier.cols());
            max_inlier_num = inlier_num;
            max_inlier_idx.block(0,0,inlier.rows(),inlier.cols()) = inlier;
            max_inlier_idx.block(inlier.rows(),0,idx.rows(),1) = idx;
        }
        
    }
    
    if(max_inlier_num >0){
        Eigen::Matrix<float,-1,-1> NewLines;
        Eigen::Matrix<float,-1,-1> NewEndpoints;
        Eigen::Matrix<float,-1,-1> New_R;
        Eigen::Matrix<int,-1,-1> EraseIdx;
	
	Eigen::MatrixXf d_est_inlier(max_inlier_idx.rows(),d_est.cols());
	Eigen::MatrixXf theta_est_inlier(max_inlier_idx.rows(),theta_est.cols());
	for(int inlierId_id = 0;inlierId_id<max_inlier_idx.rows();inlierId_id++){
	  d_est_inlier.row(inlierId_id) = d_est.row(max_inlier_idx(inlierId_id,0));
	  theta_est_inlier.row(inlierId_id) = theta_est.row(max_inlier_idx(inlierId_id,0));
	}

	
        LineSegment_Diff_method(d_est_inlier,theta_est_inlier,seg_dist_thresh,min_num,Cov_Laser,NewLines,NewEndpoints,New_R,EraseIdx);

	//std::cout << "finish segment" <<std::endl;
	Eigen::MatrixXi idx_sort = sort_matrix_int(EraseIdx);
        for(int row_id=idx_sort.rows()-1;row_id>0;row_id--){
	    int remove_idx =  idx_sort(row_id);
            removeRow(max_inlier_idx, remove_idx);
	}
        AddMatrixToCol(line,NewLines);
        AddMatrixToCol(endpoints,NewEndpoints);
        AddMatrixToRow(R_mat,New_R);
        // clear inliers
	Eigen::MatrixXi idx_sort_2 = sort_matrix_int(max_inlier_idx);
        for(int row_id=idx_sort_2.rows()-1;row_id>0;row_id--){
	    int RemoveLineIdx = idx_sort_2(row_id);
            removeRow(d_est, RemoveLineIdx);
            removeRow(theta_est, RemoveLineIdx);
        }
    }
    iter = iter + 1;
}

}


void LineExtraction_RANSAC::LineSegment_Diff_method(Eigen::MatrixXf d_inlier, Eigen::MatrixXf theta_inlier,float seg_dist_thresh,int min_num,Eigen::Matrix2f Cov_Laser,Eigen::Matrix<float,-1,-1>& NewLines,Eigen::Matrix<float,-1,-1>& NewEndpoints,Eigen::Matrix<float,-1,-1>& R_total,Eigen::Matrix<int,-1,-1>& EraseIdx){
    /* This code is for line segmentation
% This code is slightly different from Wu Fei's code
% The differences: 
%   1. To determine whether to split line segmements, I use (x,y)
%   coordinate. If the two neighbor points distance are greater than seg_dist_thresh,
%   they are from two different segments.(The advantage is 
%   we can use the seg_dist_thresh to be as the robot's safe pass distance)
%   2. The begin_idx used in Wu Fei's code is confusing. begin_idx should
%   be updated whenever there are segments
%   3. For the grouped points, whose total number is too small to form a
%   line, I didn't re-use it in the next RANSAC line extraction.
% Input:
% d_inlier & theta_inlier: measurement value from 2D laser scanner inlier
% seg_dist_thresh: maxmum distance between neightbor points that are in the
%                  same line segment
% min_mum: minimum number of points on a line
%
% Output: 
% NewLines: 2xN matrix, each column represents 2 parameters (d,theta) of a line in
% 2D space (polar coordinate). 
% NewEndpoints: 4*N matrix, each column represents two end points for a
% line(theta from small to large).
% R_total: 2N*2 matrix, each 2x2 matrix represents the covariance of the
% current line parameters d_star(the first one), theta_star
% EraseIdx: whole index.
*/

// 1. sort theta from small to large

//Eigen::MatrixXi I = sort_indexes(theta_inlier);
Eigen::MatrixXi I(theta_inlier.rows(),1);
for(int inlier_idx=0;inlier_idx<theta_inlier.rows();inlier_idx++)
  I(inlier_idx,0) = inlier_idx;
Eigen::MatrixXf x_pos = d_inlier.array()*theta_inlier.array().cos();
Eigen::MatrixXf y_pos = d_inlier.array()*theta_inlier.array().sin();
//float mean_theta = abs(theta_inlier.mean());
int begin_idx = 0;

for(int t = 0;t<d_inlier.rows();t++){
    /* 2. compute difference of theta between two neighbors. 
    % If difference is larger than threshold: seg_dist_thresh,
    % segment the line, and treat the front part as a temporary new line
    */
    if(t<d_inlier.rows()-1){
        
        float dist_neighbor = sqrt(pow((x_pos(I(t))-x_pos(I(t+1))),2)+pow((y_pos(I(t))-y_pos(I(t+1))),2));
        // dist_theta = abs(theta_inlier(I(t)) - theta_inlier(I(t+1)));
        // if dist_theta > theta_thresh*mean_theta
        if (dist_neighbor > seg_dist_thresh){
        /*
        % 3. If the front part satisfied minimum-points-condition, add    
        % two end points and recompute the line parameters
        % try to include another criteria: length of the line segment
        */    
            //if ((t - begin_idx + 1) >= min_num) && sqrt((x_pos(I(begin_idx))-x_pos(I(t)))^2+(y_pos(I(begin_idx))-y_pos(I(t)))^2)>300
            if ((t - begin_idx + 1) >= min_num){
                Eigen::Matrix<float,4,1> NewBreakPoint;
		NewBreakPoint << d_inlier(I(begin_idx,0),0),theta_inlier(I(begin_idx,0),0),d_inlier(I(t,0),0),theta_inlier(I(t,0),0);
                AddMatrixToCol(NewEndpoints,NewBreakPoint);

                // recompute the line parameter by Linear Least Square
                Eigen::MatrixXf u_stack(t-begin_idx+1,2);
		Eigen::MatrixXi I_select_id = I.block(begin_idx,0,t - begin_idx + 1,1);
		Eigen::MatrixXf d_inlier_I(I_select_id.rows(),1);
		Eigen::MatrixXf theta_inlier_I(I_select_id.rows(),1);
		for(int I_select_id_iter=0;I_select_id_iter<I_select_id.rows();I_select_id_iter++){
		  d_inlier_I(I_select_id_iter,0) = d_inlier(I_select_id(I_select_id_iter),0);
		  theta_inlier_I(I_select_id_iter,0) = theta_inlier(I_select_id(I_select_id_iter),0);
		}
                u_stack.col(0) = d_inlier_I.array()*theta_inlier_I.array().cos();
                u_stack.col(1) = d_inlier_I.array()*theta_inlier_I.array().sin();
                Eigen::MatrixXf u_vec = u_stack.colwise().mean().transpose(); //2x1 vector
                Eigen::MatrixXf ones = Eigen::MatrixXf::Ones(1,t-begin_idx+1);
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

                /*
                % 4. compute the optimal solution theta_est, d_est in order
                % to calculate the covariance of the estimated theta_est,
                % d_est
                */
                float theta_est = atan2(coef_est(1),coef_est(0));
                float d_est = coef_est(2);

                /*
                % setup the sigma_d and sigma_theta
                % according to the material of 5551, distance accuracy is
                % +-15mm, angular resolution 0.5deg
                %Cov_single_point = [25,0;0,1e-4];%based on the material of 5551
                */
                Eigen::Matrix2f Cov_single_point = Cov_Laser;
                Eigen::MatrixXf R_est = Ob_Line_Cov_cal(d_inlier_I,theta_inlier_I,d_est,theta_est,Cov_single_point);
                AddMatrixToRow(R_total,R_est);
                Eigen::Matrix<float,2,1> coef_est_2;
		        coef_est_2 << d_est,theta_est;
                AddMatrixToCol(NewLines,coef_est_2);
		
            }
            else{
                // otherwise, add all the points in this part into outliers.
	        Eigen::MatrixXi I_select_id = I.block(begin_idx,0,t - begin_idx + 1,1);
                AddMatrixToRow(EraseIdx,I_select_id);
            }
            // change the begin point index
            begin_idx = t+1;
	}        
    }
    else{
        /*
        % 3. If the front part satisfied minimum-points-condition, add
        % two end points and recompute the line parameters
        */
        if ((t - begin_idx + 1) >= min_num){
            Eigen::Matrix<float,4,1> NewBreakPoint;
            NewBreakPoint << d_inlier(I(begin_idx)),theta_inlier(I(begin_idx)),d_inlier(I(t)),theta_inlier(I(t));
            AddMatrixToCol(NewEndpoints,NewBreakPoint);

            // recompute the line parameter by Linear Least Square
            Eigen::MatrixXf u_stack(t-begin_idx+1,2);
	    Eigen::MatrixXi I_select_id = I.block(begin_idx,0,t - begin_idx + 1,1);
	    Eigen::MatrixXf d_inlier_I(I_select_id.rows(),1);
	    Eigen::MatrixXf theta_inlier_I(I_select_id.rows(),1);
	    for(int I_select_id_iter=0;I_select_id_iter<I_select_id.rows();I_select_id_iter++){
		d_inlier_I(I_select_id_iter,0) = d_inlier(I_select_id(I_select_id_iter),0);
		theta_inlier_I(I_select_id_iter,0) = theta_inlier(I_select_id(I_select_id_iter),0);
	    }
            u_stack.col(0) = d_inlier_I.array()*theta_inlier_I.array().cos();
            u_stack.col(1) = d_inlier_I.array()*theta_inlier_I.array().sin();
            Eigen::MatrixXf u_vec = u_stack.colwise().mean().transpose(); //2x1 vector
	    Eigen::MatrixXf ones = Eigen::MatrixXf::Ones(1,t-begin_idx+1);
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
            Eigen::Matrix2f Cov_single_point;
            Cov_single_point = Cov_Laser;//based on the material of 5551
            float theta_est = atan2(coef_est(1),coef_est(0));
            float d_est = coef_est(2);
            Eigen::MatrixXf R_est = Ob_Line_Cov_cal(d_inlier_I,theta_inlier_I,d_est,theta_est,Cov_single_point);
            AddMatrixToRow(R_total,R_est);
            Eigen::Matrix<float,2,1> coef_est_2;
	    coef_est_2 << d_est,theta_est;
            AddMatrixToCol(NewLines,coef_est_2);
            }
        else{
            // otherwise, add all the points in this part into outliers.
	    Eigen::MatrixXi I_select_id = I.block(begin_idx,0,t - begin_idx + 1,1);
            AddMatrixToRow(EraseIdx,I_select_id); 

        }

    }
  }
}


void LineExtraction_RANSAC::line_para_cal(Eigen::MatrixXf d_est,Eigen::MatrixXf theta_est,float& d_H,float& theta_H){
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



Eigen::MatrixXf LineExtraction_RANSAC::Ob_Line_Cov_cal(Eigen::MatrixXf d_all, Eigen::MatrixXf theta_all,float d_star, float theta_star,Eigen::MatrixXf Cov_single_point){
/*
% This function is to calculate the new observed line covariance matrix
% Input: d_all, nx1 vector, polar distance to the points
%        theta_all, nx1 vector, polar angle of the points
%        d_star, scalar, optimal linear solution for the line segment
%        theta_star, scalar, optimal linear solution for the line segment
%        Cov_single_point, 2x2 matrix, covariance of the polar measurement
%        for one single point. Assume all points are same, and
%        uncorrelated.
% Output: R_est, 2x2 matrix, covariance of the estimated line parameter
%         d_star(the first one), theta_star
*/
    Eigen::MatrixXf x_all = d_all.array()*theta_all.array().cos();
    Eigen::MatrixXf y_all = d_all.array()*theta_all.array().sin();
    float x_mean = x_all.mean();
    float y_mean = y_all.mean();
    float S_x_sq = (x_all.array()-x_mean).array().pow(2).sum();
    float S_y_sq = (y_all.array()-y_mean).array().pow(2).sum();
    float S_x_y = ((x_all.array()-x_mean).array()*(y_all.array()-y_mean).array()).sum();
    int n = d_all.rows();

    Eigen::MatrixXf R_est = Eigen::MatrixXf::Zero(2,2);
    for(int i =0;i<n;i++){
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(2,2);
        A(1,0) = ((y_mean-y_all(i))*(S_y_sq-S_x_sq)+2*S_x_y*(x_mean-x_all(i)))/(pow((S_y_sq-S_x_sq),2)+4*pow(S_x_y,2));
        A(1,1) = ((x_mean-x_all(i))*(S_y_sq-S_x_sq)-2*S_x_y*(y_mean-y_all(i)))/(pow((S_y_sq-S_x_sq),2)+4*pow(S_x_y,2));
        A(0,0) = cos(theta_star)/n-x_mean*sin(theta_star)*A(1,0)+y_mean*cos(theta_star)*A(1,0);
        A(0,1) = sin(theta_star)/n-x_mean*sin(theta_star)*A(1,1)+y_mean*cos(theta_star)*A(1,1);
        Eigen::Matrix2f B;
        B << cos(theta_all(i)), -1*d_all(i)*sin(theta_all(i)),sin(theta_all(i)), d_all(i)*cos(theta_all(i));
        Eigen::MatrixXf J = A*B;
        R_est = R_est + J*Cov_single_point*J.transpose();
    }
    int aaaa  = 1;
    return R_est;
}

Eigen::MatrixXi LineExtraction_RANSAC::sort_matrix_int(Eigen::MatrixXi v){
    std::vector<int> v1;
    for(int i=0;i<v.rows();i++)
        v1.push_back(v(i,0));
    std::sort (v1.begin(), v1.begin()+v1.size());
    for(int i=0;i<v.rows();i++)
        v(i,0) = v1[i];
    return v;
}


Eigen::MatrixXi LineExtraction_RANSAC::sort_indexes(Eigen::MatrixXf v){
    std::vector<float> v1;
    std::map<float,int> map_v;
    Eigen::MatrixXi idx(v.rows(),1);
    for(int i=0;i<v.rows();i++){
        v1.push_back(v(i,0));
        map_v[v(i,0)] = i;
    }
    std::sort (v1.begin(), v1.begin()+v1.size());
    std::map<float,int>::iterator it;
    for(int i=0;i<v.rows();i++){
       it = map_v.find(v1[i]);
       if (it != map_v.end()){
	  //std::cout << v1[i] <<std::endl;
	  //std::cout << it->first <<std::endl;
	  //std::cout << it->second <<std::endl;
	  idx(i,0) = it->second;
      }
    }
   
    return idx;
}


void LineExtraction_RANSAC::removeRow(Eigen::MatrixXf& matrix, int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void LineExtraction_RANSAC::removeRow(Eigen::MatrixXi& matrix, int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void LineExtraction_RANSAC::removeColumn(Eigen::MatrixXf& matrix, int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void LineExtraction_RANSAC::AddMatrixToRow(Eigen::MatrixXf& M,Eigen::MatrixXf m){
    if(M.rows() == 0){
      M = m;
    }
    else{
      Eigen::MatrixXf f(M.rows()+m.rows(),M.cols());
      f.block(0,0,M.rows(),M.cols()) = M;
      f.block(M.rows(),0,m.rows(),m.cols()) = m;
      M = f;}
}

void LineExtraction_RANSAC::AddMatrixToRow(Eigen::MatrixXi& M,Eigen::MatrixXi m){
    if(M.rows() == 0){
      M = m;
    }
    else{   
    //M.block(M.rows(),0,m.rows(),m.cols()) = m;}
    Eigen::MatrixXi f(M.rows()+m.rows(),M.cols());
    f.block(0,0,M.rows(),M.cols()) = M;
    f.block(M.rows(),0,m.rows(),m.cols()) = m;
    M = f;}
}

void LineExtraction_RANSAC::AddMatrixToCol(Eigen::MatrixXf& M,Eigen::MatrixXf m){
    if(M.rows() == 0){
      M = m;
    }
    else{
    Eigen::MatrixXf f(M.rows(),M.cols()+m.cols());
    f.block(0,0,M.rows(),M.cols()) = M;
    f.block(0,M.cols(),m.rows(),m.cols()) = m;
    M = f;}
}
void LineExtraction_RANSAC::AddRow(Eigen::MatrixXi& M,int m){
    if(M.rows() == 0){
      M.resize(1,1);
      M(0,0) = m;}
    else{
      Eigen::MatrixXi f(M.rows()+1,1);
      f.block(0,0,M.rows(),1) = M;
      f(M.rows(),0) = m;
      M.resize(M.rows()+1,1);
      M = f;
    }
}

Eigen::MatrixXi LineExtraction_RANSAC::GenerateUniqueInt(int tam, int number){
        std::vector<int> generatedValues;
        Eigen::MatrixXi generatedValuesMatrix(number,1);
        for (int i=0;i<number;i++){
             int num = rand() % tam;
             while(!Contains(generatedValues, num)){ //You'll need a function to check whether the num is in this vector
                  num = rand() % tam;  
             }
             generatedValues.push_back(num);
             generatedValuesMatrix(i,0) = generatedValues[i];
        }        
       return generatedValuesMatrix;
   }
   
bool LineExtraction_RANSAC::Contains(std::vector<int> generatedValues,int num){
        int flag = 1;
	if(generatedValues.size() == 0)
	  return true;
        for(int i=0;i<generatedValues.size();i++)
	  if(generatedValues[i] == num)
	    flag = 0;
        if(flag == 0)
	  return false;
	else
	  return true;
   }
