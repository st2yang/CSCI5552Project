function [line,endpoints,R_mat] = LineExtraction_RANSAC_2(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh,Cov_Laser)
% This code is for line extraction by RANSAC
%
% Input:
% d_est & theta_est: measurement value from 2D laser scanner
% thresh: maxmum distance from points to lines
% min_mum: minimum number of points on a line
% maxIter: maximum number of iterations
% min_rest: minimum number of points left after iterations
%
% Output: 
% line: 2xN matrix, each column represents 2 parameters (d,theta) of a line in
% 2D space (polar coordinate). 
% endpoints: 4*N matrix, each column represents two end points for a
% line(theta from small to large).
% R_mat: 2N*2 matrix, each 2x2 matrix, represent the corresponding line
% parameters' covariance matrix, d_star(the first one), theta_star


line = [];
% choice_inlier = [];
iter = 0;
endpoints = [];
R_mat = [];
outIter = 5;
%d_est_group = []
while (iter < outIter) && (length(d_est) > min_rest)
    
    max_inlier_num = 0;
    for inner = 1:maxIter
    %while (inner < maxIter) && (length)
        
        % 1. randomly select two points, not just two, a small set
        idx = randperm(size(d_est,1),min_num);

        
        % 2.compute the line parameter by Linear Least Square
        [d_H, theta_H] = line_para_cal(d_est(idx),theta_est(idx));
        
        
        % 3. compute distance from other points to this line
        inlier = [];
        inlier_num = 0;
        d_test = d_est; % use only un-selected points
        d_test(idx) = [];
        theta_test = theta_est;
        theta_test(idx) = [];
        
        for i = 1:length(d_test)
            %dist = abs(cos(theta_H - theta_est(i)) - d_H/d_est(i));
            dist = abs(d_test(i)*cos(theta_H - theta_test(i)) - d_H);
            %dist_total(i) = dist;
            if dist < thresh
                inlier = [inlier;i];
                inlier_num = inlier_num + 1;
            end
        end
        
        % 4. find the best line for the current remaining data points
        if inlier_num >= min_num && inlier_num >=max_inlier_num
            max_inlier_num = inlier_num;
            max_inlier_idx = [inlier;idx'];
        end
        
    end
    
    if max_inlier_num >0
        [NewLines,NewEndpoints,New_R,EraseIdx] = LineSegment_Diff_method(d_est(max_inlier_idx),theta_est(max_inlier_idx),seg_dist_thresh,min_num,Cov_Laser);
        max_inlier_idx(EraseIdx) = [];
        line = [line NewLines];
        endpoints = [endpoints NewEndpoints];
        R_mat = [R_mat;New_R];
        % clear inliers
        d_est(max_inlier_idx) = [];
        theta_est(max_inlier_idx) = [];
    end
    iter = iter + 1;
end

    
    
    
        