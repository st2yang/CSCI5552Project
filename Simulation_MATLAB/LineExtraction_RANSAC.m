function [line,endpoints,R_mat] = LineExtraction_RANSAC(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh)
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
% line: N*3 matrix, each row represents 3 parameters (cos(theta),sin(theta),d) of a line in
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
while (iter < maxIter) && (length(d_est) > min_rest)
    % 1. randomly select two points
    idx = randperm(size(d_est,1),1);
    [idx1_2, dist1_2] = knnsearch([d_est theta_est],[d_est(idx) theta_est(idx)], 'k', 5,'Distance',@DISTFUN);
    idx = [idx;idx1_2(2)];
    
    x_est = d_est(idx).*cos(theta_est(idx));
    y_est = d_est(idx).*sin(theta_est(idx));
    
    % 2. compute line equation
    theta_H = atan2(-(x_est(2) - x_est(1)),y_est(2) - y_est(1));
    d_H = x_est(1)*cos(theta_H) + y_est(1)*sin(theta_H);
%     abs(cos(theta_H - theta_est(idx(1))) - d_H/d_est(idx(1)))
%     abs(cos(theta_H - theta_est(idx(2))) - d_H/d_est(idx(2)))

    % 3. compute distance from other points to this line
    inlier = [];
    inlier_num = 0;
    dist_total = zeros(size(d_est));
    for i = 1:length(d_est)
        dist = abs(cos(theta_H - theta_est(i)) - d_H/d_est(i));
        dist_total(i) = dist;
        if dist < thresh
            inlier = [inlier;i];
            inlier_num = inlier_num + 1;
        end
    end
    
    % 4. check if the the number of points on this line satisfies minimum
    % number of points. If so, add this line and remove all inliers.
    if inlier_num >= min_num
        % find corret segement and end points on the line (by theta)
        %[NewLines,NewEndpoints,EraseIdx] = LineSegment(d_est(inlier),theta_est(inlier),theta_thresh,min_num);
        %seg_dist_thresh = 10.;
        [NewLines,NewEndpoints,New_R,EraseIdx] = LineSegment_Diff_method(d_est(inlier),theta_est(inlier),seg_dist_thresh,min_num,theta_thresh);
        inlier(EraseIdx) = [];
        line = [line;NewLines];
        endpoints = [endpoints NewEndpoints];
        
        R_mat = [R_mat;New_R];
        % clear inliers
        d_est(inlier) = [];
        theta_est(inlier) = [];
    end
    iter = iter + 1;
end
end

function D2 = DISTFUN(ZI, ZJ)

D2 = zeros(size(ZJ,1),1);
for i = 1:size(ZJ,1)
    D2(i) = abs(cos(ZI(2) - ZJ(i,2)) - ZI(1)/ZJ(i,1));
end
end