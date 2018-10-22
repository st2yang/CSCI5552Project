function [NewLines,NewEndpoints,EraseIdx] = LineSegment(d_inlier,theta_inlier,theta_thresh,min_num)
% This code is for line segmentation
%
% Input:
% d_inlier & theta_inlier: measurement value from 2D laser scanner inlier
% theta_thresh: ratio of maxmum distance between points (by theta) divided
% by average distance.
% min_mum: minimum number of points on a line
%
% Output: 
% NewLines: N*3 matrix, each row represents 3 parameters (cos(theta),sin(theta),d) of a line in
% 2D space (polar coordinate). 
% NewEndpoints: 4*N matrix, each column represents two end points for a
% line(theta from small to large).
% EraseIdx: outlier index.

% 1. sort theta from small to large
[~,I] = sort(theta_inlier);
mean_theta = abs(mean(theta_inlier));
begin_idx = 1;
NewLines = [];
NewEndpoints = [];
EraseIdx = [];
for t = 1:length(d_inlier)-1
        % 2. compute difference of theta between two neighbors. If difference is larger than threshold, which is determined by
        % certain times of average difference, segment this line and treat
        % the front part as a temporary new line
        dist_theta = abs(theta_inlier(I(t)) - theta_inlier(I(t+1)));
        if dist_theta > theta_thresh*mean_theta
            
           % 3. If the front part satisfied minimum-points-condition, add
           % two end points and recompute the line parameters
           if (t - begin_idx + 1) >= min_num  
              NewEndpoints = [NewEndpoints d_inlier(I(begin_idx)) theta_inlier(I(begin_idx)) d_inlier(I(t)) theta_inlier(I(t))];
    
              % recompute the line parameter by Linear Least Square
              A = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))) d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t))) -ones(t - begin_idx + 1,1)];
              [~,~,V] = svd(A);
              coef_est = V(:,end);
              coef_est = coef_est/(coef_est(1)^2 + coef_est(2)^2)^0.5;
              NewLines = [NewLines;coef_est'];
              
              % change the begin point index
              begin_idx = t+1;
           else
               
              % otherwise, add all the points in this part into outliers.
              EraseIdx = [EraseIdx;I(begin_idx:t)]; 
           end
        end
end

% Deal with the rest points (process is similar to previous)
if (length(d_inlier) - begin_idx) >= min_num  
   NewEndpoints = [NewEndpoints d_inlier(I(begin_idx)) theta_inlier(I(begin_idx)) d_inlier(I(t)) theta_inlier(I(t))];              
   
   % recompute the line parameter by Linear Least Square
   A = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))) d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t))) -ones(t - begin_idx + 1,1)];
   [~,~,V] = svd(A);
   coef_est = V(:,end);
   coef_est = coef_est/(coef_est(1)^2 + coef_est(2)^2)^0.5;
   NewLines = [NewLines;coef_est'];
end

NewEndpoints = reshape(NewEndpoints,4,length(NewEndpoints)/4);