
function [NewLines,NewEndpoints,R_total,EraseIdx] = LineSegment_Diff_method(d_inlier,theta_inlier,seg_dist_thresh,min_num,Cov_Laser)
% This code is for line segmentation
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


% 1. sort theta from small to large
[~,I] = sort(theta_inlier);
x_pos = d_inlier.*cos(theta_inlier);
y_pos = d_inlier.*sin(theta_inlier);
mean_theta = abs(mean(theta_inlier));
begin_idx = 1;
%end_idx = 1;
NewLines = [];
R_total = [];
NewEndpoints = [];
EraseIdx = []; % Erase all inlier points
for t = 1:length(d_inlier)
    % 2. compute difference of theta between two neighbors. 
    % If difference is larger than threshold: seg_dist_thresh,
    % segment the line, and treat the front part as a temporary new line
    if t<length(d_inlier)
        
        dist_neighbor = sqrt((x_pos(I(t))-x_pos(I(t+1)))^2+(y_pos(I(t))-y_pos(I(t+1)))^2);
        
        %dist_theta = abs(theta_inlier(I(t)) - theta_inlier(I(t+1)));
        %if dist_theta > theta_thresh*mean_theta
        if (dist_neighbor > seg_dist_thresh)
        
        % 3. If the front part satisfied minimum-points-condition, add    
        % two end points and recompute the line parameters
        % try to include another criteria: length of the line segment
            
            %if ((t - begin_idx + 1) >= min_num) && sqrt((x_pos(I(begin_idx))-x_pos(I(t)))^2+(y_pos(I(begin_idx))-y_pos(I(t)))^2)>300
            if ((t - begin_idx + 1) >= min_num)
                NewEndpoints = [NewEndpoints d_inlier(I(begin_idx)) theta_inlier(I(begin_idx)) d_inlier(I(t)) theta_inlier(I(t))];

                % recompute the line parameter by Linear Least Square
                u_stack = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))), d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t)))];
                u_vec = mean(u_stack,1)'; %2x1 vector
                dif_u_p = u_vec*ones(1,t-begin_idx+1)-u_stack';
                A = dif_u_p*dif_u_p';
                [V,D] = eig(A);
                [~,idx_min] = min(diag(D));
                V_min = V(:,idx_min);
                cos_theta= V_min(1);
                sin_theta= V_min(2);
                d_theta = u_vec'*V_min;
                %A = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))) d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t))) -ones(t - begin_idx + 1,1)];
                %[~,~,V] = svd(A);
                coef_est = [cos_theta,sin_theta,d_theta]';
                coef_est = sign(d_theta)*coef_est;
                
                % 4. compute the optimal solution theta_est, d_est in order
                % to calculate the covariance of the estimated theta_est,
                % d_est
                theta_est = atan2(coef_est(2),coef_est(1));
                d_est = coef_est(3);
                % setup the sigma_d and sigma_theta
                % according to the material of 5551, distance accuracy is
                % +-15mm, angular resolution 0.5deg
                %Cov_single_point = [25,0;0,1e-4];%based on the material of 5551
                Cov_single_point = Cov_Laser;
                R_est = Ob_Line_Cov_cal(d_inlier(I(begin_idx:t)),theta_inlier(I(begin_idx:t)),d_est,theta_est,Cov_single_point);
                R_total = [R_total;R_est];
                NewLines = [NewLines [d_est,theta_est]'];
            else
                % otherwise, add all the points in this part into outliers.
                EraseIdx = [EraseIdx;I(begin_idx:t)]; 
            end
        % change the begin point index
            begin_idx = t+1;
        end
    else
        
        % 3. If the front part satisfied minimum-points-condition, add
        % two end points and recompute the line parameters
        if ((t - begin_idx + 1) >= min_num) 
            NewEndpoints = [NewEndpoints d_inlier(I(begin_idx)) theta_inlier(I(begin_idx)) d_inlier(I(t)) theta_inlier(I(t))];

            % recompute the line parameter by Linear Least Square
            u_stack = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))), d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t)))];
            u_vec = mean(u_stack,1)'; %2x1 vector
            dif_u_p = u_vec*ones(1,t-begin_idx+1)-u_stack';
            A = dif_u_p*dif_u_p';
            [V,D] = eig(A);
            [~,idx_min] = min(diag(D));
            V_min = V(:,idx_min);
            cos_theta= V_min(1);
            sin_theta= V_min(2);
            d_theta = u_vec'*V_min;
            %A = [d_inlier(I(begin_idx:t)).*cos(theta_inlier(I(begin_idx:t))) d_inlier(I(begin_idx:t)).*sin(theta_inlier(I(begin_idx:t))) -ones(t - begin_idx + 1,1)];
            %[~,~,V] = svd(A);
            coef_est = [cos_theta,sin_theta,d_theta]';
            coef_est = sign(d_theta)*coef_est;
            %coef_est = V(:,end);
            %coef_est = coef_est/(coef_est(1)^2 + coef_est(2)^2)^0.5;
            %Cov_single_point = [60,0;0,1e-4];%based on the material of 5551
            Cov_single_point = Cov_Laser;
            theta_est = atan2(coef_est(2),coef_est(1));
            d_est = coef_est(3);
            R_est = Ob_Line_Cov_cal(d_inlier(I(begin_idx:t)),theta_inlier(I(begin_idx:t)),d_est,theta_est,Cov_single_point);
            R_total = [R_total;R_est];
            NewLines = [NewLines [d_est, theta_est]'];
        else
            % otherwise, add all the points in this part into outliers.
            EraseIdx = [EraseIdx;I(begin_idx:t)]; 

        end

    end
end
NewEndpoints = reshape(NewEndpoints,4,length(NewEndpoints)/4);   
  
        
        
        
        
        
        
