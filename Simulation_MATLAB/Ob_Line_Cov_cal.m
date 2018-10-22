
function R_est = Ob_Line_Cov_cal(d_all, theta_all,d_star,theta_star,Cov_single_point)

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

x_all = d_all.*cos(theta_all);
y_all = d_all.*sin(theta_all);
x_mean = mean(x_all);
y_mean = mean(y_all);
S_x_sq = sum((x_all-x_mean).^2);
S_y_sq = sum((y_all-y_mean).^2);
S_x_y = sum((x_all-x_mean).*(y_all-y_mean));
n = length(d_all);

R_est = zeros(2);
for i = 1:n
    A = zeros(2);
    %B = zeros(2);
    A(2,1) = ((y_mean-y_all(i))*(S_y_sq-S_x_sq)+2*S_x_y*(x_mean-x_all(i)))/((S_y_sq-S_x_sq)^2+4*S_x_y^2);
    A(2,2) = ((x_mean-x_all(i))*(S_y_sq-S_x_sq)-2*S_x_y*(y_mean-y_all(i)))/((S_y_sq-S_x_sq)^2+4*S_x_y^2);
    A(1,1) = cos(theta_star)/n-x_mean*sin(theta_star)*A(2,1)+y_mean*cos(theta_star)*A(2,1);
    A(1,2) = sin(theta_star)/n-x_mean*sin(theta_star)*A(2,2)+y_mean*cos(theta_star)*A(2,2);
    B = [cos(theta_all(i)), -1*d_all(i)*sin(theta_all(i));sin(theta_all(i)), d_all(i)*cos(theta_all(i))];
    J = A*B;
    R_est = R_est + J*Cov_single_point*J';
end
aaaa  = 1;
