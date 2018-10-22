function [d,theta] = gen_sim_laser_data(p_r,heading_r,sigma_d)
% generate distance d corresponding to the angle theta
% input: pose (position and heading) of the robot and noise in distance
% note p_r should line in the hall way

% noise in distance measurement
%sigma_d = 0.0;

x_shift = 1;
% lines equations
x1 = 0+x_shift;
width = 3;
x2 = x1+width;
% y1 is actually hidden
y1 = 12;
y2 = y1+width;

% visualize
figure();
line([x1,x1],[0,12],'Color','black');
hold on;
line([x2,x2],[0,y2],'Color','black');
hold on;
line([x2,-width],[y2,y2],'Color','black');

% robot parameters
max_range = 10;

% robot initial pose
%p_r = [x1+width/2;width];
%heading_r = pi/6;
%heading_r = pi/2;

alpha = atan((y2-p_r(2,1))/(x2-p_r(1,1)));
beta = atan((y1-p_r(2,1))/(p_r(1,1)-x1));

d = [];

for i = 0:1:180
    
    theta = (i-90)/180*pi+heading_r;
    
    if theta < alpha
        
        d1 = (x2-p_r(1,1))/cos(theta);
        d1 = d1+randn(1)*sigma_d;
        % comment this to abandon the line > max_range
%         if d1>max_range
%             d1 = max_range;
%         end
        line([p_r(1,1),p_r(1,1)+d1*cos(theta)],[p_r(2,1),p_r(2,1)+d1*sin(theta)],'Color','red');
        hold on;
        d = [d;d1];
        
    elseif theta >= alpha && theta < (pi-beta)
        
        d2 = (y2-p_r(2,1))/sin(theta);
        d2 = d2+randn(1)*sigma_d;
%         if d2>max_range
%             d2 = max_range;
%         end
        line([p_r(1,1),p_r(1,1)+d2*cos(theta)],[p_r(2,1),p_r(2,1)+d2*sin(theta)],'Color','green');
        hold on;
        d = [d;d2];
        
    elseif theta >= (pi-beta)
        
        d3 = (p_r(1,1)-x1)/cos(pi-theta);
        d3 = d3+randn(1)*sigma_d;
%         if d3>max_range
%             d3 = max_range;
%         end
        line([p_r(1,1),p_r(1,1)+d3*cos(theta)],[p_r(2,1),p_r(2,1)+d3*sin(theta)],'Color','blue');
        hold on;
        d = [d;d3];
        
    else
        
        disp('wrong, check the code');
        
    end
    
end

theta = (0:1:180)/180*pi;