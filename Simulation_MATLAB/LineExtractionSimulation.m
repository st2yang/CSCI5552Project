clear all;
close all
clc;

%% create 2D environment
NumOfLines = 10;
NumOfPoint = 10;

% 1. create lines and points
% coef = rand(3,NumOfLines);
% x = randi([1 100],NumOfPoint,NumOfLines);
% y = zeros(size(x));
% for i = 1:NumOfLines
%     y(:,i) = (-coef(3,i) - coef(1,i)*x(:,i))/coef(2,i);
% end
x_num = randi(NumOfLines,1,1); y_num = NumOfLines - x_num;
x = zeros(NumOfLines*NumOfPoint,1);
y = zeros(NumOfLines*NumOfPoint,1);
x(1:x_num*NumOfPoint) = reshape(repmat(10*randperm(10,x_num),NumOfPoint,1),x_num*NumOfPoint,1);
for i = 1:x_num
    y(NumOfPoint*(i-1)+1:NumOfPoint*i) = 10*randn(NumOfPoint,1);
end
y(x_num*NumOfPoint+1:end) = reshape(repmat(10*randperm(10,y_num),NumOfPoint,1),y_num*NumOfPoint,1);
for i = 1:y_num
    x(x_num*NumOfPoint+NumOfPoint*(i-1)+1:x_num*NumOfPoint+NumOfPoint*i) = 10*randn(NumOfPoint,1);
end


plot(x,y,'r*','MarkerSize',10);hold on


% 2. translate points into polar coordinate
d = (x.^2 + y.^2).^0.5;
theta = atan2(y,x);
% d_true = abs(coef(3,:))./(coef(1,:).^2 + coef(2,:).^2).^0.5;
% theta_true = atan2(coef(2,:),coef(1,:));


% 3. add noise 
sigma_d = 1; sigma_t = 0.1;
d_est = reshape(d + rand()*sigma_d^2,[size(d,1)*size(d,2),1]);
theta_est = reshape(theta + rand()*sigma_t^2,[size(d,1)*size(d,2),1]);

% dist = (coef'*[d_est.*cos(theta_est) d_est.*sin(theta_est) ones(length(d_est),1)]')./((coef(1,:).^2 + coef(2,:).^2).^0.5 .* coef(3,:))';
%% line extraction
maxIter = 1000;
min_num = 6;
thresh = 0.01;
min_rest = 6;
theta_thresh = 3;
seg_dist_thresh = 10;
% 1. RANSAC
[line_RANSAC,endpoints] = LineExtraction_RANSAC(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh);

% 2. Incremental Algorithm
% L = 2;
% thresh = 0.1;
% line_Incremental = LineExtraction_Incremental(d_est,theta_est,L,thresh,min_num,maxIter,min_rest);

% check results
% abs(coef'*[d_H*cos(theta_H);d_H*sin(theta_H);1])
% abs(coef(2)/coef(1) - tan(theta_H))


%% display
x_est = d_est.*cos(theta_est);
y_est = d_est.*sin(theta_est);
% plot(x,y,'b-');hold on
plot(x_est,y_est,'go','MarkerSize',10);hold on 
x_draw = x(:,1);

% Since all the lines are either horizontal or vertical in this simulation, we temporarily
% draw them in this way. In practice, we need to draw them in standered
% way.
for i = 1:size(line_RANSAC,1)
%     y_draw = (line(i,1) - x_draw *cos(line(i,2)))/sin(line(i,2));
    if abs(line_RANSAC(i,1)) < 1e-2
        x_draw = endpoints(3,i)*cos(endpoints(4,i)):1:endpoints(1,i)*cos(endpoints(2,i));
        y_draw = repmat(line_RANSAC(i,3)/line_RANSAC(i,2),1,length(x_draw));
    elseif abs(line_RANSAC(i,2)) < 1e-2
        y_draw = endpoints(1,i)*sin(endpoints(2,i)):1:endpoints(3,i)*sin(endpoints(4,i));
        x_draw = repmat(line_RANSAC(i,3)/line_RANSAC(i,1),1,length(y_draw));
    end
    plot(x_draw,y_draw,'r-');hold on
end


% for i = 1:size(line_Incremental,1)
% %     y_draw = (line(i,1) - x_draw *cos(line(i,2)))/sin(line(i,2));
%     if abs(line_Incremental(i,1)) < 1e-3
%         x_draw = -10:1:10;
%         y_draw = repmat(-1 /line_Incremental(i,2),1,length(x_draw));
%     elseif abs(line_Incremental(i,2)) < 1e-3
%         y_draw = -10:1:10;
%         x_draw = repmat(-1 /line_Incremental(i,1),1,length(y_draw));
%     end
%     plot(x_draw,y_draw,'r-');hold on
% end


