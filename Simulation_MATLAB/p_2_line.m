function [state_end_1,state_end_2,measure_end_1,measure_end_2] = p_2_line(line,endpoint_ref_1,endpoint_ref_2,end1_test_gl,end2_test_gl)

% This function is to calculate the projected points in line for each
% point

% Input:
%       line:2xm line parameter matrix, each column is [d,theta]'
%       endpoint_ref_1: 2xm ending points matrix for the state line
%       endpoint_ref_2: another 2xm ending points matrix for the state line
%       end1_test_gl: 2xm ending points matrix for the measurement line
%       end2_test_gl: another 2xm ending points matrix for measurement line

% Output:
%       updated ending points that are projected to the lines
state_end_1 = [];
state_end_2 = [];
measure_end_1 = [];
measure_end_2 = [];
for i = 1:size(line,2)
    theta_ref = line(2,i);
    d_ref = line(1,i);
    % projected points calculation
    s_1 = zeros(2,1);
    s_1(1,1) = sin(theta_ref)*(sin(theta_ref)*endpoint_ref_1(1,i)-cos(theta_ref)*endpoint_ref_1(2,1))...
            -cos(theta_ref)*(-1*d_ref);
    s_1(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_ref_1(1,1)+cos(theta_ref)*endpoint_ref_1(2,1))...
            -sin(theta_ref)*(-1*d_ref);  
    state_end_1 = [state_end_1,s_1];
    
    s_2 = zeros(2,1);
    s_2(1,1) = sin(theta_ref)*(sin(theta_ref)*endpoint_ref_2(1,i)-cos(theta_ref)*endpoint_ref_2(2,1))...
            -cos(theta_ref)*(-1*d_ref);
    s_2(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*endpoint_ref_2(1,1)+cos(theta_ref)*endpoint_ref_2(2,1))...
            -sin(theta_ref)*(-1*d_ref);  
    state_end_2 = [state_end_2,s_2];
    
    m_1 = zeros(2,1);
    m_1(1,1) = sin(theta_ref)*(sin(theta_ref)*end1_test_gl(1,i)-cos(theta_ref)*end1_test_gl(2,1))...
            -cos(theta_ref)*(-1*d_ref);
    m_1(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*end1_test_gl(1,1)+cos(theta_ref)*end1_test_gl(2,1))...
            -sin(theta_ref)*(-1*d_ref);  
    measure_end_1 = [measure_end_1,m_1];    

    m_2 = zeros(2,1);
    m_2(1,1) = sin(theta_ref)*(sin(theta_ref)*end2_test_gl(1,i)-cos(theta_ref)*end2_test_gl(2,1))...
            -cos(theta_ref)*(-1*d_ref);
    m_2(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*end2_test_gl(1,1)+cos(theta_ref)*end2_test_gl(2,1))...
            -sin(theta_ref)*(-1*d_ref);  
    measure_end_2 = [measure_end_2,m_2];   
    
end
    
    
