function endpoints_in_state_new = update_mached_line_endpoint(state_all,endpoints_in_state,line_idx,endpoints_line_measure)

% This function is to update the matched lines' ending points in state
% vector. This should be done after measurement update.

% Input: 
%   state_all: (3+2N)x1 vector. State vector after measurement update
%   endpoints_in_state: 4xN matrix, each column is one line's ending points
%                       in the order x1,x2,y1,y2. This is from state vector's line
%   line_idx: mx1 vector, each element is the matched line in state vector's starting
%                 index.
%   endpoints_line_measure: 4Xm matrix, each column is mached measurement line's ending points 

% Output:
%   endpoints_in_state_new: 4xN matrix, each column is one line's updated
%                           ending points. The line is from state vector.

%% 1. Extrac line parameter from state vector
mached_num = length(line_idx);
line_state = zeros(2,mached_num);
%mached_end_in_state = zeros(2,mached_num);
endpoint_ref_1 = [];
endpoint_ref_2 = [];
for i = 1:mached_num
    idx = (line_idx(i)-4)/2+1;
    line_state(:,i) = state_all(line_idx(i):line_idx(i)+1,1);
    mached_end_in_state = endpoints_in_state(:,idx);
    endpoint_ref_1 = [endpoint_ref_1 [mached_end_in_state(1);mached_end_in_state(3)]];
    endpoint_ref_2 = [endpoint_ref_2 [mached_end_in_state(2);mached_end_in_state(4)]];
end

%% 2. Change the measured lines ending points to global coordinate
end1_test_gl = [];% ending point 1
end2_test_gl = [];% ending point 2
robot_vector = state_all(1:3,1);
for i = 1:mached_num
    x_part = endpoints_line_measure(1:2,i);
    y_part = endpoints_line_measure(3:4,i);    
    
    [end_point_1,end_point_2,~] = end_point_2_global(robot_vector,x_part,y_part);
    end1_test_gl = [end1_test_gl,end_point_1];
    end2_test_gl = [end2_test_gl,end_point_2];
end
%% 2. Project all the ending points to the current lines due to the measurement update

[state_end_1,state_end_2,measure_end_1,measure_end_2] = p_2_line(line_state,endpoint_ref_1,endpoint_ref_2,end1_test_gl,end2_test_gl);

%% 3. Update the ending points
endpoints_in_state_new = endpoints_in_state;
for i = 1:mached_num
    
    % get the line index that needs to update the end points
    idx = (line_idx(i)-4)/2+1;
    % reference line already in the state
    line_single = line_state(:,i);
    theta_ref = line_single(2);
    d_ref = line_single(1);

    % compute the relative angle w.r.t. the line parameter theta

    end1_test_theta = atan2(measure_end_1(2,i),measure_end_1(1,i))-theta_ref;
    end1_test_theta = atan2(sin(end1_test_theta),cos(end1_test_theta)); 
    
    end2_test_theta = atan2(measure_end_2(2,i),measure_end_2(1,i))-theta_ref;
    end2_test_theta = atan2(sin(end2_test_theta),cos(end2_test_theta));
    
    end1_ref_theta = atan2(state_end_1(2,i),state_end_1(1,i))-theta_ref;
    end1_ref_theta = atan2(sin(end1_ref_theta),cos(end1_ref_theta));
    
    end2_ref_theta = atan2(state_end_2(2,i),state_end_2(1,i))-theta_ref;
    end2_ref_theta = atan2(sin(end2_ref_theta),cos(end2_ref_theta));
    
    if end1_test_theta > end2_test_theta
        temp = end1_test_theta;
        end1_test_theta = end2_test_theta;
        end2_test_theta = temp;
        temp = measure_end_1(:,i);
        measure_end_1(:,i) = measure_end_2(:,i);
        measure_end_2(:,i) = temp;
    end
    
    if end1_ref_theta > end2_ref_theta
        temp = end1_ref_theta;
        end1_ref_theta = end2_ref_theta;
        end2_ref_theta = temp;
        temp = state_end_1(:,i);
        state_end_1(:,i) = state_end_2(:,i);
        state_end_2(:,i) = temp;
        
    end
    
    if end1_test_theta >= end1_ref_theta && end2_test_theta<=end2_ref_theta % test line segment within the reference segment
        %didn't update the ending points
        aa =1;
    elseif end1_test_theta>=end1_ref_theta && end1_test_theta<=end2_ref_theta % test seg starting point within reference seg, while end not
        state_end_2(:,i) = measure_end_2(:,i);
    elseif end2_test_theta>=end1_ref_theta && end2_test_theta<=end2_ref_theta % test seg end point within reference seg, while start not
        state_end_1(:,i) = measure_end_1(:,i);
    elseif end1_ref_theta>= end1_test_theta && end2_ref_theta<= end2_test_theta && end2_ref_theta >= end1_test_theta % test seg encompass the reference seg
        state_end_1(:,i) = measure_end_1(:,i);
        state_end_2(:,i) = measure_end_2(:,i);
    end
    endpoints_in_state_new(:,idx) = [state_end_1(1,i),state_end_2(1,i),state_end_1(2,i),state_end_2(2,i)];

end