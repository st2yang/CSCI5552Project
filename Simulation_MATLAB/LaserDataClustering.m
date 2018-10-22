

% This is written by Jianjun Yuan
clear all;
close all;

%fileID = fopen('LaserData_forward.txt','r');
%fileID = fopen('LaserData_left.txt','r');
fileID = fopen('LaserData_right.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
dataNum = 361;
A = reshape(A,2,dataNum);
A = A';

filterOutLength = 3000;
filterOutAxis = 2000;
j=0;
for i = 1:length(A)
    temp = A(i,:);
    
    if norm(temp)<filterOutLength && temp(1,1)<filterOutAxis &&temp(1,2)<filterOutAxis
        j = j+1;
        B(j,:) = temp;
    end
end

figure;

scatter(B(:,1),B(:,2));
xlabel('X position(mm)');
ylabel('Y position(mm)');
title('Static Data Samples Distribution');
clusterNum = 15;
[IDX, C] = kmeans(B, clusterNum);
hold on;
scatter(C(:,1),C(:,2),'M','filled');
hold off;
[IDX, D] = kmeans(C, 6);
hold on;
scatter(D(:,1),D(:,2),'T','filled');
hold off;
aaa= 1;
