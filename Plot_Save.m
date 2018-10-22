
%% This is used to plot the current progress of the SLAM

% Note:::
%         The unit is m instead of mm

function [] = Plot_Save()

clear all
close all

savevideo=1;

cleanupObj = onCleanup(@() cleanMeUp(savevideo)); %handle cntrl-c event

Odomrun=[];
Odomrunbytesize=0;
Odomrunlength=0;

Scan=[];
NewScanPlot=plot(0,0,'MarkerSize', 0.00001);

Scanbytesize=0;
Scanlength=0;

Features=[];
Featuresbytesize=0;
Featureslength=0;

NewEndPlot=plot(0,0,'MarkerSize', 0.00001);
End_point_bytesize = 0;
End_point_length = 0;

%knownFeatures=[];
%knownFeaturesbytesize=0;
%knownFeatureslength=0;
%knownFeaturesPlot=plot(0,0,'MarkerSize', 0.00001);


Pxy=[];
Pxybytesize=0;
Pxylength=0;
Ellipseplot=plot(0,0,'MarkerSize', 0.00001);
currXY=[0;0];

fighandle=figure(1);
hold on;
dt=.1;

if savevideo
V=VideoWriter('video.avi');
V.FrameRate=1/dt;
open(V);
fighandle.Resize='off';
%fighandle.Position=[0 0 1920 1080];
end


while(1)
    %% Robot x y plot
    Odomrunfinfo=dir('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotpose.txt');
    %Odomrunfinfo=dir('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/robotpose.txt');
    if size(Odomrunfinfo,1)>0&&~isempty(Odomrunfinfo.bytes)
    if Odomrunfinfo.bytes>Odomrunbytesize
        newA=dlmread('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotpose.txt','',Odomrunlength,0);
        %newA=dlmread('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/robotpose.txt','',Odomrunlength,0);
        Odomrunbytesize=Odomrunfinfo.bytes;
        plot(newA(:,1)/1000.,newA(:,2)/1000.,'b*');
        %plot(newA(:,1),newA(:,2),'b*');
        
        Odomrunlength=Odomrunlength+size(newA,1);
        currXY=[newA(end,1)/1000;newA(end,2)/1000];
        %currXY=[newA(end,1);newA(end,2)];
    end
    end

    Pxyinfo=dir('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotcov.txt');
    %Pxyinfo=dir('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/robotcov.txt');
    if  size(Pxyinfo,1)>0&&~isempty(Pxyinfo.bytes)
    if Pxyinfo.bytes>Pxybytesize
        NewPxy=dlmread('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotcov.txt','',Pxylength,0);
        %NewPxy=dlmread('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/robotcov.txt','',Pxylength,0);
        Pxy=[NewPxy(end,1)/1e6,NewPxy(end,2)/1e6;NewPxy(end,3)/1e6,NewPxy(end,4)/1e6];
        %Pxy=[NewPxy(end,1),NewPxy(end,2);NewPxy(end,3),NewPxy(end,4)];
        
        Pxybytesize=Pxyinfo.bytes;
        

        Pxylength=Pxylength+size(NewPxy,1);
        delete( Ellipseplot)
       [EllipseX,EllipseY]=  plot_error_ellipse_plotting(currXY,Pxy);
       Ellipseplot=plot(EllipseX,EllipseY,'g');
    end
    end
    
    %% Laser scan plot
%      Scanifo=dir('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/scan_update.txt');
%     if size(Scaninfo,1)>0&&~isempty(Scaninfo.bytes)
%      if Scaninfo.bytes>Scanbytesize
%         newScan=dlmread('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/scan_update.txt','',Scanlength,0);
%         Scanbytesize=Scaninfo.bytes;
%         delete(NewScanPlot);
%         plot(newScan(:,1),newScan(:,2),'k*');% old scans are black, current one is red
%         NewScanPlot=plot(newScan(:,1),newScan(:,2),'r*');
%         %Scan=[Scan;newScan];
%         %Scanlength=size(Scanlength,1);
%         Scanlength=Scanlength+size(newScan,1);
%     end
%     end
    
    %% measured lines plot
    Featuresinfo=dir('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints.txt');
    %Featuresinfo=dir('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/endpoints.txt');
    if size(Featuresinfo,1)>0&&~isempty(Featuresinfo.bytes)
    if Featuresinfo.bytes>Featuresbytesize
        newFeature=dlmread('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints.txt','',Featureslength,0);
        %newFeature=dlmread('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/endpoints.txt','',Featureslength,0);
        Featuresbytesize=Featuresinfo.bytes;
        for i = 1:size(newFeature,1)
            %end_point_gl = endpoints_in_state(:,i);
            end_point_gl = [newFeature(i,1),newFeature(i,2),newFeature(i,3),newFeature(i,4)];
            plot([end_point_gl(1),end_point_gl(2)]'/1000,[end_point_gl(3),end_point_gl(4)]'/1000,'g-');
            %plot(newFeature(:,1),newFeature(:,2),'g*');
        end
        %Features=[Features;newFeature];
        %Featureslength=size(Featureslength,1);
        Featureslength=Featureslength+size(newFeature,1);
    end
    end
    
    %% newest end-points plot
    update_feature=dir('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints_update.txt'); 
    %update_feature=dir('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/endpoints_update.txt'); 
    if size(update_feature,1)>0&&~isempty(update_feature.bytes)
        if update_feature.bytes > End_point_bytesize
            
            update_end_all=dlmread('/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints_update.txt','',End_point_length,0);
            %update_end=dlmread('/Users/yuanjianjun/Documents/Courses_MN/CSCI5552EstimationinRobotics/5552PioneerProjectRepo/data/endpoints_update.txt','',0,0);
            update_end = update_end_all(end,:);
            line_num = length(update_end)/4;
            update_reshape = reshape(update_end,4,line_num);
            delete(NewEndPlot);
            NewEndPlot = plot(update_reshape(1:2,:)/1000.,update_reshape(3:4,:)/1000.,'r-.','LineWidth',3);
            End_point_bytesize=update_feature.bytes;
            End_point_length=End_point_length+size(update_end_all,1);
            
        end
    end
        
  
    
%      knownFeaturesinfo=dir('../data/features/knownfeaturesRun.txt');
%     if size(knownFeaturesinfo,1)>0&& ~isempty(knownFeaturesinfo.bytes)
%     if knownFeaturesinfo.bytes>knownFeaturesbytesize
%         knownnewFeature=dlmread('../data/features/knownfeaturesRun.txt','',knownFeatureslength,0);
%         knownFeaturesbytesize=knownFeaturesinfo.bytes;
%         delete(knownFeaturesPlot);
%         knownFeaturesPlot=plot(knownnewFeature(:,1),knownnewFeature(:,2),'g*');
%         %Features=[Features;newFeature];
%         %Featureslength=size(Featureslength,1);
%         knownFeatureslength=knownFeatureslength+size(knownnewFeature,1);
%     end
%     end   
    
    pause(dt)
    if savevideo
    currFrame = getframe;
    writeVideo(V,currFrame)
    end
end

function cleanMeUp(f)
    % saves data to file (or could save to workspace)
    fprintf('saving stuff to file\n');
    if f
        close(V)
    end
end

end
