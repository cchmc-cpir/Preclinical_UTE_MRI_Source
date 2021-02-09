function [EchoTraj,EchoFID] = Separate_Echoes(traj,FID,NEcho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A function to take data acquired in a multi-echo UTE sequence and separate
%the k-space data into a series of trajectory and FID matrices that will
%create one image per echo.
%
%Arguments
%         traj -  trajectory matrix that is 3xNFIDPtsxNPro
%         FID -  FID matrix that is NFIDPtsxNPro
%         NEcho -  The Number of echoes
%Returns 
%    EchoTraj1 - 3xNPtsxNPro matrix of trajectories for the first flyout
%    EchoFID1 - NPtsxNPro matrix of FID points for the first flyout
%    EchoTraj - 3x2NPtsxNPro x NEcho-1 matrix of trajectories for flyback
%    and flyout for all echoes after the first
%    EchoFID 2NPtsxNProxNecho-1 matrix of FID points for flyback and flyout
%    of all echoes after the first.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First, find the radii of the trajectory points
trajrad = squeeze(sqrt(traj(1,:,:).^2 + traj(2,:,:).^2 + traj(3,:,:).^2));
FID = FID/max(max(abs(FID)));
%Now, find indices of maximum radii for each echo - theoretically it should
%be the same for each echo, but I doubt that it is in practice
% First, need the number of points per flyout and back
NptsperEcho = floor(size(FID,1)/NEcho);

MaxR = zeros(1,NEcho+1); %index of maximum R
for i = 1:NEcho
    [maxr,indx] = max(trajrad((NptsperEcho*(i-1)+1):NptsperEcho*i,1));
    MaxR(i+1) = indx+NptsperEcho*(i-1); %Make it so that the first index to work with is 0
end

EchoTraj = zeros(3,floor(size(FID,1)/NEcho+50),size(FID,2),NEcho);
EchoFID = zeros(floor(size(FID,1)/NEcho+50),size(FID,2),NEcho); %Make these matrices much larger than it needs to be. I'll just delete the zeros before reconning
xfig = figure('Name','X Splitting');
yfig = figure('Name','Y Splitting');
zfig = figure('Name','Z Splitting');
colors = {'k','b','r','g','y'};
for i = 1:NEcho
    if i ==1
        indextraj = MaxR(i+1)-MaxR(i);
        EchoTraj(:,1:indextraj,:,i) = traj(:,(MaxR(i)+1):MaxR(i+1),:);
        EchoFID(1:indextraj,:,i) = FID((MaxR(i)+1):MaxR(i+1),:);
        
        figure(xfig)
        hold on
        plot(1:MaxR(i+1),traj(1,1:MaxR(i+1),100),colors{mod(i,5)},1:MaxR(i+1),abs(FID(1:MaxR(i+1),1)),colors{mod(i,5)})
        hold off
        figure(yfig)
        hold on
        plot(1:MaxR(i+1),traj(2,1:MaxR(i+1),100),colors{mod(i,5)},1:MaxR(i+1),abs(FID(1:MaxR(i+1),100)),colors{mod(i,5)})
        hold off
        figure(zfig)
        hold on
        plot(1:MaxR(i+1),traj(3,1:MaxR(i+1),100),colors{mod(i,5)},1:MaxR(i+1),abs(FID(1:MaxR(i+1),100)),colors{mod(i,5)})
        hold off
    else
        indextraj = MaxR(i+1) - MaxR(i)+1;
        EchoTraj(:,1:indextraj,:,i) = traj(:,MaxR(i):MaxR(i+1),:);
        EchoFID(1:indextraj,:,i) = FID(MaxR(i):MaxR(i+1),:);
        
        ind = mod(i,5);
        if ind == 0
            ind = 5;
        end
        
        figure(xfig)
        hold on
        plot(MaxR(i):MaxR(i+1),traj(1,MaxR(i):MaxR(i+1),100),colors{ind},MaxR(i):MaxR(i+1),abs(FID(MaxR(i):MaxR(i+1),100)),colors{ind})
        hold off
        figure(yfig)
        hold on
        plot(MaxR(i):MaxR(i+1),traj(2,MaxR(i):MaxR(i+1),100),colors{ind},MaxR(i):MaxR(i+1),abs(FID(MaxR(i):MaxR(i+1),100)),colors{ind})
        hold off
        figure(zfig)
        hold on
        plot(MaxR(i):MaxR(i+1),traj(3,MaxR(i):MaxR(i+1),100),colors{ind},MaxR(i):MaxR(i+1),abs(FID(MaxR(i):MaxR(i+1),100)),colors{ind})
        hold off
    end
end



