function recon_ute(path)

%% If no path is passed, get path
if nargin == 0
    path = uigetdir();
end

cd(path)

%% Load Method Parameters and Trajectories
[traj,Method_Params] = read_method_ute(path);

%% Sanity Check - View first 100 trajectory points (Can be commented out)
figure('Name','Trajectory Sanity Check')
hold on
for i = 1:1000
    plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i))
end
hold off

%% Trajectory Delay Correction
disp('Perform Trajectory Delay Correction')
%Feel free to implement this trajectory delay in a more user-friendly way -
%This is quick and easy
traj_delay = [0 0 0];%[0 0 1];
traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);

%% Read in FIDs
disp('Read FID File')
FIDs = Bruker_Load('fid');
%Reshape FIDs
FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NumTEs);
FIDMat = zeros(size(FIDs,1),Method_Params.NPro,Method_Params.NumTEs);
for i = 1:Method_Params.NumTEs
    FIDMat(:,:,i) = FIDs(:,i:Method_Params.NumTEs:end);
end
FIDs = FIDMat;

%FIDs are "zero-filled" - That is, a bunch of zeros are added to the end of
%every FID to get the number of points to a nice power of 2. We need to get
%rid of those.
findzeros = find(FIDs(:,1,1)==0);
FIDs(findzeros,:,:) = [];

FIDs(1:Method_Params.AcqShift,:,:) = [];
traj(:,1:Method_Params.AcqShift,:) = [];

%% FID sanity Check
figure('Name','FID Sanity Check')
subplot(1,2,1)
plot(abs(FIDs(:,1,1)))
subplot(1,2,2)
imagesc(abs(FIDs(:,:,1)))
colormap(gray)
axis off

for i = 1:Method_Params.NumTEs
    %% Retrospective Gating
    %Check if retrogating has been performed previously - if not do it.
    %If it has, then just load in the retrogating file
    if isfile(['RetroGating' num2str(i) '.mat'])
        load(['RetroGating' num2str(i) '.mat'])
    else
        gen_retrogating_file_basic(squeeze(FIDs(:,:,i)),Method_Params,i);
        load(['RetroGating' num2str(i) '.mat'])
    end
    if ~isnan(Insp_indx)
        ExpFID = squeeze(FIDs(:,Exp_indx,i));
        ExpTraj = squeeze(traj(:,:,Exp_indx));
        InspFID = FIDs(:,Insp_indx,i);
        InspTraj = traj(:,:,Insp_indx);

        Exp_FID_1 = reshape(ExpFID,1,[])';
        ExpTraj_1x = reshape(ExpTraj(1,:,:),1,[])';
        ExpTraj_1y = reshape(ExpTraj(2,:,:),1,[])';
        ExpTraj_1z = reshape(ExpTraj(3,:,:),1,[])';
        ExpTraj_1 = [ExpTraj_1x ExpTraj_1y ExpTraj_1z];

        rad_1 = sqrt(ExpTraj_1(:,1).^2+ExpTraj_1(:,2).^2+ExpTraj_1(:,3).^2);
        del_pts = find(rad_1 > 0.5);
        Exp_FID_1(del_pts) = [];
        ExpTraj_1(del_pts,:) = [];

        InspFID_2 = reshape(InspFID,1,[])';
        InspTraj_2x = reshape(InspTraj(1,:,:),1,[])';
        InspTraj_2y = reshape(InspTraj(2,:,:),1,[])';
        InspTraj_2z = reshape(InspTraj(3,:,:),1,[])';
        InspTraj_2 = [InspTraj_2x InspTraj_2y InspTraj_2z];

        % Remove Trajectories and FID Points with radius > 0.5
        rad_2 = sqrt(InspTraj_2(:,1).^2+InspTraj_2(:,2).^2+InspTraj_2(:,3).^2);
        del_pts = find(rad_2 > 0.5);
        InspFID_2(del_pts) = [];
        InspTraj_2(del_pts,:) = [];
        ImSize = Method_Params.MatrixSize;

        [Exp_Image(:,:,:,i),~,~] = base_noncart_recon(ImSize,Exp_FID_1,ExpTraj_1);
        [Insp_Image(:,:,:,i),~,~] = base_noncart_recon(ImSize,InspFID_2,InspTraj_2);
        
       % imslice(abs(Insp_Image),['Inspiration Image Echo ' num2str(i)])
       % imslice(abs(Exp_Image),['Expiration Image Echo ' num2str(i)])
        niftiwrite(abs(Exp_Image(:,:,:,i)),['Expiration Image Echo ' num2str(i)])
        niftiwrite(abs(Insp_Image(:,:,:,i)),['Inspiration Image Echo ' num2str(i)])
        
        save(['Recon Workspace ' num2str(i) '.mat'],'FIDs','traj','Exp_Image','Insp_Image');
    else
        %% Reshape for reconstruction
        trajx = reshape(traj(1,:,:),1,[])';
        trajy = reshape(traj(2,:,:),1,[])';
        trajz = reshape(traj(3,:,:),1,[])';

        recon_traj = [trajx trajy trajz];

        recon_fid = reshape(FIDs(:,:,i),1,[])';

        %trajectories bigger than 0.5 cause issues
        rad = sqrt(recon_traj(:,1).^2+recon_traj(:,2).^2+recon_traj(:,3).^2);
        toobig = find(rad>0.5);
     %   disp(['Discarded Points from End of Readout: ' num2str(nnz(toobig))]);
        recon_traj(toobig,:) = [];
        recon_fid(toobig) = [];

        %% Reconstruct
        ImSize = Method_Params.MatrixSize;
        %I'm not bothering to write out gridded k-space or recon parameters, but
        %feel free to do so if you want them.
        [Image(:,:,:,i),~,~] = base_noncart_recon(ImSize,recon_fid,recon_traj);
      %  imslice(abs(Image),['Echo Image ' num2str(i)])
        niftiwrite(abs(Image(:,:,:,i)),['Echo Image ' num2str(i)])
        save(['Recon Workspace ' num2str(i) '.mat'],'FIDs','traj','Image');
    end
end
load(['RetroGating' num2str(i) '.mat'])
if ~isnan(Insp_indx)
    imslice(abs(Insp_Image),'Inspiration Image')
    imslice(abs(Exp_Image),'Expiration Image')
else
    imslice(abs(Image),'Echo Image')
end
