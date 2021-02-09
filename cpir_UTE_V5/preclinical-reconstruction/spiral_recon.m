function [Dir_Name,recon_def] = spiral_recon(path,traj,Method_Params,Config_Params,fast_recon)

cd(path)

Version = 'V1';
recon_def.function = ['spiral_recon.m version ' Version];

%I'll need to calibrate this to be proper values for proton and xenon
traj_delay = Config_Params.traj_delay;
if strcmp(Method_Params.Traj_Type,'measured') && strcmp('<129Xe>',Method_Params.Nucleus)
    traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);%2.75);
elseif strcmp(Method_Params.Traj_Type,'measured')
    traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);
elseif strcmp(Method_Params.Traj_Type,'theoretical')
    traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);
end

if fast_recon
    Dir_Name = ['Fast Recon ' datestr(now,30)];
else
    Dir_Name = ['Recon ' datestr(now,30)];
end
mkdir(Dir_Name);
%% Read in FIDs
disp('Reading FID File')
FIDs = Bruker_Load('fid');

if contains(Method_Params.Dims,'2D')
    if strcmp(Method_Params.DiffusionYN,'No')
        FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NSlices*Method_Params.Repetitions);
    else
        FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NSlices*Method_Params.Nbvalue*Method_Params.Repetitions);
    end
else
    if strcmp(Method_Params.DiffusionYN,'No')
        FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.Repetitions);
    else
        FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.Nbvalue*Method_Params.Repetitions);
    end
end
numZeroFill = find(FIDs(:,1)==0);
numZeroFill(numZeroFill <10)=[]; %Sometimes the first couple of points are 0, but we don't want to delete those
cutoffindx = [];
for i =length(numZeroFill):(-1):2
    if (numZeroFill(i)-numZeroFill(i-1))>1
        cutoffindx = i;
        break
    end
end
if isempty(cutoffindx)
    cutoffindx = 1;
end
numZeroFill(numZeroFill<numZeroFill(cutoffindx)) = [];
FIDs(numZeroFill,:)=[];

%Now, FIDs should be similarly shaped to trajectories and both still have
%extra points at the beginning - delete those here
FIDs(1:Method_Params.AcqShift,:) = [];
traj(:,1:Method_Params.AcqShift,:) = [];
%Make sure that the trajectories are within acceptable bounds

k0mag = abs(FIDs(1,:));
figure('Name','k0 Plot')
plot(1:length(k0mag),k0mag,'k-*')
%set the first point of each trajectory to true k0 so that deapodization
%works
for i = 1:size(traj,3)
    for j = 1:3
        traj(j,:,i) = traj(j,:,i)-traj(j,1,i);
    end
end
maxtraj = max(traj(:));
%check = traj/maxtraj/2.223;
%traj = check;
figure('Name','Trajectories')
xlim([-.5 .5])
ylim([-.5 .5])
zlim([-.5 .5])
hold on
for i = 1:Method_Params.NPro
    plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i))
   % pause(0.2)
end
hold off
%Temporary display of FIDs to make sure things are okay...
figure('Name','FIDs')
imagesc(abs(FIDs))

%Throw out projections based on info from configuration file
%I need to do this slightly differently to make this work correctly
if Config_Params.DiscardFromEnd ~= 0
    FIDs(:,(end - Config_Params.DiscardFromEnd):end) = NaN;% = FIDs(:,1:Config_Params.UseProj);
    %traj = traj(:,:,1:Config_Params.UseProj);
end
if Config_Params.DiscardFromBeginning ~= 0
    FIDs(:,1:Config_Params.DiscardFromBeginning) = NaN;% = FIDs(:,1:Config_Params.UseProj);
    %traj = traj(:,:,1:Config_Params.UseProj);
end

%If we're imaging xenon, I want to see k0
%if strcmp('<129Xe>',Method_Params.Nucleus)
    k0mag = abs(FIDs(1,:));
    figure('Name','k0 Plot')
    plot(1:length(k0mag),k0mag,'k-*')
%end

if strcmp(Method_Params.DiffusionYN,'Yes')
    FIDMat = zeros(Method_Params.NPts-Method_Params.AcqShift,Method_Params.NPro,Method_Params.NSlices*Method_Params.Nbvalue,Method_Params.Repetitions);
else
    FIDMat = zeros(Method_Params.NPts-Method_Params.AcqShift,Method_Params.NPro,Method_Params.NSlices,Method_Params.Repetitions);
end

for j = 1:Method_Params.Repetitions
    if strcmp(Method_Params.DiffusionYN,'Yes')
        SliceIndx = Method_Params.NSlices*Method_Params.Nbvalue;
    else
        SliceIndx = Method_Params.NSlices;
    end
    for i = 1:SliceIndx
        sliceIndex = i:SliceIndx:(size(FIDs,2)/Method_Params.Repetitions);
        sliceIndex = sliceIndex+(j-1)*(size(FIDs,2)/Method_Params.Repetitions);
        FIDMat(:,:,i,j) = FIDs(:,sliceIndex);
    end
end
if Method_Params.Repetitions>1
    sum = zeros(size(FIDMat,1),size(FIDMat,2),size(FIDMat,3));
    for i = 1:Method_Params.Repetitions
        sum = sum + FIDMat(:,:,:,i);
    end
    FIDMat(:,:,:,Method_Params.Repetitions+1) = sum/Method_Params.Repetitions;
end
trajx = reshape(traj(1,:,:),1,[])';
trajy = reshape(traj(2,:,:),1,[])';
trajz = reshape(traj(3,:,:),1,[])';
trajC = [trajx trajy trajz];

ImSize = Method_Params.MatrixSize;
if contains(Method_Params.Dims,'3D')
    Images = zeros(ImSize(1),ImSize(2),ImSize(3),size(FIDMat,3),size(FIDMat,4));
    Kspace = zeros(2*ImSize(1),2*ImSize(2),2*ImSize(3),size(FIDMat,3),size(FIDMat,4));
    if Config_Params.VentBPM ~=0 && contains(Method_Params.Dims,'3D') && ~fast_recon
        Norm_Breath_Images = zeros(ImSize(1),ImSize(2),ImSize(3),size(FIDMat,3),size(FIDMat,4));
        Norm_k0_Images = zeros(ImSize(1),ImSize(2),ImSize(3),size(FIDMat,3),size(FIDMat,4));
        T1_corr_Images = zeros(ImSize(1),ImSize(2),ImSize(3),size(FIDMat,3),size(FIDMat,4));
    else
        Norm_Breath_Images = 0;
        Norm_k0_Images = 0;
        T1_corr_Images = 0;
    end
else
    ImSize(3) = 1;
    Images = zeros(ImSize(1),ImSize(2),1,size(FIDMat,3),size(FIDMat,4));
    Kspace = zeros(2*ImSize(1),2*ImSize(2),2,size(FIDMat,3),size(FIDMat,4));
    Norm_Breath_Images = 0;
    Norm_k0_Images = 0;
    T1_corr_Images = 0;
end

%Preallocate to make sure these get set
M0 = 0;
FA_Map = 0;
KeyImages = 0;
%% Deal with Retrospective Gating - Need to make this better - get code from Alex
if strcmp(Method_Params.Trigger,'Off') && strcmp(Method_Params.Nucleus,'<1H>') && contains(Method_Params.Dims,'3D')

    if strcmp(Method_Params.DiffusionYN,'Yes')
        loop_end = Method_Params.Nbvalue;
    else
        loop_end = 1;
    end
    Images = zeros(ImSize(1),ImSize(2),ImSize(3),2*loop_end);

    for i = 1:loop_end
        %Check if retrogating has been performed previously - if not do it.
        %If it has, then just load in the retrogating file
        if isfile(['RetroGating' num2str(i) '.mat'])
            load(['RetroGating' num2str(i) '.mat'])
        else
            gen_retrogating_file(squeeze(FIDMat(:,:,i)),i,Method_Params);
            load(['RetroGating' num2str(i) '.mat'])
        end
        ExpFID = squeeze(FIDMat(:,Exp_indx,i));
        ExpTraj = squeeze(traj(:,:,Exp_indx));
        InspFID = FIDMat(:,Insp_indx,i);
        InspTraj = traj(:,:,Insp_indx);
        
        Exp_FID_1 = reshape(ExpFID,1,[])';
        ExpTraj_1x = reshape(ExpTraj(1,:,:),1,[])';
        ExpTraj_1y = reshape(ExpTraj(2,:,:),1,[])';
        ExpTraj_1z = reshape(ExpTraj(3,:,:),1,[])';
        ExpTraj_1 = [ExpTraj_1x ExpTraj_1y ExpTraj_1z];
        
        % Remove Trajectories and FID Points with radius > 0.5
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
        
        if ~isempty(ExpTraj_1)
            [Images(:,:,:,1+2*(i-1)),Gridded_KSpace(:,:,:,1+2*(i-1)),Recon_Params] = base_noncart_recon(ImSize,Exp_FID_1,ExpTraj_1);
        else
            Images(:,:,:,1+2*(i-1)) = zeros(ImSize);
        end
        if ~isempty(InspTraj_2)
            [Images(:,:,:,2+2*(i-1)),Gridded_KSpace(:,:,:,2+2*(i-1)),Recon_Params] = base_noncart_recon(ImSize,InspFID_2,InspTraj_2);
        else
            Images(:,:,:,2+2*(i-1)) = zeros(ImSize);
        end
        FileName = ['Image ' num2str(i) ' Gated_1'];
        FileName2 = ['Image ' num2str(i) ' Gated_2'];
        Complete_FileName = fullfile(pwd,Dir_Name,FileName);
        Complete_FileName2 = fullfile(pwd,Dir_Name,FileName2);
        DicomWrite(Complete_FileName,squeeze(Images(:,:,:,1+2*(i-1))),Method_Params.Resolution)
        DicomWrite(Complete_FileName2,squeeze(Images(:,:,:,2+2*(i-1))),Method_Params.Resolution)
        %Write out Niftis that are formatted properly for autosegmentation
        niftiwrite(abs(flip(rot90(squeeze(Images(:,:,:,1+2*(i-1))),3))),Complete_FileName)
        niftiwrite(abs(flip(rot90(squeeze(Images(:,:,:,2+2*(i-1))),3))),Complete_FileName2)
        
        imslice(abs(Images),'Reconstructed Image')
        KeyImages = 0;
        nbk0_KeyImages = 0;
        T1Corr_KeyImages = 0;
        T1n_fid = 0;
    end
else
    if strcmp(Method_Params.Nucleus,'<129Xe>') && strcmp(Method_Params.VFAYN,'No') && ~fast_recon
        Keyhole = true;
        KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
        nbk0_KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
        T1Corr_KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
    else
        Keyhole = false;
        KeyImages = 0;
        nbk0_KeyImages = 0;
        T1Corr_KeyImages = 0;
        T1n_fid = 0;
    end
    for j = 1:size(FIDMat,4)
        for i = 1:size(FIDMat,3)
            traj_recon = trajC;
            tmp_FID = reshape(FIDMat(:,:,i,j),1,[])';
            if (Config_Params.DiscardFromBeginning ~= 0) || (Config_Params.DiscardFromEnd ~= 0)
                find_zeros = find(isnan(tmp_FID));
                tmp_FID(find_zeros) = [];
                traj_recon(find_zeros,:) = [];
            end
            rad = sqrt(traj_recon(:,1).^2+traj_recon(:,2).^2+traj_recon(:,3).^2);
            del_pts = find(rad > 0.5);
            tmp_FID(del_pts) = [];
            traj_recon(del_pts,:) = [];
            [Images(:,:,:,i,j),Kspace(:,:,:,i,j),Recon_Params] = base_noncart_recon(ImSize,tmp_FID,traj_recon);
            %Write out Niftis that are formatted properly for autosegmentation
            Complete_FileName = fullfile(pwd,Dir_Name,['Image ' num2str(i) num2str(j)]);
            if contains(Method_Params.Dims,'3D')
                niftiwrite(abs(flip(rot90(squeeze(Images(:,:,:,i,j)),3))),Complete_FileName)
            end
            if Config_Params.VentBPM ~=0 && contains(Method_Params.Dims,'3D') && ~fast_recon
                %Here's where I need to do al that stuff... Start
                %with the easy cases:
                %tmp_FID = reshape(FIDMat(:,:,i,j),1,[])';
                tmp_FID = squeeze(FIDMat(:,:,i,j));
                tmp_fid = tmp_FID;
                %Normalized first projection of each breath
%                 nb_fid = tmp_fid;
%                 for ii = 1:Method_Params.ProjPerTrig:Method_Params.NPro
%                     scale_factor = abs(tmp_fid(1,1))/abs(tmp_fid(1,ii));
%                     indices = ii:(ii+Method_Params.ProjPerTrig-1);
%                     indices(indices>Method_Params.NPro) = [];
%                     nb_fid(:,indices) = tmp_fid(:,indices)*scale_factor;
%                 end
%                 %normalized all projections
%                 n_fid = tmp_fid;
%                 for ii = 1:Method_Params.NPro
%                     scale_factor = abs(tmp_fid(1,1))/abs(tmp_fid(1,ii));
%                     n_fid(:,ii) = tmp_fid(:,ii)*scale_factor;
%                 end
                meank0 = mean(abs(tmp_FID(1,:)));
                npts = size(tmp_FID,1);
                n_fid = tmp_FID;
                for ii = 1:Method_Params.NPro
                    scale_factor = meank0/abs(tmp_FID(1,ii));
                    n_fid(round(npts/2):end,ii) = tmp_FID(round(npts/2):end,ii)*scale_factor;
                end
                
                %Fitting to T1 decay:
                pts = abs(tmp_FID(1,1:Method_Params.ProjPerTrig:Method_Params.NPro));
                %Using Moller, Cleveland, Driehuys:
                %Approximate surface area of bag:
                A = 328; %Surface area in cm^2 - This is probably right for us
                %Config file will have breaths per minute (#/min), tidal
                %volume (mL), O2 Percent (%)
                %Change in volume in mL per min
                Vdot = (Config_Params.VentTV - Config_Params.VentTV * Config_Params.VentO2Pct/100) * Config_Params.VentBPM;
                %In this case, the Breath TR is 1/BPM
                BTR = 1/Config_Params.VentBPM;
                Bt = 0:BTR:(BTR*(length(pts)-1));

                pts = pts/pts(1);
                %Let's throw out points greater than 10% higher
                %than the first point
                bad_pts = find(pts>1.1);
                good_pts = pts;
                Bt_good = Bt;
                Bt_good(bad_pts) = [];
                good_pts(bad_pts) = [];
                fiteq = fittype('a*exp(-x/b)*((c-n*x)/c)^(m*d/n)','problem',{'n','m'},'independent','x','coefficients',{'a','b','c','d'});
                %In this equation, 
                %a = Mz(0)
                %b = T1_reservoir
                %c = V(0)
                %d = kappa
                %n = Vdot
                %m = A
                %x = BTR
                try
                    T1fit = fit(Bt_good',good_pts',fiteq,'problem',{Vdot,A},'StartPoint',[1,120,300,0.5],'Upper',[Inf,500,500,3],'Lower',[0,0,0,0]);
                    goodfit = true;
                catch
                    try 
                        warning('Fit Failed - Trying fit on smoothed points')
                        smoothpts = smooth(good_pts',15);
                        T1fit = fit(Bt_good',smoothpts,fiteq,'problem',{Vdot,A},'StartPoint',[1,120,300,0.5],'Upper',[Inf,300,500,3],'Lower',[0,0,0,0]);
                        goodfit = true;
                    catch
                        warning('Fit Failed. Skipping T1 fitting and correction')
                        goodfit = false;
                    end
                end
                %Temporary Plot to see how fitting went - It works
                %- don't need to bother displaying these anymore
%                         figure('Name','Normalized first proj')
%                         plot(abs(nb_fid(1,:)),'-*');
%                         figure('Name','Normalized all')
%                         plot(abs(n_fid(1,:)),'-*');
%                         figure('Name','Check Fit')
%                         plot(T1fit,Bt,pts);
                if goodfit
                    T1_fid = tmp_fid/abs(tmp_fid(1,1));
                    T1_pts = T1fit(Bt);
                    diff = abs(pts-T1_pts')./T1_pts';
                    %Toss Points that are more than 20% away from the fit line 
                    throw_away = find(diff > 0.2);
                    Check_Fit_Fig = figure('Name','Examine Reservoir T1 Fit and Discarded Points');
                    plot(Bt,pts,'ko','MarkerFaceColor','k')
                    hold on
                    plot(Bt,T1_pts,'r','LineWidth',2)
                    plot(Bt(throw_away),pts(throw_away),'ro','MarkerFaceColor','r')
                    hold off    
                    saveas(Check_Fit_Fig,'Check_T1_Figure'); 
                    counter = 1;
                    for ii = 1:Method_Params.ProjPerTrig:Method_Params.NPro
                        scale_factor = 1/abs(T1_pts(counter));
                        indices = ii:(ii+Method_Params.ProjPerTrig-1);
                        indices(indices>Method_Params.NPro) = [];
                        T1_fid(:,indices) = T1_fid(:,indices)*scale_factor;
                        counter = counter + 1;
                    end
                    
                    T1n_fid = T1_fid;
                    for ii = 1:Method_Params.NPro
                        scale_factor = meank0/abs(T1_fid(1,ii));
                        T1n_fid(round(npts/2):end,ii) = T1_fid(round(npts/2):end,ii)*scale_factor;
                    end
                    
                    throw_away_pts = [];
                    for ii = 1:length(throw_away)
                        throw_away_pts_hold = (throw_away(ii)-1) * Method_Params.ProjPerTrig + (1:Method_Params.ProjPerTrig);
                        throw_away_pts = [throw_away_pts throw_away_pts_hold];
                    end
                    T1n_fid(:,throw_away_pts) = [];
                    figure('Name','Check k0 after fitting, scaling, and discarding')
                    plot(abs(T1n_fid(1,:)),'*k')

                    %Get ready for reconstruction
                    T1_traj = traj;
                    T1_traj(:,:,throw_away_pts) = [];
                    T1_trajx = reshape(T1_traj(1,:,:),1,[])';
                    T1_trajy = reshape(T1_traj(2,:,:),1,[])';
                    T1_trajz = reshape(T1_traj(3,:,:),1,[])';
                    T1_trajC = [T1_trajx,T1_trajy,T1_trajz];
                    T1_fidC = reshape(T1n_fid,1,[])';
                    T1_trajC1 = T1_trajC;
                else
                    T1n_fid = 0;
                end

                nfid_C = reshape(n_fid,1,[])';
                %nbfid_C = reshape(nb_fid,1,[])';

                traj_recon_1 = traj_recon;
                rad = sqrt(traj_recon_1(:,1).^2+traj_recon_1(:,2).^2+traj_recon_1(:,3).^2);
                del_pts = find(rad > 0.5);
                nfid_C(del_pts) = [];
                %nbfid_C(del_pts) = [];
                traj_recon_1(del_pts,:) = [];

                %[Norm_Breath_Images(:,:,:,i,j),~,~] = base_noncart_recon(ImSize,nbfid_C,traj_recon_1);
                [Norm_k0_Images(:,:,:,i,j),~,~] = base_noncart_recon(ImSize,nfid_C,traj_recon_1);
                if goodfit
                    rad1 = sqrt(T1_trajC(:,1).^2+T1_trajC(:,2).^2+T1_trajC(:,3).^2);
                    del_pts = find(rad1 > 0.5);
                    T1_fidC(del_pts) = [];
                    T1_trajC(del_pts,:) = [];
                    [T1_corr_Images(:,:,:,i,j),~,~] = base_noncart_recon(ImSize,T1_fidC,T1_trajC);
                else
                    T1_corr_Images = 0;
                end

                 if strcmp(Method_Params.VFAYN,'No') && Method_Params.ProjPerTrig>1 
                    nbk0_KeyImages = preclinical_keyhole(n_fid,traj_recon,Method_Params);
                    if goodfit
                        T1Corr_KeyImages = preclinical_keyhole(T1n_fid,T1_trajC1,Method_Params);
                    else
                        T1Corr_KeyImages = 0;
                    end
                 end
            end
        end
    end
    imslice(abs(squeeze(Images)))
    if contains(Method_Params.Nucleus,'129Xe') && strcmp(Method_Params.VFAYN,'No') && contains(Method_Params.Dims,'3D') && Method_Params.ProjPerTrig>1 && ~fast_recon 
        PPT = Method_Params.ProjPerTrig;
        %I think that we're generally only going to be concerned with the
        %first FID -
        %Get FIDs and trajectories - need to account for both basic UTE and
        %Flyback cases
        KeyFID = squeeze(FIDMat(:,:,1));
        trajx = reshape(traj(1,:,:),1,[])';
        trajy = reshape(traj(2,:,:),1,[])';
        trajz = reshape(traj(3,:,:),1,[])';
        trajC = [trajx trajy trajz];
        KeyTraj = trajC;
        %Get Key Radius - Make it have a minimum radius of 1/4 of the
        %original number of points along the projection
        KeyRadLim = 0.25;
        %Generate Keyhole
        Keyhole = KeyFID;
        rad1 = squeeze(sqrt(traj(1,:,:).^2 + traj(2,:,:).^2 + traj(3,:,:).^2));
        
        KeyRad = KeyRadLim *max(rad1(:));
        KeyKSpace = zeros(size(KeyFID,1),size(KeyFID,2),PPT);
        for ii = 1:PPT
            KeyKSpace(:,:,ii) = Keyhole;
            for kk = 1:size(FIDMat,2)
                KeyKSpace(rad1(:,kk)<KeyRad,kk,ii) = NaN;
                if mod(kk,PPT) == ii
                    KeyKSpace(rad1(:,kk)<KeyRad,kk,ii) = KeyFID(rad1(:,kk)<KeyRad,kk);
                elseif mod(kk,PPT) == 0
                    KeyKSpace(rad1(:,kk)<KeyRad,kk,ii) = KeyFID(rad1(:,kk)<KeyRad,kk);
                end
            end
            keyk0 = abs(KeyFID(1,ii:PPT:size(KeyFID,2)));
            MeanKeyk0 = mean(keyk0);
            for jj = 1:size(KeyFID,2)
                KeyKSpace(:,jj,ii) = KeyKSpace(:,jj,ii)*MeanKeyk0/abs(KeyFID(1,jj));
            end
            fid = reshape(squeeze(KeyKSpace(:,:,ii)),1,[])';
            findzeros = find(isnan(fid));
            ThisKeyTraj = KeyTraj;
            ThisKeyTraj(findzeros,:) = [];
            fid(findzeros) = [];
            
            %Remove Points that are beyond a trajectory radius of 0.5
            traj_recon_1 = ThisKeyTraj;
            rad = sqrt(traj_recon_1(:,1).^2+traj_recon_1(:,2).^2+traj_recon_1(:,3).^2);
            del_pts = find(rad > 0.5);
            fid(del_pts) = [];
            traj_recon_1(del_pts,:) = [];
            
            [KeyImages(:,:,:,ii),~,~] = base_noncart_recon(ImSize,fid,traj_recon_1);
        end
        imslice(abs(KeyImages),'Reconstructed Key Images')
        [~,Mask] = erode_dilate(squeeze(KeyImages(:,:,:,1)),1,5);
        [M0,FA_Map] = key_image_fitting(KeyImages,Mask);
        imslice(abs(M0),'Magnetization Based on Key Fitting')
        imslice(FA_Map,'Flip Angle Map')
    else
        M0 = 0;
        FA_Map = 0;
    end
end

Recon.Keys = KeyImages;
Recon.FA_Map = FA_Map;
Recon.M0 = M0;
Recon.Images = Images;
Recon.normk0_Images = Norm_k0_Images;
%Recon.normbreath_Images = Norm_Breath_Images;
Recon.T1_corr_Images = T1_corr_Images;
Recon.KSpace = Kspace;
Raw.FIDs = FIDs;
Raw.Traj = traj;
Recon.nbk0_Keys = nbk0_KeyImages;
Recon.T1Corr_Keys = T1Corr_KeyImages;

Recon_Params.traj_delay = traj_delay; %Make sure to write out what trajectory delay was used for recon
recon_def.Params = Recon_Params;

fclose('all');
if fast_recon
    save(fullfile(pwd,Dir_Name,'Fast_Recon_WorkSpace.mat'),'Recon','Raw','Method_Params','Config_Params')
else
    save(fullfile(pwd,Dir_Name,'Recon_WorkSpace.mat'),'Recon','Raw','Method_Params','Config_Params')
end

