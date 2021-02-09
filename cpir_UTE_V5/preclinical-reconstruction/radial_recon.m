function [Dir_Name,recon_def] = radial_recon(path,traj,Method_Params,Config_Params,fast_recon)

%Total reconstruction for "clean UTE" in Bruker - currently named
%pn_cleanUTE - will be renamed for the final version
%Make this very modular -
%Steps:
%----Validate that we are in the proper path and have the correct sequence for recon 
%----Generate/Read Configuration File
%----Read Method (maybe this first so that we only read the method once)
%----Get trajectories - I think I'll do this from the method file rather
%than the traj file
%----Read FIDs
%----Apply Trajectory Delay (if desired)
%----Reshape FIDs and Trajectories
%----Reconstruct Images
%----Perform Keyhole Analysis/T2* Analysis if desired
%----Save images, raw data, trajectories, analysis (Do in a version
%controlled way so that each reconstruction has its own folder rather than
%overwriting previous recons)

cd(path)

Version = 'V1';
recon_def.function = ['radial_recon.m version ' Version];
traj_delay = Config_Params.traj_delay;
%% Trajectory Delay - would be really nice to have a gui sort of thing for this
if strcmp(Method_Params.Traj_Type,'measured')
    traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);
elseif strcmp(Method_Params.Traj_Type,'theoretical')
    traj = traj_delay_correction(traj,Method_Params.Dwell,Method_Params.Dwell*traj_delay);
end
%% Read in FIDs
disp('Reading FID File')
FIDs = Bruker_Load('fid');
%Reshape FIDs accounting for NPro, NumTEs, and Number of b values
if strcmp(Method_Params.FlybackYN,'Yes') 
    FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.Repetitions);%,Method_Params.Repetitions);
elseif strcmp(Method_Params.DiffusionYN,'Yes')
    FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NumTEs*Method_Params.Nbvalue*Method_Params.Repetitions);%,Method_Params.Repetitions);
else
    FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NumTEs*Method_Params.Repetitions);%,Method_Params.Repetitions);
end
%If more than one repetition, need to further reshape matrices
if Method_Params.Repetitions > 1
    nProPerRep = size(FIDs,2)/Method_Params.Repetitions;
    for i = 1:Method_Params.Repetitions
        NewFID(:,:,i) = FIDs(:,(1:nProPerRep)+nProPerRep*(i-1));
    end
    FIDs = NewFID;
end 

numZeroFill = find(squeeze(FIDs(:,1,1))==0);
numZeroFill(numZeroFill <10)=[]; %Sometimes the first couple of points are 0, but we don't want to delete those
FIDs(numZeroFill,:,:)=[];

%Now, FIDs should be similarly shaped to trajectories and both still have
%extra points at the beginning - delete those here
FIDs(1:Method_Params.AcqShift,:,:) = [];
traj(:,1:Method_Params.AcqShift,:) = [];

%When doing xenon, the last point on the projection is high intensity and
%gives a target-like artifact on the reconstructed image - kill that here
if contains(Method_Params.Nucleus,'129Xe')
    FIDs(end,:) = [];
    traj(:,end,:) = [];
end

%set the first point of each trajectory to true k0 so that deapodization
%works
for i = 1:size(traj,3)
    for j = 1:3
        traj(j,:,i) = traj(j,:,i)-traj(j,1,i);
    end
end

%Throw out projections based on info from configuration file
if Config_Params.DiscardFromEnd ~= 0
    FIDs(:,(end - Config_Params.DiscardFromEnd):end,:) = [];
    traj(:,:,(end - Config_Params.DiscardFromEnd):end) = [];
end
if Config_Params.DiscardFromBeginning ~=0
    FIDs(:,1:Config_Params.DiscardFromBeginning,:) = [];
    traj(:,:,1:Config_Params.DiscardFromBeginning) = [];
end

Method_Params.NPro = size(traj,3);

if size(traj,2) > size(FIDs,1)
    traj = traj(:,1:size(FIDs,1),:);
end

%If we're imaging xenon, I want to see k0
if strcmp('<129Xe>',Method_Params.Nucleus)
    k0mag = abs(FIDs(1,:));
    figure('Name','k0 Plot')
    plot(0:Method_Params.TR:(Method_Params.TR*(length(k0mag)-1)),k0mag,'k-*')
end

for i = 1:Method_Params.Repetitions
    %I probably need to update this for diffusion combined with flyback
    if strcmp(Method_Params.FlybackYN,'Yes') 
        [EchoTraj,EchoFID] = Separate_Echoes(traj,FIDs,Method_Params.NumTEs);
    elseif strcmp(Method_Params.DiffusionYN,'Yes')
        FID_Mat(:,:,:,i) = splitFID(FIDs(:,:,i),Method_Params.NumTEs*Method_Params.Nbvalue);
    elseif Method_Params.NumTEs > 1
        FID_Mat(:,:,:,i) = splitFID(FIDs(:,:,i),Method_Params.NumTEs);
    else
        FID_Mat = FIDs;
    end
end
%Sum together FIDs for signal averaging the repetitions
if strcmp(Method_Params.FlybackYN,'No') 
    if strcmp(Method_Params.DiffusionYN,'Yes') || Method_Params.NumTEs > 1
        sum = zeros(size(FID_Mat,1),size(FID_Mat,2),size(FID_Mat,3));
        for i = 1:Method_Params.Repetitions
            sum = sum + FID_Mat(:,:,:,i);
        end
        sum = sum/Method_Params.Repetitions;
        FID_Mat(:,:,:,Method_Params.Repetitions+1)=sum;
    elseif Method_Params.Repetitions > 1
        sum = zeros(size(FID_Mat,1),size(FID_Mat,2));
        for i = 1:Method_Params.Repetitions
            sum = sum + FID_Mat(:,:,i);
        end
        sum = sum/Method_Params.Repetitions;
        FID_Mat(:,:,Method_Params.Repetitions+1)=sum;
    end
end
%Define T1_fid, M0, and FA Map so that it is properly populated
M0 = 0;
FA_Map = 0;
T1_fid = 0;

%% Preallocate Memory
%Preallocate memory for Images
ImSize = Method_Params.MatrixSize;
if strcmp(Method_Params.DiffusionYN,'Yes')
    EchoImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs*Method_Params.Nbvalue,Method_Params.Repetitions);
    if Config_Params.VentBPM ~=0 && ~fast_recon
       % Norm_Breath_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs*Method_Params.Nbvalue,Method_Params.Repetitions);
        Scale_HF_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs*Method_Params.Nbvalue,Method_Params.Repetitions);
        T1_corr_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs*Method_Params.Nbvalue,Method_Params.Repetitions);
    else
%         Norm_Breath_Images = 0;
        Scale_HF_Images = 0;
        T1_corr_Images = 0;
    end
else
    EchoImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs,Method_Params.Repetitions);
    if Config_Params.VentBPM ~=0 && ~fast_recon
%         Norm_Breath_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs,Method_Params.Repetitions);
        Scale_HF_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs,Method_Params.Repetitions);
        T1_corr_Images = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.NumTEs,Method_Params.Repetitions);
    else
%         Norm_Breath_Images = 0;
        Scale_HF_Images = 0;
        T1_corr_Images = 0;
    end
end
Gridded_KSpace = zeros(ImSize(1)*2,ImSize(2)*2,ImSize(3)*2,Method_Params.NumTEs,Method_Params.Repetitions);

%% Before Reconstruction, Create Folder where Date will be saved
% Create a Reconstruction Folder with a timestamp so that I don't overwrite
% anything
if fast_recon
    Dir_Name = ['Fast Recon ' datestr(now,30)];
else
    Dir_Name = ['Recon ' datestr(now,30)];
end
mkdir(Dir_Name);

%% Deal with Retrospective Gating
if strcmp(Method_Params.FlybackYN,'No') && strcmp(Method_Params.Trigger,'Off') && strcmp(Method_Params.Nucleus,'<1H>')
    %If we're imaging proton in a mouse and it's not triggered, we need to
    %do retrospective gating - The config file is used to crop points, so I
    %don't need to bother with that here.
    if strcmp(Method_Params.DiffusionYN,'Yes')
        loop_end = Method_Params.NumTEs*Method_Params.Nbvalue;
    else
        loop_end = Method_Params.NumTEs;
    end
    EchoImages = zeros(ImSize(1),ImSize(2),ImSize(3),2*loop_end);
    %Loop through echoes, Check if things have been retrogated, and load in
    %retrogating files
    for i = 1:loop_end
        %Check if retrogating has been performed previously - if not do it.
        %If it has, then just load in the retrogating file
        if isfile(['RetroGating' num2str(i) '.mat'])
            load(['RetroGating' num2str(i) '.mat'])
        else
            gen_retrogating_file(squeeze(FID_Mat(:,:,i)),i,Method_Params);
            load(['RetroGating' num2str(i) '.mat'])
        end
        ExpFID = squeeze(FID_Mat(:,Exp_indx,i));
        ExpTraj = squeeze(traj(:,:,Exp_indx));
        InspFID = FID_Mat(:,Insp_indx,i);
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
        
        [EchoImages(:,:,:,1+2*(i-1)),Gridded_KSpace(:,:,:,1+2*(i-1)),Recon_Params] = base_noncart_recon(ImSize,Exp_FID_1,ExpTraj_1);
        [EchoImages(:,:,:,2+2*(i-1)),Gridded_KSpace(:,:,:,2+2*(i-1)),Recon_Params] = base_noncart_recon(ImSize,InspFID_2,InspTraj_2);
        
        FileName = ['Image ' num2str(i) ' Expiration'];
        FileName2 = ['Image ' num2str(i) ' Inspiration'];
        Complete_FileName = fullfile(pwd,Dir_Name,FileName);
        Complete_FileName2 = fullfile(pwd,Dir_Name,FileName2);
        DicomWrite(Complete_FileName,squeeze(EchoImages(:,:,:,1+2*(i-1))),Method_Params.Resolution)
        DicomWrite(Complete_FileName2,squeeze(EchoImages(:,:,:,2+2*(i-1))),Method_Params.Resolution)
        %Write out Niftis that are formatted properly for autosegmentation
        niftiwrite(abs(flip(rot90(squeeze(EchoImages(:,:,:,1+2*(i-1))),3))),Complete_FileName)
        niftiwrite(abs(flip(rot90(squeeze(EchoImages(:,:,:,2+2*(i-1))),3))),Complete_FileName2)
        
        KeyImages = 0;
        scale_KeyImages = 0;
        T1Corr_KeyImages = 0;
        T1n_fid =0;
    end
   % [ExpFID,ExpTraj,InspFID,InspTraj,GatingFig]= RetroGating_Comp(traj,FIDs,500,Method_Params.TR,70);
   % [InspImage,ExpImage] = retro_gated_recon(ExpFID,ExpTraj,InspFID,InspTraj,ImSize);
   % EchoImages(:,:,:,1) = InspImage;
   % EchoImages(:,:,:,2) = ExpImage;
%Deal with the case of wash out imaging
elseif Config_Params.WashinOut
    EchoImages = washout_recon(FID_Mat,traj,Method_Params,Config_Params);
    KeyImages = 0;
    Gridded_KSpace = 0;
    M0 = 0;
    FA_Map = 0;
else    
    %% Reconstruct Images
    disp('Reconstructing Images')
    if strcmp(Method_Params.Nucleus,'<129Xe>') && strcmp(Method_Params.VFAYN,'No')
        Keyhole = true;
        KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
        if ~fast_recon
            scale_KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
            T1Corr_KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),Method_Params.ProjPerTrig);
        else
            scale_KeyImages = 0;
            T1Corr_KeyImages = 0;
            T1n_fid = 0;
        end
    else
        Keyhole = false;
        KeyImages = 0;
        scale_KeyImages = 0;
        T1Corr_KeyImages = 0;
        T1n_fid = 0;
    end

    if strcmp(Method_Params.DiffusionYN,'Yes')
        loop_end = Method_Params.NumTEs*Method_Params.Nbvalue;
    else
        loop_end = Method_Params.NumTEs;
    end
    trajx = reshape(traj(1,:,:),1,[])';
    trajy = reshape(traj(2,:,:),1,[])';
    trajz = reshape(traj(3,:,:),1,[])';
    traj_recon = [trajx trajy trajz];
    if Method_Params.Repetitions > 1
        reps = Method_Params.Repetitions+1;
    else
        reps = Method_Params.Repetitions;
    end
    for j = 1:reps
        for i = 1:loop_end
            if strcmp(Method_Params.Trigger,'On') || strcmp(Method_Params.Nucleus,'<129Xe>')
                if strcmp(Method_Params.FlybackYN,'No') 
                    %I think I would like to have some sort of way to
                    %reconstruct xenon images with k0 corrections - fit to
                    %reservoir decay, normalizing first projection of each
                    %breath, and normalizing all projections
                    
                    %Update 6/3 - After looking at things a little closer,
                    %I think the better way is to only scale high frequency
                    %data to the mean of k0 - In point spread simulations,
                    %this totally removes imaging artifacts, and it
                    %preserves signal intensity
                    if Config_Params.VentBPM ~=0 && ~fast_recon
                        %Here's where I need to do al that stuff... Start
                        %with the easy cases:
                         if ndims(FID_Mat)>3
                            tmp_fid = squeeze(FID_Mat(:,:,i,j));
                        elseif loop_end>1
                            tmp_fid = squeeze(FID_Mat(:,:,i));
                        elseif reps>1
                            tmp_fid = squeeze(FID_Mat(:,:,i));
                         else
                            tmp_fid = FID_Mat;
                        end
%                         %Normalized first projection of each breath
%                         nb_fid = tmp_fid;
%                         for ii = 1:Method_Params.ProjPerTrig:Method_Params.NPro
%                             scale_factor = abs(tmp_fid(1,1))/abs(tmp_fid(1,ii));
%                             indices = ii:(ii+Method_Params.ProjPerTrig-1);
%                             indices(indices>Method_Params.NPro) = [];
%                             nb_fid(:,indices) = tmp_fid(:,indices)*scale_factor;
%                         end
%                         %normalized all projections
%                         n_fid = tmp_fid;
%                         for ii = 1:Method_Params.NPro
%                             scale_factor = abs(tmp_fid(1,1))/abs(tmp_fid(1,ii));
%                             n_fid(:,ii) = tmp_fid(:,ii)*scale_factor;
%                         end
                        meank0 = mean(abs(tmp_fid(1,:)));
                        npts = size(tmp_fid,1);
                        n_fid = tmp_fid;
                        for ii = 1:Method_Params.NPro
                            scale_factor = meank0/abs(tmp_fid(1,ii));
                            n_fid(round(npts/2):end,ii) = tmp_fid(round(npts/2):end,ii)*scale_factor;
                        end
                        %Fitting to T1 decay:
                        pts = abs(tmp_fid(1,1:Method_Params.ProjPerTrig:Method_Params.NPro));
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
                                smoothpts = smooth(good_pts',20);
                                T1fit = fit(Bt_good',smoothpts,fiteq,'problem',{Vdot,A},'StartPoint',[1,120,300,0.5],'Upper',[Inf,300,500,3],'Lower',[0,0,0,0]);
                                goodfit = true;
                            catch
                                warning('Fit did not work. Skipping T1 fitting and correction')
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
                            
                            meank0 = mean(abs(T1_fid(1,:)));
                            npts = size(T1_fid,1);
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
                        [Scale_HF_Images(:,:,:,i,j),~,~] = base_noncart_recon(ImSize,nfid_C,traj_recon_1);
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
                            scale_KeyImages = preclinical_keyhole(n_fid,traj_recon,Method_Params);
                            if goodfit
                                T1Corr_KeyImages = preclinical_keyhole(T1n_fid,T1_trajC1,Method_Params);
                            else
                                T1Corr_KeyImages = 0;
                            end
                         end
                    end
                    if ndims(FID_Mat)>3
                        tmp_fid = squeeze(FID_Mat(:,:,i,j));
                        fid = reshape(tmp_fid,1,[])';
                    elseif loop_end>1
                        fid = reshape(FID_Mat(:,:,i),1,[])';
                    elseif reps>1
                        fid = reshape(FID_Mat(:,:,j),1,[])';
                    else
                        fid = reshape(FID_Mat,1,[])';
                    end
                    %Remove Points that are beyond a trajectory radius of 0.5
                    traj_recon_1 = traj_recon;
                    rad = sqrt(traj_recon_1(:,1).^2+traj_recon_1(:,2).^2+traj_recon_1(:,3).^2);
                    del_pts = find(rad > 0.5);
                    fid(del_pts) = [];
                    traj_recon_1(del_pts,:) = [];
                    
                    [EchoImages(:,:,:,i,j),Gridded_KSpace(:,:,:,i,j),Recon_Params] = base_noncart_recon(ImSize,fid,traj_recon_1);
                    %Loop order is b values on the inner, TEs on the outer
                    if strcmp(Method_Params.DiffusionYN,'Yes')
                        TE_index = ceil(i/Method_Params.Nbvalue);
                        b_index = mod(i,Method_Params.Nbvalue);% + (TE_index-1)*Method_Params.Nbvalue;
                        if mod(i,Method_Params.Nbvalue) == 0
                            b_index = b_index+Method_Params.Nbvalue;
                        end
                        FileName = ['Image for TE ' num2str(Method_Params.TE(TE_index)) 'and b value ' num2str(Method_Params.bvalues(b_index)) ' Repetition ' num2str(j)];
                    else
                        FileName = ['Image for TE ' num2str(Method_Params.TE(i)) ' Repetition ' num2str(j)];
                    end

                    Complete_FileName = fullfile(pwd,Dir_Name,FileName);
                    DicomWrite(Complete_FileName,squeeze(EchoImages(:,:,:,i,j)),Method_Params.Resolution)
                    %Write out Niftis that are formatted properly for autosegmentation
                    niftiwrite(abs(flip(rot90(squeeze(EchoImages(:,:,:,i,j)),3))),Complete_FileName)
                else
                    %Create Temporary variables to get the multi echo data and delete extra
                    %zeros before reshaping and reconning
                    EchoTrajTmp = squeeze(EchoTraj(:,:,:,i));
                    EchoFIDTmp = EchoFID(:,:,i);
                    indx0 = find(EchoFIDTmp(:,1) == 0); 
                    EchoTrajTmp(:,indx0,:) = [];
                    EchoFIDTmp(indx0,:) =[];
                    EchoTrajx = reshape(squeeze(EchoTrajTmp(1,:,:)),[],1);
                    EchoTrajy = reshape(squeeze(EchoTrajTmp(2,:,:)),[],1);
                    EchoTrajz = reshape(squeeze(EchoTrajTmp(3,:,:)),[],1);
                    EchoTrajC = [EchoTrajx,EchoTrajy,EchoTrajz];
                    EchoFIDC = reshape(EchoFIDTmp,[],1);
                    %Remove Points that are beyond a trajectory radius of 0.5
                    traj_recon_1 = EchoTrajC;
                    rad = sqrt(traj_recon_1(:,1).^2+traj_recon_1(:,2).^2+traj_recon_1(:,3).^2);
                    del_pts = find(rad > 0.5);
                    EchoFIDC(del_pts) = [];
                    traj_recon_1(del_pts,:) = [];
                    if ~isempty(EchoFIDC)
                        [EchoImages(:,:,:,i),Gridded_KSpace(:,:,:,i),Recon_Params] = base_noncart_recon(ImSize,EchoFIDC,traj_recon_1);   
                        FileName = ['Image for TE ' num2str(Method_Params.TE(i))];
                        Complete_FileName = fullfile(pwd,Dir_Name,FileName);
                        DicomWrite(Complete_FileName,squeeze(EchoImages(:,:,:,i)),Method_Params.Resolution)
                    end
                end
            end
        end  %Need to add keyholing to this, but for now, this should be fine.
    end
    %Add keyhole image reconstruction for xenon imaging sequences when not
    %doing VFA
    if contains(Method_Params.Nucleus,'129Xe') && strcmp(Method_Params.VFAYN,'No') && contains(Method_Params.Dims,'3D') && Method_Params.ProjPerTrig>1 && (~(strcmp(Method_Params.FlybackYN,'No') && Method_Params.NumTEs > 1)) && ~fast_recon
        PPT = Method_Params.ProjPerTrig;
        %I think that we're generally only going to be concerned with the
        %first FID -
        %Get FIDs and trajectories - need to account for both basic UTE and
        %Flyback cases
        if strcmp(Method_Params.FlybackYN,'No') 
            if ndims(FID_Mat)>3
                error('Apparently there are situations where you have to worry about FID_Mat having more than 3 Dimensions')
            end
            KeyFID = squeeze(FID_Mat(:,:,1));
            KeyTraj = traj_recon;
        else
            KeyFID = squeeze(EchoFID(:,:,1));
            KeyTraj = squeeze(EchoTraj(:,:,:,1));
            indx0 = find(EchoFIDTmp(:,1) == 0);
            KeyFID(:,indx0,:) = [];
            KeyTraj(indx0,:) =[];
            KeyTrajx = reshape(squeeze(KeyTraj(1,:,:)),[],1);
            KeyTrajy = reshape(squeeze(KeyTraj(2,:,:)),[],1);
            KeyTrajz = reshape(squeeze(KeyTraj(3,:,:)),[],1);
            KeyTraj = [KeyTrajx,KeyTrajy,KeyTrajz];
        end
        %Get Key Radius - Make it have a minimum radius of 1/4 of the
        %original number of points along the projection
        KeyRad = max(floor(size(KeyFID,1)/4),floor(sqrt(size(KeyFID,2)/4/pi/PPT)));
        %Generate Keyhole
        Keyhole = KeyFID;
        Keyhole(1:KeyRad,:) = NaN;
        KeyKSpace = zeros(size(KeyFID,1),size(KeyFID,2),PPT);
        for ii = 1:PPT
            KeyKSpace(:,:,ii) = Keyhole;
            keyk0 = abs(KeyFID(1,ii:PPT:size(KeyFID,2)));
            MeanKeyk0 = mean(keyk0);
            for jj = 1:size(KeyFID,2)
                KeyKSpace(:,jj,ii) = KeyKSpace(:,jj,ii)*MeanKeyk0/abs(KeyFID(1,jj));
            end
            KeyKSpace(1:KeyRad,ii:PPT:size(KeyFID,2),ii) = KeyFID(1:KeyRad,ii:PPT:size(KeyFID,2));
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
        %[M0,FA_Map] = key_image_fitting(KeyImages,Mask);
        %imslice(abs(M0),'Magnetization Based on Key Fitting')
        %imslice(FA_Map,'Flip Angle Map')
    else
        M0 = 0;
        FA_Map = 0;
    end
end
imslice(abs(EchoImages),'Reconstructed Image');


Recon.Images = EchoImages;
Recon.Scale_HF_Images = Scale_HF_Images;
%Recon.normbreath_Images = Norm_Breath_Images;
Recon.T1_corr_Images = T1_corr_Images;
Recon.Keys = KeyImages;
Recon.KSpace = Gridded_KSpace;
Recon.M0 = M0;
Recon.FA_Map = FA_Map;
Recon.Scale_HF_Keys = scale_KeyImages;
Recon.T1Corr_Keys = T1Corr_KeyImages;

Raw.FIDs = FIDs;
Raw.Traj = traj;
Raw.T1_fit_Corrected_FIDs = T1n_fid;


Recon_Params.traj_delay = traj_delay; %Make sure to write out what trajectory delay was used for recon
recon_def.Params = Recon_Params;
if fast_recon
    save(fullfile(pwd,Dir_Name,'Fast_Recon_WorkSpace.mat'),'Recon','Raw','Method_Params','Config_Params')
    save(fullfile(pwd,Dir_Name,'Images.mat'),'EchoImages')
else
    save(fullfile(pwd,Dir_Name,'Recon_WorkSpace.mat'),'Recon','Raw','Method_Params','Config_Params')
    save(fullfile(pwd,Dir_Name,'Images.mat'),'EchoImages')
end