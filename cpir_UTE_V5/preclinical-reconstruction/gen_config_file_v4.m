function gen_config_file_v4(path)

%Function to write out a configuration file for xenon images of mice -
%write to a text file things like date of birth, scan type, mouse type,
%etc...

%% Start by loading in FIDs so that I can see how many rays that need to be kept
%% Make sure we have the path we need
if nargin == 0
    [path]=uigetdir('C:\','Select Folder in Which Method File is Located');
    cd(path)
else
    cd(path)
end

%I want to pull info from old config files if they already exist - This is
%pretty much only useful for me, right now... but still
ConfigPresent = exist(fullfile(path,'ConfigFile_V3.txt'),'file');
if ConfigPresent
    Config_Params = read_config_file_v3(path);
    ScanType = Config_Params.Subject;
    Mouse_Strain = Config_Params.Strain;
    Mouse_Sex = Config_Params.Sex;
	Mouse_Mass = Config_Params.Mass;
	Mouse_DOB = Config_Params.DOB;
	BreathType = Config_Params.BreathType;
	VentBPM = num2str(Config_Params.VentBPM);
	VentTV = num2str(Config_Params.VentTV*1000);
	VentInDur = num2str(Config_Params.VentInDur);
	VentHoldDur = num2str(Config_Params.VentHoldDur);
	VentO2Pct = num2str(Config_Params.VentO2Pct);
    From_Begin = num2str(Config_Params.DiscardFromBeginning);
	From_End = num2str(Config_Params.DiscardFromEnd);
    Traj_Delay = num2str(Config_Params.traj_delay);
    if ~strcmp(BreathType,'N/A')
        WIO = questdlg('Was This a Washin/Washout Image?','Wash in Wash out','Yes','No','No');
        if strcmp(WIO,'Yes')
            prompt = {'Wash In Breaths','Wash Out Breaths'};
            dims = [1 100];
            definput = {'5','5'};
            WIO_Params = inputdlg(prompt,name,dims,definput);
            Washin = VentParams{1,1};
            Washout = VentParams{2,1};
        else
            Washin = '0';
            Washout = '0';
        end
    else
        WIO = 'No';
        Washin = '0';
        Washout = '0';
    end
else
    
    %% Get Method File/Files
    methodfiles = dir('*ethod*');
    % If there are multiple method files, we need to separate into measured and
    % theoretical
    if length(methodfiles)>1
        if methodfiles(1).bytes > methodfiles(2).bytes
            theo_method = methodfiles(2).name;
            meas_method = methodfiles(1).name;
        elseif methodfiles(2).bytes > methodfiles(1).bytes
            theo_method = methodfiles(1).name;
            meas_method = methodfiles(2).name;
        else
            error('Unable to differentiate measured and theoretical method file')
        end
    else
        method_file = methodfiles.name;
        theo_method = method_file;
        meas_method = method_file;
    end

    fid=fopen(char(theo_method));
    methodRead=textscan(fid,'%s','delimiter','\n');
    methodRead=methodRead{1};
    for index=1:size(methodRead,1)
        testStr=char(methodRead{index});
        if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
            AcqShift=str2num(testStr(13:end));
        end
        if contains(testStr,'##$NPro') %Number of Projections
            NPro=str2num(testStr(9:end));
        end
        if contains(testStr,'##$NumTEs=') %NumTEs
            NumTEs=str2num(testStr(11:end));
        end
        if contains(testStr,'##$FlybackYN') %Flyback Yes or No
            FlybackYN = testStr(14:end);
        end
        if contains(testStr,'##$Method')
            Method=(testStr(11:end));
        end
        if contains(testStr,'##$DiffusionYN=') %Diffusion Encoding YN
            DiffusionYN = testStr(16:end);
        end    
        if contains(testStr,'##$Nbvalue=') %Diffusion Encoding YN
            Nbvalue = str2num(testStr(12:end));
        end  
        if contains(testStr,'##$PVM_SPackArrNSlices=') %Undersampling Parameters
            NSlices = str2num(methodRead{index+1});
        end
        if contains(testStr,'##$PVM_NRepetitions=') %Repetitions
            Repetitions = str2num(testStr(21:end));
        end
        if contains(testStr,'##$PVM_Nucleus1Enum=') %Nucleus
            Nucleus = testStr(21:end);
        end
    end

    if contains(Method,'ute','IgnoreCase',true)
        Seq = 'Radial';
    elseif contains(Method,'spiral','IgnoreCase',true)
        Seq = 'Spiral';
    else
        error('Cannot Determine Sequence Type')
    end


    FIDs = Bruker_Load('fid');
    if strcmp(Seq,'Radial')
        if strcmp(FlybackYN,'Yes')
            FIDs = reshape(FIDs,[],NPro*Repetitions);
            Nk0 = NPro*Repetitions;
        elseif strcmp(DiffusionYN,'Yes')
            FIDs = reshape(FIDs,[],NPro*NumTEs*Nbvalue*Repetitions);
            Nk0 = NPro*NumTEs*Nbvalue*Repetitions;
        else
            FIDs = reshape(FIDs,[],NPro*NumTEs*Repetitions);
            Nk0 = NPro*NumTEs*Repetitions;
        end
    else
       % if strcmp(DiffusionYN,'No')
            FIDs = reshape(FIDs,[],NPro*NSlices*Repetitions);
            Nk0 = NPro*NSlices*Repetitions;
     %   else
     %       FIDs = reshape(FIDs,[],NPro*NSlices*Nbvalue*Repetitions);
     %       Nk0 = NPro*NSlices*Repetitions*Nbvalue;
     %   end
    end

    %Find the number of zeros used in zerofilling
    numZeroFill = find(FIDs(:,1)==0);
    numZeroFill(numZeroFill < 10) = [];
    %Kill the zero-filled FID points
    FIDs(numZeroFill,:)=[];
    %Delete those points from the FIDs
    FIDs(1:AcqShift,:) = [];

    %% Get parameters from the scan
    ScanType = questdlg('Did you image a mouse or phantom?','Scan Type','Mouse','Phantom','Mouse'); 
    if strcmp(ScanType,'Mouse')
        %get sex
        Mouse_Sex = questdlg('Was the mouse male or female?','Mouse Sex','Male','Female','Male'); 
        %get DOB
        Mouse_DOB = datestr(uigetdate);
        Mouse_DOB = Mouse_DOB(1:11);
        close;
        prompt = {'Mouse Strain:','Mouse Mass'};
        name = 'Mouse Details';
        dims = [1 100];
        definput = {'C57BL-6J','25'};
        MouseDet = inputdlg(prompt,name,dims,definput);
        Mouse_Strain = MouseDet{1,1};
        Mouse_Mass = MouseDet{2,1};
    else
        Mouse_Sex = 'N/A';
        Mouse_DOB = 'N/A';
        Mouse_Strain = 'N/A';
        Mouse_Mass = 'N/A';
    end

    %If doing Xenon, get ventilator settings
    if contains(Nucleus,'129Xe')
        name = 'Ventilator Settings';
        prompt = {'Imaging Location (i.e. FRC)','Ventilator BPM','Ventilator Tidal Volume (uL)','Ventilator Inhalation Dur (ms)','Ventilator Breathhold Dur (ms)','O2 Percentage'};
        dims = [1 100];
        definput = {'End Inspiration','80','250','200','200','21'};
        VentParams = inputdlg(prompt,name,dims,definput);
        BreathType = VentParams{1,1};
        VentBPM = VentParams{2,1};
        VentTV = VentParams{3,1};
        VentInDur = VentParams{4,1};
        VentHoldDur = VentParams{5,1};
        VentO2Pct = VentParams{6,1};
        WIO = questdlg('Was This a Washin/Washout Image?','Wash in Wash out','Yes','No','No');
        if strcmp(WIO,'Yes')
            prompt = {'Wash In Breaths','Wash Out Breaths'};
            dims = [1 100];
            definput = {'5','5'};
            WIO_Params = inputdlg(prompt,name,dims,definput);
            Washin = WIO_Params{1,1};
            Washout = WIO_Params{2,1};
        else
            Washin = '0';
            Washout = '0';
        end
    else
        BreathType = 'N/A';
        VentBPM = '0';
        VentTV = '0';
        VentInDur = '0';
        VentHoldDur = '0';
        VentO2Pct = '0';
        WIO = 'No';
        Washin = '0';
        Washout = '0';
    end

    %Now, display k0 so that the user can throw away points if they so desire
    k0figA = figure('Name','First Point of FIDs','units','normalized','outerposition',[0 0.04 1 .96]);
    plot(1:Nk0,squeeze(abs(FIDs(1,:))),'k-*')

    Throw_Away = questdlg('Do any points need to be discarded?','Discard Points','Yes','No','No'); 
    if strcmp(Throw_Away,'Yes')
        prompt = {'Discard from Beginning','Discard From End'};
        name = 'Discarding Points';
        dims = [1 100];
        definput = {'0','0'};
        Discard_Pts = inputdlg(prompt,name,dims,definput);
        From_Begin = Discard_Pts{1,1};
        From_End = Discard_Pts{2,1};
    else
        From_Begin = '0';
        From_End = '0';
    end

    prompt = {'Trajectory Delay (In Multiples of Dwell Time)'};
    name = 'Trajectory Delay';
    dims = [1 100];
    definput = {'0'};
    traj_delay = inputdlg(prompt,name,dims,definput);
    Traj_Delay = traj_delay{1,1};

    close;
end
%% Write out parameters to text file
fileID = fopen('ConfigFile_V4.txt','w');
fprintf(fileID,['##$Subject=' ScanType '\n']);
fprintf(fileID,['##$MouseStrain=' Mouse_Strain '\n']);
fprintf(fileID,['##$MouseSex=' Mouse_Sex '\n']);
fprintf(fileID,['##$MouseMass=' Mouse_Mass '\n']);
fprintf(fileID,['##$MouseDOB=' Mouse_DOB '\n']);
fprintf(fileID,['##$BreathType=' BreathType '\n']);
fprintf(fileID,['##$VentBPM=' VentBPM '\n']);
fprintf(fileID,['##$VentTV=' VentTV '\n']);
fprintf(fileID,['##$VentInDur=' VentInDur '\n']);
fprintf(fileID,['##$VentHoldDur=' VentHoldDur '\n']);
fprintf(fileID,['##$VentO2Pct=' VentO2Pct '\n']);
fprintf(fileID,['##$Washinout=' WIO '\n']);
fprintf(fileID,['##$NWashin=' Washin '\n']);
fprintf(fileID,['##$NWashout=' Washout '\n']);
fprintf(fileID,['##$DiscardFromBeginning=' From_Begin '\n']);
fprintf(fileID,['##$DiscardFromEnd=' From_End '\n']);
fprintf(fileID,['##$TrajDelay=' Traj_Delay '\n']);
fclose(fileID);

