function recon(path,fast_recon)

%set fast_recon to "true" (or "1") to skip keyholing and other slow steps so that you can quickly get a reconstructed image 
if nargin == 0
    [path]=uigetdir('C:\','Select Folder in Which Data is Stored');
    fast_recon = 0;
    cd(path)
elseif nargin == 1
    fast_recon = 0;
    cd(path)
else
    cd(path)
end

%Make sure that folder contains a method file and a fid file
%Sometimes I rename method file things like measured_Method, so make sure
%to accomodate this sort of thing
Method_Files = dir('*ethod*');
fid_File = dir('fid*');
if isempty(Method_Files)
    error('Unable to find a method file in the specified path')
elseif isempty(fid_File)
    error('Unable to find a fid file in the specified path')
elseif length(fid_File) > 1
    error('Supply a path containing only one fid file')
elseif length(Method_Files) > 2
    error('Supply a path containing a maximum of 2 method files')
end

%% Read in Configuration File
%If it doesn't exist, create one
ConfigPresent = exist(fullfile(path,'ConfigFile_V4.txt'),'file');
if ConfigPresent == 0
    gen_config_file_v4(path)
end
Config_Params = read_config_file_v4(path);

%% Read Method File
disp('Reading Method File')
[traj,Method_Params] = read_method(path);

%If we are doing a phantom, want to pretend that we are triggered so that
%we don't try to retrogate (which would be pointless)
if strcmpi(Config_Params.Subject,'Phantom') %non case sensitive string compare
    Method_Params.Trigger = 'On';
end

%% Pass to required Reconstruction code
%If - else should be fine here, because read_method will break if it can't
%determine whether we are spiral or radial.
if strcmp(Method_Params.Sequence,'Radial')
    [Dir_Name,recon_def] = radial_recon(path,traj,Method_Params,Config_Params,fast_recon);
else
    [Dir_Name,recon_def] = spiral_recon(path,traj,Method_Params,Config_Params,fast_recon);
end

%% Write out Recon Details
fileID = fopen(fullfile(pwd,Dir_Name,'Recon_Details.txt'),'w');

if isunix() 
    UserName = getenv('USER'); 
else 
    UserName = getenv('username'); 
end

fprintf(fileID,[Config_Params.Subject ' imaged on ' Method_Params.ScanDate ' at ' Method_Params.ScanTime '\n']);
fprintf(fileID,['Imaged using the sequence ' Method_Params.SequenceName '\n']);
fprintf(fileID,['Reconstructed by ' UserName ' on ' datestr(now) ' using ' recon_def.function '\n']);
fprintf(fileID,'\n');
fprintf(fileID,['*****************************************************************' '\n']);
fprintf(fileID,'\n');
fprintf(fileID,['***Significant General Imaging Parameters***' '\n']);

if Method_Params.NSlices > 1 || contains(Method_Params.Dims,'2D')
    MatSizeStr = [num2str(Method_Params.MatrixSize(1)) 'x' num2str(Method_Params.MatrixSize(2)) 'x' num2str(Method_Params.NSlices)];
    FOVStr = [num2str(Method_Params.FOV(1)) 'x' num2str(Method_Params.FOV(2)) 'x' num2str(Method_Params.SliceThick*Method_Params.NSlices)];
    ResStr = [num2str(Method_Params.Resolution(1)) 'x' num2str(Method_Params.Resolution(2)) 'x' num2str(Method_Params.SliceThick)];
else
    MatSizeStr = [num2str(Method_Params.MatrixSize(1)) 'x' num2str(Method_Params.MatrixSize(2)) 'x' num2str(Method_Params.MatrixSize(3))];
    FOVStr = [num2str(Method_Params.FOV(1)) 'x' num2str(Method_Params.FOV(2)) 'x' num2str(Method_Params.FOV(3))];
    ResStr = [num2str(Method_Params.Resolution(1)) 'x' num2str(Method_Params.Resolution(2)) 'x' num2str(Method_Params.Resolution(3))];
end
fprintf(fileID,['Nucleus = ' Method_Params.Nucleus '\n']);
fprintf(fileID,['Image Size = ' MatSizeStr '\n']);
fprintf(fileID,['FOV = ' FOVStr '\n']);
fprintf(fileID,['Resolution = ' ResStr '\n']);
fprintf(fileID,['TR = ' num2str(Method_Params.TR) '\n']);
fprintf(fileID,['TE = ' num2str(Method_Params.TE) '\n']);
fprintf(fileID,['Acquisition Time = ' num2str(Method_Params.AcqTime) '\n']);
fprintf(fileID,['Bandwidth = ' num2str(Method_Params.Bandwidth) '\n']);
fprintf(fileID,['Number of Acquisition Shift Points = ' num2str(Method_Params.AcqShift) '\n']);
fprintf(fileID,['Averages = ' num2str(Method_Params.Averages) '\n']);
fprintf(fileID,['Repetitions = ' num2str(Method_Params.Repetitions) '\n']);
fprintf(fileID,['Excitation Pulse Parameters = ' Method_Params.PulseParams '\n']);
fprintf(fileID,['Slab Selection = ' Method_Params.SlabYN '\n']);
if strcmp(Method_Params.SlabYN,'Yes')
    fprintf(fileID,['Slab Thickness = ' num2str(Method_Params.SlabThickness) '\n']);
end
fprintf(fileID,['Diffusion Encoding = ' Method_Params.DiffusionYN '\n']);
if strcmp(Method_Params.DiffusionYN,'Yes')
    fprintf(fileID,['Number of b values = ' num2str(Method_Params.Nbvalue) '\n']);
    fprintf(fileID,['b values = ' num2str(Method_Params.bvalues) '\n']);
    fprintf(fileID,['delta = ' num2str(Method_Params.delta) '\n']);
    fprintf(fileID,['Delta = ' num2str(Method_Params.Delta) '\n']);
end
fprintf(fileID,'\n');
fprintf(fileID,['*****************************************************************' '\n']);
fprintf(fileID,'\n');
fprintf(fileID,['***Method Specific Imaging Parameters***' '\n']);
if strcmp(Method_Params.Sequence,'Radial')
    fprintf(fileID,['Measure EchoTimes with Flyback = ' Method_Params.FlybackYN '\n']);
    fprintf(fileID,['Number of Projections = ' num2str(Method_Params.NPro) '\n']);
    fprintf(fileID,['Sampling Percentage = ' num2str((1/Method_Params.USamp)*100) '%%' '\n']);
    fprintf(fileID,['Play out Rewinder = ' Method_Params.Rewinder '\n']);
else
    fprintf(fileID,['Number of Spiral Projections = ' num2str(Method_Params.NPro) '\n']);
    fprintf(fileID,['Spiral Interleaves Parameter = ' num2str(Method_Params.Interleaves) '\n']);
    fprintf(fileID,['Spiral Readout Time = ' num2str(Method_Params.SpROTime) '\n']);
    if Method_Params.NSlices == 1
        fprintf(fileID,['Number of Hubs = ' num2str(Method_Params.Hubs) '\n']);
    else
        fprintf(fileID,['Number of Spiral Slices = ' num2str(Method_Params.NSlices) '\n']);
    end
end
fprintf(fileID,'\n');
fprintf(fileID,['*****************************************************************' '\n']);
fprintf(fileID,'\n');
fprintf(fileID,['***Triggering Parameters***' '\n']);
fprintf(fileID,['Trigger Module = ' Method_Params.Trigger '\n']);
if strcmp(Method_Params.Trigger,'On') && strcmp(Method_Params.Nucleus,'<129Xe>')
    if contains(Method_Params.Sequence,'Radial')
        fprintf(fileID,['Projections Per Trigger = ' num2str(Method_Params.ProjPerTrig) '\n']);
        fprintf(fileID,['Variable Flip Angle = ' Method_Params.VFAYN '\n']);
    end
end
fprintf(fileID,'\n');
fprintf(fileID,['*****************************************************************' '\n']);
fprintf(fileID,'\n');
fprintf(fileID,['***Reconstruction Parameters***' '\n']);
fprintf(fileID,['Trajectory Delay = ' num2str(recon_def.Params.traj_delay) ' (*Dwell)' '\n']);
fprintf(fileID,['Recon Function = ' recon_def.Params.Recon_Func '\n']);
fprintf(fileID,['Kernel Sharpness = ' num2str(recon_def.Params.KernelSharpness) '\n']);
fprintf(fileID,['Kernel Extent = ' num2str(recon_def.Params.KernelExtent) '\n']);
fprintf(fileID,['Overgridding = ' num2str(recon_def.Params.OverGridding) '\n']);
fprintf(fileID,['Number of Density Compensation Iterations = ' num2str(recon_def.Params.DCIterations) '\n']);
fprintf(fileID,['Number of Threads = ' num2str(recon_def.Params.nThreads) '\n']);
if recon_def.Params.Deapodize
    Deap = 'True';
else
    Deap = 'False';
end
fprintf(fileID,['Deapodize Image = ' Deap '\n']);
fclose(fileID);
%% Close all open files
fclose('all');