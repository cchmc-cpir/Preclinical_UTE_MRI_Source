function [traj,Method_Params] = read_method_ute(path)
%% Function to Read a Bruker Method File and return any parameters of interest
%Return a structure containing all the data that we care about

%% Make sure we have the path we need
if nargin == 0
    [path]=uigetdir('C:\','Select Folder in Which Method File is Located');
    cd(path)
else
    cd(path)
end

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

%% Read Method File to get measured trajectories
locx = 0;
locy = 0;
locz = 0;
locend = 0;
%Read measured method file to get trajectories
fid=fopen(char(meas_method));
methodRead=textscan(fid,'%s','delimiter','\n');
methodRead=methodRead{1};
for index=1:size(methodRead,1)
    testStr=char(methodRead{index});
    if contains(testStr,'##$Method')
        Method=(testStr(11:end));
        Method_Params.Sequence = Method;
        Method_Params.SequenceName = Method;
    end
end
for index=1:size(methodRead,1)
    testStr=char(methodRead{index});
    if contains(testStr,'##$PVM_TrajKx')
        locx = index; %Index where X Trajectory begins 
        Length_X = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajKy')
        locy = index; %Index where Y Trajectory begins (and X ends)
        Length_Y = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajKz')
        locz = index; %Index where Z Trajectory begins (and Y ends)
        Length_Z = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajBx')
        loc_end = index; %Index where Z Trajectory ends
    end
    if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
        AcqShift=str2num(testStr(13:end));
    end
    if contains(testStr,'##$PVM_SpatDimEnum=') %Acquisition Time
        Dims = testStr(20:end);
    end
end
%Check length of Trajectories - if longer than 1, we have measured
%trajectories and should read them here. Otherwise move on
if Length_X > 1
    %Measured Trajectories should work the same for both Spiral and Radial,
    %so only need minimal changes here.
    Method_Params.Traj_Type = 'measured';
    %Read out Trajectory Measurements
    trajx = [];
    trajy = [];
    trajz = [];
    for index = (locx+1):(locy-1)
        readOutStr = methodRead{index};
        kx = str2num(readOutStr(1:end));
        trajx = [trajx kx];
    end
    for index = (locy+1):(locz-1)
        readOutStr = methodRead{index};
        ky = str2num(readOutStr(1:end));
        trajy = [trajy ky];
    end
    for index = (locz+1):(loc_end-1)
        readOutStr = methodRead{index};
        kz = str2num(readOutStr(1:end));
        trajz = [trajz kz];
    end
    trajx_shift_val = mean(trajx(1:AcqShift));
    trajy_shift_val = mean(trajy(1:AcqShift));
    trajz_shift_val = mean(trajz(1:AcqShift));
    
    trajx = trajx - trajx_shift_val; %These have a non-zero value at the beginning - get rid of that here
    trajy = trajy - trajy_shift_val;
    trajz = trajz - trajz_shift_val;
else
    Method_Params.Traj_Type = 'theoretical';
    fid=fopen(char(theo_method));
    methodRead=textscan(fid,'%s','delimiter','\n');
    methodRead=methodRead{1};
    ReadGrad = 0;
    PhaseGrad = 0;
    SliceGrad = 0;
    for index=1:size(methodRead,1)
        testStr=char(methodRead{index});
        if contains(testStr,'##$GradShape') %Theoretical Gradient Shape
            GradShapeStart = index;
            GradShapePts = str2num(testStr(14:end));
        end
        if contains(testStr,'##$GradRes') %Gradient Shape Resolution
            GradRes=str2num(testStr(12:end));
            Method_Params.GradRes = GradRes;
            GradShapeEnd = index;
        end
        if contains(testStr,'##$PVM_DigNp') %Number of points before Ramp Compensation
            NumPts=str2num(testStr(14:end));
        end
        if contains(testStr,'##$RampPoints') %Number of points for Ramp Compensation
            RampPoints = str2num(testStr(15:end));
        end
        if contains(testStr,'##$PVM_DigDw') %Dwell Time
            Dwell = str2num(testStr(14:end));
        end
        if contains(testStr,'##$ReadGrad')
            ReadGrad = str2num(testStr(13:end));
        end
        if contains(testStr,'##$PhaseGrad')
            PhaseGrad = str2num(testStr(14:end));
        end
        if contains(testStr,'##$SliceGrad')
            SliceGrad = str2num(testStr(14:end));
        end
        if contains(testStr,'##$PVM_GradCalConst')
            GradCalConst = str2num(testStr(21:end));
        end
    end
    GradShape = [];
    %PJN - Edit this to accomodate Flyback - need to add Flyback Pts to
    %method
    %Get the actual number of points
    Tot_Num_Pts = NumPts;
    %Get the gradient shape
    for index = (GradShapeStart+1):(GradShapeEnd-1)
        readOutStr = methodRead{index};
        gradval = str2num(readOutStr(1:end));
        GradShape = [GradShape gradval];
    end
    %Integrate the gradient
    GradTime = 0:GradRes:(GradRes*(length(GradShape)-1));
    TrajShape = cumtrapz(GradTime,GradShape);

    SampleTime = 0:Dwell:(Dwell*(Tot_Num_Pts-1));
    %Get these in Physical Units if the correct, updated version is
    %used
    if ReadGrad ~= 0
        TrajShapex = GradCalConst * TrajShape * ReadGrad/100 /1000;
        TrajShapey = GradCalConst * TrajShape * PhaseGrad/100 /1000;
        TrajShapez = GradCalConst * TrajShape * SliceGrad/100 /1000;
        trajx = interp1(GradTime,TrajShapex,SampleTime);
        trajy = interp1(GradTime,TrajShapey,SampleTime);
        trajz = interp1(GradTime,TrajShapez,SampleTime);
    else    
        %Get the trajectory points at actual sampling points
        traj = interp1(GradTime,TrajShape,SampleTime);
        trajx = traj;
        trajy = traj;
        trajz = traj;
    end
end

%Save the shapes in the method params structure
Method_Params.Base_Shape_x = trajx;
Method_Params.Base_Shape_y = trajy;
Method_Params.Base_Shape_z = trajz;

%% Read Parameters from Method File
fid=fopen(char(theo_method));
methodRead=textscan(fid,'%s','delimiter','\n');methodRead=methodRead{1};

Method_Params.ReadGrad = 0; %Line to make sure this gets set In versions before V3 of the code, this wasn't written to method file
for index=1:size(methodRead,1)
   	testStr=char(methodRead{index});
    if contains(testStr,'##OWNER=') %Scan Date, Time, Location
        DateTime = methodRead{index+1};
        ScanDate = DateTime(4:13);
        ScanTime = DateTime(15:22);
        Method_Params.ScanDate = ScanDate;
        Method_Params.ScanTime = ScanTime;
        FileLoc = methodRead{index+2};
        Method_Params.File = FileLoc;
    end
    if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
        AcqShift=str2num(testStr(13:end));
        EchoTimesEnd = index; %AcqShift comes after EchoTimes, so for a scan with lots of echo times, we want the end index
        Method_Params.AcqShift = AcqShift;
    end
    if contains(testStr,'##$NPro') %Number of Projections
        NPro=str2num(testStr(9:end));
        Method_Params.NPro = NPro;
    end
    if contains(testStr,'##$PVM_EchoTime=') %Echo Time
        TE=str2num(testStr(17:end));
        Method_Params.TE = TE;
    end
    if contains(testStr,'##$PVM_RepetitionTime=') %Repetition Time
        TR = str2num(testStr(23:end));
        Method_Params.TR = TR;
    end
    if contains(testStr,'##$ProUndersampling=') %Undersampling
        USamp = str2num(testStr(21:end));
        Method_Params.USamp = USamp;
    end
    if contains(testStr,'##$PVM_NAverages=') %Averages
        Averages = str2num(testStr(18:end));
        Method_Params.Averages = Averages;
    end
    if contains(testStr,'##$PVM_NRepetitions=') %Repetitions
        Repetitions = str2num(testStr(21:end));
        Method_Params.Repetitions = Repetitions;
    end
    if contains(testStr,'##$PVM_DigNp') %Number of points before Ramp Compensation
        NumPts=str2num(testStr(14:end));
        Method_Params.NPts = NumPts;
    end
    if contains(testStr,'##$RampPoints') %Number of points for Ramp Compensation
        RampPoints = str2num(testStr(15:end));
        Method_Params.RampPoints = RampPoints;
    end
    if contains(testStr,'##$PVM_DigDw') %Dwell Time
        Dwell = str2num(testStr(14:end));
        Method_Params.Dwell = Dwell;
    end
    if contains(testStr,'##$ExcPulse1Enum') %Pulse Shape
        PulseShape = testStr(18:end);
        Method_Params.PulseShape = PulseShape;
    end
    if contains(testStr,'##$PVM_Nucleus1Enum=') %Nucleus
        Nucleus = testStr(21:end);
        Method_Params.Nucleus = Nucleus;
    end
    if contains(testStr,'##$RefPowCh1=') %Reference Power
        RefPow = str2num(testStr(14:end));
        Method_Params.RefPow = RefPow;
    end
    if contains(testStr,'##$PVM_FrqWork=') %Working Frequency
        Freq = str2num(char(methodRead{index+1}));
        Freq = Freq(1);
        Method_Params.Frequency = Freq;
    end
    if contains(testStr,'##$PVM_Matrix=') %Matrix Size
        Matrix = str2num(char(methodRead{index+1}));
        Method_Params.MatrixSize = Matrix;
    end
    if contains(testStr,'##$PVM_EffSWh=') %Bandwidth
        Bandwidth = str2num(testStr(15:end));
        Method_Params.Bandwidth = Bandwidth;
    end
    if contains(testStr,'##$PVM_AcquisitionTime=') %Acquisition Time
        AcqTime = str2num(testStr(24:end));
        Method_Params.AcqTime = AcqTime;
    end
    if contains(testStr,'##$SpoilAmp=') %Spoiler Amplitude
        SpoilerAmp = str2num(testStr(15:end));
        Method_Params.SpoilerAmp = SpoilerAmp;
    end
    if contains(testStr,'##$Spoiling=') %Amount of spoiling
        Spoiling = str2num(testStr(16:end));
        Method_Params.Spoiling = Spoiling;
    end
    if contains(testStr,'##$SpoilDur=') %Spoiler Duration
        SpoilDur = str2num(testStr(13:end));
        Method_Params.SpoilDur = SpoilDur;
    end
    if contains(testStr,'##$RampTime=') %Ramp Time
        RampTime = str2num(testStr(13:end));
        Method_Params.RampTime = RampTime;
    end
    if contains(testStr,'##$DummyScans=') %Dummy Scans
        NumDummies = str2num(testStr(19:end));
        Method_Params.Dummies = NumDummies;
    end
    if contains(testStr,'##$PVM_TriggerModule=') %Trigger on or Off
        Trigger = testStr(22:end);
        Method_Params.Trigger = Trigger;
    end
    if contains(testStr,'##$ExcPulse1=') %Excitation Pulse Parameters
        PulseParams = testStr(14:end);
        PulseIndex = 1;
        while ~contains(char(methodRead{index+PulseIndex}),'##$ExcPulse1Ampl=')
            PulseParams = [PulseParams, char(methodRead{index+PulseIndex})];
            PulseIndex = PulseIndex+1;
        end
        Method_Params.PulseParams = PulseParams;
    end
    if contains(testStr,'##$PVM_SpatResol') %Resolution
        Resolution = str2num(methodRead{index+1});
        Method_Params.Resolution = Resolution;
    end
    if contains(testStr,'##$PVM_Fov=') %Field of View
        FOV = str2num(methodRead{index+1});
        Method_Params.FOV = FOV;
    end    
    if contains(testStr,'##$GradRes') %Gradient Shape Resolution
        GradRes=str2num(testStr(12:end));
        Method_Params.GradRes = GradRes;
    end
    if contains(testStr,'##$PVM_SpatDimEnum=') %Number of Dimensions
        Dims = testStr(20:end);
        Method_Params.Dims = Dims;
    end
    if contains(testStr,'##$ReadGrad')
        Method_Params.ReadGrad = str2num(testStr(13:end));
    end
    if contains(testStr,'##$PhaseGrad')
        Method_Params.PhaseGrad = str2num(testStr(14:end));
    end
    if contains(testStr,'##$SliceGrad')
        Method_Params.SliceGrad = str2num(testStr(14:end));
    end
    if contains(testStr,'##$NumTEs')
        Method_Params.NumTEs = str2num(testStr(11:end));
    end
    if contains(testStr,'##$EchoTimes=') %Echo Time
        EchoStartind = index;
    end
    if contains(testStr,'$$ @vis=') %Echo Time
        EchoEndind = index;
    end
end
%% Get EchoTime values
TE = [];
for i = (EchoStartind+1):(EchoEndind-1)
    TEstr=char(methodRead{i});
    TE = [TE str2num(TEstr)];
end
Method_Params.TE = TE;
%% Rotate Projections according to golden means acquisition
phi1 = 0.46557123;
phi2 = 0.6823278;
gs = 1;
gr = 1;
gp = 1;

r = zeros(1,NPro);
p = zeros(1,NPro);
s = zeros(1,NPro);
%Rotation code from Jinbang's UTE sequence
for i = 0:(NPro-1)
    kz = (i*phi1-floor(i*phi1))*2-1;
    ts = kz*gs;
    alpha = (i*phi2-floor(i*phi2))*2*pi;
    tr = sqrt(1-kz*kz)*cos(alpha)*gr;
    tp = sqrt(1-kz*kz)*sin(alpha)*gp;
    r(i+1) = tr;
    p(i+1) = tp;
    s(i+1) = -ts;
end

trajx_rot = zeros(length(trajx),NPro);
trajy_rot = trajx_rot;
trajz_rot = trajx_rot;
for i = 1:NPro
    trajx_rot(:,i) = r(i)*trajx';
    trajy_rot(:,i) = p(i)*trajy';
    trajz_rot(:,i) = s(i)*trajz';
end

traj = cat(3,trajx_rot,trajy_rot,trajz_rot);
traj = permute(traj,[3,1,2]);
%Trajectories are in physical Units at this point: Convert to Pixel
%Units
if Method_Params.ReadGrad == 0 && strcmp(Method_Params.Traj_Type, 'theoretical')
    radius = squeeze(sqrt(traj(1,:,:).^2+traj(2,:,:).^2+traj(3,:,:).^2));   
    maxradius = max(max(radius));
    traj = traj/maxradius/2;
    %traj(3,:,:) = traj(3,:,:) /(FOV(3)/FOV(1)); %Need to make the slice field of view behave itself
else
    kFOV_desired = 1./Method_Params.Resolution;
    kMax_desired = kFOV_desired/2;
    max_k = max(kMax_desired);
    traj = traj/max_k/2;

    traj(3,:,:) = traj(3,:,:) * (FOV(3)/FOV(1)); %PJN - Add to get slice field of view correct... seems to have issues.
end


