function Config_Params = read_config_file_v3(path)

cd(path)
%Function to read in previously generated configuration file and pass back
%the values
%Pass the path where the file is located.
fileID = fopen('ConfigFile_V3.txt');

configRead=textscan(fileID,'%s','delimiter','\n');
configRead=configRead{1};

%% Since I am the one writing this file, I know where things are located 
%Just need to read them in and grab the values
Subject = configRead{1};
Subject(1:11) = [];
Config_Params.Subject = Subject;

MouseType = configRead{2};
MouseType(1:15) = [];
Config_Params.Strain = MouseType;

MouseSex = configRead{3};
MouseSex(1:12) = [];
Config_Params.Sex = MouseSex;

%Even though the mass is a number, I never need to use it as a number, so
%just keep it as a string
MouseMass = configRead{4};
MouseMass(1:13) = [];
Config_Params.Mass = MouseMass;

MouseDOB = configRead{5};
MouseDOB(1:12) = [];
Config_Params.DOB = MouseDOB;

BreathType = configRead{6};
BreathType(1:14) = [];
Config_Params.BreathType = BreathType;

VentBPM = configRead{7};
VentBPM(1:11) = [];
VentBPM = str2num(VentBPM);
Config_Params.VentBPM = VentBPM;

VentTV = configRead{8};
VentTV(1:10) = [];
VentTV = str2num(VentTV);
Config_Params.VentTV = VentTV/1000; %User gives TV in uL, so convert to mL here

VentInDur = configRead{9};
VentInDur(1:13) = [];
VentInDur = str2num(VentInDur);
Config_Params.VentInDur = VentInDur;

VentHoldDur = configRead{10};
VentHoldDur(1:15) = [];
VentHoldDur = str2num(VentHoldDur);
Config_Params.VentHoldDur = VentHoldDur;

VentO2Pct = configRead{11};
VentO2Pct(1:13) = [];
VentO2Pct = str2num(VentO2Pct);
Config_Params.VentO2Pct = VentO2Pct;

DiscardFromBeginning = configRead{12};
DiscardFromBeginning(1:24) = [];
DiscardFromBeginning = str2num(DiscardFromBeginning);
Config_Params.DiscardFromBeginning = DiscardFromBeginning;

DiscardFromEnd = configRead{13};
DiscardFromEnd(1:18) = [];
DiscardFromEnd = str2num(DiscardFromEnd);
Config_Params.DiscardFromEnd = DiscardFromEnd;

traj_delay = configRead{14};
traj_delay(1:13) = [];
traj_delay = str2num(traj_delay);
Config_Params.traj_delay = traj_delay;