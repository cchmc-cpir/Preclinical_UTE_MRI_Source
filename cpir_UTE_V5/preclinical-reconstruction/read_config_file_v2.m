function Config_Params = read_config_file_v2(path)

cd(path)
%Function to read in previously generated configuration file and pass back
%the values
%Pass the path where the file is located.
fileID = fopen('ConfigFile_V2.txt');

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
Config_Params.TriggerInfo = BreathType;

KeepRays = configRead{7};
KeepRays(1:12) = [];
KeepRays = str2num(KeepRays);
Config_Params.UseProj = KeepRays;

traj_delay = configRead{8};
traj_delay(1:13) = [];
traj_delay = str2num(traj_delay);
Config_Params.traj_delay = traj_delay;