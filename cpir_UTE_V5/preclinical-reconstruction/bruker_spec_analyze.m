function bruker_spec_analyze(path)

if nargin == 0
    path = uigetdir;
end
cd(path)
%Call function to load in Bruker FID file - make sure the method file is in
%the same directory
raw = nspect_load('rawdata.job0');

%Look in method file for the pulse powers, pulse length, working frequency and bandwidth 
methodfile = 'method';
fid=fopen(char(methodfile));
methodRead=textscan(fid,'%s','delimiter','\n');methodRead=methodRead{1};
for index=1:size(methodRead,1)
    testStr=char(methodRead{index});
    if length(testStr)>13
		if strcmp(testStr(1:14),'##$PVM_SpecSWH')==1
            readOutStr=methodRead{index+1};
            Method_Params.BW = str2num(readOutStr);
		end
    end
    if length(testStr)>11
		if strcmp(testStr(1:12),'##$PVM_DigDw')==1
			readOutStr=methodRead{index};
			Method_Params.Dwell=str2num(readOutStr(14:length(readOutStr)));
		end
    end
    if contains(testStr,'##$ExcPulse1Enum') %Pulse Shape
        PulseShape = testStr(18:end);
        Method_Params.PulseShape = PulseShape;
    end
    if contains(testStr,'##$PVM_FrqWork=') %Working Frequency
        Freq = str2num(char(methodRead{index+1}));
        Freq = Freq(1);
        Method_Params.Frequency = Freq;
    end
    if contains(testStr,'##$PVM_RepetitionTime=') %Repetition Time
        TR = str2num(testStr(23:end));
        Method_Params.TR = TR;
    end
    if contains(testStr,'##$PVM_NRepetitions=') %Repetitions
        Repetitions = str2num(testStr(21:end));
        Method_Params.Repetitions = Repetitions;
    end
    if contains(testStr,'##$PVM_SpecMatrix=') %Repetitions
        Num_Spec_Pts=str2num(methodRead{index+1});
        Method_Params.NPts = Num_Spec_Pts;
    end
end

%Separate FIDs into individual powers
%Powers = str2num(methodReadPowers);
%CentFreq = str2num(Method_Params.Frequency);
Bandwidth = Method_Params.BW;
Dwell = Method_Params.Dwell;

%Get the number of points per spectrum
%NPts = length(raw)/Method_Params.Repetitions;

time_domain = 0:Dwell:(Method_Params.NPts-1)*Dwell;
time_domain = time_domain/1000; %get this in units that work

Method_Params.Repetitions = length(raw)/Method_Params.NPts;

%Preallocate memory
FIDs = zeros(Method_Params.NPts,[]);
FIDs = complex(FIDs,0);
FFT_spec = FIDs;
x = 1:Method_Params.NPts;
IntInt = zeros(1,Method_Params.Repetitions);
Peak = zeros(1,Method_Params.Repetitions);
Time_Axis = zeros(1,Method_Params.NPts);
%Loop through FIDs and get FFT
for i = 1:Method_Params.Repetitions
    FIDs(:,i) = raw((Method_Params.NPts*(i-1)+1):(Method_Params.NPts*i));
    FFT_spec(:,i)=abs(fftshift(fft(FIDs(:,i))));
    Time_Axis(i,:) = ones(1,Method_Params.NPts)*(i-1)*Method_Params.TR;
end
%Normalize and Display Spectra
FFT_spec = FFT_spec/max(max(FFT_spec));
freq_domain = (linspace(-Bandwidth,Bandwidth,Method_Params.NPts))';
freq_domain = freq_domain*1e-6 + Method_Params.Frequency;
figure('Name','Spectra')
hold on
for i = 1:Method_Params.Repetitions
    plot3(freq_domain,Time_Axis(i,:),FFT_spec(:,i))
end
set(gca,'Xdir','reverse')
view(21,16);
xlabel('Frequency (Hz)')
ylabel('Time (ms)')
zlabel('Signal Intensity (Arb. Units)')

save('Spectra Analysis Workspace.mat')

    
    