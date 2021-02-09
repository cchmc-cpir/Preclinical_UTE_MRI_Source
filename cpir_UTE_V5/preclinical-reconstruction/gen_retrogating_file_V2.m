%function gen_retrogating_file_V2(fid,number,Method_Params)

%Function to retrospectively gate images based on Jinbang's and Ian's code
%Special Thanks to 
%Jinbang Guo 
%Ian Stecker

%Pass the FID file - No direct output - it just writes out a file
%containing a vector with the end expiration indices and the end
%inspiration indices

%Also pass a number to name the file for the case of multiple TE or b value
%retrogating
if nargin == 1
    number = 1;
end

%% display k0 points starting halfway through the scan
dispstart = round(size(fid,2)/2);

%It's useful to have different limits for insp and expiration - we can be
%more discerning in the expiration points we select since there's more of
%them
%Set some initial values to test
eLL_val = -30;
eUL_val = 30;

iLL_val = -30;
iUL_val = 30;

exp_thres = 0.5;
insp_thres = 2;

NPro = size(fid,2);
%We're going to check with the user to see if they like the gating
happy_gate = 0;

PorM = 'm';

%% Repeat this process until the user approves
while happy_gate == 0
    if strcmp(PorM,'m')
        k0 = squeeze(abs(fid(1,:)));
    else
        k0 = squeeze(angle(fid(1,:)));
    end
    smk0 = smooth(k0,30);
    %% Start by displaying k0 magnitude and the smoothed version
    h = figure('Name','Retrospective Gating');
    set(h,'Units','Normalized','Position',[.05 .05 .9 .9],'Color','w')
    subplot(2,3,1);
    plot(k0, '.k', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'k', 'linewidth', 2);
    
    
    
    if strcmp(PorM,'m')
        title('Gating off Magnitude')
    else
        title('Gating off Phase')
    end
    xlim([dispstart dispstart+1000]);
    ylim([(min(k0) - min(k0)*0.001) (max(k0) + max(k0)*0.004)]);
    % xlabel('Radial Views', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('Signal Intensity', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');

    subplot(2,3,2);
    plot(smk0, '-', 'linewidth', 3, 'color',...
        [105/256, 105/256, 105/256]);

    ylabel('Signal Intensity', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');
    xlim([dispstart dispstart+1000]);
    ylim([(min(k0) - min(k0)*0.001) (max(k0) + max(k0)*0.004)]);

    %% Plot the first derivative
    Magnitude1stDer = diff(smk0, 1);
    Magnitude1stDerSmooth = smooth(Magnitude1stDer, 30);

    subplot(2,3,3);
%     rectangle('position', [dispstart LL_val 1000 -LL_val+UL_val], 'facecolor',...
%            [245/256 230/256 245/256], 'edgecolor', 'w'); hold on;
       
    rectangle('position', [dispstart iLL_val 1000 -iLL_val+iUL_val], 'facecolor',...
           'b', 'edgecolor', 'w'); hold on;

    rectangle('position', [dispstart eLL_val 1000 -eLL_val+eUL_val], 'facecolor',...
           'r', 'edgecolor', 'w'); hold on;

    plot(Magnitude1stDerSmooth, '-k', 'linewidth', 3);

%     line([0 NPro], [iLL_val iLL_val], 'color', [175/256, 60/256, 155/256],...
%         'linewidth', 4, 'linestyle', '--');
% 
%     line([0 NPro], [iUL_val iUL_val], 'color', [175/256, 60/256, 155/256],...
%         'linewidth', 4, 'linestyle', '--'); hold off; 
    
    line([0 NPro], [iLL_val iLL_val], 'color', 'b',...
        'linewidth', 4, 'linestyle', '--');

    line([0 NPro], [iUL_val iUL_val], 'color', 'b',...
        'linewidth', 4, 'linestyle', '--'); hold off; 
    
    line([0 NPro], [eLL_val eLL_val], 'color', 'r',...
        'linewidth', 4, 'linestyle', '--');

    line([0 NPro], [eUL_val eUL_val], 'color', 'r',...
        'linewidth', 4, 'linestyle', '--'); hold off; 

    ylabel('Signal Intensity / Projection', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');
    xlim([dispstart dispstart+1000]);
%     ylim([-(round(max(Magnitude1stDerSmooth(50:(end-50))) + max(Magnitude1stDerSmooth(50:(end-50)))*.13,-1))...
%         (round(max(Magnitude1stDerSmooth(50:(end-50))) + max(Magnitude1stDerSmooth(50:(end-50)))*.13,-1))]);
    
    ylim([-((max(Magnitude1stDerSmooth(50:(end-50)))))*1.1 ...
        ((max(Magnitude1stDerSmooth(50:(end-50)))))*1.1 ]);

    %I'm going to force this always to get both end inspiration and end
    %expiration always.

    Magnitude2ndDer = diff(Magnitude1stDerSmooth, 1);
    Magnitude2ndDerSmooth = smooth(Magnitude2ndDer, 30);
%% Plot second derivative
    subplot(2,3,4);
    rectangle('position', [dispstart -max(Magnitude2ndDerSmooth)/2 1000....
        (max(Magnitude2ndDerSmooth)/2 + exp_thres)],...
        'facecolor', [243/255, 151/255, 151/255], 'edgecolor', 'w'); hold on;

    rectangle('position', [dispstart insp_thres 1000 max(Magnitude2ndDerSmooth)/2],...
        'facecolor', [171/255, 209/255, 255/255], 'edgecolor', 'w');

    plot(Magnitude2ndDerSmooth, '-k', 'linewidth', 3) 
    %'color', [105/256, 105/256, 105/256]);

    line([0 NPro], [exp_thres exp_thres], 'color', 'r',...
        'linewidth', 4, 'linestyle', '--');

    line([0 NPro], [insp_thres insp_thres], 'color', 'b',...
        'linewidth', 4, 'linestyle', '--');
    
    ylabel('Signal Intensity / Projection^2', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');
    xlabel('Radial Projection', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');
    
    xlim([dispstart dispstart+1000]);
%     ylim([-(round(max(Magnitude2ndDerSmooth(50:(end-50))) + max(Magnitude2ndDerSmooth(50:(end-50)))*.2, 1))...
%         (round(max(Magnitude2ndDerSmooth(50:(end-50))) + max(Magnitude2ndDerSmooth(50:(end-50)))*.2, 1))]);
    
    ylim([-((max(Magnitude2ndDerSmooth(50:(end-50)))))*1.1...
        ((max(Magnitude2ndDerSmooth(50:(end-50)))))*1.1]);
    
    %Need to fill the derivatives to have the same length as k0
    Magnitude2ndDerSmooth(NPro - 1) = Magnitude2ndDerSmooth(NPro - 2);
    Magnitude2ndDerSmooth(NPro)     = Magnitude2ndDerSmooth(NPro - 2);

    Magnitude1stDerSmooth(NPro) = Magnitude1stDerSmooth(NPro - 1);
    
    StartGateSmoothExp  = smk0;
    StartGateSmoothInsp = smk0;
    
    % Selecting projections at end-expiration
% This makes the "gating window" bars essentially
    StartGateSmoothExp(Magnitude1stDerSmooth > eLL_val &...
    Magnitude1stDerSmooth < eUL_val &...
    Magnitude2ndDerSmooth < exp_thres) = max(smk0(:));
    % This variable contains locations of k0 points at expiration
    SelectVectorExp = Magnitude1stDerSmooth > eLL_val &...
        Magnitude1stDerSmooth < eUL_val &...
        Magnitude2ndDerSmooth < exp_thres;

    % This makes the "gating window" bars essentially
% Selecting projections at end-inspiration
    StartGateSmoothInsp(Magnitude1stDerSmooth > iLL_val &...
    Magnitude1stDerSmooth < iUL_val &...
    Magnitude2ndDerSmooth > insp_thres) =  min(smk0(:));

% This variable contains locations of k0 points at inspiration
    SelectVectorInsp = Magnitude1stDerSmooth > iLL_val &...
        Magnitude1stDerSmooth < iUL_val &...
        Magnitude2ndDerSmooth > insp_thres;
    disp('Finished Gating at End-Inspiration');
    
    % This is strictly for plotting purposes
    StartGateSmoothExp(StartGateSmoothExp   ~= max(StartGateSmoothExp)) = nan;
    StartGateSmoothInsp(StartGateSmoothInsp ~= min(StartGateSmoothInsp)) = nan;

    SelectedProjExp  = SelectVectorExp.*smk0; SelectedProjExp(SelectedProjExp == 0) = nan;
    SelectedProjInsp = SelectVectorInsp.*smk0; SelectedProjInsp(SelectedProjInsp == 0) = nan;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,5);

    plot(StartGateSmoothExp, '-k', 'linewidth', 3); hold on;

    plot(smk0, '-', 'color', [105/256 105/256 105/256], 'linewidth', 3)

    plot(SelectedProjExp, '-r', 'linewidth', 3);

    plot(StartGateSmoothInsp, '-k', 'linewidth', 3);

    plot(SelectedProjInsp, '-b', 'linewidth', 3); hold off;

    ylabel('Signal Intensity', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');
    xlabel('Radial Projection', 'FontSize', 12, 'FontWeight', 'bold', 'Color','k');

    xlim([dispstart dispstart+1000]);
    ylim([(min(k0) - min(k0)*0.001) (max(k0) + max(k0)*0.004)]);
    
    subplot(2,3,6);
    FFT_k0 = abs(fftshift(fft(smk0)));
    
    SampF = 1/(Method_Params.TR/1000);
    %Get x frequency space
    Freq_Axis = SampF*((-(NPro/2)+1):(NPro/2))/NPro;
    plot(Freq_Axis,FFT_k0,'k')
    xlim([0.3 10]);
    
    peak = FFT_k0((NPro/2+20):end);
    freq = Freq_Axis((NPro/2+20):end);
    
    [pmax,ind] = max(peak);
    
    br = freq(ind);
    br_act = 60*br;
    text(1.01*br,1.01*pmax,['Breath Rate: ' num2str(br_act) ' bpm']);
    
    answer = questdlg('Is this Retrospective Gating okay?');       
    if ~strcmp(answer,'Yes') 
        prompt = {['Insp Lower Threshold (Was ' num2str(iLL_val) ')'],['Insp Upper Threshold (Was ' num2str(iUL_val) ')'],['Exp Lower Threshold (Was ' num2str(eLL_val) ')'],['Exp Upper Threshold (Was ' num2str(eUL_val) ')'],['Expiration Threshold (Was ' num2str(exp_thres) ')'],['Inspiration Threshold (Was ' num2str(insp_thres) ')'],'Gate off of (p for phase, m for magnitude'};
        dlgtitle = 'New Retrospective Gating Parameters';
        dims = [1 50];
        definput = {num2str(iLL_val),num2str(iUL_val),num2str(eLL_val),num2str(eUL_val),num2str(exp_thres),num2str(insp_thres),char(PorM)};
        user_input = inputdlg(prompt,dlgtitle,dims,definput);
        iLL_val = str2num(user_input{1});
        iUL_val = str2num(user_input{2});
        eLL_val = str2num(user_input{3});
        eUL_val = str2num(user_input{4});
        exp_thres = str2num(user_input{5});
        insp_thres = str2num(user_input{6});
        PorM = user_input{7};
    else
        happy_gate = 1;
        saveas(h,'Respiratory RetroGating Summary.png')
    end
    close
end

Insp_indx = SelectVectorInsp;
Exp_indx = SelectVectorExp;
%% Cardiac Gating Section
answer = questdlg('Would you like to cardiac gate as well?');   
if ~strcmp(answer,'Yes') 
    save(['V2_RetroGating' num2str(number) '.mat'],'Insp_indx','Exp_indx')
else
    %% Start by displaying k0 magnitude and the smoothed version
    happy_gate = 0;
    D1UL = 0.25;
    D1LL = -0.25;
    D2UL = -0.05;
    D2LL = 0.05;
    sm_wind = 15;
    while happy_gate == 0
        if strcmp(PorM,'m')
            k0 = squeeze(abs(fid(1,:)));
        else
            k0 = squeeze(angle(fid(1,:)));
        end
        smk0 = smooth(k0,30);
        
        h = figure('Name','Cardiac Retrospective Gating');
        set(h,'Units','Normalized','Position',[.05 .05 .9 .8],'Color','w')
        subplot(2,5,1);
     
        SampF = 1/(Method_Params.TR/1000);
        %Get x frequency space
        Freq_Axis = SampF*((-(NPro/2)+1):(NPro/2))/NPro;
        FFT_k0 = abs(fftshift(fft(k0)));
        Time_Ax = (0:(NPro-1))*Method_Params.TR/1000;
        
        dispstart2 = Time_Ax(dispstart);
        dispend = Time_Ax(dispstart+500);
        plot(Time_Ax,k0, '.k', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', 'k', 'linewidth', 2);
        xlim([dispstart2 dispend]);
        if strcmp(PorM,'m')
            title('Gating off Magnitude')
        else
            title('Gating off Phase')
        end
        
        %FFT of Raw k0
        subplot(2,5,2);
        plot(Freq_Axis,FFT_k0);
        xlim([0.3 10]);
        %Find the dominant Frequency
        peak = FFT_k0((NPro/2+20):(end-NPro/4));
        freq = Freq_Axis((NPro/2+20):(end-NPro/4));
    
        [pmax,ind] = max(peak);
    
        dom_freq = freq(ind);
        br_act = 60*dom_freq;
        text(1.05*br,0.85*pmax,['Resp. Rate: ' num2str(br_act) ' bpm']);
        
        %Need to try to filter out harmonics - Let's try to murder out to
        %the 10th harmonic
        %We're going to have a crapload of points, so our respiration peak
        %should be super narrow. And we don't want to nuke the cardiac...
        %Let's try 0.1 Hz
        pass_width = 0.1;
        filt_k0 = k0;
        for i = 1:10
            filt_k0 = bandstop(filt_k0,[dom_freq*i-pass_width/2 dom_freq*i+pass_width/2],SampF);
%             subplot(2,3,4)
%             plot(Time_Ax,filt_k0,'.k')
%             xlim([dispstart2 dispend]);
%             title('Time Domain - Filtered k0')
            
          %  FFT_Filt = abs(fftshift(fft(filt_k0)));
%             subplot(2,3,5)
%             plot(Freq_Axis,FFT_Filt)
%             xlim([0.3 10]);
%             title('Freq. Domain - Filtered k0')
        end
        
        %Now that I've filtered out the first 10 harmonics, let's nuke
        %any higher frequencies (higher than 20 Hz) with a lowpass filter
        filt_k0 = lowpass(filt_k0,20,SampF);
        FFT_Filt = abs(fftshift(fft(filt_k0)));
        
        %Now, display heart rate to screen
        hr_peak = FFT_Filt((NPro/2+20):(end-NPro/4));
        hr_freq = Freq_Axis((NPro/2+20):(end-NPro/4));
        
        %heart rate obviously won't be less than 1 or greater than 10
        getrid = find(hr_freq < 1);
        hr_peak(getrid) = [];
        hr_freq(getrid) = [];
        
        getrid2 = find(hr_freq > 10);
        hr_peak(getrid2) = [];
        hr_freq(getrid2) = [];
        
        [pmax,ind] = max(hr_peak);
        dom_hr_freq = hr_freq(ind);
        
        %Now that I have the heart rate, do a final band-pass filter to try
        %to fully isolate the heart beat frequency - Use a nice broad
        %filter (2 Hz)
        
        HR_passwidth = 2;
        filt_k0 = bandpass(filt_k0,[1 40],SampF);
        FFT_Filt = abs(fftshift(fft(filt_k0)));
        
        subplot(2,5,3)
        plot(Time_Ax,filt_k0,'.k')
        xlim([dispstart2 dispend]);
        title('Time Domain - Filtered k0')

        subplot(2,5,4)
        plot(Freq_Axis,FFT_Filt)
        xlim([0.3 10]);
        title('Freq. Domain - Filtered k0')
        
        
        text(1.05*dom_hr_freq,0.85*pmax,['Heart Rate: ' num2str(dom_hr_freq*60) ' bpm']);
        %All this filtering totally screws with the first and last
        %points... Kill the first and last 3000 points (I'm assuming we'll
        %have plenty of projections for doing such things
        CutPts = 3000;
        filt_k0(1:CutPts) = [];
        filt_k0((end-CutPts):end) = [];
        Time_Ax2 = Time_Ax((CutPts+1):(end-CutPts-1));

        NPro2 = length(Time_Ax2);
        %Now, I want to do a double smoothing using a window that is 0.5
        %times as wide as 1 whole cardiac cycle
        card_cy_t = 1/dom_hr_freq;
        card_cy_pts = card_cy_t/(Method_Params.TR/1000);
        sm_wind = round(0.5*card_cy_pts);
        
        %Filter the first time
        sm_filt_k0 = smooth(filt_k0,sm_wind);
        sm_filt_k0 = smooth(sm_filt_k0,sm_wind);
        %Filter the second time with a narrower smooth window to iron out
        %the last few troubles
        sm_filt_k0 = smooth(sm_filt_k0,round(sm_wind/2));
        %And one more time for good measure
        sm_filt_k0 = smooth(sm_filt_k0,round(sm_wind/4));
        
        subplot(2,5,6)
        plot(Time_Ax2,sm_filt_k0,'.-k')
        xlim([dispstart2 dispend]);
        title('Time Domain - Smooth Filtered k0')
        
        %Measure first derivative
        diff_k0 = diff(sm_filt_k0,1);
        diff_k0(NPro2) = diff_k0(NPro2-1);
        
        sm_diff_k0 = smooth(diff_k0,sm_wind);
        
        subplot(2,5,7)
        plot(Time_Ax2,sm_diff_k0,'.-k')
        xlim([dispstart2 dispend]);
        title('1st Derivative')
        hold on
        plot([0 max(Time_Ax2)],[D1UL D1UL],'r--','LineWidth',2);
        plot([0 max(Time_Ax2)],[D1LL D1LL],'r--','LineWidth',2);
         
        %Measure second derivative
        diff2_k0 = diff(sm_diff_k0,1);
        diff2_k0(NPro2) = diff2_k0(NPro2-1);
        sm_diff2_k0 = smooth(diff2_k0,sm_wind);
        
        subplot(2,5,8)
        plot(Time_Ax2,sm_diff2_k0,'.-k')
        xlim([dispstart2 dispend]);
        title('2nd Derivative')
        hold on
        plot([0 max(Time_Ax2)],[D2UL D2UL],'r--','LineWidth',2);
        plot([0 max(Time_Ax2)],[D2LL D2LL],'r--','LineWidth',2);
        
         %Pick stuff based on first derivative
        inflec1 = (sm_diff_k0>D1LL & sm_diff_k0<D1UL & sm_diff2_k0<D2UL);
        inflec2 = (sm_diff_k0>D1LL & sm_diff_k0<D1UL & sm_diff2_k0>D2LL);
        
        subplot(2,5,9)
        plot(Time_Ax2,sm_filt_k0,'.-k')
        xlim([dispstart2 dispend]);
        hold on
        plot(Time_Ax2(inflec1),sm_filt_k0(inflec1),'.r')
        plot(Time_Ax2(inflec2),sm_filt_k0(inflec2),'.b')
        title('Selected Points')
        
        subplot(2,5,10)
        %Finally, combine cardiac gating with respiratory gating
        card_high = inflec1 & Exp_indx;
        card_low = inflec2 & Exp_indx;
       
        plot(Time_Ax,smk0,'.-k')
        %Make my dimensions twice as big as before
        xlim([dispstart2 (dispend+(dispend-dispstart2))]);
        hold on
        plot(Time_Ax(card_high),smk0(card_high),'.r','MarkerSize',10)
        plot(Time_Ax(card_low),smk0(card_low),'.b','MarkerSize',10)
        title('Selected Points_final')
        
        answer = questdlg('Are you happy with this cardiac gating?');       
        if ~strcmp(answer,'Yes') 
            prompt = {'Smoothing Window','1st Derivative Lower Limit','1st Derivative Upper Limit','Second Derivative Lower Limit (Should be positive)','Second Derivative Upper Limit (Should be Negative)','Gate off of (p for phase, m for magnitude'};
            dlgtitle = 'New Cardiac Retrospective Gating Parameters';
            dims = [1 50];
            definput = {num2str(sm_wind),num2str(D1LL),num2str(D1UL),num2str(D2LL),num2str(D2UL),char(PorM)};
            user_input = inputdlg(prompt,dlgtitle,dims,definput);
            sm_wind = str2num(user_input{1});
            D1LL = str2num(user_input{2});
            D1UL = str2num(user_input{3});
            D2LL = str2num(user_input{4});
            D2UL = str2num(user_input{5});
            PorM = user_input{6};
            close;
        else
            happy_gate = 1;
            saveas(h,'Cardiac RetroGating Summary.png')
        end
        
    end
    save(['V2_RetroGating' num2str(number) '.mat'],'Insp_indx','Exp_indx','card_high','card_low')
end


%save(['V2_RetroGating' num2str(number) '.mat'],'Insp_indx','Exp_indx')




