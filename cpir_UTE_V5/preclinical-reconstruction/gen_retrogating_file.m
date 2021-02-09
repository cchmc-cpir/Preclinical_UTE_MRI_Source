function gen_retrogating_file(fid,number,Method_Params)

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

%display k0 points starting halfway through the scan
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

%Repeat this process until the user approves
while happy_gate == 0
    if strcmp(PorM,'m')
        k0 = squeeze(abs(fid(1,:)));
    else
        k0 = squeeze(angle(fid(1,:)));
    end
    smk0 = smooth(k0,30);
    %% Start by displaying k0 magnitude and the smoothed version
    h = figure('Name','Retrospective Gating');
    set(h,'Units','Normalized','Position',[.05 .05 .9 .8],'Color','w')

    % Display 
    subplot(2,3,1);
    plot(k0, '.', 'MarkerEdgeColor', 'k',...
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
    rectangle('position', [dispstart -100 1000....
        100 + exp_thres],...
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
        saveas(h,'RetroGating Summary.png')
    end
    close
end

Insp_indx = SelectVectorInsp;
Exp_indx = SelectVectorExp;

%WHen this function is called, we're already in the correct folder, so we
%can just save with no path definition

save(['RetroGating' num2str(number) '.mat'],'Insp_indx','Exp_indx')




