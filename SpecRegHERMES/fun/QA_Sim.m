function out = QA_Sim(MRS_struct, showPlots, simInd)
% f_results is vector of Navg freq offset estimates or minus the correction factors (in Hz)
% ph_results is vector of Navg phase offset estimates (in deg)

% Set some initial variables
f_obvs = MRS_struct.out.f_results(simInd,:)';
ph_obvs = MRS_struct.out.ph_results(simInd,:)';

% Load freq and phase offsets added to simulated data (ground truth)
load(fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat'));
f_truth = fs';
f_truth = f_truth(:);
ph_truth = phs';
ph_truth = ph_truth(:);

% Calculate QA metrics
out.error.mean_f = mean(f_obvs - f_truth);
out.error.rms_f = sqrt(mean((f_obvs - f_truth).^2));
out.error.std_f = std(abs(f_obvs - f_truth));

out.error.mean_ph = mean(ph_obvs - ph_truth);
out.error.rms_ph = sqrt(mean((ph_obvs - ph_truth).^2));
out.error.std_ph = std(abs(ph_obvs - ph_truth));

% Normalize outcomes relative to doing nothing
% 1 = perfect; 0 = did nothing; < 0 = worse than doing nothing
rms_f0 = sqrt(mean((zeros(320,1)-f_truth).^2));
rms_ph0 = sqrt(mean((zeros(320,1)-ph_truth).^2));
out.quality.f = 1 - out.error.rms_f/rms_f0;
out.quality.ph = 1 - out.error.rms_ph/rms_ph0;

out.quality.overall = (out.quality.f + out.quality.ph)/2;

% Plots
if showPlots
    
    w = 0.65;
    h = 0.5;
    l = (1-w)/2;
    b = 0.3;
    
    h1 = figure(222);
    set(h1, 'Units', 'normalized', 'OuterPosition', [l b w h], ...
        'Name', 'QA of Simulated Data', 'NumberTitle', 'off');
    
    % Freq offsets
    subplot(1,2,1);
    h2 = plot(1:numel(f_truth), f_truth, 'r', 1:numel(f_obvs), f_obvs, 'b');
    ylabel('freq (Hz)');
    set(gca, 'tickdir', 'out', 'xlim', [1 numel(f_obvs)], 'xtick', [1 50:50:numel(f_obvs)]);
    legend(h2, {'TRUE','ESTIMATED'}, 'location', 'best');
    legend boxoff;
    
    % Phase offsets
    subplot(1,2,2);
    plot(1:numel(ph_truth), ph_truth, 'r', 1:numel(ph_obvs), ph_obvs, 'b');
    ylabel('phase (deg)');
    set(gca, 'tickdir', 'out', 'xlim', [1 numel(ph_obvs)], 'xtick', [1 50:50:numel(ph_obvs)]);
    
%     set(gcf,'NextPlot','add');
%     axes;
%     t = title(sprintf('PRESS SPACE TO CONTINUE\n'));
%     set(gca,'Visible','off');
%     set(t,'Visible','on');
    
%     pause;
    
%     close all;
    
end


