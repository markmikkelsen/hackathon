function out = QA_InVivo(MRS_struct, showPlots, simInd)
% showPlots = false;
% Set some initial variables
n = 1:numel(MRS_struct.gabafile); % spectra to run QA on
n(simInd) = []; % exclude simulated data
freq = MRS_struct.spec.freq;
GABA_diff_noalign = real(MRS_struct.spec.vox1.GABAGlx.diff_noalign);
GABA_diff = real(MRS_struct.spec.vox1.GABAGlx.diff);
GSH_diff_noalign = real(MRS_struct.spec.vox1.GSH.diff_noalign);
GSH_diff = real(MRS_struct.spec.vox1.GSH.diff);

sum_prealign = real(MRS_struct.spec.subspec.prealign.sum);
sum_postalign = real(MRS_struct.spec.subspec.postalign.sum);
resid_prealign = real(MRS_struct.spec.subspec.prealign.resid);
resid_postalign = real(MRS_struct.spec.subspec.postalign.resid);

GABA_diff_noalign(1:10,:) = GABA_diff_noalign(1:10,:)*1.13;
GABA_diff_noalign(11:20,:) = GABA_diff_noalign(11:20,:)*74;
GABA_diff_noalign(21:30,:) = GABA_diff_noalign(21:30,:)*2.8;
GABA_diff_noalign(31:40,:) = GABA_diff_noalign(31:40,:)*0.9;
GABA_diff(1:10,:) = GABA_diff(1:10,:)*1.13;
GABA_diff(11:20,:) = GABA_diff(11:20,:)*74;
GABA_diff(21:30,:) = GABA_diff(21:30,:)*2.8;
GABA_diff(31:40,:) = GABA_diff(31:40,:)*0.9;

GSH_diff_noalign(1:10,:) = GSH_diff_noalign(1:10,:)*1.13;
GSH_diff_noalign(11:20,:) = GSH_diff_noalign(11:20,:)*74;
GSH_diff_noalign(21:30,:) = GSH_diff_noalign(21:30,:)*2.8;
GSH_diff_noalign(31:40,:) = GSH_diff_noalign(31:40,:)*0.9;
GSH_diff(1:10,:) = GSH_diff(1:10,:)*1.13;
GSH_diff(11:20,:) = GSH_diff(11:20,:)*74;
GSH_diff(21:30,:) = GSH_diff(21:30,:)*2.8;
GSH_diff(31:40,:) = GSH_diff(31:40,:)*0.9;


grey = [0.6 0.6 0.6];
shading = 0.3;
w = 0.3;
h = 0.65;
l = (1-w)/2;
b = 0.3;

for ii = n
    
    % Set frequency ranges
    lb = find(freq(ii,:) >= 3.185-0.05+0.02);
    ub = find(freq(ii,:) <= 3.185+0.05+0.02);
    ChoRange = intersect(lb,ub);
    
    lb = find(freq(ii,:) >= 10);
    ub = find(freq(ii,:) <= 11);
    noiseRange = intersect(lb,ub);
    
    lb = find(freq(ii,:) >= 3.2-0.75);
    ub = find(freq(ii,:) <= 3.2+0.75);
    plotRange = intersect(lb,ub);
    
    lb = find(freq(ii,:) >= 0.5);
    ub = find(freq(ii,:) <= 4.25);
    metabRange = intersect(lb,ub);
    
%     df = abs(freq(ii,1)-freq(ii,2));
    
    % Calculate QA metrics
    out.SA.GABA.prealign_std(ii) = std(GABA_diff_noalign(ii,ChoRange));
    out.SA.GABA.prealign_max(ii) = max(GABA_diff_noalign(ii,ChoRange));
    out.SA.GABA.postalign_std(ii) = std(GABA_diff(ii,ChoRange));
    out.SA.GABA.postalign_max(ii) = max(GABA_diff(ii,ChoRange));
    out.noise.GABA.prealign(ii) = std(detrend(GABA_diff_noalign(ii,noiseRange)));
    out.noise.GABA.postalign(ii) = std(detrend(GABA_diff(ii,noiseRange)));
    out.SA.GABA.prealign_SNR(ii) = out.SA.GABA.prealign_max(ii)/out.noise.GABA.prealign(ii);
    out.SA.GABA.postalign_SNR(ii) = out.SA.GABA.postalign_max(ii)/out.noise.GABA.postalign(ii);
        
    out.SA.GSH.prealign_std(ii) = std(GSH_diff_noalign(ii,ChoRange));
    out.SA.GSH.prealign_max(ii) = max(GSH_diff_noalign(ii,ChoRange));
    out.SA.GSH.postalign_std(ii) = std(GSH_diff(ii,ChoRange));
    out.SA.GSH.postalign_max(ii) = max(GSH_diff(ii,ChoRange));
    out.noise.GSH.prealign(ii) = std(detrend(GSH_diff_noalign(ii,noiseRange)));
    out.noise.GSH.postalign(ii) = std(detrend(GSH_diff(ii,noiseRange)));
    out.SA.GSH.prealign_SNR(ii) = out.SA.GSH.prealign_max(ii)/out.noise.GSH.prealign(ii);
    out.SA.GSH.postalign_SNR(ii) = out.SA.GSH.postalign_max(ii)/out.noise.GSH.postalign(ii);
    
    % Calculate power (with rms)
    out.power.sum_prealign(ii) = rms(sum_prealign(ii,metabRange)) / rms(sum_prealign(ii,noiseRange));
    out.power.sum_postalign(ii) = rms(sum_postalign(ii,metabRange)) / rms(sum_postalign(ii,noiseRange));
    out.power.resid_prealign(ii) = rms(resid_prealign(ii,metabRange)) / rms(resid_prealign(ii,noiseRange));
    out.power.resid_postalign(ii) = rms(resid_postalign(ii,metabRange)) / rms(resid_postalign(ii,noiseRange));
    
    out.power.sum_prctChange(ii) = 100 * ((out.power.sum_postalign(ii) - out.power.sum_prealign(ii)) / ...
                                   out.power.sum_prealign(ii));
    out.power.resid_prctChange(ii) = 100 * ((out.power.resid_postalign(ii) - out.power.resid_prealign(ii)) / ...
                                     out.power.resid_prealign(ii));
    
    % Calculate power (by integration)
    out.power2.sum_prealign(ii) = (1/length(sum_prealign(ii,metabRange))) * sum(abs(sum_prealign(ii,metabRange)).^2);
    out.power2.sum_postalign(ii) = (1/length(sum_postalign(ii,metabRange))) * sum(abs(sum_postalign(ii,metabRange)).^2);
    out.power2.resid_prealign(ii) = (1/length(resid_prealign(ii,metabRange))) * sum(abs(resid_prealign(ii,metabRange) - mean(resid_prealign(ii,metabRange))).^2);
    out.power2.resid_postalign(ii) = (1/length(resid_postalign(ii,metabRange))) * sum(abs(resid_postalign(ii,metabRange) - mean(resid_postalign(ii,metabRange))).^2);
    
    out.power2.sum_prctChange(ii) = 100 * ((out.power2.sum_postalign(ii) - out.power2.sum_prealign(ii)) / ...
                                   out.power2.sum_prealign(ii));
    out.power2.resid_prctChange(ii) = 100 * ((out.power2.resid_postalign(ii) - out.power2.resid_prealign(ii)) / ...
                                     out.power2.resid_prealign(ii));
    
end

% MEDIAN_PRE_GABA = [repmat(median(out.SA.GABA.prealign_std(1:10)),[1 10]) ...
%                    repmat(median(out.SA.GABA.prealign_std(11:20)),[1 10]) ...
%                    repmat(median(out.SA.GABA.prealign_std(21:30)),[1 10]) ...
%                    repmat(median(out.SA.GABA.prealign_std(31:40)),[1 10])];
% MEDIAN_PRE_GSH = [repmat(median(out.SA.GSH.prealign_std(1:10)),[1 10]) ...
%                   repmat(median(out.SA.GSH.prealign_std(11:20)),[1 10]) ...
%                   repmat(median(out.SA.GSH.prealign_std(21:30)),[1 10]) ...
%                   repmat(median(out.SA.GSH.prealign_std(31:40)),[1 10])];
% MEDIAN_PRE_GABA = ones(1,length(n));
% MEDIAN_PRE_GSH = ones(1,length(n));

%     out.SA.GABA.prealign_std(1:10)=out.SA.GABA.prealign_std(1:10)*1.13;
%    out.SA.GABA.prealign_std(11:20)=out.SA.GABA.prealign_std(11:20)*74;
%    out.SA.GABA.prealign_std(21:30)=out.SA.GABA.prealign_std(21:30)*2.8;
%    out.SA.GABA.prealign_std(31:40)=out.SA.GABA.prealign_std(31:40)*0.9;
%    
%    out.SA.GSH.prealign_std(1:10)=out.SA.GSH.prealign_std(1:10)*1.13;
%    out.SA.GSH.prealign_std(11:20)=out.SA.GSH.prealign_std(11:20)*74;
%    out.SA.GSH.prealign_std(21:30)=out.SA.GSH.prealign_std(21:30)*2.8;
%    out.SA.GSH.prealign_std(31:40)=out.SA.GSH.prealign_std(31:40)*0.9;
%    
%    
%    out.SA.GABA.postalign_std(1:10)=out.SA.GABA.postalign_std(1:10)*1.13;
%    out.SA.GABA.postalign_std(11:20)=out.SA.GABA.postalign_std(11:20)*74;
%    out.SA.GABA.postalign_std(21:30)=out.SA.GABA.postalign_std(21:30)*2.8;
%    out.SA.GABA.postalign_std(31:40)=out.SA.GABA.postalign_std(31:40)*0.9;
%    
%    out.SA.GSH.postalign_std(1:10)=out.SA.GSH.postalign_std(1:10)*1.13;
%    out.SA.GSH.postalign_std(11:20)=out.SA.GSH.postalign_std(11:20)*74;
%    out.SA.GSH.postalign_std(21:30)=out.SA.GSH.postalign_std(21:30)*2.8;
%    out.SA.GSH.postalign_std(31:40)=out.SA.GSH.postalign_std(31:40)*0.9;
   
   MEDIAN_PRE_GABA = median(out.SA.GABA.prealign_std);
   MEDIAN_PRE_GSH = median(out.SA.GSH.prealign_std);
   
%     MEDIAN_PRE_GABA_post = median(out.SA.GABA.postlign_std);
%    MEDIAN_PRE_GSH_post = median(out.SA.GSH.postlign_std);

for ii = n
    
    % Normalize outcomes relative to pre-aligned data
    % 1 = perfect; 0 = did nothing; < 0 = worse than no alignment
    out.SA.GABA.quality_pre(ii) = 1 - out.SA.GABA.prealign_std(ii)./MEDIAN_PRE_GABA;
    out.SA.GSH.quality_pre(ii) = 1 - out.SA.GSH.prealign_std(ii)./MEDIAN_PRE_GSH;
    out.SA.GABA.quality(ii) = 1 - out.SA.GABA.postalign_std(ii)./MEDIAN_PRE_GABA;
    out.SA.GSH.quality(ii) = 1 - out.SA.GSH.postalign_std(ii)./MEDIAN_PRE_GSH;
    
    out.SA.GABA.quality_prag(ii) = (max([out.SA.GABA.quality(ii) ; out.SA.GABA.quality_pre(ii)]));
    out.SA.GSH.quality_prag(ii) = (max([out.SA.GSH.quality(ii) ; out.SA.GSH.quality_pre(ii)]));
    
    % Plots
    if showPlots
        
        h1 = figure(111);
        clf;
        set(h1, 'Units', 'normalized', 'OuterPosition', [l b w h], ...
            'Name', 'QA of In Vivo Data', 'NumberTitle', 'off');
        % GABA DIFF spectra
        subplot(2,1,1);
        hold on;
        patch([freq(ii,ChoRange), fliplr(freq(ii,ChoRange))], ...
            [repmat(1e6, size(freq(ii,ChoRange))), fliplr(repmat(-1e6, size(freq(ii,ChoRange))))], 1, ...
            'facecolor', grey+(1-grey)*(1-shading), 'edgecolor', 'none');
        h2 = plot(freq(ii,:), GABA_diff_noalign(ii,:), 'r', freq(ii,:), GABA_diff(ii,:), 'b');
        hold off;
        set(gca, 'xdir', 'reverse', 'tickdir', 'out', 'xlim', [3.2-0.75 3.2+0.75], ...
            'ylim', [min([GABA_diff_noalign(ii,plotRange), GABA_diff(ii,plotRange)]), max([GABA_diff_noalign(ii,plotRange), GABA_diff(ii,plotRange)])]);
        text(0.5, 0.85, 'Cho SA', 'Rotation', 90, 'HorizontalAlignment', 'right', 'Units', 'normalized');
        title(sprintf('PRESS SPACE TO CONTINUE\n\nGABA'));
        legend(h2, {'PRE','POST'}, 'location', 'northeast');
        legend boxoff;
        
        % GSH DIFF spectra
        subplot(2,1,2);
        hold on;
        patch([freq(ii,ChoRange), fliplr(freq(ii,ChoRange))], ...
            [repmat(1e6, size(freq(ii,ChoRange))), fliplr(repmat(-1e6, size(freq(ii,ChoRange))))], 1, ...
            'facecolor', grey+(1-grey)*(1-shading), 'edgecolor', 'none');
        plot(freq(ii,:), GSH_diff_noalign(ii,:), 'r', freq(ii,:), GSH_diff(ii,:), 'b');
        hold off;
        title('GSH');
        set(gca, 'xdir', 'reverse', 'tickdir', 'out', 'xlim', [3.2-0.75 3.2+0.75], ...
            'ylim', [min([GSH_diff_noalign(ii,plotRange), GSH_diff(ii,plotRange)]), max([GSH_diff_noalign(ii,plotRange), GSH_diff(ii,plotRange)])]);
                
%         % Power spectra
%         h1 = figure(112);
%         clf;
%         set(h1, 'Units', 'normalized', 'OuterPosition', [l b w h], ...
%             'Name', 'Sum & "Residual" Spectra', 'NumberTitle', 'off');
%         subplot(2,1,1);
%         h2 = plot(freq(ii,metabRange), sum_prealign(ii,metabRange), ...
%             freq(ii,metabRange), sum_postalign(ii,metabRange));
%         set(gca, 'xdir', 'reverse', 'tickdir', 'out');
%         legend(h2, {'PRE','POST'}, 'location', 'northeast');
%         legend boxoff;
%         
%         subplot(2,1,2);
%         plot(freq(ii,metabRange), resid_prealign(ii,metabRange) - mean(resid_prealign(ii,metabRange)), ...
%             freq(ii,metabRange), resid_postalign(ii,metabRange) - mean(resid_postalign(ii,metabRange)));
%         set(gca, 'xdir', 'reverse', 'tickdir', 'out');
        
%         % Power spectra (2)
%         h1 = figure(113);
%         clf;
%         set(h1, 'Units', 'normalized', 'OuterPosition', [l b w h], ...
%             'Name', 'Sum & "Residual" Spectra', 'NumberTitle', 'off');
%         subplot(2,1,1);
%         h2 = plot(freq(ii,metabRange), (1/length(sum_prealign(ii,metabRange))) * abs(sum_prealign(ii,metabRange)).^2, ...
%             freq(ii,metabRange), (1/length(sum_postalign(ii,metabRange))) * abs(sum_postalign(ii,metabRange)).^2);
%         set(gca, 'xdir', 'reverse', 'tickdir', 'out');
%         legend(h2, {'PRE','POST'}, 'location', 'northeast');
%         legend boxoff;
%         
%         subplot(2,1,2);
%         plot(freq(ii,metabRange), (1/length(resid_prealign(ii,metabRange))) * abs(resid_prealign(ii,metabRange) - mean(resid_prealign(ii,metabRange))).^2, ...
%             freq(ii,metabRange), (1/length(resid_postalign(ii,metabRange))) * abs(resid_postalign(ii,metabRange) - mean(resid_postalign(ii,metabRange))).^2);
%         set(gca, 'xdir', 'reverse', 'tickdir', 'out');
        
        fprintf(['\n' MRS_struct.gabafile{ii}]);
        pause;
        
    end
    
end

out.SA.GABA.overall_quality_pre = median(out.SA.GABA.quality_pre);
out.SA.GSH.overall_quality_pre = median(out.SA.GSH.quality_pre);
out.SA.GABA.overall_quality = median(out.SA.GABA.quality);
out.SA.GSH.overall_quality = median(out.SA.GSH.quality);


% out.SA.overall_quality = (out.SA.GABA.overall_quality + out.SA.GSH.overall_quality)/2;
% out.SA.overall_quality_prag = (out.SA.GABA.overall_quality_prag + out.SA.GSH.overall_quality_prag)/2;

out.power.sum_prealign_overall = mean(out.power.sum_prealign);
out.power.sum_postalign_overall = mean(out.power.sum_postalign);
out.power.resid_prealign_overall = mean(out.power.resid_prealign);
out.power.resid_postalign_overall = mean(out.power.resid_postalign);

out.power.sum_prctChange_overal = mean(out.power.sum_prctChange);
out.power.resid_prctChange_overal = mean(out.power.resid_prctChange);

% Power (2)
out.power2.sum_prealign_overall = mean(out.power2.sum_prealign);
out.power2.sum_postalign_overall = mean(out.power2.sum_postalign);
out.power2.resid_prealign_overall = mean(out.power2.resid_prealign);
out.power2.resid_postalign_overall = mean(out.power2.resid_postalign);

out.power2.sum_prctChange_overal = mean(out.power2.sum_prctChange);
out.power2.resid_prctChange_overal = mean(out.power2.resid_prctChange);

close all;

