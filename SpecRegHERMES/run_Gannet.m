%%%%%%%%%%%%% Batch-run Gannet3.0 (hackathon version) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script *must* be run from within the 'hackathon' directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Add directories to user's search path (ensures 'Gannet3.0' and 'data' are always at the top of the search path)
addpath(fullfile(pwd,'Gannet3.0'), fullfile(pwd,'data'));

% Create 'output' directory
if ~exist('output','dir')
    mkdir('output');
end

% Create list of SDAT files to process
exp = {'GO','IW','KC','MS'};
subj = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};
% exp = {'IW'};
% subj = {'S02'};

metab = cell(numel(exp)*numel(subj),1);
water = metab;
c = 1;
for ii = 1:length(exp)
    for jj = 1:length(subj)
        metab{c} = fullfile(pwd, 'data', 'in_vivo', exp{ii}, subj{jj}, [exp{ii} '_' subj{jj} '_HERMES_act.sdat']);
        water{c} = fullfile(pwd, 'data', 'in_vivo', exp{ii}, subj{jj}, [exp{ii} '_' subj{jj} '_HERMES_ref.sdat']);
        c = c + 1;
    end
end
metab{end+1} = fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat');
water{end+1} = fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat');
% metab{1} = fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat');
% water{1} = fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat');

% Run Gannet
MRS = GannetLoad(metab, water);
close all;

% Save GannetLoad data with unique filename (*_run#.mat)
newFilename = false;
runNum = 1;
filename = ['hackathon_data_run' num2str(runNum) '.mat'];
filename = fullfile(pwd, 'output', filename);
while ~newFilename
    if exist(filename,'file')
        runNum = runNum + 1;
        filename = fullfile(pwd, 'output', ['hackathon_data_run' num2str(runNum) '.mat']);
    else
        newFilename = true;
    end
end
save(filename,'MRS');

