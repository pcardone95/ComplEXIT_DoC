%% Ketamine pipeline
% Most of the scripts will be taken from the Mohawk, with some small
% modification to take out the bad epochs

%% Import data
% Organizing paths
mohawkpath = ...
'C:\Users\Paolo\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\MOHAWK\Volumes\bigdisk\Work\MOHAWK';

filepath = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';
cd(filepath)
% Select file
[filename,filepath] = uigetfile('*.vhdr','MOHAWK - Select file to process',filepath);
if filename == 0
    return
end
[filename,ext] = strtok(filename,'.');


% Specify a name for the final file
answers = inputdlg2({'Specify dataset base name:'},'MOHAWK Dataset',1,{filename});

if isempty(answers)
    return
elseif isempty(answers{1})
    error('Basename cannot be empty.');
else
    basename = answers{1};
end

% fid = fopen(loadpathsloc,'w');
% fprintf(fid,'filepath=''%s'';',filepath);
% fclose(fid);

if ~exist(sprintf('%s%sfigures',filepath,filesep),'dir')
    mkdir(sprintf('%s%sfigures',filepath,filesep));
end

cur_wd = pwd;
cd(mohawkpath);

%% Data import
% From: dataimport(filename,basename,'VHDR');


fullfile = [filepath filename '.vhdr'];
EEG = pop_fileio(fullfile);
EEG.chanlocs(47).labels = 'FCz';
fprintf('Loading default channel locations.');
EEG = pop_chanedit(EEG, 'lookup', which('standard-10-5-cap385.elp'));
fprintf('Checking for additional channels.');
EEG = pop_select(EEG,'nochannel',{'ECG'});


EEG = eeg_checkset(EEG);

%%%% PREPROCESSING

%REDUCE SAMPLING RATE TO 250HZ
if EEG.srate > 250
    fprintf('Downsampling to 250Hz.\n');
    EEG = pop_resample(EEG,250);
elseif EEG.srate < 250
    error('Sampling rate too low!');
end

%Filter
hpfreq = 0.5;
lpfreq = 45;
fprintf('Low-pass filtering below %.1fHz...\n',lpfreq);
EEG = pop_eegfiltnew(EEG, 0, lpfreq);
fprintf('High-pass filtering above %.1fHz...\n',hpfreq);
EEG = pop_eegfiltnew(EEG, hpfreq, 0);

%Remove line noise. Change line noise frequency below if needed.
fprintf('Removing line noise.\n');
EEG = rmlinenoisemt(EEG, 50);

EEG.setname = sprintf('%s_orig',basename);
EEG.filename = sprintf('%s_orig.set',basename);
EEG.filepath = filepath;

EEG = eeg_checkset(EEG);

fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);

%% Reject data during behavioural assessment
EEG = pop_loadset('filepath',filepath,'filename',[basename '_orig.set']);

% original_event = eeg_eventtable(EEG); % Save original info about events

% pop_eegplot( EEG, 1, 1, 1)

% %% Automatically do it:
% % Need for consistency across recordings: for the moment is like this:
% % 1) Begin of the recording
% % 2) Dose 0.15
% % 3) Dose 0.3
% % 4) Dose 0.45
% % 5) SECONDs 30
% % 6) End SECONDs 30
% % 7) Dose 0.6
% % 8) Dose 0.75
% % 9) SECONDs 60
% % 10) End SECONDs 60
% % 9) SECONDs 90
% % 10) End SECONDs 90
% EEG = eeg_eegrej( EEG, [EEG.urevent(5).latency EEG.urevent(6).latency; ...
%     EEG.urevent(9).latency EEG.urevent(10).latency] );

% 
% %%%K1 SESS1
% EEG = eeg_eegrej( EEG, [(EEG.urevent(3).latency + EEG.srate*60*10) EEG.urevent(4).latency; ...
%    EEG.urevent(6).latency EEG.pnts] );


% %%%K1 SESS2
% EEG = eeg_eegrej( EEG, [EEG.urevent(1).latency EEG.urevent(2).latency; ...
%     EEG.urevent(5).latency EEG.urevent(6).latency; ...
%     (EEG.urevent(8).latency + EEG.srate*60*10 +1) EEG.urevent(9).latency;...
%     EEG.urevent(11).latency EEG.urevent(12).latency] );

% % %%%K2 SESS1
% EEG = eeg_eegrej( EEG, [(EEG.urevent(4).latency + EEG.srate*60*10) EEG.urevent(5).latency + 1; ...
%    EEG.urevent(10).latency-1 EEG.pnts] );

% %%%K2 SESS2
% EEG = eeg_eegrej( EEG, [(EEG.urevent(2).latency + EEG.srate*60*10) EEG.urevent(5).latency-1; ... % Noise
%    (EEG.urevent(9).latency + EEG.srate*60*10 +1) EEG.urevent(11).latency-1;... % SECOND 30
%    (EEG.urevent(13).latency + EEG.srate*60*10 +1) EEG.urevent(15).latency+1;... % SECOND 60
%    (EEG.urevent(16).latency + EEG.srate*60*15 +1) EEG.pnts... % SECOND 90
%    ] );


% %%%K3 SESS1
% EEG = eeg_eegrej( EEG, [(EEG.urevent(2).latency + EEG.srate*60*10) EEG.urevent(5).latency-1; ... % Noise
%    EEG.urevent(7).latency-1  EEG.urevent(8).latency+1;... % SECOND 30
%    EEG.urevent(13).latency-1 EEG.urevent(14).latency+1 ... % SECOND 60
%    ] );

%%%K3 SESS2
EEG = eeg_eegrej( EEG, [EEG.urevent(1).latency+1  EEG.urevent(2).latency-1; ... % Before 0.15
   EEG.urevent(11).latency-1  EEG.urevent(12).latency-1;... % SECOND 30
   EEG.urevent(21).latency-1 EEG.urevent(22).latency+1; ... % SECOND 60
   EEG.urevent(33).latency EEG.pnts... % Spasticity and end
   ] );


%%%%
EEG.setname = sprintf('%s_orig_rej',basename);
EEG.filename = sprintf('%s_orig_rej.set',basename);
% EEG.filepath = filepath;

EEG = eeg_checkset(EEG);
pop_saveset(EEG,'filename', EEG.filename, 'filepath', filepath)
%% Run rest of the analysis
% OLD Maybe to take out
% % Epoch data into 10 sec epochs
% fprintf('\n*** EPOCHING DATA ***\n');
% EEG = pop_loadset('filepath',filepath,'filename',[basename '_orig_rej.set']);
% 
% % Length of each epoch in seconds
% epochlength = 10;
% 
% events = (0:epochlength:EEG.xmax)';
% events = cat(2,repmat({'EVNT'},length(events),1),num2cell(events));
% %assignin('base','events',events);
% 
% EEG = pop_importevent(EEG,'event',events,'fields',{'type','latency'});
% clear events
% EEG = eeg_checkset(EEG,'makeur');
% EEG = eeg_checkset(EEG,'eventconsistency');
% 
% fprintf('\nSegmenting into %d sec epochs.\n',epochlength);
% EEG = pop_epoch(EEG,{'EVNT'},[0 epochlength]);
% 
% EEG = pop_rmbase(EEG,[]);
% 
% EEG = eeg_checkset(EEG);
% 
% EEG.setname = [basename '_epochs'];
% EEG.filename = [basename '_epochs.set'];
% fprintf('Saving %s%s.\n',EEG.filepath,EEG.filename);
% pop_saveset(EEG,'filename', EEG.filename, 'filepath', filepath);

epochdata(basename);
% Parse data into doses - important for later.
EEG = pop_loadset('filename',[basename '_epochs.set'],'filepath',filepath);

EEG.epoch(1).dose =[];
% Assign "dose" to epochs
% index of the event of the dose i.e. 0.15
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};

for dosei = 1:length(dose_l)
    dose_s = dose_l{dosei};
    index_dose = find(strcmp({EEG.event.type}, dose_s));
    
    % index of the event of the event when dose beceome i.e. 0.15
    index_trial = find(cellfun(@(x) sum(x==index_dose), {EEG.epoch.event}, 'UniformOutput', 1));
    
    EEG.epoch(index_trial).dose =[];
    temp = repmat({dose_s}, 1, EEG.trials-index_trial);
    [EEG.epoch(index_trial+1:end).dose] =deal(temp{:});
end

% NOW FOR SECONDS
sec_l = {{'SEC30B', 'SEC30E'}, {'SEC60B', 'SEC60E'}, {'SEC90B', 'SEC90E'}};

for seci = 1:length(sec_l)
    [sec_beg, sec_end] = sec_l{seci}{:};
    
    %Beginning
    index_sec_b = find(strcmp({EEG.event.type}, sec_beg));
    index_trial_b = find(cellfun(@(x) sum(x==index_sec_b), {EEG.epoch.event}, 'UniformOutput', 1));
    if isempty(index_sec_b) % Not in the recording for any reasons - maybe should never happen
        warning('%s not found', sec_beg)
        continue;
    end
    
    %End
    index_sec_e = find(strcmp({EEG.event.type}, sec_end));
    if isempty(index_sec_e) % after the end of the recording
        index_trial_e = EEG.trials;
    else
        index_trial_e = find(cellfun(@(x) sum(x==index_sec_e), {EEG.epoch.event}, 'UniformOutput', 1));
    end
    
    temp = repmat({'SEC'}, 1, index_trial_e-index_trial_b+1);
    [EEG.epoch(index_trial_b:index_trial_e).dose] =deal(temp{:});
end
pop_saveset(EEG,'filename',[basename '_epochs.set'],'filepath',filepath);

%% MANUAL STEP
% First pass of quasi-automatic rejection of noisy channels and epochs based on variance
% thresholding
fprintf('\n*** SELECT BAD CHANNELS AND TRIALS ***\n');
rejartifacts([basename '_epochs'], 1, 4, 1);
%%

fprintf('\n*** COMPUTING IC DECOMPOSITION ***\n');
% Run ICA decomposition with optional PCA pre-processing
computeic([basename '_epochs'])

%% MANUAL
fprintf('\n*** SELECT BAD ICs ***\n');
% Visually identify and reject noisy ICs, e.g., eye movements, muscle
% activity, etc.
rejectic(basename, 'prompt', 'off')

%% REJECT DATA DURING SECONDs
fprintf('\nDeleting trials during SECONDs...\n');
EEG = pop_loadset('filename',[basename '_clean.set'],'filepath',filepath);
% Find bad trials
badtrials = find(strcmp({EEG.epoch.dose}, 'SEC'));

EEG = pop_select(EEG, 'notrial', badtrials);
if isfield(EEG,'rejepoch')
    EEG.rejepoch = [EEG.rejepoch badtrials];
else
    EEG.rejepoch = badtrials;
end
pop_saveset(EEG,'filename',[basename '_clean.set'],'filepath',filepath);

%% MANUAL
% Second and final pass of quasi-automatic rejection of noisy channels and epochs based on variance
% thresholding, to remove any remaining noisy data.
fprintf('\n*** SELECT ANY REMAINING BAD CHANNELS AND TRIALS AND INTERPOLATE ***\n');
rejartifacts([basename '_clean'], 2, 4, 0, [], 500, 250);
%%

fprintf('\n*** REFERENCING DATA TO COMMON AVERAGE ***\n');
% re-reference data to common average for connectivity estimation.
rereference(basename, 1);

%% Full length LZC - TO CHECK
% Only for here, make things easier
EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);

results = struct();

% LZC calculation all over
data_binar = zeros(size(EEG.data(:,:,:)));

% Preallocate
strings = cell(size( EEG.data,1),size( EEG.data,3)); % electrode x trial
C = zeros(size(EEG.data,1),size( EEG.data,3));
H = cell(size(EEG.data,1),size( EEG.data,3));

% Calcule - from Figure 3 (Paper 2)
for electi = 1:size( EEG.data,1)
    for triali = 1:size( EEG.data,3)
        data_binar(electi,:,triali) = abs(hilbert( EEG.data(electi,:,triali)));
        data_binar(electi,:,triali) = ...
            data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
        strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
        [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
    end
end

results(1).C = C;
results(1).H = H;
results(1).C_mean = mean(C,2);
results(1).concentration = 'all';
results(1).ntrials = EEG.trials;

for dosei = 1:length(dose_l) % -1 because last index is the end
    dose_s = dose_l{dosei};
    dose_trials = find(strcmp({EEG.epoch.dose},dose_s));
    
    % Take only the data with the same concentration
    EEG_temp = pop_select(EEG, 'trial', dose_trials);
    pop_saveset(EEG_temp,'filename',sprintf('%s_%s.set', basename, dose_s),...
        'filepath',filepath);
  
    results(dosei+1).C = C(:,dose_trials);
    results(dosei+1).H = H(:,dose_trials);
    results(dosei+1).C_mean = mean(C(:,dose_trials),2);
    results(dosei+1).concentration = dose_s;
    results(dosei+1).ntrials = EEG_temp.trials;
end

% Check for transition - necessary?


% Save file
cd(filepath)
save([basename '_LZC.mat'], 'results')

% Write results in table
results_table = table(sort(repmat([0.15, 0.30, 0.45,0.60,0.75]',127,1)), repmat({EEG.chanlocs.labels}', 4, 1),...
    [results(1).C_mean; results(2).C_mean; results(3).C_mean; results(4).C_mean;...
    results(5).C_mean], 'VariableNames', {'Concentration', 'Electrode', 'LZC'});
writetable(results_table, [basename '_results_LZC'])

%% Divide per concentration
% fprintf('\n*** RETAINING 10 MINUTES (60 EPOCHS) OF DATA ***\n');
% % optionally fix number of epochs contributing to connectivity estimation.
% % 60 epochs below will effectively use 10 minutes of clean data.
% checktrials(basename,60,''); % Prob not need for this
for dosei = 1:length(dose_l)
    dose_s = dose_l{dosei};
    
    fprintf('\n*** CALCULATING MULTITAPER SPECTRUM ***\n');
    % calculate power spectrum using the multi-taper method
    calcftspec(sprintf('%s_%s.set', basename, dose_s));
    
    fprintf('\n*** PLOTTING SPECTRUM ***\n');
    % visualise and save the power spectrum of all channels
    plotftspec(sprintf('%s_%s.set', basename, dose_s));
    
    fprintf('\n*** CALCULATING CONNECTIVITY ***\n');
    % estimate dwPLI connectivity between pairs of channels
    ftcoherence(sprintf('%s_%s.set', basename, dose_s));
    
    fprintf('\n*** CALCULATING GRAPH-THEORETIC NETWORK METRICS ***\n');
    % calculate graph theory metrics
    calcgraph(sprintf('%s_%s.set', basename, dose_s));
    
    % %% plot 3D connectivity topographs in the delta, theta and alpha bands.
    plothead(sprintf('%s_%s.set', basename, dose_s),1);
    plothead(sprintf('%s_%s.set', basename, dose_s),2);
    plothead(sprintf('%s_%s.set', basename, dose_s),3);
    %
    % %%
    % fprintf('\n*** PLOTTING METRICS ***\n');
    % plotmetric(basename,'participation coefficient',3,'ylabel','Network centrality')
    % plotbands(basename,'participation coefficient','title','Network centrality');
    %
end
%% Calculate LZC
pop_loadset([basename, '_clean.set'], filepath)
rereference(basename, 1, '','LZC');

% PARSE DATA
% Check from the even when somethign happens
real_events  = not(or(strcmp({EEG.urevent.type},'EVNT'), strcmp({EEG.urevent.type},'boundary')));
real_events_i = find(real_events);

% Trasform index of events in index of *trial*
% Basically, real_events_i gives use info about the number of the event we
% are looking at. However, we work here as a function of trial!
% trial_cond = [0 EEG.urevent(real_events_i(2:end)-1).init_index size(EEG.data,3)+1];
% trial_cond = [0 EEG.urevent(real_events_i([2,3,6,7])-1).init_index size(EEG.data,3)+1]

% % K2 S1
% trial_cond = [EEG.urevent(real_events_i([2,3,4, 5, 6 ])-1).init_index size(EEG.data,3)+1]
% 
% K2 S2
% % trial_cond = [EEG.urevent(real_events_i([2, 5, 9, 11, 13])-1).init_index size(EEG.data,3)+1]
% trial_cond = sort([EEG.urevent(real_events_i([2, 9, 13])-1).init_index ...
%     size(EEG.data,3)+1 EEG.urevent(real_events_i([5, 11])-2).init_index])

% K3 S1
trial_cond = sort([EEG.urevent(real_events_i([2, 4, 5, 9, 12])-1).init_index ...
    size(EEG.data,3)+1 EEG.urevent(real_events_i([4, 9])-2).init_index]);


% % K3 S2
% trial_cond = sort([0 EEG.urevent(real_events_i([2, 7, 14])-1).init_index ...
%     size(EEG.data,3)+1 EEG.urevent(real_events_i([4, 11])-2).init_index]);



% Putting -1, we just take the index of the trial where the dose change 
% happens! If the dose changes in the middle of trial 22, then  we will
% have the index 22. Note that because of the way we implement the code,
% the first integer will be 0 and the last numb of trial +1 


% temp = struct2table(EEG.chanlocs); % convert the struct array to a table
% temp = cellfun(@(x) x(2:end), temp.labels, 'UniformOutput', false);
% temp = sprintf('%s,', temp{:});
% temp = sscanf(temp, '%g,');
% [~ , index] = sort(temp);
% %EEG.chanlocs = EEG.chanlocs(index);
% %EEG.data = EEG.data(index);

for dosei = 1:length(trial_cond)-1 % -1 because last index is the end
    % Take only the data with the same concentration
    data_temp = EEG.data(:,:,trial_cond(dosei)+1:trial_cond(dosei+1)-1); % channel x time x trial
    data_binar = zeros(size(data_temp));
    
    % Preallocate
    strings = cell(size(data_temp,1),size(data_temp,3)); % electrode x trial
    C = zeros(size(data_temp,1),size(data_temp,3));
    H = cell(size(data_temp,1),size(data_temp,3));
    
    % Calcule - from Figure 3 (Paper 2)
    for electi = 1:size(data_temp,1)
        for triali = 1:size(data_temp,3)
            data_binar(electi,:,triali) = abs(hilbert(data_temp(electi,:,triali)));
            data_binar(electi,:,triali) = ...
                data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
            strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
            [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
        end
    end
    
    results(dosei).C = C;
    results(dosei).H = H;
    results(dosei).C_mean = mean(C,2);
end

cd(filepath)
save([basename '_LZC.mat'], 'results')

% Write results in table
% results_table = table(sort(repmat((0.15:0.15:0.75)',127,1)), repmat({EEG.chanlocs.labels}', 5, 1),...
%     [results(1).C_mean; results(2).C_mean; results(3).C_mean; results(4).C_mean;...
%     results(5).C_mean], 'VariableNames', {'Dose', 'Electrode', 'LZC'});
results_table = table(sort(repmat([0.15,0.45,0.60,0.75]',127,1)), repmat({EEG.chanlocs.labels}', 4, 1),...
    [results(1).C_mean; results(2).C_mean; results(3).C_mean; results(4).C_mean;...
    ], 'VariableNames', {'Dose', 'Electrode', 'LZC'});
writetable(results_table, [basename 'results_LZC'])
%% Quick figure
figure;
hold on
for dosei = 1:4
    plot(dosei, results(dosei).C_mean, '.')
end
xlim([0 6])
ylabel('LZC')
xlabel('Dose')
xticklabels({'', '0.15', '0.30', '0.45', '0.60', '0.75', ''})
title(strrep(basename, '_', ' '))
ax = gca; 
ax.FontSize = 16; 
hold off
savefig([basename '_LZC.fig'])