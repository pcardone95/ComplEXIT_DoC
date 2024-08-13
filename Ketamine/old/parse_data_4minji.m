%% Parse data for Minji
%% Set up things
data_wd = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';% pwd();
minji_wd = fullfile(data_wd,'minji');


mohawkpath = ...
'C:\Users\Paolo\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\MOHAWK\Volumes\bigdisk\Work\MOHAWK';
loadpathsloc = fullfile(mohawkpath,'loadpaths.m');
cd(data_wd)

drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02", "KET03", "KET03"],...
    "drug", ["Placebo", "Ketamine", "Ketamine", "Placebo", "Placebo", "Ketamine"]);
basename_names = {'Pilot_KET01_exp', 'PilotKET02_sess2_exp',...
    'PK02_S1_exp','PK02_S2_exp', 'PK3_S1_exp', 'PK3_S2_exp'};
doses = {'0.15', '0.30', '0.45', '0.60', '0.75'};

% 60 Channels
load('eximia_addchn.mat') % Structure with good additional channels
channel_name = {'Fp1', 'Fpz', 'Fp2', 'AF1', 'AFz', 'AF2', ...
'F7' ,'F3', 'F1', 'Fz','F2', 'F4', 'F8', 'FT9', ...
'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', ...
'T7','C5', 'C3', 'C1',...
'Cz',...% <-- check this
'C2','C4','C6','T8', 'TP9', ...% 
'TP7', 'CP5',  'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...%
'P9', 'P7', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P8', 'P10', 'PO3', 'POz',...
'PO4','O1','Oz','O2', 'Iz'};

% Missing ones
% {'Fpz','AF1','AF2','Cz','Iz'}

%% Parse data
% Loop throught the data
fold_pat = ls("Pilot*");
for pati = 1%1:size(fold_pat,1) % Microsoft
    if ~exist(fullfile(minji_wd, fold_pat(pati,:)),'dir')
        mkdir(fullfile(minji_wd, fold_pat(pati,:)));
    end
    cd(fullfile(data_wd, fold_pat(pati,:)))
    sess = ls('Sess*');
    %cd(mohawkpath)
    
    for sessi = 1:2
        %% Import data
        basename = basename_names{(pati-1)*2 + sessi};
        filepath = fullfile(data_wd, fold_pat(pati,:),sess(sessi,:));
        
        fid = fopen(loadpathsloc,'w');
        fprintf(fid,'filepath=''%s'';',filepath);
        fclose(fid);
        %cd(filepath)
        if ~exist(fullfile(minji_wd, fold_pat(pati,:),sess(sessi,:)),'dir')
            mkdir(fullfile(minji_wd, fold_pat(pati,:),sess(sessi,:)));
        end
        EEG = pop_loadset('filepath',filepath, 'filename', [basename, '_clean.set']);
        EEG.filepath = fullfile(minji_wd, fold_pat(pati,:),sess(sessi,:));
        rereference(basename, 1, '','_minji2parse');

        %% Select epochs
        basename = [basename '_minji2parse'];
        EEG = pop_loadset([basename, '.set']);
        real_events  = not(or(strcmp({EEG.urevent.type},'EVNT'), strcmp({EEG.urevent.type},'boundary')));
        real_events_i = find(real_events);
        switch (pati-1)*2 + sessi
            case 1
                trial_cond = [0 EEG.urevent(real_events_i(2:end)-1).init_index size(EEG.data,3)+1];
            case 2
                trial_cond = [0 EEG.urevent(real_events_i([2,3,6,7])-1).init_index size(EEG.data,3)+1];
            case 3
                trial_cond = [EEG.urevent(real_events_i([2,3,4, 5, 6 ])-1).init_index size(EEG.data,3)+1];
            case 4
                trial_cond = sort([EEG.urevent(real_events_i([2, 9, 13])-1).init_index ...
                    size(EEG.data,3)+1 EEG.urevent(real_events_i([5, 11])-2).init_index]);
            case 5
                  trial_cond = sort([EEG.urevent(real_events_i([2, 12])-1).init_index ...
                      EEG.urevent(real_events_i(5)-2).init_index... %Previous is empty
                    size(EEG.data,3)+1 EEG.urevent(real_events_i([4, 9])-1).init_index]);
            case 6
                trial_cond = sort([EEG.urevent(real_events_i([7, 14])-1).init_index ...
                      EEG.urevent(real_events_i(2)-2).init_index... %Previous is empty
                    size(EEG.data,3)+1 EEG.urevent(real_events_i(11)-2).init_index...
                    EEG.urevent(real_events_i(4)-1).init_index]);
        end
        % trial_cond = [trial_cond 0 0]; % Two zeros at the end for next loop
        %% Separation for dose
        % Controlled rejected trials
            fid = fopen(loadpathsloc,'w');
            fprintf(fid,'filepath=''%s'';',[fullfile(minji_wd, fold_pat(pati,:)),sess(sessi,:), '\']);
            fclose(fid);
            cd(fullfile(fullfile(minji_wd, fold_pat(pati,:)),sess(sessi,:)))
            [~, channels_ECI]= ismember(channel_name, {EEG.chanlocs.labels});
            channels_ECI = nonzeros(channels_ECI);
            for dosei = 1%:5%[1, 3, 5] %0.15-0.30; 0.45-0.6; 0.75
                
                EEG_temp = pop_select(EEG,'trial', ...
                    [trial_cond(dosei)+1:trial_cond(dosei+1)-1]);
                pop_saveset(EEG_temp,'filepath',fullfile(minji_wd, fold_pat(pati,:),sess(sessi,:)),...
                    'filename',[sprintf('%s_%s',basename,doses{dosei}) '.set']);
                
                % 60 Eximia channels
                EEG_temp.chanlocs = [EEG_temp.chanlocs, eximia_addchn];
                
                % Create 0 data
                EEG_temp.data = [EEG_temp.data(:,:,:); zeros(length(eximia_addchn),2500,size(EEG_temp.data(:,:,:),3))];
                EEG_temp.nbchan = EEG_temp.nbchan + length(eximia_addchn);
                %
                % Interpolate
                EEG_temp = eeg_interp(EEG_temp, 128:EEG_temp.nbchan);
                
                EEG_temp = pop_select(EEG_temp, 'channel', [channels_ECI 128:132]);
                pop_saveset(EEG_temp,'filepath',fullfile(minji_wd, fold_pat(pati,:),sess(sessi,:)),...
                    'filename',[sprintf('%s_%s',basename,doses{dosei}) '_60ch.set']);         
            end
                
    end
end
cd(fullfile(data_wd, 'results'))



