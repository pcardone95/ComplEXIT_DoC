function parse_data_4ECI_temp(basename, filepath)
%% Parse data for Minji
%% Set up things
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};

% 60 Channels
load('C:\Users\Paolo\OneDrive - Universite de Liege\Documents\GitHub\ComplEXIT_DoC\Ketamine\eximia_addchn.mat',...
    'eximia_addchn') % Structure with good additional channels
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


for dosei = 1:length(dose_l) 
    EEG = pop_loadset('filename',sprintf('%s_%s.set', basename,dose_l{dosei}),'filepath',filepath);
    
    % Control where are the good electrodes
    [~, channels_ECI]= ismember(channel_name, {EEG.chanlocs.labels});
    channels_ECI = nonzeros(channels_ECI);
    
    % 60 Eximia channels
    EEG.chanlocs = [EEG.chanlocs, eximia_addchn];
    
    % Create 0 data
    EEG.data = [EEG.data(:,:,:); zeros(length(eximia_addchn),2500,size(EEG.data(:,:,:),3))];
    EEG.nbchan = EEG.nbchan + length(eximia_addchn);
    %
    % Interpolate
    EEG = eeg_interp(EEG, 128:EEG.nbchan);
    
    EEG = pop_select(EEG, 'channel', [channels_ECI' 128:132]);
    pop_saveset(EEG,'filepath', filepath,...
        'filename',sprintf('%s_%s_60ch.set', basename,dose_l{dosei}));
end

