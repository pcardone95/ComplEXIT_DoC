%% Template homogenise events
% Code to homogenise the ensure that all the events are properly saved across different recording. 
% One code for patient
% % Need for consistency across recordings. for the moment is like this:
% 1) Begin of the recording
% 2) Concentration 0.15
% 3) Concentration 0.3
% 4) Concentration 0.45
% 5) SECONDs 30
% 6) End SECONDs 30
% 7) Concentration 0.6
% 8) Concentration 0.75
% 9) SECONDs 60
% 10) End SECONDs 60
% 11) SECONDs 90
% 12) End SECONDs 90

% Set up things as Mohawk
mohawkpath = ...
    'C:\Users\Paolo\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\MOHAWK\Volumes\bigdisk\Work\MOHAWK';
loadpathsloc = fullfile(mohawkpath,'loadpaths.m');

% if exist(loadpathsloc,'file')
%     run(loadpathsloc);
% else
%     filepath = userpath;
% end

%
data_wd = pwd;
sess = ls('Sess*');
filepath = data_wd;
for sessi = 1:2
    cd(fullfile(data_wd,sess(sessi,:)))
    [filename,filepath] = uigetfile('*.vhdr','Select file to homogenize',filepath);
    if filename == 0
        return
    end
    [filename,ext] = strtok(filename,'.');
    
    basename = sprintf('PK[NUMBER_PATIENT]_S%i_exp', sessi);
    
    fid = fopen(loadpathsloc,'w');
    fprintf(fid,'filepath=''%s'';',filepath);
    fclose(fid);
    
    %% Import data
    dataimport(filename,basename, 'VHDR');
    
    EEG = pop_loadset('filepath',filepath,'filename',[basename '_orig.set']);
    
    switch sessi
        case 1
            % Rename
            % For example, when Concentration has been increased to 0.15, it was written "Concentration to 0.15"
            EEG.event(2).type = '0.15'; EEG.urevent(2).type = '0.15';
            
            % Add
	    % If some events are forgotten, then here we put them. I.e., we forgot to put the event for the end of SECONDs at 30,
            % right before we increase Concentration to 0.6
            num = size(EEG.event,2); % Counter for new events
            
            num = num +1;% New event, number to be changed 
            temp = struct('type',{'SEC30E'},'value',{'Comment'},'latency',...
                {(EEG.event(5).latency -1)},'duration',{1});
            EEG.urevent(num) = temp; temp.urevent = num; temp.duration= 0.5; EEG.event(num) = temp;
	    % Note that this is the most conservative method, by saying it finished in the timepoint before concentration was increased
            
            % Remove
            % To make things cleaner. For example EEG.event(5) is a duplicate of EEG.event(4) 
	    EEG.event(5) =[]; EEG.urevent(5) =[]; 
        case 2
	    % Same as before
    end
    
    pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
end

%% Sandbox
