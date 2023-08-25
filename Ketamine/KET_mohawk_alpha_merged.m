%% Ketamine pipeline
% Most of the scripts will be taken from the Mohawk, with some small
% modification to take out the bad epochs

%% Setting up things
% Path
mohawkpath = ...
'C:\Users\Paolo\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\MOHAWK\Volumes\bigdisk\Work\MOHAWK';
filepath = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';

cur_wd = pwd();
cd(filepath) %cd ..
data_wd = filepath;% pwd(); 
% Variables
drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02", "KET03", "KET03"],...
    "drug", ["Placebo", "Ketamine", "Ketamine", "Placebo", "Placebo", "Ketamine"]);
basename_names = {'Pilot_KET01_exp', 'PilotKET02_sess2_exp',...
    'PK02_S1_exp','PK02_S2_exp', 'PK3_S1_exp', 'PK3_S2_exp'};
results_merged_table = table();
elecs = 127; % numb electrodes
loadpathsloc = fullfile(mohawkpath,'loadpaths.m');
dose_ranges = {'0.15-0.30', '0.45-0.60', '0.75'};

%%% TO CHECK %%%
trange = ones(1,37); trange(1:4) = 0; trange = logical(trange);
%%% END CHECK %%%

% load previous results
load('C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine\results\pilotKET_alpha_merge.mat')
results_merged_table.Properties.VariableNames = {'Var1','Var2','Var3','value_alphai'};

% Loop throught the data
fold_pat = ls("Pilot*");
for pati = 3%1:size(fold_pat,1) % Microsoft
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
        cd(filepath)
        if ~exist(fullfile(filepath,'alpha_merged'),'dir')
            mkdir(fullfile(filepath,'alpha_merged'));
            mkdir(fullfile(filepath,'alpha_merged','figures'))
        end
        EEG = pop_loadset('filepath',filepath, 'filename', [basename, '_clean.set']);
        rereference(basename, 1, '','_merged');

        %% Select epochs
        basename = [basename '_merged'];
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
        trial_cond = [trial_cond 0 0]; % Two zeros at the end for next loop
        %% Separation for dose
        % Controlled rejected trials
        
        for dosei = [1, 3, 5] %0.15-0.30; 0.45-0.6; 0.75
            fid = fopen(loadpathsloc,'w');
            fprintf(fid,'filepath=''%s'';',[fullfile(filepath,'alpha_merged'), '\']);
            fclose(fid);
            cd(fullfile(filepath,'alpha_merged'))
            EEG_temp = pop_select(EEG,'trial', ...
                [trial_cond(dosei)+1:trial_cond(dosei+1)-1 ...
                trial_cond(dosei+1)+1:trial_cond(dosei+2)-1]);
            % % New way
%             EEG_temp = pop_select(EEG,'trial', ...
%                 [trial_cond(dosei)+1:trial_cond(dosei+1)-1 ...
%                 trial_cond(dosei+1)+1:trial_cond(dosei+2)-1]);
            pop_saveset(EEG_temp,'filepath',fullfile(filepath,'alpha_merged'),'filename',[sprintf('%s_range%i',basename,dosei) '.set']);

            
            %% Analisis
            
            fprintf('\n*** CALCULATING MULTITAPER SPECTRUM ***\n');
            % calculate power spectrum using the multi-taper method
            calcftspec(sprintf('%s_range%i',basename,dosei));
            
            fprintf('\n*** PLOTTING SPECTRUM ***\n');
            % visualise and save the power spectrum of all channels
            plotftspec(sprintf('%s_range%i',basename,dosei));
            
            fprintf('\n*** CALCULATING CONNECTIVITY ***\n');
            % estimate dwPLI connectivity between pairs of channels
            ftcoherence(sprintf('%s_range%i',basename,dosei));
            
            fprintf('\n*** CALCULATING GRAPH-THEORETIC NETWORK METRICS ***\n');
            % calculate graph theory metrics
            calcgraph(sprintf('%s_range%i',basename,dosei));
            
            % %% plot 3D connectivity topographs in the delta, theta and alpha bands.
            % plothead(basename,1);
            % plothead(basename,2);
            % plothead(basename,3);
            
            %
            fprintf('\n*** PLOTTING METRICS ***\n');
            plotmetric(sprintf('%s_range%i',basename,dosei),'participation coefficient',3,'ylabel','Network centrality')
            plotbands(sprintf('%s_range%i',basename,dosei),'participation coefficient','title','Network centrality');
            
            load(sprintf('%s_range%i_mohawk.mat',basename,dosei))
            value_alphai = squeeze(mean(std(graphdata{8,2}(3,trange,:),[],3),2));
            results_merged_table = [results_merged_table; ...
            table(sprintf("KET0%d", pati), ... %Sub
            drug_sess.drug((pati-1)*2 + sessi),... %Condition
            dose_ranges(ceil(dosei/2)),...  % Dose
            value_alphai)] % alpha
        % pause()
        end
                
    end
end
cd(fullfile(data_wd, 'results'))
results_merged_table.Properties.VariableNames = {'Subject', 'Condition','Dose', 'Alpha_centrality'}%% End&Save
save('pilotKET_alpha_merge.mat_v2','results_merged_table')
writetable(results_merged_table, 'pilotKET_merged_R_alpha_v2.txt')
cd(cur_wd)