data_wd = '...';

alpha_dose_l = {{'0.15', '0.30'}, {'0.45', '0.60'}, {'0.75', '0.75'}};

for pati = 1:3
    for sessi = 1:2
        filepath = fullfile(data_wd, sprintf('PilotKET0%i\\Sess%i\\',pati, sessi));
        basename = sprintf('PK%i_S%i_exp',pati, sessi);
        
        EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);
        
        for dosei = 1:length(alpha_dose_l)
            [dose_beg, dose_end] = alpha_dose_l{dosei}{:};
            dose_beg_t = find(strcmp({EEG.epoch.dose},dose_beg),1,'first');
            dose_end_t = find(strcmp({EEG.epoch.dose},dose_end),1,'last');
            
            % Take only the data with the same concentration
            EEG_temp = pop_select(EEG, 'trial', dose_beg_t:dose_end_t);
            
            % Change basename
            pop_saveset(EEG_temp,'filename',sprintf('%s_alpha_%s-%s.set', basename, dose_beg, dose_end),...
                'filepath',filepath);
            
            basename =  sprintf('%s_alpha_%s-%s', basename, dose_beg, dose_end);
            
            fprintf('\n*** CALCULATING MULTITAPER SPECTRUM ***\n');
            % calculate power spectrum using the multi-taper method
            calcftspec(basename,filepath);
            
            fprintf('\n*** PLOTTING SPECTRUM ***\n');
            % visualise and save the power spectrum of all channels
            plotftspec(basename,[], filepath);
            
            fprintf('\n*** CALCULATING CONNECTIVITY ***\n');
            % estimate dwPLI connectivity between pairs of channels
            ftcoherence(basename,filepath);
            
            fprintf('\n*** CALCULATING GRAPH-THEORETIC NETWORK METRICS ***\n');
            % calculate graph theory metrics
            calcgraph(basename, 'filepath', filepath);
            
            % %% plot 3D connectivity topographs in the delta, theta and alpha bands.
            plothead(basename,1, 'filepath',filepath);
            plothead(basename,2, 'filepath',filepath);
            plothead(basename,3, 'filepath',filepath)
        end
    end
end

%% To add - Alpha label
trange = ones(1,37); trange(1:4) = 0; trange = logical(trange);
res_alpha_merged_table = table();
res_alpha_merged_table.Properties.VariableNames = {'Var1', 'Var2','Var3','value_alphai'};
for pati = 1:3
    for sessi = 1:2
        filepath = fullfile(data_wd, sprintf('PilotKET0%i\\Sess%i\\Baseline\\',pati, sessi));
        %basename = sprintf('PK%i_S%i_exp',pati, sessi);
        cd(filepath)
        temp_nam = ls(fullfile(filepath,'*mohawk.mat'))
        pause(1)
        load(temp_nam)
        value_alphai = squeeze(mean(std(graphdata{8,2}(3,trange,:),[],3),2));
        res_alpha_merged_table = [res_alpha_merged_table; ...
                table(sprintf("KET0%d", pati), ... %Sub
                drug_sess.drug((pati-1)*2 + sessi),... %Condition
                {'0'},...  % Concentration
                value_alphai)] % alpha
         pause()
    end
end

save('alpha_KET.mat','res_alpha_merged_table')
writetable(res_alpha_merged_table, 'pilotKET_merged_R_alpha.txt')


