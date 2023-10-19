%% merge_KET_LZC
% This scripts looks for the data in the LZC files and merge them in bunch
% of range of doses: we'll gather 0.15-0.3, 0.45 - 0.6, 0.75.
% This is done to allow the same number of points with the mohawk

% Prepare info
results_merged = table();
cur_wd = pwd();
data_wd = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';% pwd();
cd(data_wd)
drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02", "KET03", "KET03"],...
    "drug", ["Placebo", "Ketamine", "Ketamine", "Placebo", "Placebo", "Ketamine"]);
results_merged_table = table();
elecs = 127; % numb electrodes

% % Loop throught the data
% fold_pat = ls("Pilot*");
% for pati = 1:size(fold_pat,1) % Microsoft
%     cd(fullfile(data_wd, fold_pat(pati,:)))
%     sess = ls('Sess*');
%     for sessi = 1:2
%          cd(fullfile(data_wd, fold_pat(pati,:), sess(sessi,:)))
%          load(ls('*LZC.mat'))
%          
% %          results_merged_table = [results_merged_table; ...
% %              table(repmat(sprintf('KET0%d', pati),127*3,1), ... %Sub
% %              repmat(drug_sess.drug((pati-1)*2 + sessi),127*3,1),... %Condition
% %              sort(repmat({'0.15-0.30'; '0.75'; '0.40-0.65'}, 127, 1)),...  % Dose
% %              repmat((1:127)',3,1),... % Electrodes
% %              [mean([results(1).C, results(2).C],2); ...
% %              mean([results(3).C, results(4).C],2);  results(5).C_mean])] % LZC
% 
%          results_merged_table = [results_merged_table; ...
%              table(repmat(sprintf('KET0%d', pati),127*5,1), ... %Sub
%              repmat(drug_sess.drug((pati-1)*2 + sessi),127*5,1),... %Condition
%              sort(repmat({'0.15'; '0.30'; '0.45'; '0.60'; '0.75'}, 127, 1)),...  % Dose
%              repmat((1:127)',5,1),... % Electrodes
%              [results(1).C_mean; results(2).C_mean; results(3).C_mean;...
%              results(4).C_mean; results(5).C_mean])]; % LZC
%          % pause()
%          
%          % Baseline
%          cd(fullfile(data_wd, fold_pat(pati,:), sess(sessi,:),'Baseline'))
%          load(ls('*LZC.mat'))
%     end
% end
% 
% results_merged_table.Properties.VariableNames = {'Subject', 'Condition',...
%     'Dose', 'Electrode', 'LZC'};
% 
% %% Save data
% cd(fullfile(data_wd,'results'))
% writetable(results_merged_table, 'pilotKET_merged_R.txt')
% save('pilotKET_merged.mat', 'results_merged_table')

%% Version 2
data_wd = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';% pwd();
cd(data_wd)
drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02", "KET03", "KET03"],...
    "drug", ["Placebo", "Ketamine", "Ketamine", "Placebo", "Placebo", "Ketamine"]);
results_merged_table = table();
elecs = 127; % numb electrodes

temp_all_merged = table();
temp_PerConc_merged = table();
for pati = 1:3
    for sessi = 1:2
        cd(fullfile(data_wd, sprintf('PilotKET0%i\\Sess%i\\',pati, sessi)));
        temp_all = readtable(sprintf('PK%i_S%i_exp_LZC_all',pati, sessi));
        temp_PerConc = readtable(sprintf('PK%i_S%i_exp_LZC_PerConcentration',pati, sessi));
        
        temp_all.Condition = repmat(drug_sess.drug((pati-1)*2 + sessi),127,1);
        temp_PerConc.Condition = repmat(drug_sess.drug((pati-1)*2 + sessi),127*5,1);
        
        temp_all_merged = [temp_all_merged; temp_all]; 
        temp_PerConc_merged = [temp_PerConc_merged; temp_PerConc];
    end
end

cd('C:/Users/Paolo/OneDrive - Universite de Liege/Bureau/OldComputer/D/Complexit_doc/Ketamine/results')
writetable(temp_all_merged, 'pilotKET_merged_R_V3_all.txt')
writetable(temp_PerConc_merged, 'pilotKET_merged_R_V3_PerConc.txt')
save('pilotKET_merged_V3.mat', 'temp_all_merged', 'temp_PerConc_merged')