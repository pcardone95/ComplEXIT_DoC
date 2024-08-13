%% merge_KET_LZC
% This scripts looks for the data in the LZC files and merge them in bunch
% of range of doses: we'll gather 0.15-0.3, 0.45 - 0.6, 0.75.
% This is done to allow the same number of points with the mohawk

% Prepare info
results_merged = table();
cur_wd = pwd();
cd ..
data_wd = pwd();
drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02"], "drug",...
    ["Placebo", "Ketamine", "Ketamine", "Placebo"]);
results_merged_table = table();
elecs = 127; % numb electrodes

% Loop throught the data
fold_pat = ls("Pilot*");
for pati = 1:size(fold_pat,1) % Microsoft
    cd(fullfile(data_wd, fold_pat(pati,:)))
    sess = ls('Sess*');
    for sessi = 1:2
         cd(fullfile(data_wd, fold_pat(pati,:), sess(sessi,:)))
         load(ls('*LZC.mat'))
         
         results_merged_table = [results_merged_table; ...
             table(repmat(sprintf('KET0%d', pati),127*3,1), ... %Sub
             repmat(drug_sess.drug((pati-1)*2 + sessi),127*3,1),... %Condition
             sort(repmat({'0.15-0.30'; '0.75'; '0.40-0.65'}, 127, 1)),...  % Dose
             repmat((1:127)',3,1),... % Electrodes
             [mean([results(1).C, results(2).C],2); ...
             mean([results(3).C, results(4).C],2);  results(5).C_mean])] % LZC
         % pause()
    end
end

results_merged_table.Properties.VariableNames = {'Subject', 'Condition',...
    'Dose', 'Electrode', 'LZC'};

%% Save data
cd(fullfile(data_wd,'results'))
writetable(results_merged_table, 'pilotKET_merged_R.txt')
save('pilotKET_merged.mat', 'results_merged_table')

%% Adding third pilot

% Loop throught the data
path_data_pilot3 = 'D:\Complexit_doc\Ketamine\PilotKET03';
cd(path_data_pilot3)
sess = ls('Sess*');
pati = 3;
drug_sess = struct("sub", ["KET01", "KET01", "KET02", "KET02", "KET03", "KET03"],...
    "drug", ["Placebo", "Ketamine", "Ketamine", "Placebo", "Placebo", "Ketamine"]);
for sessi = 1:2
    cd(fullfile(path_data_pilot3, sess(sessi,:)))
    load(ls('*LZC.mat'))
    switch sessi
        case 1
            temp = table(repmat(sprintf('KET0%d', pati),127*3,1), ... %Sub
                repmat(drug_sess.drug((pati-1)*2 + sessi),127*3,1),... %Condition
                sort(repmat({'0.15-0.30'; '0.75'; '0.40-0.65'}, 127, 1)),...  % Dose
                repmat((1:127)',3,1),... % Electrodes
                [results(1).C_mean; ...
                mean([results(2).C, results(3).C],2);  results(4).C_mean]);% LZC
            temp.Properties.VariableNames = {'Subject', 'Condition',...
    'Dose', 'Electrode', 'LZC'};
            results_merged_table = [results_merged_table; temp]
            % pause()
        case 2
            temp = table(repmat(sprintf('KET0%d', pati),127*3,1), ... %Sub
             repmat(drug_sess.drug((pati-1)*2 + sessi),127*3,1),... %Condition
             sort(repmat({'0.15-0.30'; '0.75'; '0.40-0.65'}, 127, 1)),...  % Dose
             repmat((1:127)',3,1),... % Electrodes
             [mean([results(1).C, results(2).C],2); ...
             mean([results(3).C, results(4).C],2);  results(5).C_mean]);% LZC
            temp.Properties.VariableNames = {'Subject', 'Condition',...
    'Dose', 'Electrode', 'LZC'};
            results_merged_table = [results_merged_table; temp]
    end
end
cd('D:\Complexit_doc\Ketamine\results')
writetable(results_merged_table, 'pilotKET_merged_R_complete.txt')
save('pilotKET_merged_complete.mat', 'results_merged_table')


