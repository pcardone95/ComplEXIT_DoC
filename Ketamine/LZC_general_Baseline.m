%% Script
% Go to the data folder
pathwork = '...';
cd(pathwork)
to_check = {};

results = struct();
name_fold = ls("Pilot*"); % If you have multiple sessions
for foldi = 1:size(name_fold,1)
    cd(fullfile(pathwork,name_fold(foldi,:)))
    sessl = ls("Sess*");
    for sessi = 1:size(sessl,1)
        cd(fullfile(pathwork,name_fold(foldi,:), sessl(sessi,:), 'Baseline'))
        % Import data - BASED ON MOHAWK
        files_list = ls("*.set");
        EEG = pop_loadset(strrep(files_list(1,:),' ','')); %import data
        
        if not(strcmp('averef',EEG.ref))
            to_check = [to_check, {[name_fold(foldi,:), '_', sessl(sessi,:)]}];
            continue
        end
        % Re-order the data
        % After interpolation, the order of the channels is changed - here we
        % put them again in order
        temp = struct2table(EEG.chanlocs); % convert the struct array to a table

        [~ , index] = sort(lower(temp.labels)); % Lower used in case of strange problems with 
%         EEG.chanlocs = EEG.chanlocs(index);
%         EEG.data = EEG.data(index);
        
        % LZC CODE FROM HERE
        data = EEG.data(index,:,:); % channel x time x trial
        data_binar = zeros(size(data));
        
        % Preallocate
        strings = cell(size(data,1),size(data,3)); % electrode x trial
        C = zeros(size(data,1),size(data,3));
        H = cell(size(data,1),size(data,3));
        
        % Calculate
        for electi = 1:size(data,1)
            for triali = 1:size(data,3)
                data_binar(electi,:,triali) = abs(hilbert(data(electi,:,triali)));
                data_binar(electi,:,triali) = ...
                    data_binar(electi,:,triali) > mean(data_binar(electi,:,triali));
                strings{electi, triali} = binary_seq_to_string(data_binar(electi,:,triali));
                [C(electi, triali), H{electi, triali}] = calc_lz_complexity(data_binar(electi,:,triali), "exhaustive", 1);
            end
        end
        i = 2*(foldi-1)+sessi;
        results(i).C = C; results(i).H = H; %results(folderi).strings = strings;
        results(i).C_mean = mean(results(i).C,2);
        results(i).label = [name_fold(foldi,:), '_', sessl(sessi,:)];    
    end
end
cd(pathwork)
save('LZC_baseline.mat', 'results')


%% Convert in table 
% Could have done it before, but whatever
EEG.chanlocs = EEG.chanlocs(index);

LZC_table = table();
for i = 1:6
    labels = strsplit(results(i).label,'_');
    LZC_table = [LZC_table; ...
        table(repmat(labels{1}(6:end), 127,1), repmat(labels{2}(end), 127,1),repmat('0', 127,1),...
        {EEG.chanlocs.labels}', results(i).C_mean)];
end
LZC_table.Properties.VariableNames = {'Subject', 'Session', 'Concentration', 'Electrode', 'LZC'};
cd(pathwork)
writetable(LZC_table, 'LZC_KET_baseline')
save('LZC_table.mat', 'LZC_table')
