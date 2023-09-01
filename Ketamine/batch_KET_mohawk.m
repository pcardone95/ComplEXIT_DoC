function batch_KET_mohawk(patients_number, typenalysis)
%% Batch KET Mohawk
% This code aims to perform the KET mohawk for several patients at the same
% time. This is done as well for several sessions.

if nargin <1
    warning('Please put an array of intergers for the patients')
end

typenalysis_list = {'fullanalysis','UptoTrialRejection','ICA',...
    'ICArej_toReref', 'LZC_Alpha', 'ECI'};

if ~exist('typenalysis','var') || any(strcmp('a',typenalysis_list))
    [idx, tf]=listdlg('ListString', typenalysis_list,...
        'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');
    
    if tf
        typenalysis = typenalysis_list{idx};
    else
        warning('Please select one analysis');
        return
    end
end


data_wd = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';
for pati = patients_number
    for sessi = 1:2
        filepath = fullfile(data_wd, sprintf('PilotKET0%i\\Sess%i\\',pati, sessi));
        basename = sprintf('PK%i_S%i_exp',pati, sessi);
        
        switch typenalysis
            case 'fullanalysis'
                KET_mohawk(basename, filepath);
            case 'UptoTrialRejection'
                % Epoch data
                fprintf('\n*** EPOCHING DATA ***\n');
                epochdata(basename, filepath);
                
                % Parse data
                fprintf('\n*** PARSE DATA IN DIFFERENT CONCENTRATION ***\n');
                parsedata_KET(basename, filepath);
                
                %% REJECT DATA DURING SECONDs
                fprintf('\n DELETE TRIALS DURING SECONDs...\n');
                rejectSECONDs_KET([basename '_epochs.set'], filepath)
                
                % Trial rejection
                fprintf('\n*** SELECT BAD CHANNELS AND TRIALS ***\n');
                rejartifacts([basename '_epochs'], 1, 4, 1, [], [], [], filepath);
            case 'ICA'
                fprintf('\n*** COMPUTING IC DECOMPOSITION ***\n');
                % Run ICA decomposition with optional PCA pre-processing
                computeic([basename '_epochs'], [], [], filepath)
            case 'ICArej_toReref'
                fprintf('\n*** SELECT BAD ICs ***\n');
                % Visually identify and reject noisy ICs, e.g., eye movements, muscle
                % activity, etc.
                rejectic(basename, 'prompt', 'off', 'filepath', filepath)
               
                
                %% MANUAL
                % Second and final pass of quasi-automatic rejection of noisy channels and epochs based on variance
                % thresholding, to remove any remaining noisy data.
                fprintf('\n*** SELECT ANY REMAINING BAD CHANNELS AND TRIALS AND INTERPOLATE ***\n');
                rejartifacts([basename '_clean'], 2, 4, 0, [], 500, 250, filepath);
                %% Manual
                fprintf('\n*** REFERENCING DATA TO COMMON AVERAGE ***\n');
                % re-reference data to common average for connectivity estimation.
                rereference(basename, 1, [],[], filepath);
            case 'LZC_Alpha'
                dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};
                fprintf('\n*** Measuring complexity and dividing in sessions ***\n');
                calc_LZC_KET(basename, filepath)
                for dosei = 1:length(dose_l)
                    dose_s = dose_l{dosei};
                    
                    fprintf('\n*** CALCULATING MULTITAPER SPECTRUM ***\n');
                    % calculate power spectrum using the multi-taper method
                    calcftspec(sprintf('%s_%s', basename, dose_s),filepath);
                    
                    fprintf('\n*** PLOTTING SPECTRUM ***\n');
                    % visualise and save the power spectrum of all channels
                    plotftspec(sprintf('%s_%s', basename, dose_s),[], filepath);
                    
                    fprintf('\n*** CALCULATING CONNECTIVITY ***\n');
                    % estimate dwPLI connectivity between pairs of channels
                    ftcoherence(sprintf('%s_%s', basename, dose_s),filepath);
                    
                    fprintf('\n*** CALCULATING GRAPH-THEORETIC NETWORK METRICS ***\n');
                    % calculate graph theory metrics
                    calcgraph(sprintf('%s_%s', basename, dose_s), 'filepath', filepath);
                    
                    % %% plot 3D connectivity topographs in the delta, theta and alpha bands.
                    plothead(sprintf('%s_%s', basename, dose_s),1, 'filepath',filepath);
                    plothead(sprintf('%s_%s', basename, dose_s),2, 'filepath',filepath);
                    plothead(sprintf('%s_%s', basename, dose_s),3, 'filepath',filepath);
                    %
                    
                end
            case 'ECI'
                parse_data_4ECI_temp(basename, filepath)
        end
    end
end
