% Scripts for revision round 1 to construct power changes over trials
% Based on calcftspec by Srivas Chennu, srivas@gmail.com


% REMINDER: The chonological order of the acquisition is MCS- (1st
% patient), UWS (2nd), and MCS+ (3rd). To show an order, we display UWS,
% MCS-, and MCS+. In the review process, we have been asked to discuss them
% as Case 1, 2, and 3. As such, the chronological order is indees Case 2,
% 1, and 3
%% Measuring things
load freqlist.mat


% Conf info
cfg = [];
cfg.output     = 'pow';
cfg.method     = 'mtmfft';% 'wavelet';%'mtmfft';
cfg.foilim        = [0.5 45];
cfg.taper = 'dpss';
cfg.tapsmofrq = 0.3;
cfg.pad='nextpow2';
cfg.keeptrials  = 'yes';
data_wd = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';
filepath_save = fullfile(data_wd, 'r1_power');
y = [0 45]

% Import data
dose_l = {'0.15', '0.30', '0.45', '0.60', '0.75'};
                


for pati = 1:3
    for sessi = 1:2
        filepath = fullfile(data_wd, sprintf('PilotKET0%i\\Sess%i\\',pati, sessi));
        for dosei = 1:length(dose_l)
            basename = sprintf('PK%i_S%i_exp_%s.set',pati, sessi, dose_l{dosei});
            EEG = pop_loadset(fullfile(filepath, basename));
            %changes_trial = find(cellfun(@isempty, {EEG.epoch.dose}));
            
            EEG = convertoft(EEG);
            EEG = ft_freqanalysis(cfg,EEG);
            
            spectra = EEG.powspctrm;
            freqs = EEG.freq;
            
            
            
            mean_spectra = squeeze(mean(spectra, 2));
%             figure
%             imagesc(1:size(spectra,1), freqs, (mean_spectra'))
%             set(gca,'YDir','normal')
%             colorbar(); title(basename)
%             % Plot changes
%             y = [0 45];
%             hold on 
%             for xi = 1:3
%                 x = [changes_trial(xi+1) changes_trial(xi+1)]; % first is 1
%                 line(x,y, 'Color','red', 'LineWidth', 1.5)
%             end
%             hold off
%             
            mean_spectra_normal = mean_spectra ./ sum(mean_spectra,2);
%             figure
%             imagesc(1:size(spectra,1), freqs, (mean_spectra_normal'))
%             set(gca,'YDir','normal')
%             colorbar(); title([basename, ' norm'])
%             % Plot changes
%             y = [0 45];
%             hold on 
%             for xi = 1:3
%                 x = [changes_trial(xi+1) changes_trial(xi+1)]; % first is 1
%                 line(x,y, 'Color','red', 'LineWidth', 1.5)
%             end
%             hold off
            
            savefile = fullfile(filepath_save,[basename(1:end-4) '.mat']);
            save(savefile, 'freqs', 'mean_spectra', 'mean_spectra_normal', 'spectra');
        end
        
    end
end
       

%% Make figures
cd(filepath_save)

% W/o annotations
% Ketamine
ket = {{1 2}, {2 1}, {3 2}};

figure;
for keti = 1:3
    subplot(3, 1, keti)
    pati = ket{keti}{1}; sessi = ket{keti}{2};
    mat_sess = []; % matrix with all the datapoints
    n_trials = zeros(1, 5);
    for dosei = 1:length(dose_l)
         load(sprintf('PK%i_S%i_exp_%s.mat',pati, sessi, dose_l{dosei}))
         n_trials(dosei) = size(mean_spectra, 1);
         mat_sess = [mat_sess; mean_spectra_normal];
    end
    imagesc(1:size(mat_sess,1), freqs, (mat_sess'));
    set(gca,'YDir','normal')
    colorbar();
    title(sprintf('Case %i', pati))
    caxis([0 0.05])
    
    n_trials = cumsum(n_trials);
    for xi = 1:4
        x = [n_trials(xi)+1 n_trials(xi)+1]; % first is 1
        line(x,y, 'Color','red', 'LineWidth', 1.5)
    end
end
sgtitle('Ketamine normal')
saveas(gcf,'KetNorm.png'); saveas(gcf,'KetNorm.fig')


% Placebo
placebo = {{1 1}, {2 2}, {3 1}};

figure;
for plci = 1:3
    subplot(3, 1, plci)
    pati = placebo{plci}{1}; sessi = placebo{plci}{2};
    mat_sess = []; % matrix with all the datapoints
    n_trials = zeros(1, 5);
    for dosei = 1:length(dose_l)
         load(sprintf('PK%i_S%i_exp_%s.mat',pati, sessi, dose_l{dosei}))
         n_trials(dosei) = size(mean_spectra, 1);
         mat_sess = [mat_sess; mean_spectra_normal];
    end
    imagesc(1:size(mat_sess,1), freqs, (mat_sess'));
    set(gca,'YDir','normal')
    colorbar();
    title(sprintf('Case %i', pati))
    caxis([0 0.05])
    
      n_trials = cumsum(n_trials);
    for xi = 1:4
        x = [n_trials(xi)+1 n_trials(xi)+1]; % first is 1
        line(x,y, 'Color','red', 'LineWidth', 1.5)
    end
end
sgtitle('Placebo normal')
saveas(gcf,'PlcNorm.png'); saveas(gcf,'PlcNorm.fig')

%% W/ annotations
% Set up
y = [0 45];
% Ketamine
ket = {{2 1}, {1 2},  {3 2}};

figure('units','normalized','outerposition',[0 0 0.6 0.6])
for keti = 1:3
    subplot(3, 1, keti)
    pati = ket{keti}{1}; sessi = ket{keti}{2};
    mat_sess = []; % matrix with all the datapoints
    n_trials = zeros(1, 5);
    for dosei = 1:length(dose_l)
         load(sprintf('PK%i_S%i_exp_%s.mat',pati, sessi, dose_l{dosei}))
         n_trials(dosei) = size(mean_spectra, 1);
         mat_sess = [mat_sess; mean_spectra_normal];
    end
    imagesc(1:size(mat_sess,1), freqs, (mat_sess'));
    set(gca,'YDir','normal')
    colorbar();
    title(sprintf('Case %i', keti))
    caxis([0 0.05])
    
    % Annotate text
    text(round(n_trials(1)/2), 40, dose_l{1}, 'Color', 'Red', 'FontSize', 15, ...
    'HorizontalAlignment', 'Center')
    for dosei = 1:length(dose_l)
                text(sum(n_trials(1:dosei-1)) + round(n_trials(dosei)/2),...
                    40, dose_l{dosei}, 'Color', 'Red', 'FontSize', 15, ...
    'HorizontalAlignment', 'Center')
    end
    
    n_trials = cumsum(n_trials);
    for xi = 1:4
        x = [n_trials(xi)+1 n_trials(xi)+1]; % first is 1
        line(x,y, 'Color','red', 'LineWidth', 1.5)
    end
    xlabel('Trial number'); ylabel('Frequency (Hz)')
end
sgtitle('Ketamine - Normalized')
saveas(gcf,'KetNorm_anno.png'); saveas(gcf,'KetNorm_anno.fig')


% Placebo
placebo = {{2 2}, {1 1},  {3 1}};

figure('units','normalized','outerposition',[0 0 0.6 0.6])
for plci = 1:3
    subplot(3, 1, plci)
    pati = placebo{plci}{1}; sessi = placebo{plci}{2};
    mat_sess = []; % matrix with all the datapoints
    n_trials = zeros(1, 5);
    for dosei = 1:length(dose_l)
         load(sprintf('PK%i_S%i_exp_%s.mat',pati, sessi, dose_l{dosei}))
         n_trials(dosei) = size(mean_spectra, 1);
         mat_sess = [mat_sess; mean_spectra_normal];
    end
    imagesc(1:size(mat_sess,1), freqs, (mat_sess'));
    set(gca,'YDir','normal')
    colorbar();
    title(sprintf('Case %i', plci))
    caxis([0 0.03])
    
    % Annotate text
    text(round(n_trials(1)/2), 40, dose_l{1}, 'Color', 'Red', 'FontSize', 15, ...
    'HorizontalAlignment', 'Center')
    for dosei = 2:length(dose_l)
                text(sum(n_trials(1:dosei-1)) + round(n_trials(dosei)/2),...
                    40, dose_l{dosei}, 'Color', 'Red', 'FontSize', 15, ...
    'HorizontalAlignment', 'Center')
    end
    
    n_trials = cumsum(n_trials);
    for xi = 1:4
        x = [n_trials(xi)+1 n_trials(xi)+1]; % first is 1
        line(x,y, 'Color','red', 'LineWidth', 1.5)
    end
    xlabel('Trial number'); ylabel('Frequency (Hz)')
end
sgtitle('Placebo - Normalized')
saveas(gcf,'PlcNorm_anno.png'); saveas(gcf,'PlcNorm_anno.fig')

%% Power spectra
% Setting stuff
fontname = 'Helvetica';
fontsize = 16;
xlims = [0.01 40];
ylims = [-30 10];
bands = {
    'delta'
    'theta'
    'alpha'
    'beta'
    'gamma'
    };

data_moh = 'C:\Users\Paolo\OneDrive - Universite de Liege\Bureau\OldComputer\D\Complexit_doc\Ketamine';
Colmat = jet(5); %redbluecmap(5);% hsv();

% Ketamine
ket = {{2 1}, {1 2},  {3 2}};
placebo = {{2 2}, {1 1},  {3 1}};

figure('units','normalized','outerposition',[0 0 0.6 0.6]);

for keti = 1:3
    subplot(3, 2, 2*keti)
    hold on 
    pati = ket{keti}{1}; sessi = ket{keti}{2};
    mat_sess = []; % matrix with all the datapoints
    for dosei = 1:length(dose_l)
        load(fullfile(data_moh, sprintf('PilotKET0%i',pati), sprintf('Sess%i',sessi),...
            sprintf('PK%i_S%i_exp_%s_mohawk.mat',pati, sessi, dose_l{dosei})))
        plot(freqs,10*log10(mean(spectra)),'LineWidth',2, 'Color', Colmat(dosei,:));
    end
    hold off
    set(gca,'XLim',xlims,'YLim',ylims,'FontSize',fontsize,'FontName',fontname);
    xlabel('Frequency (Hz)','FontSize',fontsize,'FontName',fontname);
    ylabel('Power (dB)','FontSize',fontsize,'FontName',fontname);
    ylimits = ylims;
    for f = 1:5
        line([freqlist(f,1) freqlist(f,1)],ylims,'LineWidth',1,...
            'LineStyle','-.','Color','black');
        line([freqlist(f,2) freqlist(f,2)],ylims,'LineWidth',1,...
            'LineStyle','-.','Color','black');
         text(freqlist(f,1)+2,ylimits(2)+3,...
            sprintf('\\%s',bands{f}),'FontName',fontname,'FontSize',fontsize);
    end
    %title(sprintf('Case %i', pati))
    box off
end


for plci = 1:3
    subplot(3, 2, 2*plci-1)
    hold on 
    pati = placebo{plci}{1}; sessi = placebo{plci}{2};
    mat_sess = []; % matrix with all the datapoints
    for dosei = 1:length(dose_l)
        load(fullfile(data_moh, sprintf('PilotKET0%i',pati), sprintf('Sess%i',sessi),...
            sprintf('PK%i_S%i_exp_%s_mohawk.mat',pati, sessi, dose_l{dosei})))
        plot(freqs,10*log10(mean(spectra)),'LineWidth',2, 'Color', Colmat(dosei,:));
    end
    hold off
    set(gca,'XLim',xlims,'YLim',ylims,'FontSize',fontsize,'FontName',fontname);
    xlabel('Frequency (Hz)','FontSize',fontsize,'FontName',fontname);
    ylabel(sprintf('\\bfCase %i\\rm  \n Power (dB)', plci),'FontSize',fontsize,'FontName',fontname);
    ylimits = ylims;
    for f = 1:5
        line([freqlist(f,1) freqlist(f,1)],ylims,'LineWidth',1,...
            'LineStyle','-.','Color','black');
        line([freqlist(f,2) freqlist(f,2)],ylims,'LineWidth',1,...
            'LineStyle','-.','Color','black');
        text(freqlist(f,1)+2,ylimits(2)+3,...
            sprintf('\\%s',bands{f}),'FontName',fontname,'FontSize',fontsize);
    end
    %title(sprintf('Case %i', pati))
    box off
end

lgd = legend('0.15', '0.30', '0.45', '0.60', '0.75',...
    'Orientation', 'horizontal', 'NumColumns',1);
lgd.Title.String = 'Concentration';
legend('boxoff')
% change position Manually

% Change title of columns + rows
hLF1 = subplot(3,2,1);
text((max(hLF1.XLim)-min(hLF1.XLim))/2+min(hLF1.XLim),max(hLF1.YLim)+3,'Placebo','EdgeColor','none',...
    'FontSize',fontsize+5,'FontName',fontname,'HorizontalAlignment', 'center','VerticalAlignment','Bottom', 'FontWeight', 'bold')

hLF2 = subplot(3,2,2);
text((max(hLF2.XLim)-min(hLF2.XLim))/2+min(hLF2.XLim),max(hLF2.YLim)+3,'Ketamine','EdgeColor','none',...
    'FontSize',fontsize+5,'FontName',fontname,'HorizontalAlignment', 'center','VerticalAlignment','Bottom', 'FontWeight', 'bold')

% h1 = annotation('textbox',[0.25 0.95 0.2 0.05],...
%     'string','Placebo', 'FontSize',fontsize,'FontName',fontname, ...
%     'FontWeight', 'bold');
% h2 = annotation('textbox',[0.65 0.95 0.2 0.05],...
%     'string','Ketamine', 'FontSize',fontsize,'FontName',fontname, ...
%     'FontWeight', 'bold','HorizontalAlignment', 'center');
% set([h1 h2], 'fitboxtotext','on',...
%     'edgecolor','none')

% 
% title('a')
% h1 = text(-0.25, 0.5,'Case 1');
% set(h1, 'rotation', 90)
% text(0.35,1.2,'column 1');
% 
% subplot(3,2,2)
% text(1,1.2,'Ketamine');
% 
% 
% h2 = text(-0.25, 0.5, 'Case');
% set(h2, 'rotation', 90)
% subplot(2,2,4)
% title('d')


%% Changing limits
cd(filepath_save)
saveas(gcf,'PowerSpecto.png'); saveas(gcf,'PowerSpect.fig')
% % Save for the paper
% saveas(gcf,'PowerSpecto.png', 'Resolution',300);

