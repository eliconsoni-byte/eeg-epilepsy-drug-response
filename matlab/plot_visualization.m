% Requires the workspace variable psd_exp: a Brainstorm PSD export (fields TF and Freqs) of one patient.
% The second section plots the numbers of significant features, typed in by hand from the output of
% statistical_analysis_2.m.
% Visualization of PSD and its transforms
PSD = squeeze(psd_exp.TF);
Freqs = psd_exp.Freqs;


% Log and magnitude transforms definition
PSD_log = 10 * log10(PSD);   
PSD_mag = sqrt(PSD);          

% Subplot with the three transformations
figure('Color', 'w', 'Position', [100 100 700 600]);

subplot(3,1,1);
plot(Freqs, PSD, 'k', 'LineWidth', 1.2);
title('Raw Power Spectral Density');
ylabel('Power (\muV^2/Hz)');
xlim([0 70]);
grid on;

subplot(3,1,2);
plot(Freqs, PSD_log, 'b', 'LineWidth', 1.2);
title('Log-Transformed PSD (10·log_{10})');
ylabel('Power (dB)');
xlim([0 70]);
grid on;

subplot(3,1,3);
plot(Freqs, PSD_mag, 'r', 'LineWidth', 1.2);
title('Magnitude-Transformed PSD (\surd)');
ylabel('Amplitude (\muV/\surdHz)');
xlabel('Frequency (Hz)');
xlim([0 70]);
grid on;

sgtitle('Comparison of PSD Transformations', 'FontWeight', 'bold');


%% Histograms of significant features
clear; 
clc; 
close all;

bands      = {'\delta','\theta','\alpha','\beta','\gamma','\it spike', '\it TOT.'};
nb         = numel(bands);


% Non-normalized 
cnt_nonNorm_PSD = [0 2 0 0 2 0];
cnt_nonNorm_LOG = [4 3 0 2 3 2];
cnt_nonNorm_MAG = [3 3 0 2 3 2];

% Normalized
cnt_norm_PSD = [1 0 0 0 0 0];
cnt_norm_LOG = [0 0 0 0 0 1];
cnt_norm_MAG = [0 0 0 0 0 1];

% total column
cnt_nonNorm_PSD(end+1) = sum(cnt_nonNorm_PSD);
cnt_nonNorm_LOG(end+1) = sum(cnt_nonNorm_LOG);
cnt_nonNorm_MAG(end+1) = sum(cnt_nonNorm_MAG);

cnt_norm_PSD(end+1) = sum(cnt_norm_PSD);
cnt_norm_LOG(end+1) = sum(cnt_norm_LOG);
cnt_norm_MAG(end+1) = sum(cnt_norm_MAG);

% colours for each transformation
colPSD = [0.93 0.69 0.13];  % yellow
colLOG = [0.20 0.70 0.30];  % green
colMAG = [0.45 0.20 0.60];  % purple

% Helper function to plot one figure
plotGrouped = @(M, ttl, fname) plot_grouped_hist(M, bands, ...
                         {'PSD','LOG(PSD)','MAG(PSD)'}, ...
                         {colPSD, colLOG, colMAG}, ttl, fname);

% Figure 1: Non-normalized 
M1 = [cnt_nonNorm_PSD(:), cnt_nonNorm_LOG(:), cnt_nonNorm_MAG(:)];
plotGrouped(M1, 'Non-Normalized Significant Features', 'hist_non_normalized.png');

% --- Figure 2: Normalized ---
M2 = [cnt_norm_PSD(:), cnt_norm_LOG(:), cnt_norm_MAG(:)];
plotGrouped(M2, 'Normalized Significant Features', 'hist_normalized.png');

%% Local function
function plot_grouped_hist(M, band_labels, series_labels, series_colors, titleStr, outName)
    fig = figure('Color','w','Position',[120 120 840 420]);
    b = bar(M, 'grouped'); grid on; box on;
    for k = 1:numel(b)
        b(k).FaceColor = series_colors{k};
        b(k).EdgeColor = 'none';
        b(k).BarWidth  = 1;
    end
    set(gca, 'XTick', 1:numel(band_labels), 'XTickLabel', band_labels, ...
             'TickDir','out', 'FontSize', 12);
    xlabel('Frequency band','FontSize',13,'Interpreter','tex');
    ylabel('Number of features','FontSize',13);
    title(titleStr, 'FontWeight','bold');
    legend(series_labels, 'Location','north','Orientation','horizontal','Box','on');
    set(gca, 'FontSize',20)
    % Annotate counts above bars 
    for k = 1:numel(b)
        xk = b(k).XEndPoints;
        yk = b(k).YEndPoints;
        text(xk, yk + 0.1, string(M(:,k)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
    end

    ylim([0, max(M(:))+1]);
    exportgraphics(fig, outName, 'Resolution', 300);
end

