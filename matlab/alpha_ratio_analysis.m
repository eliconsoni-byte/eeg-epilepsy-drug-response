%%  Alpha ratios & alpha-peak ratios analysis
clear; clc; close all;

% Run from the folder containing DB_TOT_avg.csv
infile  = 'DB_TOT_avg.csv';
outdir  = fullfile(pwd,'fig_alpha_ratios_all');  
if ~exist(outdir,'dir'); mkdir(outdir); end

% Load dataset 
data = readtable(infile);

if ~ismember('outcome', data.Properties.VariableNames)
    error('Column "outcome" not found.');
end

% Split into responder/non-responder
resp_data    = data(data.outcome == 1, :);
nonresp_data = data(data.outcome == 0, :);

% Select alpha ratio columns
vnames = data.Properties.VariableNames;

is_raw_alpha_ratio      = startsWith(vnames,'spectrum_alpha_')       & contains(vnames,'_ratio');
is_raw_alpha_peak_ratio = startsWith(vnames,'spectrum_alpha_peak_')  & contains(vnames,'_ratio');

is_log_alpha_ratio      = startsWith(vnames,'spectrum_log_alpha_')       & contains(vnames,'_ratio');
is_log_alpha_peak_ratio = startsWith(vnames,'spectrum_log_alpha_peak_')  & contains(vnames,'_ratio');

is_mag_alpha_ratio      = startsWith(vnames,'spectrum_mag_alpha_')       & contains(vnames,'_ratio');
is_mag_alpha_peak_ratio = startsWith(vnames,'spectrum_mag_alpha_peak_')  & contains(vnames,'_ratio');

feat_mask = is_raw_alpha_ratio | is_raw_alpha_peak_ratio | ...
            is_log_alpha_ratio | is_log_alpha_peak_ratio | ...
            is_mag_alpha_ratio | is_mag_alpha_peak_ratio;

feat_names = vnames(feat_mask);

if isempty(feat_names)
    error('No alpha ratio / alpha-peak ratio columns found in RAW/LOG/MAG.');
end

alpha_levels = [0.001, 0.01, 0.05];                  
colors.resp    = [0.20 0.60 0.20];                   
colors.nonresp = [0.85 0.20 0.20];                   
gray = [0.50 0.50 0.50];
cat_order = {'Responders','Non-responders'};

% Collector for p-values
summary_rows = [];

for k = 1:numel(feat_names)
    feat = feat_names{k};
    x1 = resp_data.(feat);
    x2 = nonresp_data.(feat);

    % Drop NaNs
    x1 = x1(~isnan(x1));
    x2 = x2(~isnan(x2));

    if isempty(x1) || isempty(x2)
        warning('Skipping %s: one group is empty after removing NaNs.', feat);
        continue;
    end

    % Wilcoxon rank-sum test
    p = ranksum(x1, x2);

    % stars for significance levels
    stars = '';
    if p < alpha_levels(1)
        stars = '***';
    elseif p < alpha_levels(2)
        stars = '**';
    elseif p < alpha_levels(3)
        stars = '*';
    end

    % Plot 
    f = figure('Color','w','Position',[100 100 600 460]); hold on;

   
    cx1 = categorical(repmat(cat_order(1), numel(x1), 1), cat_order);
    cx2 = categorical(repmat(cat_order(2), numel(x2), 1), cat_order);

    bc1 = boxchart(cx1, x1); bc1.BoxFaceColor = colors.resp;    bc1.MarkerColor = colors.resp;
    bc2 = boxchart(cx2, x2); bc2.BoxFaceColor = colors.nonresp; bc2.MarkerColor = colors.nonresp;

    
    ax = gca; ax.FontSize = 11; ax.LineWidth = 1;
    ylabel(strrep(feat,'_',' '), 'Interpreter','none');               % y-axis label: feature name
    title(strrep(feat,'_',' '), 'Interpreter','none','FontWeight','bold');  % title: feature name
    box on; grid on;

    % Significance bar + stars 
    if ~isempty(stars)
        ymax  = max([x1; x2]);
        ymin  = min([x1; x2]);
        yspan = max(ymax - ymin, 1e-6);
        ybar  = ymax + 0.08*yspan;       
        ytxt  = ybar + 0.03*yspan;       

        plot([1 2], [ybar ybar], 'Color', gray, 'LineWidth', 1.5);
        text(1.5, ytxt, stars, 'HorizontalAlignment','center', 'FontSize', 14, 'Color', gray);
    end

    % Legend
    hResp = plot(nan, nan, 's', 'MarkerFaceColor', colors.resp,    'MarkerEdgeColor', colors.resp,    'LineStyle','none', 'DisplayName','Responders');
    hNon  = plot(nan, nan, 's', 'MarkerFaceColor', colors.nonresp, 'MarkerEdgeColor', colors.nonresp, 'LineStyle','none', 'DisplayName','Non-responders');
    legend([hResp hNon], 'Location','best');

    % Save figure
    safe_name = regexprep(feat,'[^a-zA-Z0-9]','_');
    saveas(f, fullfile(outdir, sprintf('box_%s.png', safe_name)));
    close(f);

    % collect summary
    summary_rows = [summary_rows; {feat, p, numel(x1), numel(x2)}]; %#ok<AGROW>
end

% Summary table
if ~isempty(summary_rows)
    T = cell2table(summary_rows, 'VariableNames', {'feature','p_value','n_resp','n_nonresp'});
    T = sortrows(T, 'p_value');
    writetable(T, fullfile(outdir,'alpha_ratios_wilcoxon_summary.csv'));
    fprintf('Done.\nFigures saved in: %s\n', outdir);
    fprintf('Summary CSV saved: %s\n', fullfile(outdir,'alpha_ratios_wilcoxon_summary.csv'));
else
    warning('No features processed.');
end
