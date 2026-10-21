clear all
clc
close all
% Run from the folder containing DB_TOT_avg.csv

% Load the dataset
data = readtable('DB_TOT_avg.csv');

% Extract numeric feature columns
outcome = data.outcome;
numericVars = varfun(@isnumeric, data, 'OutputFormat', 'uniform');
numericVars(strcmp(data.Properties.VariableNames, 'outcome')) = false;
features = data(:, numericVars);
featureNames = features.Properties.VariableNames;

% Define frequency bands
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'spike'};

% Split into responsive and non-responsive groups
resp_data = data(outcome == 1, :);
nonresp_data = data(outcome == 0, :);

% Initialize empty cell arrays to store the divided datasets
resp_bands = struct();
nonresp_bands = struct();

% Divide each group based on frequency band and exclude alpha ratios
for b = 1:length(bands)
    % Initialize an array to store features for the current band
    band_feats_resp = [];
    band_feats_nonresp = [];
    
    for i = 1:numel(featureNames)
        fname = featureNames{i};
        
        % Include features specific to the band (excluding the alpha ratios)
        if contains(fname, ['_' bands{b}]) && ~contains(fname, 'ratio')
            % For responders
            band_feats_resp = [band_feats_resp, resp_data.(fname)];
            % For non-responders
            band_feats_nonresp = [band_feats_nonresp, nonresp_data.(fname)];
        end
    end
    
    % Store the datasets for each band
    resp_bands.(bands{b}) = band_feats_resp;
    nonresp_bands.(bands{b}) = band_feats_nonresp;
end
%% KRUSKAL-WALLIS TEST 
fprintf('\n--- Kruskal-Wallis Test Within Groups (Across Frequency Bands) ---\n');

groups = {'Responders', 'Non-Responders'};
for g = 1:2
    if g == 1
        group_data = resp_bands;
    else
        group_data = nonresp_bands;
    end
    
    % Concatenate all features across bands
    all_vals = [];
    all_labels = {};
    
    for b = 1:length(bands)
        vals = group_data.(bands{b});
        if ~isempty(vals)
            all_vals = [all_vals; vals(:)];
            all_labels = [all_labels; repmat(bands(b), numel(vals), 1)];
        end
    end
    
    labels_cat = categorical(all_labels, bands, 'Ordinal', true);
    
    % Kruskal-Wallis test
    [p_kw, ~, stats] = kruskalwallis(all_vals, labels_cat, 'off');
    fprintf('\nKruskal-Wallis for %s:\n', groups{g});
    fprintf('  p = %.5e → %s\n', p_kw, ...
        ternary(p_kw < 0.05, ...
        'Statistically significant differences across frequency bands', ...
        'No statistically significant differences across frequency bands'));
    
    % Post-hoc comparisons if significant
    if p_kw < 0.05
        figure('Name', ['Post-hoc - ' groups{g}], 'Color', 'w');
        c = multcompare(stats, 'CType', 'bonferroni', 'Display', 'on');
        groupnames = stats.gnames;
        
        fprintf('  Post-hoc comparisons (Bonferroni corrected):\n');
        fprintf('  Group 1\tGroup 2\t\tDifference\tp-value\t\tSignificant\n');
        for i = 1:size(c,1)
            g1 = groupnames{c(i,1)};
            g2 = groupnames{c(i,2)};
            diff = c(i,4);
            pval = c(i,6);
            sig = ternary(pval < 0.05, 'YES', 'NO');
            fprintf('  %-10s\t%-10s\t%8.3f\t%8.3e\t%s\n', g1, g2, diff, pval, sig);
        end
                                                              
        [~, idxInBands] = ismember(lower(groupnames), lower(bands));
        nb = numel(bands);
        P = NaN(nb);
        for i = 1:size(c,1)
            g1 = idxInBands(c(i,1)); g2 = idxInBands(c(i,2));
            if g1==0 || g2==0, continue; end
            P(g1,g2) = c(i,6); P(g2,g1) = c(i,6);
        end

        greek_labels = {'\delta','\theta','\alpha','\beta','\gamma','\it spike'};
        sigMask    = P < 0.05;
        nonsigMask = ~isnan(P) & ~sigMask;

        % Figure
        fig1 = figure('Color','w','Name',['Post-hoc p-values - ' groups{g}]); hold on
        axis([0.5 nb+0.5 0.5 nb+0.5]); axis square
        set(gca,'XTick',1:nb,'XTickLabel',greek_labels, ...
                'YTick',1:nb,'YTickLabel',greek_labels, ...
                'TickDir','out','FontSize',11,'TickLabelInterpreter','tex');
        xtickangle(45)
        set(gca,'YDir','normal')

       
        gray = [0.8 0.8 0.8];
        for r = 1:nb
            for c2 = 1:nb
                if nonsigMask(r,c2)
                    rectangle('Position',[c2-0.5, r-0.5, 1, 1], ...
                              'FaceColor', gray, 'EdgeColor','none');
                end
            end
        end

        S = NaN(nb);
        S(sigMask) = P(sigMask);

        nGrad = 512;  
        x = linspace(0,1,nGrad).';
      
        grad = [ ...
            (0.*(1-x) + 1.*x), ...  
            (0.*(1-x) + 1.*x), ... 
            (1.*(1-x) + 0.*x)  ... 
        ];

        himg = imagesc(S, [0 0.05]);
        colormap(fig1, grad);
        set(himg, 'AlphaData', ~isnan(S));

       
        for d = 1:nb
            rectangle('Position',[d-0.5, d-0.5, 1, 1], ...
                      'FaceColor','w','EdgeColor','none');
        end

       
        cb = colorbar;
        cb.Label.String = 'p-value (celle significative)';
        cb.Ticks = [0 0.01 0.02 0.03 0.04 0.05];
        cb.Limits = [0 0.05];

        
        glabel = lower(groups{g});
        title(['Post-hoc pairwise comparisons - ' glabel]);

        % Save
        outname = ['heatmap_non_norm_posthoc_p_' lower(strrep(groups{g},'-','')) '.png'];
        exportgraphics(fig1, outname, 'Resolution', 300);


    % Boxplot across bands
    figure('Name', ['Boxplot - ' groups{g}], 'Color', 'w');
    boxplot(all_vals, labels_cat);
    title(['Frequency Band Comparison - ' groups{g}]);
    ylabel('Feature Value');
    exportgraphics(gcf, 'boxplot_nonresponders.emf', 'ContentType', 'vector');

    end
end

%% WILCOXON RANK-SUM TEST 
fprintf('\n--- Wilcoxon Rank-Sum Test Between Groups (Per Feature Per Band) ---\n');

wilcoxon_pvals = struct(); 

for b = 1:length(bands)
    vals_resp = resp_bands.(bands{b});
    vals_nonresp = nonresp_bands.(bands{b});

    % Find feature names for the current band (excluding ratios)
    band_feature_names = {};
    for i = 1:numel(featureNames)
        fname = featureNames{i};
        if contains(fname, ['_' bands{b}]) && ~contains(fname, 'ratio')
            band_feature_names{end+1} = fname; 
        end
    end

    if isempty(band_feature_names)
        continue;
    end

    fprintf('\nTested band: %s\n', bands{b});
    fprintf('Feature Name                         p-value        Significant\n');
    fprintf('---------------------------------------------------------------\n');

    % Store p-values
    wilcoxon_pvals.(bands{b}) = containers.Map();

    for j = 1:numel(band_feature_names)
        [p, ~] = ranksum(vals_resp(:, j), vals_nonresp(:, j));
        sig_msg = ternary(p < 0.05, 'YES', 'NO');
        fprintf('%-35s  %.5e    %s\n', band_feature_names{j}, p, sig_msg);

        % Save for plotting
        wilcoxon_pvals.(bands{b})(band_feature_names{j}) = p;
    end
end

%% BOXCHART VISUALIZATION WITH SIGNIFICANCE 

outDir = 'statistical_graphs';
if ~exist(outDir,'dir')
    mkdir(outDir); 
end

% Define 12 base features (prefix only)
base_features = { ...
    'spectrum_log_max_', ...
    'spectrum_log_maxfreq_', ...
    'spectrum_log_', ...
    'spectrum_log_median_', ...
    'spectrum_mag_max_', ...
    'spectrum_mag_maxfreq_', ...
    'spectrum_mag_', ...
    'spectrum_mag_median_', ...
    'spectrum_max_', ...
    'spectrum_maxfreq_', ...
    'spectrum_', ...
    'spectrum_median_'};

bands = {'delta','theta','alpha','beta','gamma','spike'};
band_labels = {'\delta', '\theta', '\alpha', '\beta', '\gamma', '\it spike'};

for f = 1:length(base_features)
    base_name = base_features{f};

    % Containers
    feature_vals_all = [];
    group_labels = {};
    band_labels_all = [];

    for b = 1:length(bands)
        band = bands{b};
        full_feature = [base_name band];

        % Check if the feature exists
        if ismember(full_feature, data.Properties.VariableNames)
            r_vals  = resp_data.(full_feature);
            nr_vals = nonresp_data.(full_feature);

            % Aggregate values
            feature_vals_all = [feature_vals_all; r_vals; nr_vals];
            group_labels = [group_labels; ...
                repmat({'Responder'},      numel(r_vals), 1); ...
                repmat({'Non-responder'}, numel(nr_vals), 1)];
            band_labels_all = [band_labels_all; ...
                repmat(b, numel(r_vals) + numel(nr_vals), 1)];
        end
    end

    if isempty(feature_vals_all)
        continue; 
    end

    % Figure
    fig = figure('Name', ['Feature: ' strrep(base_name, '_', '\_')], 'Color', 'w');
    hold on

    ymin = min(feature_vals_all);
    ymax = max(feature_vals_all);
    yrange = ymax - ymin;
    ylim([ymin - 0.1*yrange, ymax + 0.15*yrange]);

   
    for b = 1:length(bands)
        band = bands{b};
        full_feature = [base_name band];

        % Logical indexing
        idx_r  = band_labels_all == b & strcmp(group_labels, 'Responder');
        idx_nr = band_labels_all == b & strcmp(group_labels, 'Non-responder');

        x_r  = 2*b - 1;
        x_nr = 2*b;

        % Plot boxcharts
        if any(idx_r)
            boxchart(ones(sum(idx_r),1)*x_r, feature_vals_all(idx_r), ...
                'BoxFaceColor', [0.2 0.6 0.2], 'MarkerStyle','none');
        end
        if any(idx_nr)
            boxchart(ones(sum(idx_nr),1)*x_nr, feature_vals_all(idx_nr), ...
                'BoxFaceColor', [0.8 0.2 0.2], 'MarkerStyle','none');
        end

        % Significance marker (within band R vs NR)
        if isfield(wilcoxon_pvals, band) && isKey(wilcoxon_pvals.(band), full_feature)
            p = wilcoxon_pvals.(band)(full_feature);
            if p < 0.05
                y_star = ymax + 0.05*yrange;
                plot([x_r, x_nr], [y_star y_star], 'k-', 'LineWidth', 1.2);
                text(mean([x_r x_nr]), y_star, '*', ...
                    'HorizontalAlignment','center', 'FontSize', 14, 'FontWeight', 'bold');
            end
        end
    end

    xticks(1.5:2:2*length(bands));
    xticklabels(band_labels);
    xlabel('Frequency Band');

    % units of measurement
    ylabel(inferUnitsLabelFromBase(base_name), 'Interpreter','latex');

    title(prettyTitleFromBase(base_name));

    legend({'Responder','Non-responder'}, 'Location', 'best');
    grid on;

    % save
    fname = fullfile(outDir, regexprep("BOX_"+base_name+".png",'[^\w\.-]','_'));
    exportgraphics(fig,fname,'Resolution',300);
    close(fig);
end

%% Helper functions
function lbl = inferUnitsLabelFromBase(base_name)
u = inferUnitsLatexFromBase(lower(string(base_name)));
lbl = ['$[', u, ']$'];
end

function u = inferUnitsLatexFromBase(bn)
if startsWith(bn,"spectrum_log_")
    u = '\log(\mu\mathrm{V}^{2}/\mathrm{Hz})';
elseif startsWith(bn,"spectrum_mag_")
    u = '\mathrm{mag}(\mu\mathrm{V}^{2}/\mathrm{Hz})';
elseif startsWith(bn,"spectrum_")
    u = '\mu\mathrm{V}^{2}/\mathrm{Hz}';
elseif contains(bn,"freq")
    u = '\mathrm{Hz}';
elseif contains(bn,"ratio")
    u = '\mathrm{a.u.}';
else
    u = '\mathrm{a.u.}';
end
end

function ttl = prettyTitleFromBase(base_name)
s = string(base_name);
if startsWith(s,"spectrum_log_median_")
    ttl = 'Median of the band in the log(PSD)';
elseif startsWith(s,"spectrum_mag_median_")
    ttl = 'Median of the band in the mag(PSD)';
elseif startsWith(s,"spectrum_median_")
    ttl = 'Median of the band in the PSD';
elseif startsWith(s,"spectrum_log_maxfreq_")
    ttl = 'Maximum frequency of the band in the log(PSD)';
elseif startsWith(s,"spectrum_mag_maxfreq_")
    ttl = 'Maximum frequency of the band in the mag(PSD)';
elseif startsWith(s,"spectrum_maxfreq_")
    ttl = 'Maximum frequency of the band in the PSD';
elseif startsWith(s,"spectrum_log_max_")
    ttl = 'Maximum spectral power density of the band in the log(PSD)';
elseif startsWith(s,"spectrum_mag_max_")
    ttl = 'Maximum spectral power density of the band in the mag(PSD)';
elseif startsWith(s,"spectrum_max_")
    ttl = 'Maximum spectral power density of the band in the PSD';
elseif startsWith(s,"spectrum_log_")
    ttl = 'Mean of the band in the log(PSD)';
elseif startsWith(s,"spectrum_mag_")
    ttl = 'Mean of the band in the mag(PSD)';
elseif startsWith(s,"spectrum_")
    ttl = 'Mean of the band in the PSD';
else
    s = regexprep(s,'_$',''); s = strrep(s,'_',' ');
    ttl = char(upper(extractBetween(" "+s,2,2))+lower(s));
end
end
%%  WILCOXON + COHEN'S D (ALL FEATURES) 
fprintf('\n========================================\n');
fprintf('Wilcoxon Rank-Sum + Cohen''s d (All Features)\n');

all_significant_features = {};
all_p_values = [];
significant_mw = {};
cohen_d_values = [];

for i = 1:numel(featureNames)
    fname = featureNames{i};
    feat = data.(fname);
    
    % Split by responder status
    feat_resp = feat(outcome == 1);
    feat_nonresp = feat(outcome == 0);

    % Mann–Whitney U test (Wilcoxon Rank-Sum)
    [p_mw, ~] = ranksum(feat_resp, feat_nonresp);

    if p_mw < 0.05
        % Store features with p < 0.05
        all_significant_features{end+1} = fname;
        all_p_values(end+1) = p_mw;

        % Compute Cohen's d
        mean1 = mean(feat_resp, 'omitnan');
        mean2 = mean(feat_nonresp, 'omitnan');
        std1 = std(feat_resp, 'omitnan');
        std2 = std(feat_nonresp, 'omitnan');
        n1 = sum(~isnan(feat_resp));
        n2 = sum(~isnan(feat_nonresp));

        % Pooled standard deviation
        pooled_std = sqrt(((n1 - 1)*std1^2 + (n2 - 1)*std2^2) / (n1 + n2 - 2));
        d = (mean1 - mean2) / pooled_std;

        % Store if effect size is at least small
        if abs(d) > 0.2
            significant_mw{end+1} = fname;
            cohen_d_values(end+1) = d;
        end
    end
end

% Print all features with p < 0.05
fprintf('\n========================================\n');
fprintf('Features with p < 0.05 (Wilcoxon Rank-Sum): %d\n\n', numel(all_significant_features));
fprintf('%-40s\t%-10s\n', 'Feature', 'p-value');
fprintf('----------------------------------------\t-----------\n');
for j = 1:length(all_significant_features)
    fprintf('%-40s\t%10.3e\n', all_significant_features{j}, all_p_values(j));
end

% Print subset with |Cohen's d| > 0.2
fprintf('\n========================================\n');
fprintf('Significant Features (p < 0.05 and |Cohen''s d| > 0.2): %d\n\n', numel(significant_mw));
fprintf('%-40s\t%-10s\n', 'Feature', 'Cohen''s d');
fprintf('----------------------------------------\t-----------\n');
for j = 1:length(significant_mw)
    fprintf('%-40s\t%10.3f\n', significant_mw{j}, cohen_d_values(j));
end

%% TERNARY FUNCTION
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end

function cmap = divergingColormap(n)
% Simple blue-white-red diverging colormap with n entries
    if mod(n,2)==0, n = n+1; end
    m = (n-1)/2;
    up   = (0:1/m:1)';       % 0→1
    down = (1:-1/m:0)';      % 1→0
    blue  = [zeros(m,1) zeros(m,1) up(1:end-1)];   % 0,0,↑
    white = [1 1 1];
    red   = [up(1:end-1) zeros(m,1) zeros(m,1)];   % ↑,0,0
    cmap  = [blue; white; red];
end
