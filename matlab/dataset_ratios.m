clear; 
clc;

% Run this script from the data folder: one sub-folder per patient (C*), each containing the
% Brainstorm PSD export (a file with 'multiply' in its name), plus the outcome table
% 'spectrum_log_divdb.csv' (columns: Drugs, Outcome, one row per patient channel).

ID = dir('C*');
Connectivity_Measures = {'wPLI'; 'oEC'; 'PLV'; 'ciPLV'; 'COH'};

DB_TOT = dataset;
DB_TOT_norm = dataset;
FC_PSD_TOT = dataset;
PSD_TOT = dataset;

% All patient Analysis
for p = 1:length(ID)
    cd(ID(p).name);
    folder = dir('C*'); cd(folder(1).name);
    psdensity = dir('*multiply*'); psdensity = load(psdensity.name);
    
    spectrum = squeeze(psdensity.TF);
    spectrum_log = 10*log10(spectrum);
    spectrum_mag = sqrt(spectrum);
    names = {'spectrum'; 'spectrum_log'; 'spectrum_mag'};
    freq_vec = round(psdensity.Freqs);
    
    delta = [2, 4];     theta = [5, 7];     alpha = [8, 12];
    beta = [15, 29];    gamma = [30, 59];   spike = [1, 60];
    Freqs = {'delta';'theta';'alpha';'beta';'gamma';'spike'};
    DB = dataset;
    DB_norm = dataset;
    FC = dataset;
    for c = 1:length(names)
        eval(['temp =' names{c, 1} ';'])
        ID_rep = repmat(ID(p).name, size(psdensity.RowNames,2), 1);
        channel_labels = psdensity.RowNames';
        DB.ID = ID_rep; DB.channel = channel_labels;
        DB_norm.ID = ID_rep; DB_norm.channel = channel_labels;
        FC.ID = ID(p).name;
        
        for f = 1:length(Freqs)
            eval(['[~,first] = ismember(' Freqs{f} '(1), freq_vec);']);
            eval(['[~,last] = ismember(' Freqs{f} '(2), freq_vec);']);
            max_total = max(max(temp));
            max_freqband = max(max(temp(:,first:last)));
            mean_temp = nanmean(temp(:,first:last), 2);
            median_temp = nanmedian(temp(:,first:last), 2);
            [max_temp, b] = nanmax(temp(:,first:last)');
            max_temp = max_temp';
            temp_freq = freq_vec(first:last);
            max_freq = temp_freq(b)';
            
            eval(['FC.' names{c, 1} '{1,1} = temp;']);
            eval(['FC.' names{c, 1} '_max(1,1) = max_total;']);
            eval(['FC.' names{c, 1} '_max_' Freqs{f, 1} '(1,1) = max_freqband;']);
            
            eval(['DB.' names{c, 1} '_mean_' Freqs{f, 1} ' = mean_temp;']);
            eval(['DB.' names{c, 1} '_median_' Freqs{f, 1} ' = median_temp;']);
            eval(['DB.' names{c, 1} '_max_' Freqs{f, 1} ' = max_temp;']);
            eval(['DB.' names{c, 1} '_maxfreq_' Freqs{f, 1} ' = max_freq;']);
            
            eval(['DB_norm.' names{c, 1} '_mean_' Freqs{f, 1} ' = mean_temp / max(mean_temp);']);
            eval(['DB_norm.' names{c, 1} '_median_' Freqs{f, 1} ' = median_temp / max(median_temp);']);
            eval(['DB_norm.' names{c, 1} '_max_' Freqs{f, 1} ' = max_temp / max_total;']);
            eval(['DB_norm.' names{c, 1} '_maxfreq_' Freqs{f, 1} ' = max_freq / ' Freqs{f} '(2);']);
        end
        
        % Alpha-ratios definition
        bands = {'delta', 'theta', 'beta', 'gamma'};
        for ratio_idx = 1:length(bands)
            band = bands{ratio_idx};

            % Raw ratio
            eval(['mean_alpha = DB.' names{c, 1} '_mean_alpha;']);
            eval(['mean_band = DB.' names{c, 1} '_mean_' band ';']);
            alpha_ratio = mean_alpha ./ mean_band;
            eval(['DB.' names{c, 1} '_alpha_' band '_ratio = alpha_ratio;']);

            % Normalized ratio
            eval(['mean_alpha_norm = DB_norm.' names{c, 1} '_mean_alpha;']);
            eval(['mean_band_norm = DB_norm.' names{c, 1} '_mean_' band ';']);
            alpha_ratio_norm = mean_alpha_norm ./ mean_band_norm;
            eval(['DB_norm.' names{c, 1} '_alpha_' band '_ratio = alpha_ratio_norm;']);
        end

        % Alpha peak ratios (alpha max freq / other bands' max freq)
        for ratio_idx = 1:length(bands)
            band = bands{ratio_idx};

            % Raw alpha peak ratio
            eval(['alpha_maxfreq = DB.' names{c, 1} '_maxfreq_alpha;']);
            eval(['band_maxfreq = DB.' names{c, 1} '_maxfreq_' band ';']);
            alpha_peak_ratio = alpha_maxfreq ./ band_maxfreq;
            eval(['DB.' names{c, 1} '_alpha_peak_' band '_ratio = alpha_peak_ratio;']);

            % Normalized alpha peak ratio
            eval(['alpha_maxfreq_norm = DB_norm.' names{c, 1} '_maxfreq_alpha;']);
            eval(['band_maxfreq_norm = DB_norm.' names{c, 1} '_maxfreq_' band ';']);
            alpha_peak_ratio_norm = alpha_maxfreq_norm ./ band_maxfreq_norm;
            eval(['DB_norm.' names{c, 1} '_alpha_peak_' band '_ratio = alpha_peak_ratio_norm;']);
        end
    end
    
    FC_PSD_TOT = [FC_PSD_TOT; FC];
    DB_TOT = [DB_TOT; DB];
    DB_TOT_norm = [DB_TOT_norm; DB_norm];
    
    disp([ID(p).name ' ' num2str(size(DB,1)) ' ' num2str(size(DB_norm,1)) ' ' num2str(size(FC,1))])
    cd ..
    cd ..
end

%% CSV file creation
T = readtable('spectrum_log_divdb.csv');

% Extract the 'drug' and 'outcome' columns
drug_col = T.Drugs;
outcome_col = T.Outcome;

% Add to DB_TOT and DB_TOT_norm and save CSV
DB_TOT.drug = drug_col;
DB_TOT.outcome = outcome_col;

DB_TOT_norm.drug = drug_col;
DB_TOT_norm.outcome = outcome_col;

DB_TOT_table = dataset2table(DB_TOT);  
DB_TOT_table_n = dataset2table(DB_TOT_norm);  

writetable(DB_TOT_table, 'DB_TOT.csv');
writetable(DB_TOT_table_n, 'DB_TOT_norm.csv');

%% Compressed dataset (one row per patient: features averaged across the 19 channels)
% Written for both the non-normalized (DB_TOT -> DB_TOT_avg.csv) and the normalized
% (DB_TOT_norm -> DB_TOT_norm_avg.csv) tables.
for base = {'DB_TOT', 'DB_TOT_norm'}
    T = readtable([base{1} '.csv']);

    % Remove non-numeric columns
    nonFeatureCols = {'channel', 'drug', 'outcome'};

    % List of numeric feature columns
    featureCols = setdiff(T.Properties.VariableNames, [{'ID'}, nonFeatureCols]);
    T_mean = grpstats(T, 'ID', {'mean'}, 'DataVars', featureCols);

    % Extract 'drug' and 'outcome' info from original table
    [~, idx_first] = unique(T.ID, 'stable');
    T_meta = T(idx_first, {'ID', 'drug', 'outcome'});

    T_mean.GroupCount = [];
    % grpstats prefixes each variable with 'mean_'; the prefix is removed (this also shortens
    % names such as spectrum_mean_delta to spectrum_delta, which the later scripts rely on)
    T_mean.Properties.VariableNames = strrep(T_mean.Properties.VariableNames, 'mean_', '');

    % Merge mean features with metadata
    T_patient = outerjoin(T_mean, T_meta, 'Keys', 'ID', 'MergeKeys', true);

    % Save to CSV
    writetable(T_patient, [base{1} '_avg.csv']);
end
