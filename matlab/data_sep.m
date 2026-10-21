%% PSD-domain separation
% Load dataset
filename = 'DB_TOT_avg.csv'; % Adjust name of the file (e.g. DB_TOT_norm_avg.csv for the normalized features)
[~, prefix] = fileparts(filename);   % output files are named <prefix>_psd.csv, <prefix>_log_psd.csv, ...
data = readtable(filename);

% Define non-feature columns
non_feature_cols = {'ID', 'drug', 'outcome'};

% Get all variable names
var_names = data.Properties.VariableNames;

% Identify raw PSD columns (spectrum_ but not spectrum_log_ or spectrum_mag_)
raw_cols = var_names(startsWith(var_names, 'spectrum_') & ...
                     ~startsWith(var_names, 'spectrum_log_') & ...
                     ~startsWith(var_names, 'spectrum_mag_'));

% Identify log PSD columns
log_cols = var_names(startsWith(var_names, 'spectrum_log_'));

% Identify magnitude PSD columns
mag_cols = var_names(startsWith(var_names, 'spectrum_mag_'));

% Extract each dataset with non-feature columns
raw_table = data(:, [non_feature_cols, raw_cols]);
log_table = data(:, [non_feature_cols, log_cols]);
mag_table = data(:, [non_feature_cols, mag_cols]);

% Write to CSV
writetable(raw_table, [prefix '_psd.csv']);
writetable(log_table, [prefix '_log_psd.csv']);
writetable(mag_table, [prefix '_mag_psd.csv']);

disp([' - ' prefix '_psd.csv']);
disp([' - ' prefix '_log_psd.csv']);
disp([' - ' prefix '_mag_psd.csv']);
%% LEV and SCB 

filename = 'DB_TOT_avg.csv';               
SCB_set  = categorical({'CBZ','ESL','LCM','LTG','PHT'});  % SCB labels
non_feature_cols = {'ID','drug','outcome'};
M = readtable(filename);

if ~iscategorical(M.drug)
    M.drug = categorical(M.drug);
end

% Split full dataset into LEV vs SCB
isLEV = (M.drug == categorical("LEV"));
isSCB = ismember(M.drug, SCB_set);

M_LEV_full = M(isLEV, :);
M_SCB_full = M(isSCB, :);

% Save datasets
writetable(M_LEV_full, 'DB_TOT_avg_LEV.csv');
writetable(M_SCB_full, 'DB_TOT_avg_SCB.csv');


% LOG-tables only
var_names = M.Properties.VariableNames;
log_cols  = var_names(startsWith(var_names, 'spectrum_log_'));

LEV_log = M(isLEV, [non_feature_cols, log_cols]);
SCB_log = M(isSCB, [non_feature_cols, log_cols]);

% Save LOG-only dataset files 
writetable(LEV_log, 'DB_TOT_avg_log_psd_LEV.csv');
writetable(SCB_log, 'DB_TOT_avg_log_psd_SCB.csv');

disp(' - DB_TOT_avg_log_psd_LEV.csv');
disp(' - DB_TOT_avg_log_psd_SCB.csv');

%% creation balanced dataset
% Carica il dataset
T = readtable('DB_TOT_avg_log_psd.csv');

% Estrai gli indici delle due classi
idx0 = find(T.outcome == 0);
idx1 = find(T.outcome == 1);

% Conta quanti zeri ci sono
n0 = numel(idx0);

% Imposta il seed per la riproducibilità
rng(1); % puoi cambiare il numero

% Seleziona casualmente lo stesso numero di 1
sel1 = randsample(idx1, n0);

% Crea il dataset bilanciato
T_balanced = T([idx0; sel1], :);

% Mischia le righe
T_balanced = T_balanced(randperm(height(T_balanced)), :);

% Salva il nuovo dataset
writetable(T_balanced, 'DB_TOT_balanced.csv');
