% Requires the six models exported from the MATLAB Classification Learner app (workspace variables
% below). Exported models that were not trained with cross-validation are re-cross-validated here
% (crossval, 5 folds, new random partition), so fold-wise numbers can differ slightly from the app.
%% ===== Batch CV metrics (mean ± std) for 6 models + confusion matrices =====
% Your exported model variables:
models = {RUS_boosted_lev, RUS_boosted_scb, narrow_nn_lev, narrow_nn_scb, gaussian_nb_lev, gaussian_nb_scb};
modelNames = {"RUS_boosted_lev", "RUS_boosted_scb", "narrow_nn_lev", "narrow_nn_scb", "gaussian_nb_lev", "gaussian_nb_scb"};
POS_CLASS = 1;   % <-- set your positive class label here (e.g., 1 or 'R')

% ---------- Helper to fetch existing CV model (no retraining unless needed) ----------
getCVModel = @(m) localGetPartitionedModel(m, 5);  % uses existing CV if present; else 5-fold

% ---------- Output containers ----------
K_list = zeros(numel(models),1);
summaryRows = [];                 % struct array with mean/std for all metrics
perFoldTables = cell(numel(models),1);
confMats = cell(numel(models),1);

for i = 1:numel(models)
    cvMdl = getCVModel(models{i});          % ClassificationPartitionedModel
    [predLabels, scores] = kfoldPredict(cvMdl);
    trueLabels = cvMdl.Y;
    part = cvMdl.Partition;
    K = part.NumTestSets;
    K_list(i) = K;

    % Determine score column for POS_CLASS (if available)
    posIdx = [];
    if ~isempty(scores)
        classNames = cvMdl.ClassNames;
        if isstring(POS_CLASS) || ischar(POS_CLASS)
            posIdx = find(strcmp(string(classNames), string(POS_CLASS)), 1);
        else
            posIdx = find(ismember(classNames, POS_CLASS), 1);
        end
    end

    % Per-fold metrics
    acc = zeros(K,1); sens = zeros(K,1); spec = zeros(K,1);
    prec = zeros(K,1); f1 = zeros(K,1); auc = NaN(K,1); fpr = zeros(K,1);

    for k = 1:K
        testIdx = test(part, k);
        yt = trueLabels(testIdx);
        yp = predLabels(testIdx);

        % Binarize based on POS_CLASS
        if isstring(POS_CLASS) || ischar(POS_CLASS)
            isPos = (yt == string(POS_CLASS));
            isPosPred = (yp == string(POS_CLASS));
        else
            isPos = (yt == POS_CLASS);
            isPosPred = (yp == POS_CLASS);
        end

        TP = sum(isPosPred & isPos);
        TN = sum(~isPosPred & ~isPos);
        FP = sum(isPosPred & ~isPos);
        FN = sum(~isPosPred & isPos);

        acc(k)  = (TP + TN) / max(1, TP + TN + FP + FN);
        sens(k) = TP / max(1, TP + FN);                 % recall / sensitivity
        spec(k) = TN / max(1, TN + FP);
        prec(k) = TP / max(1, TP + FP);
        if (prec(k) + sens(k)) > 0
            f1(k) = 2 * prec(k) * sens(k) / (prec(k) + sens(k));
        else
            f1(k) = 0;
        end
        fpr(k) = 1 - spec(k);                            % False Positive Rate

        % AUC per fold (if scores exist)
        if ~isempty(scores) && ~isempty(posIdx)
            s = scores(testIdx, posIdx);
            if isstring(POS_CLASS) || ischar(POS_CLASS)
                ybin = yt == string(POS_CLASS);
            else
                ybin = yt == POS_CLASS;
            end
            try
                [~,~,~,auc(k)] = perfcurve(ybin, s, true);
            catch
                auc(k) = NaN;
            end
        end
    end

    % Overall confusion matrix using all CV predictions:
    confMats{i} = confusionmat(trueLabels, predLabels);

    % Store per-fold table
    perFoldTables{i} = table((1:K).', acc, sens, spec, fpr, prec, f1, auc, ...
        'VariableNames', {'Fold','Accuracy','Sensitivity','Specificity','FPR','Precision','F1','AUC'});

    % Summary (mean ± std) in [0,1]
    row = struct( ...
        'Model', string(modelNames{i}), 'KFolds', K, ...
        'Accuracy_mean', mean(acc,'omitnan'),       'Accuracy_std', std(acc,'omitnan'), ...
        'Sensitivity_mean', mean(sens,'omitnan'),   'Sensitivity_std', std(sens,'omitnan'), ...
        'Specificity_mean', mean(spec,'omitnan'),   'Specificity_std', std(spec,'omitnan'), ...
        'FPR_mean', mean(fpr,'omitnan'),            'FPR_std', std(fpr,'omitnan'), ...
        'Precision_mean', mean(prec,'omitnan'),     'Precision_std', std(prec,'omitnan'), ...
        'F1_mean', mean(f1,'omitnan'),              'F1_std', std(f1,'omitnan'), ...
        'AUC_mean', mean(auc,'omitnan'),            'AUC_std', std(auc,'omitnan') );
    summaryRows = [summaryRows; row]; %#ok<AGROW>
end

% ---------- Summary table (one row per model) ----------
summaryTbl = struct2table(summaryRows);
disp(summaryTbl);
% writetable(summaryTbl, 'cv_metrics_mean_std_6models.csv');

% ---------- Also provide a % version (mean/std ×100) ----------
summaryPct = summaryTbl;
vars = summaryTbl.Properties.VariableNames;
scaleVars = vars(contains(vars, '_mean') | contains(vars, '_std'));
for v = 1:numel(scaleVars)
    summaryPct.(scaleVars{v}) = 100 * summaryPct.(scaleVars{v});
end
disp(summaryPct);
% writetable(summaryPct, 'cv_metrics_mean_std_6models_percent.csv');

% Inspect examples:
%   perFoldTables{1}
%   confMats{1}

%% ================= Helper function =================
function cvMdl = localGetPartitionedModel(m, defaultK)
% Returns a ClassificationPartitionedModel.
% - If 'm' is already partitioned: return it.
% - If 'm' is a plain trained model: do crossval(m,'KFold',defaultK).
% - If 'm' is a struct exported by the app, try common fields to find the model,
%   and if not partitioned, wrap it with crossval.
    if isa(m, 'classreg.learning.classif.ClassificationPartitionedModel')
        cvMdl = m; return;
    end
    if isa(m, 'classreg.learning.classif.ClassificationModel')
        cvMdl = crossval(m, 'KFold', defaultK); return;
    end
    if isstruct(m)
        candidates = {'ClassificationSVM','ClassificationEnsemble','ClassificationKNN', ...
                      'ClassificationTree','ClassificationLinear','ClassificationDiscriminant', ...
                      'ClassificationNaiveBayes','ClassificationBaggedEnsemble','ClassificationECOC'};
        for c = 1:numel(candidates)
            if isfield(m, candidates{c})
                inner = m.(candidates{c});
                if isa(inner, 'classreg.learning.classif.ClassificationPartitionedModel')
                    cvMdl = inner; return;
                elseif isa(inner, 'classreg.learning.classif.ClassificationModel')
                    cvMdl = crossval(inner, 'KFold', defaultK); return;
                end
            end
        end
        fns = fieldnames(m);
        for j = 1:numel(fns)
            v = m.(fns{j});
            if isa(v, 'classreg.learning.classif.ClassificationPartitionedModel')
                cvMdl = v; return;
            elseif isa(v, 'classreg.learning.classif.ClassificationModel')
                cvMdl = crossval(v, 'KFold', defaultK); return;
            end
        end
        error('Could not locate a classification model inside the struct.');
    end
    error('Unsupported model type: %s', class(m));
end
