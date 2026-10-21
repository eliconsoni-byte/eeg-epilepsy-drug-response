clear; 
clc; 
close all;

% Load file
dataFile   = 'DB_TOT_avg.csv';
baseOutDir = 'drug_stat_analysis';
figBetween = fullfile(baseOutDir, 'between_groups');  
figLEV     = fullfile(baseOutDir, 'within_LEV');      
figSCB     = fullfile(baseOutDir, 'within_SCB');      
figQuad    = fullfile(baseOutDir, 'quad_boxplots');    
dirs = {baseOutDir, figBetween, figLEV, figSCB, figQuad};
for d = 1:numel(dirs)
    if ~exist(dirs{d},'dir')
        mkdir(dirs{d}); 
    end
end

alpha = 0.05; % Significance level
bigFS = 12;
ttlFS = 14;

% green-red colours
GREEN = [0.20 0.60 0.20];
RED   = [0.80 0.20 0.20];

% Read data
T = readtable(dataFile);
if ~iscategorical(T.drug),    T.drug    = categorical(string(T.drug)); end
if ~iscategorical(T.outcome), T.outcome = categorical(T.outcome); end

% Separate Drug group
SCB_set = categorical({'CBZ','ESL','LCM','LTG','PHT'});
T = T( T.drug=='LEV' | ismember(T.drug,SCB_set), : );
T.drugGroup = categorical( repmat("SCB",height(T),1), {'LEV','SCB'} );
T.drugGroup(T.drug=='LEV') = "LEV";
T.outcome = removecats(T.outcome);  

% Feature columns 
nonFeature = {'drug','outcome','drugGroup','ID','Channel','ChannelLabel','PatientID'};
nonFeature = intersect(nonFeature, T.Properties.VariableNames);
featCols   = setdiff(T.Properties.VariableNames, nonFeature);
isNum = varfun(@isnumeric, T, 'OutputFormat','uniform');
featCols = featCols(isNum(ismember(T.Properties.VariableNames, featCols)));
fprintf('Detected %d feature columns.\n', numel(featCols));

%% Between group analysis: LEV vs SCB
nF = numel(featCols);
p_between = nan(nF,1); d_between = nan(nF,1);

LEVmask = (T.drugGroup=='LEV'); SCBmask = (T.drugGroup=='SCB');

for k=1:nF
    f = featCols{k};
    xLEV = T{LEVmask, f}; xLEV = xLEV(~isnan(xLEV));
    xSCB = T{SCBmask, f}; xSCB = xSCB(~isnan(xSCB));
    if numel(xLEV)>1 && numel(xSCB)>1
        p_between(k) = ranksum(xLEV,xSCB);
        d_between(k) = cohen_d(xLEV,xSCB);
    end
end

Signif_between = repmat("NO",nF,1); Signif_between(p_between<alpha)="YES";
BetweenTable = table(featCols(:), p_between, d_between, Signif_between, ...
    'VariableNames', {'Feature','p_value','Cohens_d','Significant'});
BetweenSig = BetweenTable(BetweenTable.p_value<alpha,:);
if ~isempty(BetweenSig)
    BetweenSig.abs_d = abs(BetweenSig.Cohens_d);
    BetweenSig = sortrows(BetweenSig,'abs_d','descend');
end

writetable(BetweenTable, fullfile(baseOutDir,'Summary_LEV_vs_SCB_Wilcoxon.csv'));
writetable(BetweenSig,  fullfile(baseOutDir,'Summary_LEV_vs_SCB_SIGNIFICANT_sortedByAbsD.csv'));

% feature -> p_between map
pBetweenMap = containers.Map(featCols, num2cell(p_between));

%% Within group analysis
[LEV_All, LEV_Sig, pMap_LEV] = withinDrug_Wilcoxon_tables(T,'LEV',featCols,alpha);
[SCB_All, SCB_Sig, pMap_SCB] = withinDrug_Wilcoxon_tables(T,'SCB',featCols,alpha);

writetable(LEV_All, fullfile(baseOutDir,'Summary_LEV_Responders_vs_NonResponders.csv'));
writetable(LEV_Sig, fullfile(baseOutDir,'Summary_LEV_Responders_vs_NonResponders_SIGNIFICANT_sortedByAbsD.csv'));
writetable(SCB_All, fullfile(baseOutDir,'Summary_SCB_Responders_vs_NonResponders.csv'));
writetable(SCB_Sig, fullfile(baseOutDir,'Summary_SCB_Responders_vs_NonResponders_SIGNIFICANT_sortedByAbsD.csv'));

%% Between group analysis visualization 
featUnion = unique([string(LEV_Sig.Feature); string(SCB_Sig.Feature)]);
fprintf('Quad-plot for %d features (significant in LEV or SCB).\n', numel(featUnion));

for k = 1:numel(featUnion)
    f = featUnion(k);

    levR  = T{ T.drugGroup=='LEV' & T.outcome=='1', f };
    levNR = T{ T.drugGroup=='LEV' & T.outcome=='0', f };
    scbR  = T{ T.drugGroup=='SCB' & T.outcome=='1', f };
    scbNR = T{ T.drugGroup=='SCB' & T.outcome=='0', f };

    x_levR  = ones(size(levR))*1;  x_levNR = ones(size(levNR))*2;
    x_scbR  = ones(size(scbR))*4;  x_scbNR = ones(size(scbNR))*5;

    fig = figure('Color','w','Position',[220 220 980 560]); hold on; box on;
    b1 = boxchart(x_levR,  levR,  'BoxFaceColor',GREEN,'BoxFaceAlpha',0.55,'MarkerStyle','.', 'MarkerSize',5);
    b2 = boxchart(x_levNR, levNR, 'BoxFaceColor',RED,  'BoxFaceAlpha',0.55,'MarkerStyle','.', 'MarkerSize',5);
    b3 = boxchart(x_scbR,  scbR,  'BoxFaceColor',GREEN,'BoxFaceAlpha',0.55,'MarkerStyle','.', 'MarkerSize',5);
    b4 = boxchart(x_scbNR, scbNR, 'BoxFaceColor',RED,  'BoxFaceAlpha',0.55,'MarkerStyle','.', 'MarkerSize',5);
    set([b1 b2 b3 b4],'LineWidth',1.3);

    ax = gca; ax.LineWidth=1.1; ax.FontSize=bigFS;
    ax.XLim=[0 6]; ax.XTick=[1.5 4.5]; ax.XTickLabel={'LEV','SCB'};
    ax.XGrid='on'; ax.YGrid='on';

    Yall = [levR; levNR; scbR; scbNR];
    setNiceYLim(ax, Yall);
    makePlainYTicksSmart(ax);

    
    ax.YTickLabelMode = 'auto';
    expPow = computeNiceExponent(Yall);   
    ax.YAxis.Exponent = expPow;

    
    title(prettyTitleSmart(f), 'FontSize', ttlFS, 'FontWeight','bold');
    ylabel(inferUnitsLabel(f), 'Interpreter','latex');

    % legend
    hN = plot(nan,nan,'s','MarkerSize',10,'MarkerFaceColor',RED,   'MarkerEdgeColor','k');
    hR = plot(nan,nan,'s','MarkerSize',10,'MarkerFaceColor',GREEN, 'MarkerEdgeColor','k');
    hW = plot(nan,nan,'-k','LineWidth',1.6);                    
    hB = plot(nan,nan,'-','Color',[0.8 0 0],'LineWidth',1.6);   
    lgd = legend([hN hR hW hB], ...
        {'Non-responder','Responder','Within-group sig.','Between-group sig.'}, ...
        'Location','eastoutside','Orientation','vertical','Box','on', ...
        'FontSize', max(10, get(gca,'FontSize')-1));
    set(lgd,'AutoUpdate','off');

    % significance lines 
    [y_red, y_black1, y_black2] = sigHeights(ax);
    if pMap_LEV.isKey(char(f)) && pMap_LEV(char(f))<alpha
        draw_sig_line_at(1,2,ax,'k','*', y_black1);
    end
    if pMap_SCB.isKey(char(f)) && pMap_SCB(char(f))<alpha
        draw_sig_line_at(4,5,ax,'k','*', y_black2);
    end
    if pBetweenMap.isKey(char(f)) && pBetweenMap(char(f))<alpha
        draw_sig_line_at(1.5,4.5,ax,[0.8 0 0],'*', y_red);
    end

    fname = fullfile(figQuad, regexprep("QUAD_"+f+"_LEV_SCB.png",'[^\w\.-]','_'));
    exportgraphics(fig,fname,'Resolution',300);
    close(fig);
end

%% Console summary 
fprintf('\nWilcoxon summaries (alpha=%.2f)\n', alpha);
disp('Between Group Analysis: LEV vs SCB — ALL features'); disp(BetweenTable);
disp('Between Group Analysis: LEV vs SCB - Significant'); if isempty(BetweenSig), disp('None'); else, disp(BetweenSig); end
disp('Within LEV — ALL'); disp(LEV_All);  disp('Within LEV — Significant');  if isempty(LEV_Sig), disp('None'); else, disp(LEV_Sig); end
disp('Within SCB — ALL'); disp(SCB_All);  disp('Within SCB — Significant');  if isempty(SCB_Sig), disp('None'); else, disp(SCB_Sig); end

%% Summary figures
families = detectFeatureFamilies(T.Properties.VariableNames);
outDirFull = fullfile(baseOutDir, 'fullpage_band_grids_range');
if ~exist(outDirFull,'dir'), mkdir(outDirFull); end

for i = 1:numel(families)
    fam = families{i};
    fprintf('Full-page (LEV vs SCB) grid per family: %s\n', fam);
    makeFullPageBandGrid_WithSig(T, fam, outDirFull, alpha); 
end


%% Helper Functions 

function d = cohen_d(x1,x2)
x1=x1(:); x2=x2(:);
n1=sum(~isnan(x1)); n2=sum(~isnan(x2));
m1=mean(x1,'omitnan'); m2=mean(x2,'omitnan');
s1=std(x1,'omitnan');  s2=std(x2,'omitnan');
sp=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/max(n1+n2-2,1));
if sp==0 || isnan(sp), d=NaN; else, d=(m1-m2)/sp; end
end

function [AllTable,SigTable,pMap] = withinDrug_Wilcoxon_tables(T,drugName,featCols,alpha)
mask = (T.drugGroup==categorical(string(drugName))); Td=T(mask,:);
nF=numel(featCols); pvals=nan(nF,1); dvals=nan(nF,1);
mR = (Td.outcome=='1'); mN = (Td.outcome=='0');
for k=1:nF
    f=featCols{k};
    xR=Td{mR,f}; xR=xR(~isnan(xR));
    xN=Td{mN,f}; xN=xN(~isnan(xN));
    if numel(xR)>1 && numel(xN)>1
        pvals(k)=ranksum(xR,xN);
        dvals(k)=cohen_d(xR,xN);
    end
end
signif = repmat("NO",nF,1); signif(pvals<alpha)="YES";
AllTable = table(featCols(:), pvals, dvals, signif, ...
    'VariableNames',{'Feature','p_value','Cohens_d','Significant'});
AllTable = sortrows(AllTable,'p_value');
SigTable = AllTable(AllTable.p_value<alpha,:);
if ~isempty(SigTable)
    SigTable.abs_d = abs(SigTable.Cohens_d);
    SigTable = sortrows(SigTable,'abs_d','descend');
end
pMap = containers.Map(featCols, num2cell(pvals));
end

% Titles
function ttl = prettyTitleSmart(raw)
s = string(raw);
bands = ["delta","theta","alpha","beta","gamma","spike"];
greek = containers.Map( ...
    {'delta','theta','alpha','beta','gamma','spike'}, ...
    {'\delta','\theta','\alpha','\beta','\gamma','\it spike'} );

% median (PSD/log/mag)
if startsWith(s,"spectrum_median_")
    b = extract_band(s,bands); ttl = sprintf('Median of the %s band in the PSD', greekLabel(b,greek)); return;
end
if startsWith(s,"spectrum_log_median_")
    b = extract_band(s,bands); ttl = sprintf('Median of the %s band in the log(PSD)', greekLabel(b,greek)); return;
end
if startsWith(s,"spectrum_mag_median_")
    b = extract_band(s,bands); ttl = sprintf('Median of the %s band in the mag(PSD)', greekLabel(b,greek)); return;
end

% maximum spectral power density (PSD/log/mag)
if startsWith(s, "spectrum_max_")
    b = extract_band(s,bands); ttl = sprintf('Maximum spectral power density of the %s band in the PSD', greekLabel(b,greek)); return;
end
if startsWith(s, "spectrum_log_max_")
    b = extract_band(s,bands); ttl = sprintf('Maximum spectral power density of the %s band in the log(PSD)', greekLabel(b,greek)); return;
end
if startsWith(s, "spectrum_mag_max_")
    b = extract_band(s,bands); ttl = sprintf('Maximum spectral power density of the %s band in the mag(PSD)', greekLabel(b,greek)); return;
end

% maximum frequency (PSD/log/mag)
if startsWith(s, "spectrum_maxfreq_")
    b = extract_band(s,bands); ttl = sprintf('Maximum frequency of the %s band in the PSD', greekLabel(b,greek)); return;
end
if startsWith(s, "spectrum_log_maxfreq_")
    b = extract_band(s,bands); ttl = sprintf('Maximum frequency of the %s band in the log(PSD)', greekLabel(b,greek)); return;
end
if startsWith(s, "spectrum_mag_maxfreq_")
    b = extract_band(s,bands); ttl = sprintf('Maximum frequency of the %s band in the mag(PSD)', greekLabel(b,greek)); return;
end

% ratios 
low = lower(s);
if contains(low,"ratio") || contains(low,"_over_") || contains(low,"_to_")
    [b1,b2] = extract_two_bands(low,bands);
    if b1~="" && b2~=""
        if startsWith(s,"spectrum_log_")
            ttl = sprintf('Ratio of %s to %s in the log(PSD)', greekLabel(b1,greek), greekLabel(b2,greek)); return;
        elseif startsWith(s,"spectrum_mag_")
            ttl = sprintf('Ratio of %s to %s in the mag(PSD)', greekLabel(b1,greek), greekLabel(b2,greek)); return;
        else
            ttl = sprintf('Ratio of %s to %s in the PSD', greekLabel(b1,greek), greekLabel(b2,greek)); return;
        end
    else
        ttl = 'Band ratio feature'; return;
    end
end

% band means (PSD/log/mag)
if startsWith(s,"spectrum_log_")
    b = extract_band(s,bands); ttl = sprintf('Mean of the %s band in the log(PSD)', greekLabel(b,greek)); return;
end
if startsWith(s,"spectrum_mag_")
    b = extract_band(s,bands); ttl = sprintf('Mean of the %s band in the mag(PSD)', greekLabel(b,greek)); return;
end
if startsWith(s,"spectrum_")
    b = extract_band(s,bands); ttl = sprintf('Mean of the %s band in the PSD', greekLabel(b,greek)); return;
end

% fallback: clean underscores + capitalize
s = replace(s,"_"," "); ttl = char(upper(extractBetween(" "+s,2,2))+lower(s));
end

function b = extract_band(s,bands)
b = "";
for bi = bands
    if contains(lower(s),"_"+lower(bi)) || endsWith(lower(s),"_"+lower(bi))
        b = bi; return;
    end
end
end

function [b1,b2] = extract_two_bands(s,bands)
b1=""; b2="";
for bi = bands
    idxs = strfind(s, char(bi));
    if ~isempty(idxs)
        if b1=="", b1=bi; else, if b2=="", b2=bi; end, end
    end
end
end

function g = greekLabel(b,greek)
if b=="" || ~isKey(greek,char(b)), g='?'; else, g=greek(char(b)); end
end

% Units (LaTeX) 
function lbl = inferUnitsLabel(fName)
u = inferUnitsLatex(string(fName));
lbl = ['$[', u, ']$'];
end

function u = inferUnitsLatex(fn)
fn = lower(string(fn));
if startsWith(fn,"spectrum_log_")
    u = '\log(\mu\mathrm{V}^{2}/\mathrm{Hz})';
elseif startsWith(fn,"spectrum_mag_")
    u = '\mathrm{mag}(\mu\mathrm{V}^{2}/\mathrm{Hz})';
elseif startsWith(fn,"spectrum_") || contains(fn,"psd")
    u = '\mu\mathrm{V}^{2}/\mathrm{Hz}';
elseif contains(fn,"freq")
    u = '\mathrm{Hz}';
elseif contains(fn,"ratio")
    u = '\mathrm{a.u.}';
else
    u = '\mathrm{a.u.}';
end
end

% Y axis adjustments
function setNiceYLim(ax, Y)
Y = Y(~isnan(Y));
if isempty(Y)
    return; 
end
mn = min(Y); mx = max(Y);
if mn==mx
    epsv = max(1e-12, max(1e-3*abs(mn), eps));
    mn = mn - epsv; mx = mx + epsv;
end
rngv = mx - mn;
padLow  = 0.12*rngv;
padHigh = 0.22*rngv;   
ylim(ax, [mn - padLow, mx + padHigh]);
end

function makePlainYTicksSmart(ax)
ax.YAxis.Exponent = 0;
lim = ylim(ax);
rngv = max(1e-16, lim(2)-lim(1));        
step = niceStep(rngv);                  
lo = floor(lim(1)/step)*step; hi = ceil(lim(2)/step)*step;
yt = lo:step:hi;
if numel(yt) < 3, yt = linspace(lim(1), lim(2), 5); end
set(ax,'YTick', yt);


decimals = max(0, ceil(-log10(max(step, eps))) + 1);
decimals = min(decimals, 12);
labels = compose(['%.',num2str(decimals),'f'], yt);
tries = 0;
while numel(unique(labels)) < numel(labels) && decimals < 12 && tries < 6
    decimals = decimals + 1;
    labels = compose(['%.',num2str(decimals),'f'], yt);
    tries = tries + 1;
end
set(ax,'YTickLabel',labels);
end

function s = niceStep(rngv)

if rngv<=0, s=1; return; end
mag = 10^floor(log10(rngv/5));
candidates = [1 2 2.5 5 10]*mag;
[~,i] = min(abs((rngv/5) - candidates));
s = max(candidates(i), realmin);
end

% Significance line helpers 
function [y_red, y_black1, y_black2] = sigHeights(ax)
yl = ylim(ax); rng = yl(2)-yl(1);
y_red    = yl(2) - 0.06*rng;  % between (red)
y_black1 = yl(2) - 0.12*rng;  % within LEV (black)
y_black2 = yl(2) - 0.18*rng;  % within SCB (black)
end

function draw_sig_line_at(x1,x2,ax,col,starTxt,y)
% Draw line and star
plot(ax,[x1 x2],[y y],'-','Color',col,'LineWidth',1.6, ...
     'Clipping','on','HandleVisibility','off');
text(ax,mean([x1 x2]), y, starTxt, 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',18,'FontWeight','bold', ...
    'Color',col,'Clipping','on','HandleVisibility','off');
end

% Exponent helper for axis
function e = computeNiceExponent(Y)
Y = Y(~isnan(Y) & Y~=0);
if isempty(Y), e = 0; return; end
e = 3 * floor(log10(median(abs(Y))) / 3); 
end

%% FULL-PAGE GRID HELPERS (size adjustments)

function fams = detectFeatureFamilies(varNames)

bands = ["delta","theta","alpha","beta","gamma","spike"];
v = string(varNames);
candidates = strings(0,1);

for k = 1:numel(v)
    s = lower(v(k));
    hit = "";
    for b = bands
        if endsWith(s, "_" + b)
            hit = b; break;
        end
    end
    if hit ~= ""
        fam = extractBefore(s, strlength(s) - strlength(hit) + 1); 
        if ~endsWith(fam, "_"), fam = fam + "_"; end
        candidates(end+1,1) = fam; 
    end
end

whitelistStarts = [ ...
    "spectrum_median_", "spectrum_log_median_", "spectrum_mag_median_", ...
    "spectrum_max_", "spectrum_log_max_", "spectrum_mag_max_", ...
    "spectrum_maxfreq_", "spectrum_log_maxfreq_", "spectrum_mag_maxfreq_", ...
    "spectrum_log_", "spectrum_mag_", "spectrum_" ...
];
mask = false(numel(candidates),1);
for i=1:numel(candidates)
    for j=1:numel(whitelistStarts)
        if startsWith(candidates(i), whitelistStarts(j))
            mask(i)=true; break;
        end
    end
end

fams = cellstr(unique(candidates(mask)));
end
function makeFullPageBandGrid_WithSig(T, featureFamily, outDir, alpha)
    featureFamily = string(featureFamily);
    bands = ["delta","theta","alpha","beta","gamma","spike"];
    greek = {'\delta','\theta','\alpha','\beta','\gamma','\it spike'};

    % Match one column per band
    allFeat = string(T.Properties.VariableNames);
    featCols = strings(1, numel(bands));
    for i = 1:numel(bands)
        pat = featureFamily + bands(i);
        idx = find(allFeat == pat, 1);
        if ~isempty(idx), featCols(i) = allFeat{idx}; 
        end 
    end

    % Rows (LEV / SCB)
    rowDefs(1).name = 'LEV'; rowDefs(1).mask = (T.drugGroup=='LEV');
    rowDefs(2).name = 'SCB'; rowDefs(2).mask = (T.drugGroup=='SCB');

    % Row color palettes
    ROW_COLORS(1).GREEN = [0.25 0.68 0.25];
    ROW_COLORS(1).RED   = [0.85 0.28 0.28];
    ROW_COLORS(2).GREEN = [0.00 0.40 0.15];
    ROW_COLORS(2).RED   = [0.50 0.00 0.10];

   
    applyCap = isPurePSDFamily(featureFamily);  % mean/median/max della PSD (no log/mag)

    % Compute column Y-lims + exponent
    colYLims = nan(numel(bands), 2);
    colExp   = zeros(numel(bands), 1);
    for c = 1:numel(bands)
        f = featCols(c);
        if f==""; continue; end
        Ycol = [];
        for r = 1:numel(rowDefs)
            R = T{ rowDefs(r).mask & (T.outcome=='1'), f };
            N = T{ rowDefs(r).mask & (T.outcome=='0'), f };
            Ycol = [Ycol; R; N];
        end
        Ycol = Ycol(~isnan(Ycol));
        if isempty(Ycol), Ycol = 0; end

        if applyCap
            [lo, hi] = localNiceYLim_withCap(Ycol, 0.99);   
        else
            [lo, hi] = localNiceYLim(Ycol);                  
        end
        colYLims(c,:) = [lo, hi];
        colExp(c) = computeNiceExponent(Ycol);
    end

    % Full-page figure
    fig = figure('Color','w','Position',[60 100 1500 660]);
    tl  = tiledlayout(fig, 2, 6, 'Padding','tight', 'TileSpacing','loose');
    title(tl, composeTitleFromFamily(featureFamily), 'FontSize', 16, 'FontWeight','bold');

    for r = 1:numel(rowDefs)
        for c = 1:numel(bands)
            nexttile; ax = gca; box on; hold on;
            f = featCols(c);
            if f~=""
                R = T{ rowDefs(r).mask & (T.outcome=='1'), f };
                N = T{ rowDefs(r).mask & (T.outcome=='0'), f };

                boxchart(ones(size(R))*1, R, 'BoxFaceColor',ROW_COLORS(r).GREEN,'BoxFaceAlpha',0.60, ...
                         'MarkerStyle','.', 'MarkerSize',5, 'LineWidth',1.2);
                boxchart(ones(size(N))*2, N, 'BoxFaceColor',ROW_COLORS(r).RED,  'BoxFaceAlpha',0.60, ...
                         'MarkerStyle','.', 'MarkerSize',5, 'LineWidth',1.2);

                ax.XLim = [0.5 2.5];
                ax.XTick = [1 2]; ax.XTickLabel = {'R','NR'};
                ax.FontSize = 11; ax.YGrid = 'on';

              
                yl = colYLims(c,:); ylim(ax, yl);
                ax.YTickLabelMode = 'auto';
                ax.YAxis.Exponent = colExp(c);

                if r==1, title(sprintf('%s', greek{c}), 'Interpreter','tex', 'FontWeight','bold'); end
                if c==1
                    ylabel(rowDefs(r).name, 'FontWeight','bold');
                    ylabel(ax, inferUnitsLabel(f), 'Interpreter','latex');
                end

                % Line for significance
                RR = R(~isnan(R)); NN = N(~isnan(N));
                if numel(RR)>1 && numel(NN)>1
                    p = ranksum(RR, NN);
                    if p < alpha
                        y = localSigHeightLower(ax);
                        sigPad = 0.18;
                        xL = 1 + sigPad; xR = 2 - sigPad;
                        plot(ax, [xL xR], [y y], '-k', 'LineWidth', 1.2, 'Clipping','on','HandleVisibility','off');
                        text(ax, 1.5, y, '*', 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                             'FontSize',15,'FontWeight','bold','Color','k','Clipping','on','HandleVisibility','off');
                    end
                end
            else
                axis off;
                if r==1, title(sprintf('%s', greek{c}), 'Interpreter','tex', 'FontWeight','bold'); end
                if c==1, ylabel(rowDefs(r).name, 'FontWeight','bold'); end
            end
        end
    end

    % Legend 
    legAx = axes(fig,'Units','normalized','Position',[0.02 0.95 0.22 0.04]);
    axis(legAx,'off'); hold(legAx,'on');
    y=0.42; h=0.35; w=0.04; pad=0.006; gapGroups=0.30;
    x0=0.00;
    rectangle(legAx,'Position',[x0 y w h],'FaceColor',ROW_COLORS(1).GREEN,'EdgeColor','k');
    rectangle(legAx,'Position',[x0+w+pad y w h],'FaceColor',ROW_COLORS(1).RED,'EdgeColor','k');
    text(legAx,x0+2*w+2*pad+0.01, y+0.5*h,'LEV (row 1)','FontSize',10,'VerticalAlignment','middle');
    x1=x0+2*w+2*pad+gapGroups;
    rectangle(legAx,'Position',[x1 y w h],'FaceColor',ROW_COLORS(2).GREEN,'EdgeColor','k');
    rectangle(legAx,'Position',[x1+w+pad y w h],'FaceColor',ROW_COLORS(2).RED,'EdgeColor','k');
    text(legAx,x1+2*w+2*pad+0.01, y+0.5*h,'SCB (row 2)','FontSize',10,'VerticalAlignment','middle');
    xlim(legAx,[0 1]); ylim(legAx,[0 1]);

    outFigDir = fullfile(outDir, 'by_family_2x6_LEV_SCB');
    if ~exist(outFigDir,'dir'), mkdir(outFigDir); end
    outName = regexprep(['FullPage2x6_LEV_SCB_', char(featureFamily), '.png'], '[^\w\.-]', '_');
    exportgraphics(fig, fullfile(outFigDir, outName), 'Resolution', 300);
    close(fig);
end




function [mn, mx] = localNiceYLim(Y)
    Y = Y(~isnan(Y));
    if isempty(Y), mn = -1; mx = 1; return; end
    mn = min(Y); mx = max(Y);
    if mn==mx
        epsv = max(1e-12, max(1e-3*abs(mn), eps));
        mn = mn - epsv; mx = mx + epsv;
    end
    rngv = mx - mn;
    padLow  = 0.12*rngv;
    padHigh = 0.20*rngv;
    mn = mn - padLow; mx = mx + padHigh;
end

function y = localSigHeightLower(ax)
    yl = ylim(ax); rng = yl(2) - yl(1);
    y = yl(2) - 0.18*rng;  
end

function ttl = composeTitleFromFamily(fam)
fam = string(fam);

if startsWith(fam,"spectrum_median_")
    ttl = 'LEV — SCB | Median of each band in the PSD'; return;
end
if startsWith(fam,"spectrum_log_median_")
    ttl = 'LEV — SCB | Median of each band in the log(PSD)'; return;
end
if startsWith(fam,"spectrum_mag_median_")
    ttl = 'LEV — SCB | Median of each band in the mag(PSD)'; return;
end

if startsWith(fam,"spectrum_max_")
    ttl = 'LEV — SCB | Maximum spectral power density of each band in the PSD'; return;
end
if startsWith(fam,"spectrum_log_max_")
    ttl = 'LEV — SCB | Maximum spectral power density of each band in the log(PSD)'; return;
end
if startsWith(fam,"spectrum_mag_max_")
    ttl = 'LEV — SCB | Maximum spectral power density of each band in the mag(PSD)'; return;
end

if startsWith(fam,"spectrum_maxfreq_")
    ttl = 'LEV — SCB | Maximum frequency of each band in the PSD'; return;
end
if startsWith(fam,"spectrum_log_maxfreq_")
    ttl = 'LEV — SCB | Maximum frequency of each band in the log(PSD)'; return;
end
if startsWith(fam,"spectrum_mag_maxfreq_")
    ttl = 'LEV — SCB | Maximum frequency of each band in the mag(PSD)'; return;
end

if startsWith(fam,"spectrum_log_")
    ttl = 'LEV — SCB | Mean of each band in the log(PSD)'; return;
end
if startsWith(fam,"spectrum_mag_")
    ttl = 'LEV — SCB | Mean of each band in the mag(PSD)'; return;
end
if startsWith(fam,"spectrum_")
    ttl = 'LEV — SCB | Mean of each band in the PSD'; return;
end

ttl = ['LEV — SCB | ', char(fam)];
end


function tf = isPurePSDFamily(fam)
% only for selected figures
fam = string(lower(fam));
tf = (startsWith(fam,"spectrum_max_")    || ...
      startsWith(fam,"spectrum_median_") || ...
     (startsWith(fam,"spectrum_") && ~startsWith(fam,"spectrum_log_") && ~startsWith(fam,"spectrum_maxfreq_") && ~startsWith(fam,"spectrum_mag_")));
end

function [mn, mx] = localNiceYLim_withCap(Y, qTop)
Y = Y(~isnan(Y));
if isempty(Y)
    mn=-1; 
    mx=1; 
    return; 
end
if nargin<2 || isempty(qTop)
    qTop = 0.85; 
end

q1 = quantile(Y, 0.25);
q3 = quantile(Y, 0.75);
iqrY = max(q3 - q1, eps);
mn = min(Y);
mxData = max(Y);
mxCap  = quantile(Y, qTop);
mxUse  = min(mxData, mxCap + 0.10*iqrY);  

% padding
rngv = max(mxUse - mn, eps);
padLow  = 0.12*rngv;
padHigh = 0.20*rngv;
mn = mn - padLow;
mx = mxUse + padHigh;

end
