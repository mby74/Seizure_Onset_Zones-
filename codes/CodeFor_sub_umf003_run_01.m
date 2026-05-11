%% run_dynamical_geometry_ieeg.m
% Dynamical geometry analysis for invasive EEG / ECoG / SEEG time series
% Author: ChatGPT
%
% INPUT:
%   A CSV or Excel file with columns:
%     time_seconds, channel1, channel2, ...
%
% OUTPUTS:
%   1) tables of geometric/dynamical measures per channel/window
%   2) onset ranking table
%   3) interpretation text file
%   4) supporting figures
%
% MAIN GOALS:
%   - Reconstruct delay-embedded state spaces
%   - Compute geometric/dynamical measures
%   - Detect earliest transition to more ordered dynamics
%   - Rank channels by candidate onset time
%   - Summarize propagation
%
% NOTE:
%   With only ~985 rows at 500 Hz (~2 seconds), results are exploratory.
%   Longer recordings are needed for reliable preictal analysis.

clear; close all; clc;

%% ========================= USER SETTINGS ===============================
inputFile = 'sub-umf003_run-01.csv';   % <-- change to your file
timeVar   = 'time_seconds';

% Optional: if you know seizure onset time in seconds, put it here.
% If unknown, leave empty [] and the code will estimate transition times.
knownSeizureOnset = [];

% Window settings
window_sec = 0.30;      % short windows because data are short
step_sec   = 0.05;

% Embedding settings
maxLagSamples = 30;     % max delay candidate
embedDim      = 3;      % fixed embedding dimension for robustness on short data
theiler       = 5;      % temporal exclusion neighborhood for nearest neighbor search

% Recurrence settings
recurrenceQuantile = 0.10;  % threshold from pairwise distances

% Channel groups (optional, for interpretation)
groupDefs = struct();
groupDefs.Grid = regexptranslate('wildcard','G*');
groupDefs.Cing = regexptranslate('wildcard','CG*');
groupDefs.AH   = regexptranslate('wildcard','AH*');
groupDefs.PH   = regexptranslate('wildcard','PH*');
groupDefs.BO   = regexptranslate('wildcard','BO*');

outputFolder = 'dynamical_geometry_outputs';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ========================= LOAD DATA ==================================
fprintf('Loading data from %s\n', inputFile);

[~,~,ext] = fileparts(inputFile);
switch lower(ext)
    case '.csv'
        T = readtable(inputFile);
    case {'.xlsx','.xls'}
        T = readtable(inputFile);
    otherwise
        error('Unsupported file type. Use CSV or Excel.');
end

assert(any(strcmp(T.Properties.VariableNames, timeVar)), ...
    'Time column "%s" not found.', timeVar);

time = T.(timeVar);
allVars = T.Properties.VariableNames;
chanNames = setdiff(allVars, {timeVar}, 'stable');

Xraw = T{:, chanNames};
Xraw = double(Xraw);

% Remove channels that are constant or all missing
goodCh = true(1, numel(chanNames));
for k = 1:numel(chanNames)
    x = Xraw(:,k);
    if all(isnan(x)) || std(x(~isnan(x))) < eps
        goodCh(k) = false;
    end
end
chanNames = chanNames(goodCh);
Xraw = Xraw(:,goodCh);

% Sampling rate
dt = median(diff(time));
fs = 1/dt;
fprintf('Estimated sampling rate: %.3f Hz\n', fs);

%% ========================= PREPROCESS =================================
% Basic robust preprocessing:
%   - fill missing samples linearly
%   - detrend
%   - z-score each channel
X = zeros(size(Xraw));
for k = 1:size(Xraw,2)
    x = Xraw(:,k);

    if any(isnan(x))
        x = fillmissing(x, 'linear', 'EndValues', 'nearest');
    end

    x = detrend(x);
    x = (x - mean(x)) / max(std(x), eps);
    X(:,k) = x;
end

n = size(X,1);
nCh = size(X,2);

%% ========================= WINDOWING ==================================
winN  = max(round(window_sec * fs), embedDim + 20);
stepN = max(round(step_sec   * fs), 1);

winStarts = 1:stepN:(n - winN + 1);
nWins = numel(winStarts);
winCenters = zeros(nWins,1);

fprintf('Window length: %d samples (%.3f sec)\n', winN, winN/fs);
fprintf('Step length  : %d samples (%.3f sec)\n', stepN, stepN/fs);
fprintf('Number of windows: %d\n', nWins);

%% ========================= FEATURE STORAGE ============================
% Features per window x channel
featureNames = { ...
    'tau', ...
    'radius', ...
    'anisotropy', ...
    'participation_dim', ...
    'corr_dim', ...
    'lyap_proxy', ...
    'recurrence_rate', ...
    'determinism', ...
    'sample_entropy'};

F = struct();
for i = 1:numel(featureNames)
    F.(featureNames{i}) = nan(nWins, nCh);
end

%% ========================= MAIN ANALYSIS ==============================
fprintf('Computing features...\n');

for w = 1:nWins
    idx = winStarts(w):(winStarts(w)+winN-1);
    winCenters(w) = mean(time(idx));

    for c = 1:nCh
        x = X(idx,c);

        % Choose delay tau from first autocorrelation drop below exp(-1)
        tau = chooseDelayACF(x, maxLagSamples);
        tau = max(1, tau);

        % Embedding
        Y = delayEmbed(x, embedDim, tau);
        if size(Y,1) < 20
            continue;
        end

        % 1) Attractor radius
        centroid = mean(Y,1);
        r = sqrt(sum((Y - centroid).^2, 2));
        F.radius(w,c) = mean(r);

        % 2) PCA anisotropy
        C = cov(Y);
        ev = sort(eig(C), 'descend');
        if numel(ev) >= 2 && ev(2) > 0
            F.anisotropy(w,c) = ev(1) / ev(end);
        end

        % 3) Participation-ratio intrinsic dimension
        if all(ev > 0)
            F.participation_dim(w,c) = (sum(ev)^2) / sum(ev.^2);
        end

        % 4) Correlation dimension estimate
        F.corr_dim(w,c) = corrDimensionEstimate(Y);

        % 5) Largest Lyapunov proxy (Rosenstein-style short-term divergence slope)
        F.lyap_proxy(w,c) = lyapunovProxy(Y, fs, tau, theiler);

        % 6) Recurrence rate + determinism
        [RR, DET] = recurrenceMeasures(Y, recurrenceQuantile);
        F.recurrence_rate(w,c) = RR;
        F.determinism(w,c)     = DET;

        % 7) Sample entropy on original window
        F.sample_entropy(w,c) = sampen(x, 2, 0.2*std(x));

        % Store tau
        F.tau(w,c) = tau;
    end
end

fprintf('Feature computation complete.\n');

%% ========================= BUILD LONG TABLE ===========================
fprintf('Building results table...\n');

rows = [];
for w = 1:nWins
    for c = 1:nCh
        newRow = table( ...
            winCenters(w), string(chanNames{c}), ...
            F.tau(w,c), ...
            F.radius(w,c), ...
            F.anisotropy(w,c), ...
            F.participation_dim(w,c), ...
            F.corr_dim(w,c), ...
            F.lyap_proxy(w,c), ...
            F.recurrence_rate(w,c), ...
            F.determinism(w,c), ...
            F.sample_entropy(w,c), ...
            'VariableNames', {'time_center_sec','channel','tau', ...
                              'radius','anisotropy','participation_dim', ...
                              'corr_dim','lyap_proxy','recurrence_rate', ...
                              'determinism','sample_entropy'});
        rows = [rows; newRow]; %#ok<AGROW>
    end
end

writetable(rows, fullfile(outputFolder, 'window_feature_table.csv'));

%% ========================= TRANSITION / ONSET SCORE ===================
% Seizure-like transition, based on movement toward a more ordered state:
%   lower radius, lower participation_dim, lower corr_dim, lower lyap_proxy,
%   higher determinism
%
% We compute a normalized "order score" per channel over windows.
fprintf('Estimating candidate onset times...\n');

orderScore = nan(nWins, nCh);

for c = 1:nCh
    z_radius = robustZ(F.radius(:,c));
    z_pdim   = robustZ(F.participation_dim(:,c));
    z_cdim   = robustZ(F.corr_dim(:,c));
    z_lyap   = robustZ(F.lyap_proxy(:,c));
    z_det    = robustZ(F.determinism(:,c));

    % More ordered = smaller radius/dim/lyap and larger DET
    orderScore(:,c) = (-z_radius) + (-z_pdim) + (-z_cdim) + (-z_lyap) + (+z_det);
end

% Smooth
for c = 1:nCh
    orderScore(:,c) = smoothdata(orderScore(:,c), 'movmean', min(3,nWins));
end

% Candidate onset time for each channel = first large positive jump in order score
onsetTime = nan(nCh,1);
onsetStrength = nan(nCh,1);

for c = 1:nCh
    s = orderScore(:,c);
    ds = [nan; diff(s)];

    thr = nanmean(ds) + 1.0*nanstd(ds);
    idx = find(ds > thr, 1, 'first');

    if ~isempty(idx)
        onsetTime(c) = winCenters(idx);
        onsetStrength(c) = ds(idx);
    else
        [mx, idx2] = max(ds);
        onsetTime(c) = winCenters(idx2);
        onsetStrength(c) = mx;
    end
end

onsetTable = table(string(chanNames(:)), onsetTime, onsetStrength, ...
    'VariableNames', {'channel','candidate_onset_time_sec','transition_strength'});

onsetTable = sortrows(onsetTable, {'candidate_onset_time_sec','transition_strength'}, {'ascend','descend'});
writetable(onsetTable, fullfile(outputFolder, 'candidate_onset_ranking.csv'));

%% ========================= CHANNEL SUMMARY TABLE ======================
summaryTable = table(string(chanNames(:)), ...
    nanmean(F.tau,1)', ...
    nanmean(F.radius,1)', ...
    nanmean(F.anisotropy,1)', ...
    nanmean(F.participation_dim,1)', ...
    nanmean(F.corr_dim,1)', ...
    nanmean(F.lyap_proxy,1)', ...
    nanmean(F.recurrence_rate,1)', ...
    nanmean(F.determinism,1)', ...
    nanmean(F.sample_entropy,1)', ...
    onsetTime, onsetStrength, ...
    'VariableNames', {'channel','mean_tau','mean_radius','mean_anisotropy', ...
                      'mean_participation_dim','mean_corr_dim','mean_lyap_proxy', ...
                      'mean_recurrence_rate','mean_determinism','mean_sample_entropy', ...
                      'candidate_onset_time_sec','transition_strength'});

writetable(summaryTable, fullfile(outputFolder, 'channel_summary_table.csv'));

%% ========================= GROUP SUMMARY ==============================
groupNames = fieldnames(groupDefs);
groupSummary = table();

for g = 1:numel(groupNames)
    gname = groupNames{g};
    mask = false(nCh,1);
    for c = 1:nCh
        mask(c) = ~isempty(regexp(chanNames{c}, groupDefs.(gname), 'once'));
    end
    if any(mask)
        tmp = table(string(gname), ...
            mean(onsetTime(mask), 'omitnan'), ...
            mean(onsetStrength(mask), 'omitnan'), ...
            mean(summaryTable.mean_participation_dim(mask), 'omitnan'), ...
            mean(summaryTable.mean_corr_dim(mask), 'omitnan'), ...
            mean(summaryTable.mean_determinism(mask), 'omitnan'), ...
            'VariableNames', {'group','mean_onset_time_sec','mean_transition_strength', ...
                              'mean_participation_dim','mean_corr_dim','mean_determinism'});
        groupSummary = [groupSummary; tmp]; %#ok<AGROW>
    end
end

if ~isempty(groupSummary)
    writetable(groupSummary, fullfile(outputFolder, 'group_summary_table.csv'));
end

%% ========================= INTERPRETATION =============================
fprintf('Writing interpretation summary...\n');

[~, leadIdx] = sort(onsetTime, 'ascend', 'MissingPlacement', 'last');
topK = min(10, numel(leadIdx));
topChannels = chanNames(leadIdx(1:topK));

interpretationFile = fullfile(outputFolder, 'interpretation.txt');
fid = fopen(interpretationFile, 'w');

fprintf(fid, 'DYNAMICAL GEOMETRY INTERPRETATION\n');
fprintf(fid, '================================\n\n');
fprintf(fid, 'Sampling rate estimate: %.3f Hz\n', fs);
fprintf(fid, 'Recording duration: %.3f sec\n', time(end) - time(1));
fprintf(fid, 'Number of analyzed channels: %d\n', nCh);
fprintf(fid, 'Number of windows: %d\n\n', nWins);

fprintf(fid, 'Main interpretation logic:\n');
fprintf(fid, ['Channels are ranked by the earliest transition toward a more ordered state,\n' ...
              'defined here by a combination of lower geometric spread, lower effective\n' ...
              'dimension, lower divergence (Lyapunov proxy), and higher recurrence determinism.\n\n']);

fprintf(fid, 'Top candidate earliest-transition channels:\n');
for k = 1:topK
    ch = leadIdx(k);
    fprintf(fid, '%2d. %-8s  onset=%.4f sec, strength=%.4f\n', ...
        k, chanNames{ch}, onsetTime(ch), onsetStrength(ch));
end

fprintf(fid, '\nInterpretive guide:\n');
fprintf(fid, ['- Earlier candidate onset time suggests earlier transition into a seizure-like\n' ...
              '  ordered dynamical regime.\n']);
fprintf(fid, ['- Lower participation/correlation dimension suggests less complex, more constrained\n' ...
              '  state-space geometry.\n']);
fprintf(fid, ['- Lower Lyapunov proxy suggests reduced local divergence / greater order.\n']);
fprintf(fid, ['- Higher determinism suggests more repeated trajectory structure.\n']);
fprintf(fid, ['- If the same anatomical group repeatedly leads across seizures, that group is a\n' ...
              '  stronger candidate seizure onset zone.\n\n']);

if isempty(knownSeizureOnset)
    fprintf(fid, ['Known clinical seizure onset time was not provided, so the code estimates\n' ...
                  'transition times from the feature trajectories alone.\n\n']);
else
    fprintf(fid, 'Known seizure onset time used for reference: %.4f sec\n\n', knownSeizureOnset);
end

dur = time(end) - time(1);
if dur < 10
    fprintf(fid, ['WARNING: This recording is very short. Results are exploratory and mostly useful\n' ...
                  'for short-horizon onset ranking and local geometry, not for robust preictal analysis.\n']);
end

fclose(fid);

%% ========================= FIGURES ====================================
fprintf('Generating figures...\n');

% Figure 1: Raw time series for top channels
fig1 = figure('Color','w','Position',[100 100 1200 700]);
nPlot = min(6, nCh);
for i = 1:nPlot
    c = leadIdx(i);
    subplot(nPlot,1,i);
    plot(time, X(:,c), 'k', 'LineWidth', 1);
    ylabel(chanNames{c}, 'Interpreter', 'none');
    if i == 1
        title('Top candidate channels: normalized raw signals');
    end
    if ~isempty(knownSeizureOnset)
        xline(knownSeizureOnset, 'r--', 'LineWidth', 1.5);
    end
    if i < nPlot
        set(gca,'XTickLabel',[]);
    else
        xlabel('Time (s)');
    end
end
saveas(fig1, fullfile(outputFolder, 'figure1_top_channel_timeseries.png'));

% Figure 2: Delay-embedded attractor for the top channel, middle window
fig2 = figure('Color','w','Position',[100 100 1000 700]);
topCh = leadIdx(1);
midW = round(nWins/2);
idx = winStarts(midW):(winStarts(midW)+winN-1);
x = X(idx, topCh);
tau = round(nanmedian(F.tau(:,topCh)));
tau = max(1,tau);
Y = delayEmbed(x, embedDim, tau);

if size(Y,2) >= 3
    plot3(Y(:,1), Y(:,2), Y(:,3), 'k-', 'LineWidth', 1);
    grid on;
    xlabel('x(t)');
    ylabel(sprintf('x(t-%d)', tau));
    zlabel(sprintf('x(t-%d)', 2*tau));
else
    plot(Y(:,1), Y(:,2), 'k-', 'LineWidth', 1);
    grid on;
    xlabel('x(t)');
    ylabel(sprintf('x(t-%d)', tau));
end
title(sprintf('Delay-embedded attractor: %s, window centered at %.3f s', ...
    chanNames{topCh}, winCenters(midW)), 'Interpreter', 'none');
saveas(fig2, fullfile(outputFolder, 'figure2_top_channel_attractor.png'));

% Figure 3: Heatmap of order score
fig3 = figure('Color','w','Position',[100 100 1300 700]);
imagesc(winCenters, 1:nCh, orderScore');
set(gca,'YDir','normal', 'YTick', 1:nCh, 'YTickLabel', chanNames);
xlabel('Time (s)');
ylabel('Channel');
title('Order score heatmap (higher = more ordered seizure-like geometry)');
colorbar;
saveas(fig3, fullfile(outputFolder, 'figure3_order_score_heatmap.png'));

% Figure 4: Feature trajectories for top channel
fig4 = figure('Color','w','Position',[100 100 1200 800]);
subplot(3,2,1);
plot(winCenters, F.participation_dim(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Participation dim');
title(sprintf('%s', chanNames{topCh}), 'Interpreter', 'none');

subplot(3,2,2);
plot(winCenters, F.corr_dim(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Corr dimension');

subplot(3,2,3);
plot(winCenters, F.lyap_proxy(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Lyapunov proxy');

subplot(3,2,4);
plot(winCenters, F.determinism(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Determinism');

subplot(3,2,5);
plot(winCenters, F.radius(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Radius');

subplot(3,2,6);
plot(winCenters, orderScore(:,topCh), 'k-o', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Order score');
saveas(fig4, fullfile(outputFolder, 'figure4_top_channel_feature_trajectories.png'));

% Figure 5: Onset ranking
% fig5 = figure('Color','w','Position',[100 100 1200 700]);
% bar(categorical(onsetTable.channel), onsetTable.candidate_onset_time_sec);
% xlabel('Channel');
% ylabel('Candidate onset time (s)');
% title('Candidate onset ranking by channel');
% xtickangle(90);
% saveas(fig5, fullfile(outputFolder, 'figure5_onset_ranking.png'));

% Sort table (choose ascend or descend as needed)
sortedTable = sortrows(onsetTable, 'candidate_onset_time_sec', 'ascend');

% Create ordered categorical using sorted order
ch = categorical(sortedTable.channel, sortedTable.channel, 'Ordinal', true);

fig5 = figure('Color','w','Position',[100 100 1200 700]);
bar(ch, sortedTable.candidate_onset_time_sec);

xlabel('Channel');
ylabel('Candidate onset time (s)');
title('Candidate onset ranking by channel (sorted)');
xtickangle(90);

saveas(fig5, fullfile(outputFolder, 'figure5_onset_ranking.png'));

% % Figure 6: Group-level onset summary (if groups exist)
% if ~isempty(groupSummary)
%     fig6 = figure('Color','w','Position',[100 100 800 500]);
%     bar(categorical(groupSummary.group), groupSummary.mean_onset_time_sec);
%     xlabel('Group');
%     ylabel('Mean candidate onset time (s)');
%     title('Group-level onset summary');
%     saveas(fig6, fullfile(outputFolder, 'figure6_group_onset_summary.png'));
% end

% Figure 6: Group-level onset summary (sorted smallest → largest)
if ~isempty(groupSummary)

    % Sort ascending (smallest to largest)
    sortedGroup = sortrows(groupSummary, 'mean_onset_time_sec', 'ascend');

    % Preserve order in categorical
    grp = categorical(sortedGroup.group, sortedGroup.group, 'Ordinal', true);

    fig6 = figure('Color','w','Position',[100 100 800 500]);
    bar(grp, sortedGroup.mean_onset_time_sec);

    xlabel('Group');
    ylabel('Mean candidate onset time (s)');
    title('Group-level onset summary (sorted)');

    saveas(fig6, fullfile(outputFolder, 'figure6_group_onset_summary.png'));
end
fprintf('All outputs written to folder: %s\n', outputFolder);

%% ========================= DISPLAY SUMMARY ============================
disp('Top 15 candidate onset channels:');
disp(onsetTable(1:min(15,height(onsetTable)), :));

disp('Channel summary (first 15 rows):');
disp(summaryTable(1:min(15,height(summaryTable)), :));

if ~isempty(groupSummary)
    disp('Group summary:');
    disp(groupSummary);
end

fprintf('\nDone.\n');

%% ========================= LOCAL FUNCTIONS ============================
function tau = chooseDelayACF(x, maxLag)
    x = x(:) - mean(x);
    ac = xcorr(x, maxLag, 'coeff');
    ac = ac(maxLag+1:end); % nonnegative lags
    idx = find(ac < exp(-1), 1, 'first');
    if isempty(idx)
        % fallback: first local minimum
        d = diff(ac);
        idx = find(d(1:end-1) < 0 & d(2:end) > 0, 1, 'first');
        if isempty(idx)
            tau = 1;
        else
            tau = idx;
        end
    else
        tau = idx - 1;
    end
    tau = max(1, tau);
end

function Y = delayEmbed(x, m, tau)
    x = x(:);
    N = numel(x) - (m-1)*tau;
    if N <= 0
        Y = [];
        return;
    end
    Y = zeros(N, m);
    for j = 1:m
        Y(:,j) = x((1:N) + (j-1)*tau);
    end
end

function val = corrDimensionEstimate(Y)
    % crude but practical correlation dimension estimate
    val = nan;
    if size(Y,1) < 20
        return;
    end
    D = pdist(Y);
    D = D(D > 0);
    if numel(D) < 20
        return;
    end

    rmin = prctile(D, 10);
    rmax = prctile(D, 60);
    if rmin <= 0 || rmax <= rmin
        return;
    end

    r = logspace(log10(rmin), log10(rmax), 12);
    C = zeros(size(r));
    for i = 1:numel(r)
        C(i) = mean(D < r(i));
    end

    mask = C > 0 & C < 1;
    if sum(mask) < 4
        return;
    end

    p = polyfit(log(r(mask)), log(C(mask)), 1);
    val = p(1);
end

function lp = lyapunovProxy(Y, fs, tau, theiler)
    % Rosenstein-like short-term divergence slope
    lp = nan;
    N = size(Y,1);
    if N < 30
        return;
    end

    D = squareform(pdist(Y));
    D(1:N+1:end) = inf;

    % Theiler window exclusion
    for i = 1:N
        lo = max(1, i-theiler);
        hi = min(N, i+theiler);
        D(i, lo:hi) = inf;
    end

    [~, nn] = min(D, [], 2);

    maxK = min(10, floor(N/4));
    div = nan(maxK,1);

    for k = 1:maxK
        vals = [];
        for i = 1:N-k
            j = nn(i);
            if j+k <= N
                d0 = norm(Y(i,:)   - Y(j,:));
                dk = norm(Y(i+k,:) - Y(j+k,:));
                if d0 > 0 && dk > 0 && isfinite(d0) && isfinite(dk)
                    vals(end+1,1) = log(dk/d0); %#ok<AGROW>
                end
            end
        end
        if ~isempty(vals)
            div(k) = mean(vals,'omitnan');
        end
    end

    t = (1:maxK)'/fs;
    mask = isfinite(div);
    if sum(mask) >= 4
        p = polyfit(t(mask), div(mask), 1);
        lp = p(1);
    end
end

function [RR, DET] = recurrenceMeasures(Y, q)
    RR = nan; DET = nan;
    N = size(Y,1);
    if N < 20
        return;
    end

    D = squareform(pdist(Y));
    thr = quantile(D(:), q);
    R = D < thr;
    R(1:N+1:end) = 0;

    RR = mean(R(:));

    % Determinism = fraction of recurrence points belonging to diagonals of length >= 2
    minLine = 2;
    totalRec = sum(R(:));
    if totalRec == 0
        DET = 0;
        return;
    end

    diagRec = 0;
    for offset = -(N-1):(N-1)
        d = diag(R, offset);
        if isempty(d), continue; end
        runs = diff([0; d; 0]);
        starts = find(runs == 1);
        ends   = find(runs == -1) - 1;
        lens = ends - starts + 1;
        diagRec = diagRec + sum(lens(lens >= minLine));
    end

    DET = diagRec / totalRec;
end

function z = robustZ(x)
    medx = median(x, 'omitnan');
    madx = mad(x, 1);
    if madx < eps
        z = zeros(size(x));
    else
        z = (x - medx) / (1.4826*madx);
    end
end

function se = sampen(x, m, r)
    % Simple sample entropy implementation
    x = x(:);
    N = numel(x);
    if N <= m+1 || r <= 0
        se = nan;
        return;
    end

    A = 0; B = 0;
    for i = 1:(N-m)
        xi = x(i:i+m-1);
        for j = (i+1):(N-m)
            xj = x(j:j+m-1);
            if max(abs(xi - xj)) < r
                B = B + 1;
                if j <= N-m
                    if max(abs(x(i:i+m) - x(j:j+m))) < r
                        A = A + 1;
                    end
                end
            end
        end
    end

    if A == 0 || B == 0
        se = nan;
    else
        se = -log(A/B);
    end
end