% Visualise non-ideal epidemics: under-reporting and delays
% Generates Fig 1 type analysis
clearvars; clc;
close all; tic;

% Assumptions and notes
% - single true curve, no R estimation
% - view M curves under under-reporting
% - view M curves under delayed reporting

% Save data and directories
saveTrue = 0; thisDir = cd;

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(50);

%% Setup single epidemic true simulation

% Choose a scenario and serial interval
scenNo = 4; distNo = 2;

% Initialise epidemic time
tday0 = 1:151; nday0 = length(tday0);
% Number of replicates
M = 50;

% Define possible scenarios for true R and serial interval
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};
scenChoice = scenNam{scenNo}; disp(['True R scenario: ' scenChoice]);

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo}; disp(['Serial interval: ' distChoice]);

% Simulate epidemic scenarios and truncate initial 0s
Iwarn = 1; % ensure no warnings
while Iwarn
    [Iday, Lam, Rtrue, tday, Iwarn, Pomega0] = epiSimTrue(scenNo, distNo, tday0, nday0, 1);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Total number of days and cases
nday = length(tday); totcase = sum(Iday);
% Restrict Pomega
Pomega = Pomega0(1:nday);

% Times of infections of each cases (non-delayed)
tInf = zeros(1, totcase); Icheck = zeros(1, nday);
% Start and end indices for each day
caseStart = Icheck; caseEnd = Icheck;
for i = 1:nday
    % Starting case id for day i
    if i == 1
        caseStart(i) = 1;
    else
        caseStart(i) = caseEnd(i - 1) + 1;
    end
    % Ending case id for day i
    caseEnd(i) = caseStart(i) + Iday(i) - 1;
    
    % All the people infected on day i
    tInf(caseStart(i):caseEnd(i)) = i;
    
    % Check case counts
    Icheck(i) = length(find(tInf == i));
end
% Ensure breakdown is correct by reconstructing Iday
if ~all(Iday == Icheck)
    error('Assignment of case ids incorrect');
end

%% Generate M delayed versions of Iday

% Mean of delay distribution
mtau = 10; 
% NegBin delay parameters
r = 10; p = mtau/(r + mtau);

% Draws from this distribution
delay = nbinrnd(r, 1-p, [M, totcase]);
[~, vtau] = nbinstat(r, 1-p);

% PDF of delay distribution 
xdel = 0:50; Pdel = nbinpdf(xdel, r, 1-p);

% Delayed infection ids and days
tDel = zeros(M, totcase); ndaydel = zeros(1, M);
% Delayed and truncated incidence curves
Idel = cell(1, M); Itrunc = Idel;

% For every replicate add delay on infection ids
for i = 1:M
    % Delayed case ids
    tDel(i, :) = tInf + delay(i, :);
    
    % Construct an epi-curve for each delayed version
    ndaydel(i) = max(tDel(i, :));
    % Ensure maximum not < nday (e.g. due to Iday(end) = 0)
    if ndaydel(i) < nday
        ndaydel(i) = nday;
    end
    
    Itemp = zeros(1, ndaydel(i));
    for j = 1:ndaydel(i)
        % Delayed epi-curve
        Itemp(j) = length(find(tDel(i, :) == j));
    end
    % Check total cases conserved
    if sum(Itemp) ~= totcase
        error('Cases not conserved');
    else
        % Store incidence and truncated incidence
        Idel{i} = Itemp; Itrunc{i} = Itemp(1:nday);
    end
end

% Truncated incidence as a matrix
Itrunc = cell2mat(Itrunc');
% Time when saw first delayed case relative to first case
tcase0 = find(Iday, 1, 'first');
tdelay0 = zeros(1, M);
for i = 1:M
    tdelay0(i) = find(Itrunc(i, :), 1, 'first');
end

%% Generate M under-reported versions of Iday

% Mean of sampling distribution
rho = 0.6; b = 20;
% Parameters of beta distribution
fr = rho/(1 - rho); a = fr*b; 

% Under-reporting distribution
xrep = 0:0.01:1; yrep = betapdf(xrep, a, b);

% Under-reported incidence curves
Isamp = zeros(size(Itrunc));
for i = 1:M
    for j = 1:nday
        % Main downsampling
        Isamp(i, j) = binornd(Iday(j), betarnd(a, b));
    end
end

%% Visualise all noise and their distributions


figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Under-reported curves
subplot(2, 1, 1); hold on;
for i = 1:M
    stairs(tday, Isamp(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(tday, Iday, 'k', 'LineWidth', 2);
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('time $t$', 'FontSize', fnt);
ylabel('under-reported $C_t$', 'FontSize', fnt);
% Distribution for sampling
axes('Position',[0.2 0.78 0.22 0.12]);
plot(xrep, yrep/trapz(yrep), 'b', 'LineWidth', 2);
hold on; ax = gca; ax.YGrid = 'on'; 
plot([rho rho], ax.YLim, 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel(['$\rho_t$ $|$ $\bar{\rho}$ = ' num2str(rho)], 'FontSize', fnt);
title('P($\rho_t$)', 'FontSize', fnt);

% Delayed curves
subplot(2, 1, 2); hold on;
for i = 1:M
    stairs(tday, Itrunc(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(tday, Iday, 'k', 'LineWidth', 2);
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('time $t$', 'FontSize', fnt);
ylabel('delayed $C_t$', 'FontSize', fnt);
% Distribution for delays
axes('Position',[0.2 0.28 0.22 0.12]);
plot(xdel, Pdel/trapz(Pdel), 'b', 'LineWidth', 2);
hold on; ax = gca; ax.YGrid = 'on'; 
plot([mtau mtau], ax.YLim, 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel(['$\delta_{x}$ $|$ $\bar{\delta}$ = ' num2str(mtau)], 'FontSize', fnt);
title('P($\delta_x$)', 'FontSize', fnt);

% Additional axes for R
axes('Position',[0.65 0.28 0.2 0.1]);
plot(tday, Rtrue, 'Color', 'r', 'LineWidth', 2);
hold on; plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2)
grid off; box off; hold off;
xlim([tday(1) tday(end)]);
xlabel('$t$', 'FontSize', fnt);
%title('R numbers $R_t$', 'FontSize', fnt);

% % True incidence and R
% yyaxis right
% plot(tday, Rtrue, 'Color', grey1, 'LineWidth', 1);
% h = gca; h.YColor = h.XColor;
% ylim([min(Rtrue)-0.5, max(Rtrue)+0.5]);
% ylabel('reproduction no. $R_t$', 'FontSize', fnt);
% yyaxis left
% plot(tday, Iday, 'r', 'LineWidth', 2);
% h = gca; h.YColor = h.XColor;
% ylabel('incidence $I_t$', 'FontSize', fnt);
% xlabel('time $t$', 'FontSize', fnt);
% grid off; box off; xlim([tday(1) tday(end)])
% 
% subplot(2, 2, 2);
% % Noise distributions
% plot(xdel/max(xdel), Pdel, 'r', 'LineWidth', 1);
% hold on;
% plot(xrep, yrep/trapz(yrep), 'b', 'LineWidth', 1)
% hold off; grid off; box off;
% xlabel('$x$ (noise)', 'FontSize', fnt);
% 
% subplot(2, 2, 3);
% % Under-reported curves
% plot(tday, Isamp, 'Color', grey1, 'LineWidth', 2);
% hold on;
% plot(tday, Iday, 'r', 'LineWidth', 2);
% hold off; grid off; box off; xlim([tday(1) tday(end)])
% xlabel('time $t$', 'FontSize', fnt);
% ylabel('under-reported cases $C_t$', 'FontSize', fnt);
% 
% subplot(2, 2, 4);
% % Delayed curves
% plot(tday, Itrunc, 'Color', grey1, 'LineWidth', 2);
% hold on;
% plot(tday, Iday, 'r', 'LineWidth', 2);
% hold off; grid off; box off; xlim([tday(1) tday(end)])
% xlabel('time $t$', 'FontSize', fnt);
% ylabel('delayed cases $C_t$', 'FontSize', fnt);
