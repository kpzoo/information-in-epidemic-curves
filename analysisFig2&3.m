% Examine effective information in delays and under-reporting
% Generates analyses in Fig 2 and Fig 3
clearvars; clc; close all; tic;

% Assumptions and notes
% - reporting fractions via beta distribution
% - delay distributions via negative binomial
% - use effective information (theta) to compare

% Save data and directories
saveTrue = 0; thisDir = cd; saveFol = 'figures';

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(18);

%%  Under-reporting modelled with beta distributions

% Set length of time scale
T = 10000; 

% Fixed mean reporting fractions
lenr = 18; mrho = linspace(0.1, 0.9, lenr);
% Parameter b (large => more normal like)
lenb = 2000; b = logspace(-1, 2, lenb); 

% Information metric and other beta param
thetaR = zeros(lenr, lenb); a = thetaR;
% Mean and var of beta distributions draw from
mbeta = thetaR; vbeta = mbeta;

% Truncate beta distribution
pd = makedist('Beta'); pd = truncate(pd, 0.001, 1);

% For every mean get metric for each b
for j = 1:lenr
    % Fraction to get a from b
    fr = mrho(j)/(1 - mrho(j));

    % Draw sample probabilities
    rho = zeros(lenb, T);
    for i = 1:lenb
        % Other param of beta
        a(j, i) = fr*b(i);

        % Compute rho (sample rate) across time series
        rho(i, :) = betarnd(a(j, i), b(i), [1 T]);
        %pd.a = a(j, i); pd.b = b(i); rho(i, :) = random(pd, [1 T]);

        % Statistics of distribution
        [mbeta(j, i), vbeta(j, i)] = betastat(a(j, i), b(i));
        %mbeta(j, i) = mean(pd); vbeta(j, i) = var(pd);

        % Metric for this sequence
        thetaR(j, i) = geomean(rho(i, :));
    end

end

% Plot variance and effective information
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Mean of beta against theta
subplot(2, 1, 1);
hold on;
for i = 1:lenr
    plot(mbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
%plot(mbeta, thetaR, 'Color', grey1, "LineWidth", 0.5)
plot(mrho, mrho, '--', "LineWidth", 2, "Color", grey1);
hold off; grid off; box off;
xlim([mrho(1) mrho(end)]); xlabel('reporting fraction mean $\bar{\rho}$', 'FontSize', fnt);
ax = gca; ax.YGrid = 'on'; ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);
% Var of beta against theta
subplot(2, 1, 2);
hold on;
for i = 1:lenr
    plot(vbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
%plot(vbeta, thetaR, 'Color', grey1, "LineWidth", 0.5)
hold off; grid off; box off;
xlim([min(min(vbeta)) max(max(vbeta))]); xlabel('reporting fraction variance var$(\rho)$', 'FontSize', fnt);
ax = gca; ax.YGrid = 'on'; ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);

% Plot beta distribution pdfs
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Domain of distributions and indices to plot
betadom = linspace(0, 1, 200); 
idbeta = round(linspace(1, lenb, 25)); lenid = length(idbeta);
idrho = round([lenr/4 lenr/2 3*lenr/4]); lenidrho = length(idrho);
for j = 1:lenidrho
    subplot(lenidrho, 1, j);
    hold on;
    for i = 1:lenid
        pltdata = betapdf(betadom, a(idrho(j), idbeta(i)), b(idbeta(i)));
        plot(betadom, pltdata, "LineWidth", 1)
    end
    hold off; grid off; box off;
    xlim([0 1]); xlabel(['$\rho$ $|$ $\bar{\rho}$ = ' num2str(round(mrho(idrho(j)), 3))], 'FontSize', fnt);
    ax = gca; ax.XGrid = 'on'; ylabel('P$(\rho)$', 'FontSize', fnt);
end


%% Under-reporting plots at two values

% Plot variance and effective information
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Mean of beta against theta
subplot(2, 2, [1 2]);
hold on;
for i = 1:lenr
    plot(mbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
%plot(mbeta, thetaR, 'Color', grey1, "LineWidth", 0.5)
plot(mrho, mrho, '--', "LineWidth", 2, "Color", grey1);
hold off; grid off; box off;
xlim([mrho(1) mrho(end)]); xlabel('reporting fraction mean $\bar{\rho}$', 'FontSize', fnt);
ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);
% Add theta against variance for comparison
axes('Position',[0.22 0.79 0.24 0.12]);
hold on;
for i = 1:lenr
    plot(vbeta(i, :), thetaR(i, :), '-', 'MarkerSize', 10, "Color", cmap(i, :), "LineWidth", 2)
end
hold off; grid off; box off;
xlim([min(min(vbeta)) max(max(vbeta))]); xlabel('variance var$(\rho_t)$', 'FontSize', fnt);
ylabel('$\theta(C_1^{\tau})$', 'FontSize', fnt);

% Pick first value for distribution plot 
mrhos = 0.4; betadom = linspace(0, 1, 200); 
bs = logspace(-1, 2, 10); lenbs = length(bs);
% Stats and params of distributions
as = bs*mrhos/(1 - mrhos);
[~, vbetas] = betastat(as, bs);
% Actual plot of pdfs
subplot(2, 2, 3); hold on;
for i = 2:lenbs-1
    plot(betadom, betapdf(betadom, as(i), bs(i)), "Color", grey1, 'LineWidth', 1);
end
plot(betadom, betapdf(betadom, as(1), bs(1)), "Color", 'b', 'LineWidth', 1);
plot(betadom, betapdf(betadom, as(end), bs(end)), "Color", 'r', 'LineWidth', 1);
hold off; grid off; box off;
xlabel(['$\rho$ $|$ $\bar{\rho}$ = ' num2str(mrhos)], 'FontSize', fnt);
ylabel('$\theta(C_1^{\tau})$', 'FontSize', fnt);

% Draw sample probabilities
thetaRs = zeros(1, lenbs);
for i = 1:lenbs
    rhos = betarnd(as(i), bs(i), [1 T]);
    % Metric for this sequence
    thetaRs(i) = geomean(rhos);
end
% Add sample probabilities
axes('Position',[0.33 0.33 0.1 0.08]);
hold on;
plot(vbetas(2:end-1), thetaRs(2:end-1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', grey1);
plot(vbetas(1), thetaRs(1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'b');
plot(vbetas(end), thetaRs(end), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'r');
grid off; box off; hold off; ax = gca; ax.XGrid = 'on';
xlabel('var$(\rho_t)$', 'FontSize', fnt); 
ylabel('P$(\rho)$', 'FontSize', fnt);

% Pick second value for distribution plot 
mrhos = 0.8; betadom = linspace(0, 1, 200); 
bs = logspace(-1, 2, 10); lenbs = length(bs);
% Stats and params of distributions
as = bs*mrhos/(1 - mrhos);
[~, vbetas] = betastat(as, bs);
% Actual plot of pdfs
subplot(2, 2, 4); hold on;
for i = 2:lenbs-1
    plot(betadom, betapdf(betadom, as(i), bs(i)), "Color", grey1, 'LineWidth', 1);
end
plot(betadom, betapdf(betadom, as(1), bs(1)), "Color", 'b', 'LineWidth', 1);
plot(betadom, betapdf(betadom, as(end), bs(end)), "Color", 'r', 'LineWidth', 1);
hold off; grid off; box off;
xlabel(['$\rho$ $|$ $\bar{\rho}$ = ' num2str(mrhos)], 'FontSize', fnt);
ylabel('P$(\rho)$', 'FontSize', fnt);

% Draw sample probabilities
thetaRs = zeros(1, lenbs);
for i = 1:lenbs
    rhos = betarnd(as(i), bs(i), [1 T]);
    % Metric for this sequence
    thetaRs(i) = geomean(rhos);
end
% Add sample probabilities
axes('Position',[0.63 0.33 0.1 0.08]);
hold on;
plot(vbetas(2:end-1), thetaRs(2:end-1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', grey1);
plot(vbetas(1), thetaRs(1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'b');
plot(vbetas(end), thetaRs(end), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'r');
grid off; box off; hold off; ax = gca; ax.XGrid = 'on';
xlabel('var$(\rho_t)$', 'FontSize', fnt); 
title('$\theta(C_1^{\tau})$', 'FontSize', fnt);

%%  Reporting delays modelled with negbin distributions

% Set length of time scale
T = 200; xdom = 0:T-1;

% Fixed mean delays
mdel = round(linspace(5, 25, 18)); lend = length(mdel); 
% Shape param for each mean
nshape = 30; r = logspace(-1, 2, nshape); 

% Information metric and other NB param
thetaD = zeros(lend, nshape); pm1 = thetaD; pm2 = pm1;
% Mean and var of beta distributions draw from
mnegb = thetaD; vnegb = mnegb;

% For each mean examine shape
for j = 1:lend
    % NB scale parameters
    p = mdel(j)./(r + mdel(j));
    % CDF and PDF variable
    cdom = zeros(nshape, T); pdom = cdom;
   
    % CDFs for each NB distribution
    for i = 1:nshape
        % CDF and PDF values
        cdom(i, :) = nbincdf(xdom, r(i), 1-p(i));
        pdom(i, :) = nbinpdf(xdom, r(i), 1-p(i));

        % Effective information and params
        thetaD(j, i) = geomean(cdom(i, :)); 
        pm1(j, i) = r(i); pm2(j, i) = 1 - p(i);

        % Statistics of distribution
        [mnegb(j, i), vnegb(j, i)] = nbinstat(r(i), 1-p(i));
    end
end


% Plot variance and effective information
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Mean of negbin against theta
subplot(2, 1, 1);
hold on;
for i = 1:lend
    plot(mnegb(i, :), thetaD(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
hold off; grid off; box off;
xlim([mdel(1) mdel(end)]); xlabel('$\bar{\delta}$', 'FontSize', fnt);
ax = gca; ax.YGrid = 'on'; ylabel('$\theta(C_1^{\tau})$', 'FontSize', fnt);
% Var of beta against theta
subplot(2, 1, 2);
hold on;
for i = 1:lend
    plot(pm1(i, :), thetaD(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
hold off; grid off; box off;
xlim([min(r) max(r)]); xlabel('NB dispersion $k$', 'FontSize', fnt);
ax = gca; ax.XGrid = 'on'; ax.XScale = 'log';
ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);

% Plot negbin distribution pdfs
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Domain of distributions and indices to plot
nbindom = 0:100; lenbdom = length(nbindom);
idnbin = round(linspace(1, nshape, 20)); lenbin = length(idnbin);
iddel = round([lend/4 lend/2 3*lend/4]); leniddel = length(iddel);
for j = 1:leniddel
    subplot(leniddel, 1, j);
    hold on;
    for i = 1:lenbin
        pltdata = nbincdf(nbindom, pm1(iddel(j), idnbin(i)), pm2(iddel(j), idnbin(i)));
        plot(nbindom, pltdata, "LineWidth", 1);
    end
    hold off; grid off; box off;
    ylim([0 1]); xlabel(['$t$ $|$ $\bar{\delta}$ = ' num2str(round(mdel(iddel(j)), 3))], 'FontSize', fnt);
    ax = gca; ax.XGrid = 'on'; ylabel('F$_{\tau - t}$', 'FontSize', fnt);
end


%% Delay plots at 2 mean values

% Means from likely case and death delays
mdel = 3:3:30; lend = length(mdel); 
% Shape param for each mean
nshape = 100; r = logspace(-1, 2, nshape); 

% Information metric and other NB param
thetaD = zeros(lend, nshape); pm1 = thetaD; pm2 = pm1;
% Mean and var of beta distributions draw from
mnegb = thetaD; vnegb = mnegb;

% For each mean examine shape
for j = 1:lend
    % NB scale parameters
    p = mdel(j)./(r + mdel(j));
    % CDF and PDF variable
    cdom = zeros(nshape, T); pdom = cdom;
   
    % CDFs for each NB distribution
    for i = 1:nshape
        % CDF and PDF values
        cdom(i, :) = nbincdf(xdom, r(i), 1-p(i));
        pdom(i, :) = nbinpdf(xdom, r(i), 1-p(i));

        % Effective information and params
        thetaD(j, i) = geomean(cdom(i, :)); 
        pm1(j, i) = r(i); pm2(j, i) = 1 - p(i);

        % Statistics of distribution
        [mnegb(j, i), vnegb(j, i)] = nbinstat(r(i), 1-p(i));
    end
end

% Domain of distributions and indices to plot
nbindom = 0:100; lenbdom = length(nbindom);
rs = [0.1 0.5 1 10 20]; lenr = length(rs);

% Publishable figure
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Dispersion of negbin against theta
subplot(2, 2, [1 2]);
hold on;
for i = 1:lend
    plot(pm1(i, :), thetaD(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
plot(r, 0.7*ones(size(r)), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlim([min(r) max(r)]); xlabel('NB dispersion $k$', 'FontSize', fnt);
ax = gca; ax.XScale = 'log';
ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);
% Add theta against mean for comparison
axes('Position',[0.25 0.65 0.22 0.12]);
hold on;
for i = 1:lend
    plot(mnegb(i, :), thetaD(i, :), '.', 'MarkerSize', 10, "Color", cmap(i, :), "LineWidth", 2)
end
hold off; grid off; box off;
xlim([mdel(1) mdel(end)]); xlabel('mean delay $\bar{\delta}$', 'FontSize', fnt);
ax = gca; ax.YGrid = 'on'; ylabel('$\theta(C_1^{\tau})$', 'FontSize', fnt);

% Distributions for 5 and 15 day delays
subplot(2, 2, 3);
hold on; mdels = 5; gmeans = zeros(1, lenr);
for i = 1:lenr
    ps = mdels/(rs(i) + mdels);
    cdoms = nbincdf(nbindom, rs(i), 1-ps);
    pdoms = nbinpdf(nbindom, rs(i), 1-ps);
    gmeans(i) = geomean(cdoms);
    if i ~= 1 && i ~= lenr
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', grey1);
    elseif i == 1
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', 'b');
    else
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', 'r');
    end
end
hold off; grid off; box off; ax = gca; ax.XScale = 'log';
xlabel(['time $t$ $|$ $\bar{\delta}$ = ' num2str(mdels)], 'FontSize', fnt);
ax = gca; ax.XGrid = 'on'; ylabel('delay F$_{\tau - t}$', 'FontSize', fnt);
axes('Position',[0.32 0.17 0.09 0.08]);
hold on;
plot(2:lenr-1, gmeans(2:end-1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', grey1);
plot(1, gmeans(1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'b');
plot(lenr, gmeans(end), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'r');
grid off; box off; hold off;
ax = gca; ax.XTick = [1 3 5]; ax.XTickLabel = rs([1 3 5]);
xlabel('$k$', 'FontSize', fnt); ax.XTickLabelRotation = 0;
title('$\theta(C_1^{\tau})$', 'FontSize', fnt);

subplot(2, 2, 4);
hold on; mdels = 15; gmeans = zeros(1, lenr);
for i = 1:lenr
    ps = mdels/(rs(i) + mdels);
    cdoms = nbincdf(nbindom, rs(i), 1-ps);
    pdoms = nbinpdf(nbindom, rs(i), 1-ps);
    gmeans(i) = geomean(cdoms);
    if i ~= 1 && i ~= lenr
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', grey1);
    elseif i == 1
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', 'b');
    else
        stairs(nbindom, cdoms, "LineWidth", 1, 'Color', 'r');
    end
end
hold off; grid off; box off; ax = gca; ax.XScale = 'log';
xlabel(['time $t$ $|$ $\bar{\delta}$ = ' num2str(mdels)], 'FontSize', fnt);
ax = gca; ax.XGrid = 'on'; ylabel('delay F$_{\tau - t}$', 'FontSize', fnt);
axes('Position',[0.8 0.17 0.09 0.08]);
hold on;
plot(2:lenr-1, gmeans(2:end-1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', grey1);
plot(1, gmeans(1), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'b');
plot(lenr, gmeans(end), '.', 'MarkerSize', 20, "LineWidth", 1, 'Color', 'r');
grid off; box off; hold off;
ax = gca; ax.XTick = [1 3 5]; ax.XTickLabel = rs([1 3 5]);
xlabel('$k$', 'FontSize', fnt); ax.XTickLabelRotation = 0;
title('$\theta(C_1^{\tau})$', 'FontSize', fnt);
