% Examine when cases outperform death counts
% Generates analysis in Fig 4
clearvars; clc; close all; tic;

% Assumptions and notes
% - reporting fractions via beta distribution
% - comparison on geometric means vs ifr
% - delays ignored and variations on rho considered

% Save data and directories
saveTrue = 0; thisDir = cd; saveFol = 'figures';

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(18);


%%  Generate effective under-reporting case information

% Time scale and IFR for death data
T = 10000; ifr = 0.01;

% Fixed mean reporting fractions
lenr = 18; mrho = linspace(0.07, 0.38, lenr);
% Parameter b (large => more normal like)
lenb = 2000; b = logspace(-1, 2, lenb); 

% Information metric and other beta param
thetaR = zeros(lenr, lenb); a = thetaR;
% Mean and var of beta distributions draw from
mbeta = thetaR; vbeta = mbeta;


% For every mean get metric for each b
for j = 1:lenr
    % Fraction to get a from b
    fr = mrho(j)/(1 - mrho(j));

    % Draw sample probabilities
    rho = zeros(lenb, T);
    for i = 1:lenb
        % Other param of beta
        a(j, i) = fr*b(i);
        % Statistics of distribution
        [mbeta(j, i), vbeta(j, i)] = betastat(a(j, i), b(i));

        % Compute rho (sample rate) across time series
        rho(i, :) = betarnd(a(j, i), b(i), [1 T]);
        % Metric for this sequence
        thetaR(j, i) = geomean(rho(i, :));
    end

end

% Fraction of performance
dtheta = thetaR >= ifr;
perf = sum(dtheta, 2)/lenb;

% Plot variance and effective information
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Mean of beta against theta
subplot(2, 1, 1);
hold on;
for i = 1:lenr
    plot(mbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
plot(mrho, ifr*ones(size(mrho)), 'Color', 'k', "LineWidth", 2);
hold off; grid off; box off;
xlim([mrho(1) mrho(end)]); xlabel('mean reporting fraction $\bar{\rho}$', 'FontSize', fnt);
ax = gca; ax.YGrid = 'on'; ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);
% Var of beta against theta
subplot(2, 1, 2);
hold on;
for i = 1:lenr
    plot(vbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
plot(mrho, ifr*ones(size(mrho)), 'Color', 'k', "LineWidth", 2);
hold off; grid off; box off;
xlim([min(min(vbeta)) max(max(vbeta))]); xlabel('var$(\rho)$', 'FontSize', fnt);
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

%% Publishable plot

% Plot variance and effective information
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% Mean of beta against theta and comparison with IFR
subplot(2, 2, [1 2]); hold on;
for i = 1:lenr
    plot(mbeta(i, :), thetaR(i, :), '.', 'MarkerSize', 15, "Color", cmap(i, :), "LineWidth", 2)
end
plot(mrho, ifr*ones(size(mrho)), 'k--', "LineWidth", 2);
plot(mrho, mrho, '--', "LineWidth", 2, "Color", grey1);
%plot(mrho, perf, '.', "LineWidth", 2, 'MarkerSize', 20, "Color", grey2);
hold off; grid off; box off;
xlim([mrho(1) mrho(end)]); xlabel('mean reporting fraction $\bar{\rho}$', 'FontSize', fnt);
ylabel('$\theta(C_1^{\tau})$ metric', 'FontSize', fnt);
% Add plot of fraction above IFR
axes('Position',[0.19 0.77 0.2 0.12]);
plot(mrho, perf, '.', "LineWidth", 2, 'MarkerSize', 20, "Color", grey2);
hold on; plot(mrho, 0.5*ones(size(mrho)), 'k--', 'LineWidth', 1);
box off; xlim([mrho(1) mrho(end)]);
xlabel('$\bar{\rho}$', 'FontSize', fnt); title('fraction: $\theta(C_1^{\tau}) > $ ifr', 'FontSize', fnt);

% Define a space for plotting distributions
betadom = linspace(0, 1, 200); lenbdom = length(betadom);
idbeta = round(linspace(1, lenb, 10)); lenid = length(idbeta);
% Determine which params beat ifr for given mean
bset = b(idbeta); afirst = a(1, idbeta); alast = a(end, idbeta);
thetafirst = thetaR(1, idbeta) >= ifr; thetalast = thetaR(end, idbeta) >= ifr;

% At smallest and largest rho plot beta distributions
subplot(2, 2, 3); hold on;
for i = 1:lenid
    pltdata = betapdf(betadom, a(1, idbeta(i)), b(idbeta(i)));
    if thetafirst(i)
        plot(betadom, pltdata, "LineWidth", 2, 'Color', 'r');
    else
        plot(betadom, pltdata, "LineWidth", 2, 'Color', 'b');
    end
end
ax = gca; plot([mrho(1) mrho(1)], ax.YLim, 'k--', 'LineWidth', 1);
hold off; grid off; box off;
xlim([0 1]); xlabel('sample fraction $\rho_t$', 'FontSize', fnt);
ax = gca; ax.XGrid = 'on'; ylabel(['P$(\rho_t)$ $|$ $\bar{\rho}$ = ' num2str(round(mrho(1), 3))], 'FontSize', fnt);
subplot(2, 2, 4); hold on;
for i = 1:lenid
    pltdata = betapdf(betadom, a(end, idbeta(i)), b(idbeta(i)));
    if thetalast(i)
        plot(betadom, pltdata, "LineWidth", 2, 'Color', 'r');
    else
        plot(betadom, pltdata, "LineWidth", 2, 'Color', 'b');
    end
end
ax = gca; plot([mrho(end) mrho(end)], ax.YLim, 'k--', 'LineWidth', 1);
hold off; grid off; box off;
xlim([0 1]); xlabel('sample fraction $\rho_t$', 'FontSize', fnt);
ax = gca; ax.XGrid = 'on'; ylabel(['P$(\rho_t)$ $|$ $\bar{\rho}$ = ' num2str(round(mrho(end), 3))], 'FontSize', fnt);



%% Theoretical approx to geom mean of beta

% Define set of a param > 1
atheo = linspace(1.5, 150, 1000); lena = length(atheo);
Gtheo = zeros(lenr, lena); btheo = Gtheo;

% Compute metric under mean constraints
for j = 1:lenr
    Gtheo(j, :) = (atheo - 0.5)./(atheo/mrho(j) - 0.5);
    btheo(j, :) = atheo*(1 - mrho(j))/mrho(j);
end




