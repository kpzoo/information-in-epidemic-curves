% Examine when cases outperform death counts
% Generates analysis in Fig 5 (new)
clearvars; clc; close all; tic;

% Assumptions and notes
% - uses estimates from the COVID-19 literature
% - repeat with Ebola virus estimates
% - comparison on geometric means 

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Save data and directories
saveTrue = 0; thisDir = cd; saveFol = 'figures';

% Directory of some main code and plotting options
cd('main'); mainDir = cd;
cd(thisDir); addpath(mainDir);

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(18);


%%  Compute geometric mean based metrics for COVID-19

% Death delay distribution (Irons 2021)
mdel = 21; r = 1/(1+1.1); p = mdel/(r + mdel);
% Time scale and weekly bins
T = 42*4; ids = (7:7:T) - 1;
% Infection to death delay CDF at bins and mean
H = nbincdf(ids, r, 1-p); Geo_H = geomean(H);

% Reporting delay distribution (Huisman 2020)
F = geocdf(ids, 1/(1+10.8)); Geo_F = F(1);

% Under-reporting of deaths (CDC 2021)
sigbnds = [1/1.34, 1/1.29]; M = 10000;
% IFR (Meyerowitz-Katz 2020) 
ifrbnds = [0.53, 0.82]/100; 


% Sample from M trajectories of size T
Geo_sigma = zeros(1, M); Geo_ifr = Geo_sigma; 
for i = 1:M
    % Samples of reporting death rates
    psigma = sigbnds(1) + diff(sigbnds)*rand(1, T);
    Geo_sigma(i) = geomean(psigma);
    % Uncertainty on ifr
    pifr = ifrbnds(1) + diff(ifrbnds)*rand(1, T);
    Geo_ifr(i) = geomean(pifr);
end

% Under-reporting of cases 
Geo_rho = zeros(1, M); Geo_rhoErr = Geo_rho;
% Minimum case reporting fractions (Pullano 2021) 
rhobnds = [0.06, 0.08]; %rhobnds = [0.07, 0.38]; 
% Fractions from CDC or Russell 2020
%rhobnds = [1/4.7, 1/3.4]; rhobnds = [1/17.6, 1/12.4];

for i = 1:M
    % Samples of reporting case rates
    prho = rhobnds(1) + diff(rhobnds)*rand(1, T);
    Geo_rho(i) = geomean(prho);
end

% Derive ordering on cases versus deaths
caseInfoCOVID = Geo_rho*Geo_F;
deathInfoCOVID = Geo_sigma.*Geo_ifr*Geo_H;


%%  Compute geometric mean based metrics for Ebola virus

% Death delay distribution (11.4 to onset + 10 to death) 
mdel = 21.4; r = 1.5; p = mdel/(r + mdel);
% Time scale and weekly bins
T = 42*4; ids = (7:7:T) - 1;
% Infection to death delay CDF at bins
H = nbincdf(ids, r, 1-p);
% Geometric mean from delay
Geo_H = geomean(H);

% Reporting delay distribution (maximised) 
Geo_F = 1; 

% Under-reporting of deaths (maximised)
sigbnds = [1, 1]; M = 10000;
% IFR which is also a reporting rate (WHO 2014)
ifrbnds = [0.69, 0.73];
% Sample from M trajectories of size T
Geo_sigma = zeros(1, M); Geo_ifr = Geo_sigma;
for i = 1:M
    % Samples of reporting death rates
    psigma = sigbnds(1) + diff(sigbnds)*rand(1, T);
    Geo_sigma(i) = geomean(psigma);
    % Uncertainty on ifr
    pifr = ifrbnds(1) + diff(ifrbnds)*rand(1, T);
    Geo_ifr(i) = geomean(pifr);
end

% Under-reporting of cases (Dalziel 2018)
rhobnds = [0.33, 0.83]; Geo_rho = zeros(1, M); 
for i = 1:M
    % Samples of reporting case rates
    prho = rhobnds(1) + diff(rhobnds)*rand(1, T);
    Geo_rho(i) = geomean(prho);
end

% Derive ordering on cases versus deaths
caseInfoEVD = Geo_rho*Geo_F;
deathInfoEVD = Geo_sigma.*Geo_ifr*Geo_H;
% Delays vs reporting
Geo_RHS = Geo_H/Geo_F;
Geo_LHS = Geo_rho./(Geo_ifr.*Geo_sigma);


%% Publishable Figure

% Combined estimates on log scale
figure;
subplot(2, 1, 1);
h = histogram(log(caseInfoCOVID./deathInfoCOVID), 'Normalization', 'probability');
h.NumBins = 25; h.FaceAlpha = 0.2; h.EdgeAlpha = 0.1;
box off; grid off;
ylabel('COVID-19');
subplot(2, 1, 2);
h = histogram(log(caseInfoEVD./deathInfoEVD), 'Normalization', 'probability');
h.NumBins = 25; h.FaceAlpha = 0.2; h.EdgeAlpha = 0.1;
box off; grid off;
ylabel('EVD');
xlabel('$\log \theta(C_1^\tau)-\log \theta(D_1^\tau)$');


