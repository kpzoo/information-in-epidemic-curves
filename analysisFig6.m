% Examine FI of Poisson vs NB renewal model
% Reproduces the analysis of Fig 6 (appendix)
clearvars; clc; close all; 

% Assumptions and notes
% - NB has fixed r (shape)

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Define a fixed true R and shapes r
R = 1; r = logspace(-1, 3, 20); lenr = length(r);

% Possible Lam values
L = linspace(0, 100, 1000); len = length(L);

% Poisson FI
Ip = L/R;

% NB FI for r 
In = zeros(lenr, len); 

% Compute for various r
for i = 1:lenr
    In(i, :) = L/R - (L.^2)./(r(i)+ L*R);
end

% Plot results in order of L
figure;
loglog(Ip, In, 'LineWidth', 2, 'Color', grey1);
hold on;
loglog(Ip, In(1, :), 'LineWidth', 2, 'Color', 'r');
loglog(Ip, In(end, :), 'LineWidth', 2, 'Color', 'b');
hold off;
xlabel('Poisson FI');
ylabel('NB FI at various $k$');
grid off; box off;