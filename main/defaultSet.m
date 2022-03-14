% Initialisation settings
function [grey1, grey2, cmap, fnt] = defaultSet(nCol)

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
% Font sizes
fnt = 25; set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt);
% New colour choices
cmap = linspecer(nCol, 'sequential'); grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);