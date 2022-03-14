% Construct deconvolved estimate of incidence given delayed curve
function Iadj = deconvDelay(Pdelay, tD, D)

% Assumptions and notes
% - Godlstein 2013 method with known infection to event delay (Pdelay)
% - works with death or delayed incidence data D over time tD

% Lengths of time series and delays
lD = length(tD);
if length(D) ~= lD || length(Pdelay) ~= lD
    error('Input dimensions inconsistent');
end

% Find peak of delay PDF (count from 0)
[~, tpk] = max(Pdelay); tpk = tpk - 1;
% CDF of delay (assumed discrete)
Fdelay = cumsum(Pdelay);

% Initial guess at solution
x = [D(tpk+1:end) zeros(1, tpk)];

% Scaled chi-squared and iterations
chival = 0; niter = 0;
% Components of Goldstein formulation
Dn = zeros(1, lD); G = Dn;

% Main EM algorithm for deconvolution
while chival > 1 || niter >= 10000
    
    % Expected cases given delays
    for i = 1:lD
        Dn(i) = sum(x(1:i).*Pdelay(i:-1:1));
    end
    
    % Factor G based on data and expected case ratio
    facD = D./Dn;
    for j = 1:lD
        G(j) = sum(Pdelay(1:(lD-j+1)).*facD(j:lD));
    end
    
    % Richardson-Lucy deconvolution
    x = (x./flip(Fdelay)).*G;
    
    % Threshold for assessing convergence
    chival = mean(((Dn - D).^2)./Dn);
    
    % Update count of iterations
    niter = niter + 1;
    if niter == 10000
        warning('Iterations ended prematurely');
        disp(['Chi2 is ' num2str(chival)]);
    end
end

% Converged solution is adjusted cases
Iadj = round(x);
