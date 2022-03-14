% Simulate epidemic via renewal model with R and SI profiles
function [Iday, Lam, Rtrue, tday, Iwarn, Pomega] = epiSimTrue(scenNo, distNo, tday, nday, remGap)

% Assumptions and notes
% - option to remove a startup sequence of zeros
% - various R trajectories and serial intervals specified

%% Possible serial interval types available
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};

% Hyerparameters of generation time
switch(distNo)
    case 1
        % Geometric distribution - no parameter
        distvals.pm = [];
    case 2
        % Gamma distribution - shape parameter
        distvals.pm = 20;
    case 3
        % Delta distribution - odd window around mean
        distvals.pm = 7;
    case 4
        % Two Gamma distributions for flare-up
        distvals.pm = 45;
end
% Distribution type and mean
distvals.type = distNo; distvals.omega = 14.2;


% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Warning about incidence zeros
Iwarn = 0; % will be set conditionally

%% Possible true R profiles available
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};

% Functions for scenarios: R on a daily basis
switch(scenNo)
    case 1
        % Constant R
        Rtrue = 1.8*ones(1, nday);
    case 2
        % Sinusoidal R across time
        Rtrue = 1.3 + 0.8*sind(3*tday);
        %Rtrue = 1.5 + 0.5*sind(3*tday);
    case 3
        % Logistic fall in R with time
        R0 = 1.1; Rd = 0.7;
        t1 = 1; t2 = floor(nday/2 + 20);
        Rtrue = R0 + Rd*((1 + exp(-t1*t2))./(1 + exp(-t1*(t2 - tday))));
    case 4
        % Switch point halfway
        tch = ceil(nday/2);
        Rtrue = zeros(1, nday);
        Rtrue(1:tch) = 2.5;
        Rtrue(tch+1:end) = 0.5;
    case 5
        % Exponential rise and fall [0.02 0.008]
        Rtrue = zeros(size(tday)); tchange = floor(nday/4);
        trise = 1:tchange; tfall = tchange+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tchange)); Rmax = Rtrue(tchange);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.015*(tfall - tchange));
    case 6
        % Bottleneck
        Rtrue = zeros(size(tday)); tchange = floor(nday/4);
        Rtrue(1:tchange) = 2;
        Rtrue(tchange+1:2*tchange) = 1.1;
        Rtrue(2*tchange+1:end) = 0.5;
    case 7
        % Two major swings
        Rtrue = zeros(size(tday)); tchange = floor(nday/3);
        Rtrue(1:tchange) = 1;
        Rtrue(tchange+1:2*tchange) = 2;
        Rtrue(2*tchange+1:end) = 0.5;
    case 8
        % White noise with reflection
        Rtrue = abs(normrnd(0.5, 2, [1 nday]));
        % Smooth with m point averager
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
end

% Daily incidence
Iday = zeros(size(nday)); 
% Infectiousness, Poisson rate 
Lam = Iday; rate = Iday;
% Initialise epidemic
Iday(1) = 10;

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = Iday(i-1:-1:1)*Pomegat';
    % Rate for ith day incidence
    rate(i) = Lam(i)*Rtrue(i);
    % Renewal incidence
    Iday(i) = poissrnd(rate(i));
end

% Gaps between non-zero indicence values
zerogaps = diff(find(Iday ~= 0));
% Remove startup sequence of zeros if big
z1 = zerogaps(1);
if z1 > 5 && remGap
    % Update incidence and related vectors
    idz = z1+1:nday;
    % Flag zero incidence regions after startup
    if max(zerogaps(2:end)) > 5
        %warning('Zero incidences beyond startup');
        Iwarn = 1;
    end
else
    % Flag any zero incidence region
    if max(zerogaps) > 5
        %warning('Sequences of zero incidence');
        Iwarn = 1;
    end
    % Un-truncated day set
    idz = 2:nday;
end

% Adjusted vectors - including tday
clear('tday'); tday = idz;
Iday = Iday(idz); 
Rtrue = Rtrue(idz); Lam = Lam(idz);
