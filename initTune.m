%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% initTune.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function tune = initTune(tune,n)
%
% Initialize tuning parameters.
%
% Input:
%   tune   structure with explicitly set tuning parameters.
%          Missing fields are assigned default values.
%   n      problem dimension.
%
% Alternative call:
%   tune = initTune(n)
%          initializes all parameters with default values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tune = initTune(tune,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allow the call initTune(n)
if nargin == 1
    if isnumeric(tune) && isscalar(tune)
        n = tune;
        tune = struct();
    else
        error('With one input, the input must be dimension n.');
    end

elseif nargin == 2
    if isempty(tune)
        tune = struct();
    end

else
    error('initTune needs either input n or inputs tune,n.');
end

if ~isstruct(tune)
    error('The input tune must be a structure.');
end

if isempty(n) || ~isnumeric(n) || ~isscalar(n) || n <= 0
    error('The dimension n must be a positive scalar.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic MADFO tuning parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(tune,'ver')
    tune.ver = 0;
end

if ~isfield(tune,'lambda')
    tune.lambda = 4 + floor(3*log(n));
end

if ~isfield(tune,'mu')
    tune.mu = floor(tune.lambda/2);
end

if ~isfield(tune,'gammaE')
    tune.gammaE = 4;
end

if ~isfield(tune,'E')
    tune.E = inf;
end

if ~isfield(tune,'sigmainit')
    tune.sigmainit = 1;
end

if ~isfield(tune,'gamma')
    tune.gamma = 1e-12;
end

if ~isfield(tune,'sigmamin')
    tune.sigmamin = 1e-2;
end

if ~isfield(tune,'sigmamax')
    tune.sigmamax = 0.5;
end

if ~isfield(tune,'kappa')
    tune.kappa = 5;
end

if ~isfield(tune,'N')
    tune.N = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal conjugacy-based scaling parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tune.diagFormula = 2 : quadratic penalty formula
% tune.diagFormula = 3 : cubic penalty formula
% tune.diagFormula = 4 : general 2p+z penalty formula

if ~isfield(tune,'diagFormula')
    tune.diagFormula = 4;
end

if ~isfield(tune,'muConj')
    tune.muConj = 1;
end

if ~isfield(tune,'rhoConj')
    tune.rhoConj = 1e-2;
end

if ~isfield(tune,'pmin')
    tune.pmin = 1e-8;
end

if ~isfield(tune,'pmax')
    tune.pmax = 1e8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the general 2p+z penalty formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are used only when tune.diagFormula = 4.
%
% The general penalty exponent is q = 2*p + z.
% Here, p is an integer with p >= 2.
% Also, z = 0 gives an even-order penalty,
% while z = 1 gives an odd-order penalty.
%
% The cases (p,z) = (1,0) and (1,1) are handled by
% tune.diagFormula = 2 and tune.diagFormula = 3.

if ~isfield(tune,'penaltyP')
    tune.penaltyP = 20;
end

if ~isfield(tune,'penaltyZ')
    tune.penaltyZ = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Safeguards and consistency checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tune.lambda < 2 || floor(tune.lambda) ~= tune.lambda
    error('tune.lambda must be an integer with lambda >= 2.');
end

if tune.mu < 1 || tune.mu > tune.lambda || ...
        floor(tune.mu) ~= tune.mu
    error('tune.mu must be an integer with 1 <= mu <= lambda.');
end

if tune.sigmamin <= 0 || tune.sigmamax <= 0 || ...
        tune.sigmamin > tune.sigmamax
    error('Require 0 < sigmamin <= sigmamax.');
end

if tune.sigmainit <= 0
    error('The initial step size sigmainit must be positive.');
end

if tune.sigmainit > tune.sigmamax
    warning('sigmainit exceeds sigmamax; safeguard sigma.');
end

if tune.gamma <= 0
    error('tune.gamma must be positive.');
end

if tune.gammaE <= 0
    error('tune.gammaE must be positive.');
end

if tune.kappa <= 0
    error('tune.kappa must be positive.');
end

if tune.N < 1 || floor(tune.N) ~= tune.N
    error('tune.N must be a positive integer.');
end

if tune.pmin <= 0 || tune.pmax <= 0 || tune.pmin >= tune.pmax
    error('Require 0 < pmin < pmax.');
end

if tune.muConj <= 0
    error('tune.muConj must be positive.');
end

if tune.rhoConj < 0
    error('tune.rhoConj must be nonnegative.');
end

if ~ismember(tune.diagFormula,[2,3,4])
    error('tune.diagFormula must be one of 2, 3, or 4.');
end

% penaltyP and penaltyZ matter only for the general formula.
if tune.diagFormula == 4

    if tune.penaltyP < 2 || ...
            floor(tune.penaltyP) ~= tune.penaltyP
        error('For diagFormula=4, penaltyP must be integer >= 2.');
    end

    if ~ismember(tune.penaltyZ,[0,1])
        error('For diagFormula=4, penaltyZ must be 0 or 1.');
    end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
