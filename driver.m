%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for DMADFO
%
clear
clc
format compact

%%%%%%%%%%%%%%%%%%%%
% problem definition to be provided by the user
%  fun         % function handle
%  cas         % type of initial point
%               * 1: initial point as a matrix
%               * 2: initial point as a vector
%  x           % initial point
%  n,m         % problem size


%%%%%%%%%%%%%%%%%%%%%
% other information to be provided by the user
% st   % structure with start, stop, read and print information
%  .accf    % upper bound for the convergence speed
%              qf = qf(x) = (fs(x)-ftarget)/(f(x0)-ftarget).
%              * default: accf=1e-4
%              * f(x0) is the initial function value
%              * ftarget is the best function value known to us
%              * fs is the function value found by the solver s
%              * the solver s stops once a feasible point x with
%                qf(x)<=1e-4 is found
%              * ftarget=-inf disables the criterion qf(x)<=1e-4
%  .nfmax   % maximum number of function evaluations
%              * default: 5000
%  .secmax  % maximum time in seconds
%              * default: 360
%  .prt     % print level
%              *  -1: nothing
%              *   0: only improved f
%              *   1: little
%              * >=1: more and more


%%%%%%%%%%%%%%%%%%%%%
% tune    % structure with tuning parameters
%         % that are explicitly set; the others take default values
%  .ver       % type of algorithm:
%               0: defaults
%               1: no good step sizes in the mutation and
%                  recombination phases
%               2: no recombination subspace direction
%               3: no nonmonotone line search
%               4: no heuristic point
%  .lambda    % sample size
%  .mu        % number of selected search points in the populations
%  .gammaE    % parameter for expanding step size
%  .E         % maximum iterations of extrapolation step
%  .sigmainit % initial step size
%  .gamma     % parameter for line search condition
%  .sigmamin  % minimum threshold for sigma
%  .sigmamax  % maximum value for sigma
%  .kappa     % factor for the memory of the non-monotone term
%  .N         % memory for the list of function values at best points
%               to choose three vertices of the first triangle
%
% New DMADFO diagonal-scaling parameters:
%  .diagFormula = 2 : quadratic penalty formula, q = 2
%  .diagFormula = 3 : cubic penalty formula, q = 3
%  .diagFormula = 4 : general 2p+z penalty formula, q >= 4
%  .muConj       : conjugacy penalty parameter
%  .rhoConj      : conditioning/reciprocal regularization parameter
%  .pmin         : lower safeguard for diagonal entries
%  .pmax         : upper safeguard for diagonal entries
%  .penaltyP     : p in the general 2p+z penalty, used only if diagFormula=4
%  .penaltyZ     : z in the general 2p+z penalty, used only if diagFormula=4


%%%%%%%%%%%%%%%%%%%%%
%  .isOctave   % 1: run driver in Octave
%              % 0: run driver in Matlab
%
%  .det        % random generator setting
%                *  1: deterministic, for debugging
%                *  0: non-deterministic


%%%%%%%%%%%%%%%%%%%%%%
% problem definition

cas = 2; % type of starting point

switch cas
    case 1
        m = 2; n = 3;   % problem size
        x = [-1 -1 -1; -2 -2 -2];
        xx = [x(1,2); x(2,1)]; % starting point

    case 2
        m = 1; n = 2;   % problem size
        x = [-1 -1];
        xx = x;         % starting point, kept as row vector

    otherwise
        error('Unknown value of cas. Use cas=1 or cas=2.');
end

% function handle
fun = @(xx) (xx(1)-1)^2 + 100*(xx(2)-xx(1)^2)^2; % Rosenbrock function


%%%%%%%%%%%%%%%%%%%%%
% other information to be provided by the user

st.prt     = 2;
st.nfmax   = 10000;
st.secmax  = inf;
st.ftarget = 0;
st.accf    = 1e-4;
st.m       = m;
st.n       = n;
st.det     = 1;


%%%%%%%%%%%%%%%%%%%%%
% tuning parameters

tune.ver = 0;

% Choose the penalty parameters (p,z), where q = 2*p + z.
% The cases (p,z)=(1,0) and (1,1) are handled by the explicit
% quadratic and cubic formulas, respectively. Higher-order cases
% with p >= 2 are handled by the general formula.

p = 20;
z = 0;

if p == 1 && z == 0
    tune.diagFormula = 2;       % quadratic penalty, q = 2
elseif p == 1 && z == 1
    tune.diagFormula = 3;       % cubic penalty, q = 3
elseif p >= 2 && (z == 0 || z == 1)
    tune.diagFormula = 4;       % general higher-order penalty, q >= 4
    tune.penaltyP = p;
    tune.penaltyZ = z;
else
    error('Invalid penalty parameters. Use p >= 1 and z in {0,1}.');
end


%%%%%%%%%%%%%%%%%%%%%
% random generator setting

isOctave = 0;

if ~isOctave
   if st.det
      rng(0,'twister')
   else
      rng('shuffle')
   end
else
   if st.det
      rng_octave('default')
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DMADFO repeatedly calls MADFOstep

[xbest,fbest,info] = DMADFO(fun,xx,st,tune)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
