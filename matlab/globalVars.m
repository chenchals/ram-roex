
%GLOBALVARS
% Called by: P_CONST TRAPZD
%-------------------------------------
% The include file roex3.h, used in compilation, is :-
%-------------------------------------
% Translated from ROEX3.H --> #1-23
% Fortran common bloc vars are shared variables, so declare these as global
% Fortran parameters are sort of constants, so declare these as global
% Uninitialized global vars related to dataset are initialized in the
% dataset reading function...
%
% see also P_CONST TRAPZD
% diff, gain, numeric  are matlab internal functions
global el eu data calcdb gain_var diff_var xk width notch npts n_print sscount; %#ok<*NUSED>
global lin_corr numeric_var;
global shift max_shift;
global pl pu r pmax;
global cf cferb pr_erb;
global low_erb_lim up_erb_lim;

%% other vars global?
global nmax ndim;
global max_corr_data;
global c1 c2 maxexp;
% fortran parameter = constant (sort of) 
c1 = 24.673; c2 = 4.368; maxexp = 700.0;
% Maxexp should be set to something near the largest number a negative 
% exponential can have on the system being used; e.g. exp(-700.0)
% is close to 10**(-304) and exp(-182) is close to 10**(-78).
% fortran parameter = constant (sort of) 
ndim = 3; nmax = 100; max_corr_data = 1501;

%%

%Initialize these vars on startup for dataset
% nmax = max number of points, initialize to npts
% el(nmax) = 0; eu(nmax) = 0; data(nmax) = 0; calcdb(nmax) = 0;
% gain_var(nmax) = 0; diff_var(nmax) = 0; shift(nmax) = 0;
% lin_corr(max_corr_data) = 0;


