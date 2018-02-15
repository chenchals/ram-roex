function [ filterwt_out, pg, r_lin ] = filterwt( pg, r_lin )
%FILTERWT Return filter weight
% Calls: 
% Called by: NSCALC
%--------------------------------------
% Equation of the filter weighting function: pg = p * g
%-------------------------------------
% Translated from ROEX3.f90 --> #328-334
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran function:  Return only 1 value
%
% see also NSCALC

  filterwt_out = (1.0 - r_lin) * (1.0 + pg) * exp(-pg) + r_lin;
end

