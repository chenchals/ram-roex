function [ j_out ] = locate( xx, n, x, j )
%LOCATE Return index over which interpolation is to be performed
% Called by: THR_CORR
%-------------------------------------
% This subroutine picks the numbers over which the interpolation is to be
% performed.  See Numerical Recipes p. 89-90.
%-------------------------------------
% Translated from ROEXSUB.f90 --> #85-126
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
% Fortran subroutine:  Return all or few of input arguments
% Since only *j* occurs on LHS, return only that, other input vars are not
% modified (do not occur on LHS)
%
%see also THR_CORR
  jl = 0;
  ju = n+1;
  %10    if(ju - jl>1) then
  while  (ju - jl) > 1
      jm = (ju + jl) / 2;
      % eqv is to compare 2 logicals
      % a.eqv.b == true if and only if 
      % 1. both a and b are true or
      % 2. both a and b are false
      %if((xx(n)>xx(1)).eqv.(x>xx(jm)))
      if isequal((xx(n)>xx(1)),(x>xx(jm)))
          jl = jm;
      else
          ju = jm;
      end
      %goto 10
  end
  j_out = jl;
end

