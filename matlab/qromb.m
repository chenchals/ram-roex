function [ a, b, ss, p, r_lin, hz1, hz2 ] = qromb( a, b, ss, p, r_lin, hz1, hz2)
%QROMB Numerical integration by Romberg's method
% Calls: TRAPZD, POLINT
% Called by: THRESH
%-------------------------------------
% Subroutine for performing numerical integration by Romberg's method.
% See Numerical Recipes, p. 114-115.
% This routine has variables p,r_lin,hz1 and hz2 passed through it to trapzd.
%-------------------------------------
% Translated from ROEX3.f90 --> #401-424
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran subroutine:  Return all or few of input arguments
% Since only *ss* occur on LHS, return only that, other input vars are not
% modified (do not occur on LHS)
% 
% see also TRAPZD POLINT THRESH

% No call globalVars as the function *does not* include roex3.h file

  %eps = 1.0d-2 is almost certainly good enough and faster
  eps = 1.0e-5; jmax = 20; jmaxp = jmax + 1; k = 5; km = k - 1;
  s(jmaxp) = 0; h(jmaxp)= 0;
  h(1) = 1.0;
  for j = 1:jmax
      [a, b, s(j), j, p, r_lin, hz1, hz2] = trapzd(a, b, s(j), j, p, r_lin, hz1, hz2);
      if(j>=k)
          % ROEXSUB.f90 subroutine #84-126
          [h(j - km), s(j - km), k, ~, ss, dss] =  polint(h(j - km), s(j - km), k, 0.0, ss, dss);
          if(abs(dss)<eps * abs(ss))
              return
          end
          s(j + 1) = s(j);
          h(j + 1) = 0.25 * h(j);
      end
      fpause(' qromb: Too many steps');
  end
end

