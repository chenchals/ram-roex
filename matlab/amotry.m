function [ amotry_out, p, y, psum, mp, np, ndim, ihi, fac ] = amotry(p, y, psum, mp, np, ndim, funk, ihi, fac)
%AMOTRY 
% Note: funk is sscalc function handle
% Calls: 
% Called by: SIMPLEX
%--------------------------------------
% no comments in ROEXSUB.f90
%-------------------------------------
% Translated from ROEXSUB.f90 --> #207-231
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
% Fortran function:  Although a function, input args are changed by the
% call.  Since fortran is pass by reference, we will return changed input
% args- *p*, *y*, *psum* in addition to amotry_out
%
% see also SIMPLEX SSCALC

    NMAX = 20;
    fac1 = (1. - fac) / ndim;
    fac2 = fac1 - fac;
    for j = 1:ndim
        ptry(j) = psum(j) * fac1 - p(ihi, j) * fac2;
    end
    ytry = funk(ptry);
    if ytry<y(ihi)
        y(ihi) = ytry;
        for j = 1:ndim
            psum(j) = psum(j) - p(ihi, j) + ptry(j);
            p(ihi, j) = ptry(j);
        end
    end
    amotry_out = ytry;
    return
end

