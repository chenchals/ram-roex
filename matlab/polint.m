function [xa, ya, n, x, y, dy] = polint( xa, ya, n, x, y, dy ) %#ok<INUSL>
%POLINT Perform polynomial interpolation
% Called by: QROMB, THR_CORR
%-------------------------------------
% Subroutine polint performs polynomial interpolation.  See Numerical
% Recipes p. 80-83.
%-------------------------------------
% Translated from ROEXSUB.f90 --> #85-126
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
% Fortran subroutine:  Return all or few of input arguments
% Since only *y* and *dy* occur on LHS, return only those, other input vars are not
% modified (do not occur on LHS)
%
%see also QROMB THR_CORR

    nmax = 10;
    c(nmax) = 0;
    d(nmax) = 0;
    ns = 1;
    dif = abs(x - xa(1));
    for i = 1:n
        dift = abs(x - xa(i));
        if(dift<dif)
            ns = i;
            dif = dift;
        end
        c(i) = ya(i);
        d(i) = ya(i);
    end
    y = ya(ns);
    ns = ns - 1;
    for m = 1:n - 1
        for i = 1:n - m
            ho = xa(i) - x;
            hp = xa(i + m) - x;
            w = c(i + 1) - d(i);
            den = ho - hp;
            if(den==0.0), fpause( 'polint : den = 0.0 '); end
            den = w / den;
            d(i) = hp * den;
            c(i) = ho * den;
        end
        if(2 * ns<n - m)
            dy = c(ns + 1);
        else
            dy = d(ns);
            ns = ns - 1;
        end
        y = y + dy;
    end
end

