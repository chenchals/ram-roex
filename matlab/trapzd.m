function [ a, b, s, n, p, r_lin, hz1, hz2 ] = trapzd( a, b, s, n, p, r_lin, hz1, hz2 )
%TRAPZD Decrease granularity of spacing of points for numerical integration
% Called by QROMB
%-------------------------------------
% This routine decreases the fineness of spacing of points in the numerical
% integration until further changes make no significant difference.
% See Numerical Recipes, p. 111.
% This version of the routine calls function filterwt(pg,wlin) by name.
% It also calculates the threshold correction to be applied to each
% point before doing the integration.
%-------------------------------------
% Translated from ROEX3.f90 --> #425-468
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran subroutine:  Return all or few of input arguments
% Since only *s* occur on LHS, return only that, other input vars are not
% modified (do not occur on LHS)
%
% see also QROMB GLOBALVARS

    % calling globalVars as the function includes roex3.h file
    globalVars;
    if n==1
        pg = p * a;
        pgb = p * b;
        xhz = hz1;
        delhz = hz2;
        index1 = ((xhz + 5.0d0) / 10.0d0) + 1;
        index2 = ((delhz + 5.0d0) / 10.0d0) + 1;
        if(index1>max_corr_data), index1 = max_corr_data; end
        if(index2>max_corr_data), index2 = max_corr_data; end
        s = 0.5 * (b - a) * (filterwt(pg, r_lin) * lin_corr(index1)...
            + filterwt(pgb, r_lin) * lin_corr(index2));
    else
        it = 2^(n - 2);
        tnm = it;
        del = (b - a) / tnm;
        delhz = (hz2 - hz1) / tnm;
        x = a + 0.5 * del;
        xhz = hz1 + 0.5 * delhz;
        sum = 0;
        for j = 1:it
            index1 = ((xhz + 5.0d0) / 10.0d0) + 1;
            if(index1>max_corr_data), index1 = max_corr_data; end
            pg = p * x;
            sum = sum + filterwt(pg, r_lin) * lin_corr(index1);
            x = x + del;
            xhz = xhz + delhz;
        end
        s = 0.5 * (s + (b - a) * sum / tnm);
    end
end

