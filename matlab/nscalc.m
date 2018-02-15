function [ nscalc_out ] = nscalc( ns_shift )
%NSCALC Returns noise-to-signal ratio for shifted filter 
% Calls: GLOBALVARS, FILTERWT, THRESH
% Called by: SSCALC (as function handle), BRENT
%--------------------------------------
% This function calculates the noise-to-signal ratio for a filter shifted by 
% ns_shift from the on-frequency filter.
%-------------------------------------
% Translated from ROEX3.f90 --> #289-327
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran function:  Return only 1 value
%
% see also GLOBALVARS FILTERWT THRESH BRENT

    % calling globalVars as the function includes roex3.h file
    globalVars
    printDebug('in nscalc ns_shift=%1.20E\n',ns_shift)
    r_lin = 10.0^(r / 10.0);

    % The shift passed in is relative to the on-frequency filter but needs to be
    % measured from the shifted (off-frequency) filter.
    % sig_int is the weighting applied by the shifted filter at the signal frequency.

    pconst = p_const(ns_shift);
    pu_shifted = pu * pconst;
    pl_shifted = pl * pconst;
    off_shift = ns_shift / (1.0 + ns_shift);
    if off_shift==0.0
        sig_int = 1.0;
    elseif off_shift<0.0
        push = min(pu_shifted * (-off_shift),maxexp);
        %if push>maxexp, push = maxexp; end
        [sig_int, push, r_lin] = filterwt(push, r_lin);
    elseif off_shift>0.0 
        plsh =  min(pl_shifted * off_shift,maxexp);
        %if plsh>maxexp, plsh = maxexp; end
        [sig_int, plsh, r_lin] = filterwt(plsh, r_lin);
    end
    % Calulate the noise-to-signal ratio at the output of the shifted filter.
    nscalc_out = thresh(ns_shift) / sig_int;
    if nscalc_out<=0.0
        fprintf(' nscalc [%f] is le zero !! \n', nscalc_out)
        nscalc_out = exp(-maxexp);
    end
end

