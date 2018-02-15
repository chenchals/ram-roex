function [ thresh_out, ns_shift ] = thresh( ns_shift )
%THRESH Return noise passing through filter using numerical integration
% Calls: P_CONST QROMB
% Called by: SSCALC 
%-------------------------------------
% This function returns the amount of noise passing through the filter using
% either numerical integration (if a threshold correction file was specified) 
% or analytic integration (if no correction file was specified).
%-------------------------------------
% Translated from ROEX3.f90 --> #325-400
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran function:  Return only 1 value
%
% see also P_CONST QROMB SSCALC

    % calling globalVars as the function includes roex3.h file
    globalVars
    
    % init vars
    sum_low = 0;
    sum_hih = 0;
    
    r_lin = 10.0^(r / 10.0d0);
    hzl1 = (1.0 - el(notch) - width) * cf;
    if hzl1<0.0, hzl1 = 0.0; end
    hzl2 = (1.0 - el(notch)) * cf;
    hzu1 = (1.0 + eu(notch)) * cf;
    hzu2 = hzu1 + width * cf;
    
    shifted_cf = cf * (1.0 + ns_shift);
    gu = abs(hzu1 - shifted_cf) / shifted_cf;
    cu = abs(hzu2 - shifted_cf) / shifted_cf;
    gl = abs(shifted_cf - hzl2) / shifted_cf;
    cl = abs(shifted_cf - hzl1) / shifted_cf;
    
    pconst = p_const(ns_shift);
    pu_shifted = pu * pconst;
    pl_shifted = pl * pconst;
    
    if numeric_var
        %  Integrate over equally spaced 'g' values.
        [gl, cl, sum_low, pl_shifted, r_lin, hzl2, hzl1] = qromb(gl, cl, sum_low, pl_shifted, r_lin, hzl2, hzl1);
        [gu, cu, sum_hih, pu_shifted, r_lin, hzu1, hzu2] = qromb(gu, cu, sum_hih, pu_shifted, r_lin, hzu1, hzu2);
        thresh_out = sum_hih + sum_low;
        %correct for the cf of the shifted filter being used
        thresh_out = thresh_out * shifted_cf / cf;
    else
        %  Analytic integration.
        pgush = min(pu_shifted * gu,maxexp);
        pglsh = min(pl_shifted * gl,maxexp);
        %if pgush>maxexp, pgush = maxexp; end
        %if pglsh>maxexp, pglsh = maxexp; end
        pucu = min(pu_shifted * cu,maxexp);
        plcl = min(pl_shifted * cl,maxexp);
        %if plcl>maxexp, plcl = maxexp; end
        %if pucu>maxexp, pucu = maxexp; end
        thrint_l = (r_lin - 1.0d0) * (2.0d0 + plcl) * exp(-plcl) / pl_shifted + r_lin * cl;
        thresh_l = (r_lin - 1.0d0) * (2.0d0 + pglsh) * exp(-pglsh) / pl_shifted + r_lin * gl;
        thrint_l = thrint_l - thresh_l;
        thrint_u = (r_lin - 1.0d0) * (2.0d0 + pucu) * exp(-pucu) / pu_shifted + r_lin * cu;
        thresh_u = (r_lin - 1.0d0) * (2.0d0 + pgush) * exp(-pgush) / pu_shifted + r_lin * gu;
        thrint_u = thrint_u - thresh_u;
        thresh_out = thrint_l + thrint_u;
        % correct for the cf of the shifted filter being used
        thresh_out = thresh_out * shifted_cf / cf;
    end
    if thresh_out<=0.0
        fprintf(' thresh_out [%f] le 0.0\n',thresh_out)
        thresh_out = exp(-maxexp);
    end
    printDebug('   thresh called with ns_shift=%1.20E  thresh=%1.20E\n',ns_shift,thresh_out);
end

