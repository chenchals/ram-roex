function [ sscalc_out, par ] = sscalc( par )
%SSCALC Returns sum-of-squares between threshold and computed vals using
%parametes in par
% Calls: Uses NSCALC (as function handle), THRESH, BRENT
% Called by: 
%--------------------------------------
% This function returns the total sum of squares difference between the
% threshold data and the calculated values using parameters held in par (as
% determined by the starting values, or, subsequently, by the simplex subroutine).
%-------------------------------------
% Translated from ROEX3.f90 --> #237-288
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran function:  Return one value more than inputs
%
% see also NSCALC, THRESH, BRENT

    % calling globalVars as the function includes roex3.h file
    globalVars

    tol1 = 1.0e-3;
    sscount = sscount + 1;%#ok<NODEF> % global var
    pl = abs(par(1));
    pu = abs(par(2));
    r = -abs(par(3));
    sum_temp = 0.0; % sum is a fx, so use sum_temp
    ns = 0.0; %#ok<*NASGU>
    for notch = 1:npts
        printDebug('notch=         %d ****************************\n',notch);
        if el(notch)==0.0 && eu(notch)== 0.0 || max_shift<0.01
            ns = thresh(0.0);
            shift(notch) = 0.0; %#ok<AGROW> % global var
        else
            % Find the lowest noise-to-signal ratio (ns) and the shift at which it
            % occurs for each notch width.
            neg_el = max(-el(notch),-max_shift);
            %if neg_el < -max_shift, neg_el = -max_shift; end
            pos_eu = min(eu(notch),max_shift);
            %if pos_eu > max_shift,  pos_eu = max_shift; end
            printDebug('sscalc calling brent....\n')
            %[ns, neg_el, pos_eu, tol1, shift(notch)] = brent(neg_el, 0.0, pos_eu, @nscalc, tol1, shift(notch));
            [ ns ] = brent(neg_el, 0.0, pos_eu, @nscalc, tol1, shift(notch));
        end
        printDebug('notch=%d  ns=%1.20E\n',notch,ns);
        calcdb(notch) = 10.0 * log10(ns); %#ok<AGROW> % global var
        gain_var(notch) = 10.0 * log10(thresh(0.0)) - calcdb(notch); %#ok<AGROW,NASGU> % global var
        sum_temp = sum_temp - calcdb(notch) + data(notch);
    end
    xk = sum_temp / double(npts);
    sscalc_out = 0.0;
    for notch = 1:npts
        calcdb(notch) = calcdb(notch) + xk;
        diff_var(notch) = data(notch) - calcdb(notch); %#ok<AGROW> % global var
        sscalc_out = sscalc_out + diff_var(notch)^2;
    end
    % This is a safe way to set limits on the parameters.
    if (pl>pmax || pu>pmax),  sscalc_out = sscalc_out + 1000.0; end
    if (pl<0.1 || pu<0.1), sscalc_out = sscalc_out + 1000.0; end
    if n_print>0
        if mod(sscount, n_print)==0 || sscount==1
            %write(6, 9000)sscount, pl, pu, r, sscalc
            fprintf('%d pl = %7.3f  pu = %7.3f  r = %9.2f  ssq = %9.3f\n', ...
                sscount, pl, pu, r, sscalc_out);
        end
    end
end

