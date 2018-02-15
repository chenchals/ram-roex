%ROEX3
% Usage: roex3
%
% data file name? [dir-path]/TEST.DAT
% file name for results? [dir-path]/testres.txt
%see also dependent functions:
%from ROEX2.f90:
%   READINPUT, CALC_ABS, NSCALC, FILTERWT THRESH, QROMB, TRAPZD, 
%from ROEXSUB.f90:
%   THR_CORR, LOCATE, POLINT, BRENT, SIMPLEX

%--------------------------------------------------------------
% roex3.for
%** Program for deriving auditory filter shapes from notched-noise ***
%********** data assuming filters with a roex(p,r) shape ************
% The language used is FORTRAN 77
% The program, subroutines and functions are as follows:
%       Main program
%	subroutine input(title,corr_file,No)
%	subroutine calc_abs(adjust,cf_corr,corr_file)
%	real*8 function sscalc(par)
%	real*8 function nscalc(ns_shift)
%	real*8 function p_erb(shift1)
%	real*8 function filter_wt(pg,r_lin)
%	real*8 function thresh(ns_shift)
%	subroutine qromb(a,b,ss,p,r_lin,hz1,hz2)
%	subroutine trapzd(a,b,s,n,p,r_lin,hz1,hz2)
% roexsub.for contains the routines:
%	subroutine openfiles()
%	character*(*) function noblanks(text)
%	real*8 function thr_corr(freqHz,n_thr,corr_hz,corr_dB)
%	subroutine locate(xx,n,x,j)
%	subroutine polint(xa,ya,n,x,y,dy)
%	subroutine simplex(p,y,mp,np,ndim,ftol,funk,iter)
%	real*8 function brent(ax,bx,cx,f,tol,xmin)
%---------------------------------------
%  The input file should have the following format:
%  1. Title
%  2. The number of data points (npts), how often the progress of the
%     fitting should be printed, e.g every ten iterations (n_print), and
%     the maximum allowed value of pl or pu (pmax)
%  3. The noise spectrum level in dB (No), the centre frequency in Hz
%     (cf), the noise bandwidth as a proportion of cf (width), and the
%      name of the correction file (corr_file)
%  4. Start values for the parameters & maximum shift: pl,pu,r,max_shift
%  5. The threshold values in dB
%  6. The lower & upper edges of the notch for each data point, in the
%     same order as the threshold values
%  See below for an example data file
%--------------------------------------
% main program
%--------------------------------------------------------------

% calling globalVars as the function includes roex3.h file
%#ok<*SAGROW>
globalVars
nLinesPerDataSet = 10;
% Get input and Result file names:
% See ROEXSUB.f90 --> #11-25
dataFilename = input(' data file name? ', 's');
resultFilename = input(' file name for results? ', 's');

if exist(resultFilename, 'file')
   delete(resultFilename);
end

% Deviation from fortran:
% Read all lines from datafile and pass group of 10 lines for processing
rawData = readDatafile(dataFilename);
nDatasets = fix(numel(rawData)/nLinesPerDataSet);

% Translated from ROEX3.f90 --> #55-134
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
for nD = 0:nDatasets-1
    startLineNo = (nD*nLinesPerDataSet)+1;
    endLineNo = (nD+1)*nLinesPerDataSet;
    % calling globalVars as the function includes roex3.h file
    globalVars
    % init vars
    iter = 0;
    np = ndim; mp = np +1; ftol = 1e-4;
    % see ROEX3.f90 --> #46
    st_par(ndim) = 0; 
    % see ROEX3.f90 --> #48
    params(mp,np) = 0; ssq(mp)=0; temp(np)= 0; funk = @sscalc;
    % see ROEX3.f90 --> #49
    cf_corr = 0;
    
    currDataset = rawData(startLineNo:endLineNo);
    [title, corr_file, No] = inputParse(currDataset,resultFilename);
    adjust = No + 10.0d0 * log10(cf);
    cferb = c1 * (c2 * cf / 1000.0 + 1.0);
    [adjust, cf_corr, corr_file] = calc_abs(adjust, cf_corr, corr_file);
    sscheck = 5001.0;
    ssq(1) = 5000.0;
    st_diff = 1.1;
    j = 1;
    
    % see ROEX3.f90 --> #64
    % Label 30 is the return point for second and subsequent entries to simplex.
    sscount = 0;
    st_par(1) = pl;
    st_par(2) = pu;
    st_par(3) = r;
    %Decide whether to repeat simplex

    % convert if else to while loop and embed goto30 
    while (ssq(j)<sscheck - 0.1 && ssq(j)>0.5)
        sscheck = ssq(j);
        if sscheck<5000.0 && n_print>0
            fprintf(' restart ssq = %9.3f\n', ssq(j));
        end
        % Set starting vertices of simplex.
        for j = 1:mp
            for k = 1:np
                params(j, k) = st_par(k);
                if(j==k + 1), params(j, k) = params(j, k) / st_diff; end
                temp(k) = params(j, k);
            end
            % Calculate the sum of the squared deviations of the data from the
            % fitted values (ssq) for each vertex of the simplex.  See later for
            % function sscalc.
            ssq(j) = sscalc(temp);
        end
        % Perform multidimensional minimization using the downhill simplex method.
        [params, ssq, mp, np, ndim, ftol, iter] = simplex(params, ssq, mp, np, ndim, ftol, funk, iter);
        st_diff = st_diff * 2.0;
        j = 1;
        for k = 2:mp
            if(ssq(k)<ssq(j)), j = k; end
        end
        pl = abs(params(j, 1));
        pu = abs(params(j, 2));
        r = -abs(params(j, 3));
        % goto 30 --> repeat goto30 here
        sscount = 0;
        st_par(1) = pl;
        st_par(2) = pu;
        st_par(3) = r;
    end
    [temp(1), st_par] = sscalc(st_par);
    resultFid = fopen(resultFilename,'a+');
    fprintf(resultFid,'   sum of squares =%8.1f\n\n',ssq(j));
    r_lin = 10.^(r / 10.);
    erb = (1. - r_lin) / pl * (2. - (2. + low_erb_lim * pl)...
        * exp(-low_erb_lim * pl)) + low_erb_lim * r_lin;
    erb = erb + (1. - r_lin) / pu * (2. - (2. + up_erb_lim * pu)...
        * exp(-up_erb_lim * pu)) + up_erb_lim * r_lin;
    % commented out	  erb = 2./pl+2./pu
    fprintf(resultFid,' Noise spectrum level =%5.1f  cf=%8.1f Hz\n',No, cf);
    
    fprintf(resultFid,' pl =%5.1f pu =%5.1f r =%6.1fB k(dB) =%5.1f erb =%6.3f( =%6.1f Hz)\n\n',...
        pl, pu, r, xk, erb, erb * cf);
    
    fprintf(resultFid,' The erb %6.3f is calc with integration limits low side:%5.2f high side:%5.2f\n',...
        erb, low_erb_lim, up_erb_lim);
    
    fprintf(resultFid, ' \tlower\tupper\tdata\tcalc\tresid\tshift\tgain(dB)\n');
    for notch = 1:npts
        % Return data and fitted values to initial form (before correction).
        data(notch) = data(notch) + adjust + cf_corr;
        calcdb(notch) = calcdb(notch) + adjust + cf_corr;
        fprintf(resultFid,[' ' repmat('\t%10.2f',1,7),'\n'],[el(notch), eu(notch), data(notch), calcdb(notch),diff_var(notch), shift(notch), gain_var(notch)]);
    end
    fprintf(resultFid,'           -------------------------\n');
    % goto 10
    %end
    % do for other dataset lines
end
clear global
close all

%% Sub functions
% Read all lines from an ascii file as text
function [ lines ] = readDatafile(inFilename)
  lines = {};
  fid = fopen(inFilename, 'r+');
  if fid < 0
      error( ' No such file %s\n\n',inFilename);
  end
  while ~feof(fid)
      lines{end+1,1} = fgetl(fid); %#ok<AGROW>
  end
  fclose(fid);
end
