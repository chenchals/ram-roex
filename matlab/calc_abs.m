function [adjust, cf_corr, corr_file] = calc_abs(adjust, cf_corr, corr_file)
%CACL_ABS Compute correction values based on correction file
% Calls: THR_CORR
% Called by: main script
%-------------------------------------
% This subroutine uses the values in the correction file to determine
% correction values at 10-Hz spacing from 0 to 15000 Hz and stores them
% in array lin_corr as linear intensities.
%-------------------------------------
% Translated from ROEX3.f90 --> #190-236
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/REOX3.f90
% Fortran subroutine:  Return all or few of input arguments
% Updates global vars data, lin_corr
%
%see also THR_CORR

    % calling globalVars as the function includes roex3.h file
    globalVars;
    max_corr = 100;
    if ~numeric_var
        cf_corr = 0.0;
        for notch = 1:npts
            data(notch) = data(notch) - adjust;
        end
        return;
    end
    try
        fid = fopen(corr_file,'r');
        n_thr = cell2mat(textscan(fgetl(fid),'%f'));
        if n_thr > max_corr
            error(' too many points [%f] in the correction file [%s]. Max corrections [%f].\n',...
                n_thr,corr_file,max_corr);
        end
        corr_hz = cell2mat(textscan(fgetl(fid),'%f,'));
        corr_dB = cell2mat(textscan(fgetl(fid),'%f,'));
        fclose(fid);
        if numel(corr_hz)~=numel(corr_dB)~=n_thr
          error(' number of threshold corrections inconsistent in correction file [%s].\n', corr_file);   
        end
        % Calculate correction for the specified centre frequency.
        cf_corr = thr_corr(cf, n_thr, corr_hz, corr_dB);
        for notch = 1:npts
            data(notch) = data(notch) - adjust - cf_corr;
        end
        % Convert index values to frequencies in Hz, calculate threshold correction
        % at each frequency in linear intensity, and place results in array lin_corr.
        for i = 1:max_corr_data
            hz = real(i - 1) * 10.0;
            lin_corr(i) = 10.0 ^ (-thr_corr(hz, n_thr, corr_hz, corr_dB) / 10.0);
        end 
    catch me
       % on file read error stop
       error(' error opening threshold correction file [%s]\n',corr_file);
    end
end

