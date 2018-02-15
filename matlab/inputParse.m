function [ title, corr_file, No ] = inputParse( datasetLines, resultFile )
%READINPUT Reads input data one dataset at a time
% Note: Input datasetLines is the 10 lines of data that is parsed
% Note: Result file are kept open
% Note: Reads 10 lines that setup the data for the dataset. is called till
% all 10 line sets are completed
% Note: Result file is appended to with write(8,..) fortran statements
% Calls: 
% Called by: ROEX3
%-------------------------------------
% This subroutine reads the data from the input file.
%-------------------------------------
% Translated from ROEX3.f90 --> #135-189
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran subroutine:  Return all or few of input arguments
% Since *title*, *corr_file*, *No* occurs on LHS, return only those, other input vars are not
% modified (do not occur on LHS)
%
%see also ROEX3
%
%  Copyright 2018
%  Ram Lab, Department of Psychology, Vanderbilt University.
%


   % initialize global vars
   %#ok<*NASGU>
   clear global
   globalVars
   resultFid = fopen(resultFile,'a+');
    % Line 1:
    % Ref: REOX3.f90 #145-146
    title = char(datasetLines{1});
    % Line 2: nDataPoints, screenOutputFreq, maxPvalue
    % Ref: REOX3.f90 #147-150
    temp = cell2mat(textscan(datasetLines{2},'%f,'));
    npts = temp(1);
    n_print = temp(2);
    pmax = 50.0;
    if numel(temp) == 3, pmax = temp(3); end
    % Line 3: noiseDb, centralFreq, noiseBand, correction_file
    % Ref: REOX3.f90 #151-154
    temp = cell2mat(textscan(datasetLines{3},'%f,'));
    No = temp(1);
    cf = temp(2);
    width = temp(3);
    corr_file = regexp(datasetLines{3},',([a-z]*)$','tokens');
    if isempty(corr_file)
        corr_file = 'none';
    else
        corr_file = char(corr_file{1});
    end
    % Line 4: plStart, puStart, rStart, filterShift (default =
    % 0.15, set to less than 0.01 for no shift)
    % low_erb_limit defaul;t = 1.0
    % up_erb_limit default=1.0
    % Ref: REOX3.f90 #155-162
    temp = cell2mat(textscan(datasetLines{4},'%f,'));
    pl = temp(1);
    pu = temp(2);
    r = temp(3);
    % If no value is specified for max_shift, default is 0.15.  If zero shift is
    % desired, put a very small value for max_shift (0.01 or less) in data file.
    max_shift = 0.15;
    if numel(temp) >= 4, max_shift = max(temp(4),max_shift); end
    % default of 1.0 for the integration limits when calculating the erb
    low_erb_lim = 1.0;
    if numel(temp) >= 5, low_erb_lim = max(temp(5),low_erb_lim); end
    up_erb_lim = 1.0;
    if numel(temp) >= 6, up_erb_lim = max(temp(6),up_erb_lim); end
    % Read all data points
    % Ref: REOX3.f90 #163
    data = cell2mat(cellfun(@(x) cell2mat(textscan(x,'%f, ')), ...
                     datasetLines(5:7),'UniformOutput',false));
    % Read el/eu points
    % Ref: REOX3.f90 #164
    temp = cell2mat(cellfun(@(x) cell2mat(textscan(x,'%f,%f, ')), datasetLines(8:10),'UniformOutput',false));
    el = temp(:,1);
    eu = temp(:,2);
    clearvars temp
    % Ref: REOX3.f90 #165-181
    if strcmpi(corr_file,'none')
        numeric_var = false;
        fprintf(resultFid,' roex3 : ( analytical integration)\n');
        fprintf(resultFid,' No threshold correction used\n');
    else
        numeric_var = true; 
        fprintf(resultFid,' roex3 : ( numerical integration)\n');
        fprintf(resultFid,' Correction file : %s\n', corr_file);
    end
    fprintf(resultFid,' Maximum allowed shift = %6.2f\n\n', max_shift);
    fprintf(resultFid,' %s\n', title);
    % Ref: REOX3.f90 #182-188
    if n_print > 0
        fprintf(' roex3 \n');
        fprintf(' correction file : %s\n',corr_file);
        fprintf(' ''%s''\n', title);
    end
    %Initialize other vars for dataset
    % nmax = max number of points, initialize to npts
    nmax = npts;
    calcdb(nmax) = 0;
    gain_var(nmax) = 0; diff_var(nmax) = 0; shift(nmax) = 0;
    lin_corr(max_corr_data) = 0;
end




