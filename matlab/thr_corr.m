function [ thr_corr_out ] = thr_corr( freqHz, n_thr, corr_hz, corr_dB )
%THR_CORR Computes threshold correction using polynomial interpolation
% Calls: POLINT LOCATE
% Called by: CALC_ABS
%-------------------------------------
% This function returns a threshold correction in dB for a given frequency
% (freqHz) in Hz using polynomial interpolation with 3 points. See Numerical
% Recipes p. 80-83.
%-------------------------------------
% Translated from ROEXSUB.f90 --> #42-62
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
% Fortran function:  Returns thr_corr_out
%
%see also CALC_ABS POLINT LOCATE
  thr_corr_out = 0.0;
  dy = 0;
  index = 0;
  index = locate(corr_hz, n_thr, freqHz,index);
  if index == 0
      thr_corr_out = corr_dB(1);
  elseif index >= n_thr
      thr_corr_out = corr_dB(n_thr);
  elseif index == (n_thr-1)
      [thr_corr_out, ~] = polint( corr_hz(index), corr_dB(index), 2, freqHz, thr_corr_out, dy );
  else
      [thr_corr_out, ~] = polint( corr_hz(index), corr_dB(index), 3, freqHz, thr_corr_out, dy );
  end
end

