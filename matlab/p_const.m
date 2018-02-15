function [ p_const_out, ns_shift ] = p_const( ns_shift)
%P_CONST Return a factor for getting shifted pl, pu from central frequency
% Calls: GLOBALVARS
% Called by: TRAPZD NSCALC THRESH
%-------------------------------------
% This function returns a constant used to multiply the centre
% frequency pl or pl to get the shifted frequency pl or pu.
%-------------------------------------
% Translated from ROEX3.f90 --> #470-480
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEX3.f90
% Fortran function:  Return only 1 value
%
% see also GLOBALVARS TRAPZD NSCALC THRESH

    % calling globalVars as the function includes roex3.h file
    globalVars
    sh_freq = 1.0 + ns_shift;
    p_const_out = (sh_freq * cferb) / (c1 * (c2 * sh_freq * cf / 1000.0 + 1.0));
    printDebug('           p_const=%1.20E\n',p_const_out);
end



