function [  p, y, mp, np, ndim, ftol, iter ] = simplex( p, y, mp, np, ndim, ftol, funk, iter )
%SIMPLEX Summary of this function goes here
% Note: funk = SSCALC function handle
% Calls: AMOTRY
% Called by: ROEX3
% simplex(params, ssq, mp, np, ndim, ftol, sscalc, iter)
%-------------------------------------
%  Modname   simplex.f   from 2nd edition
%	Numerical Recipes: Press et al 2nd ed. p404
%-------------------------------------
% Translated from ROEXSUB.f90 --> #127-205
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
% Fortran subroutine:  Return all or few of input arguments
% Since only ????? occur on LHS, return only that, other input vars are not
% modified (do not occur on LHS)
% 
% see also ROEX3 AMOTRY SSCALC

   %p(mp, np), y(mp), funk, alpha, beta, gamma
   NMAX = 20; ITMAX = 5000; alpha = -1.0; beta = 0.5;gamma = 2.0;
   psum(NMAX) = 0.0;
   
   
   %% label 1 : Translated from ROEXSUB.f90 --> #143-149
   goto1();
   
   %% label 2 : Translated from ROEXSUB.f90 --> #151-158
   % while(true) do loop?
   while(1)
       ilo = 1;
       if y(1)>y(2)
           ihi = 1;
           inhi = 2;
       else
           ihi = 2;
           inhi = 1;
       end
       % Translated from ROEXSUB.f90 --> #159-167
       for i = 1:ndim + 1
           if y(i)<=y(ilo), ilo = i; end
           if y(i)>y(ihi)
               inhi = ihi;
               ihi = i;
           elseif y(i)>y(inhi)
               if i ~= ihi,  inhi = i; end
           end
       end
       % Translated from ROEXSUB.f90 --> #168
       rtol = 2. * abs(y(ihi) - y(ilo)) / (abs(y(ihi)) + abs(y(ilo)));
       % Translated from ROEXSUB.f90 --> #168-179
       if (rtol<ftol)
           swap = y(1);
           y(1) = y(ilo);
           y(ilo) = swap;
           for n = 1:ndim
               swap = p(1, n);
               p(1, n) = p(ilo, n);
               p(ilo, n) = swap;
           end
           return %% finally exits loop
       end
       if iter>=ITMAX, fprintf(' simplex exceeding maximum iterations\n'); end
       iter = iter + 2;
       [ytry, p, y, psum, mp, np, ndim, ihi, alpha] = amotry(p, y, psum, mp, np, ndim, funk, ihi, alpha);
       % Translated from ROEXSUB.f90 --> #183-203
       if ytry <= y(ilo)
           [ytry, p, y, psum, mp, np, ndim, ihi, gamma] = amotry(p, y, psum, mp, np, ndim, funk, ihi, gamma);
       elseif ytry >= y(inhi)
           ysave = y(ihi);
           [ytry, p, y, psum, mp, np, ndim, ihi, beta] = amotry(p, y, psum, mp, np, ndim, funk, ihi, beta);
           if(ytry >= ysave)
               for i = 1:ndim + 1
                   if i ~= ilo
                       for j = 1:ndim
                           psum(j) = 0.5 * (p(i, j) + p(ilo, j));
                           p(i, j) = psum(j);
                       end
                       y(i) = funk(psum);
                   end
               end
               iter = iter + ndim;
               goto1(); % and continue loop
           end
       else
           iter = iter - 1;
       end
       %goto 2
   end
   
    %% Nested function to access scoped variables   
    function goto1()
        for nn = 1:ndim
            sum = 0.0;
            for m = 1:ndim + 1
                sum = sum + p(m, nn);
            end
            psum(nn) = sum;
        end
    end

end

