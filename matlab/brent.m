function [ brent_out, ax, cx, tol, xmin ] = brent( ax, bx, cx, f, tol, xmin)
%BRENT Summary of this function goes here
% Calls: Uses NSCALC (argument f, as function handle)
% Called by: SSCALC
%-------------------------------------
% This is the subroutine brent which performs a one-dimensional minimization
% of the noise-to-signal ratio.  This is used for finding the optimum shift
% for each notch width. See Numerical Recipes, p. 283-286.
%-------------------------------------
% Translated from ROEXSUB.f90 --> #232-319
% See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
%  Fortran function:  Return only 1 value
%
%see also SSCALC

   printDebug('brent called with: ax=%1.20E  bx=%1.20E  cx=%1.20E  tol=%1.20E  xmin=%1.20E\n',ax,bx,cx,tol,xmin);

    itmax = 100; cgold = 0.381966; zeps = 1.0e-10;d=0;
    brent_out = 0.0;
    a = min(ax, cx);
    b = max(ax, cx);
    v = bx;
    w = v;
    x = v;
    e = 0.0;
    fx = f(x);
    fv = fx;
    fw = fx;
    for iter = 1:itmax
        xm = 0.5 * (a + b);
        tol1 = tol * abs(x) + zeps;
        tol2 = 2.0 * tol1;
        if (abs(x - xm)<=(tol2 - 0.5 * (b - a)))
            % goto 3 and end call
            printDebug('goto3@35 same as goto3@259 iter=%d\n',iter);
            goto3();
            return
        end
        if abs(e)>tol1
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0d0 * (q - r);
            if q>0.0, p = -p; end
            q = abs(q);
            etemp = e;
            e = d;
            if (abs(p)>=abs(0.5 * q * etemp) || p<=q * (a - x)...
                    || p>=q * (b - x))
                printDebug('goto1@50 same as goto1@269 iter=%d\n',iter);
                %goto1();
                goto1();
                d = cgold * e;
                goto2();
                continue
            else
                d = p / q;
                u = x + d;
                if (u - a<tol2 || b - u<tol2)
                    d = signFx(tol1, xm - x);
                    printDebug('goto2@58 same as goto1@277 iter=%d\n',iter);
                    goto2();
                    continue
                end
            end
        end
        printDebug('fall through...to label goto 1 iter=%d\n',iter);
        printDebug('a=%1.20E   b=%1.20E   e=%1.20E   d=%1.20E   x=%1.20E\n',a,b,e,d,x);

        goto1();
        d = cgold * e;
        goto2();        
        
    end
    %goto3();
    
    %% Nested functions so as to access scoped vars
    % Lines that do goto 1
    % Translated from ROEXSUB.f90 --> #275-280
    % See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
        function goto1()
            if(x>=xm)
                e = a - x;
            else
                e = b - x;
            end
        end


    % Lines that do goto 2
    % Translated from ROEXSUB.f90 --> #281-314
    % See: /Users/subravcr/IdeaProjects/ram-roex/ram-other/roex/ROEXSUB.f90
        function goto2()
            if(abs(d)>=tol1)
                u = x + d;
            else
                u = x + signFx(tol1, d);
            end
            fu = f(u);
            if(fu<=fx)
                if(u>=x)
                    a = x;
                else
                    b = x;
                end
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            else
                if(u<x)
                    a = u;
                else
                    b = u;
                end
                if(fu<=fw || w==x)
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                elseif(fu<=fv || v==x || v==w)
                    v = u;
                    fv = fu;
                end
            end
        end

    % Lines that do goto 3
        function goto3()
            xmin = x;
            brent_out = fx;
        end

end

function [val] = signFx(x,y)
    %Fortran sign function
    % CASE 1:   If y ? 0 then
    % sign(x,y) = abs(x)   ,
    % CASE 2:   If y < 0 then
    % sign(x,y) = - abs(x)   .
    if y >= 0
        val = abs(x);
    elseif y < 0
        val = -abs(x);
    end
end



