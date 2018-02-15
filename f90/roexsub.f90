! ROEXSUB.FOR. Some subroutines common to roex3, roex2, roexpr and roexp
!
! openfile
! noblanks
! thr_corr
! locate
! polint
! simplex
! brent
!---------------------------------
! Opens the data and results file; for an example data file see below.
subroutine openfile()
    implicit none
    character*40 filename

    write(6, 9000)
    9000    format(' data file name? ')
    read (5, 9010)filename
    9010    format(a)
    open (7, file = filename)
    write(6, 9020)
    9020    format(' file name for results? ')
    read (5, 9010) filename
    open (8, file = filename)
end
!-------------------------------------------
! This function removes any blanks or tabs from the string of text.
character*(*) function noblanks(text)
    implicit none
    integer i, j
    character text*(*)

    noblanks = ' '
    j = 1
    do i = 1, 40
        if(text(i : i)==' '.or.text(i : i)=='\t') goto 10
        noblanks(j : j) = text(i : i)
        j = j + 1
        10    continue
    end do
end
!----------------------------------------------
! This function returns a threshold correction in dB for a given frequency
! (freqHz) in Hz using polynomial interpolation with 3 points. See Numerical
! Recipes p. 80-83.
real(kind = 8) function thr_corr(freqHz, n_thr, corr_hz, corr_dB)
    implicit none
    integer n_thr, index
    real(kind = 8) corr_hz(*), corr_dB(*), freqHz, dy

    call locate(corr_hz, n_thr, freqHz, index)
    if(index==0) then
        thr_corr = corr_dB(1)
    else if(index>=n_thr) then
        thr_corr = corr_dB(n_thr)
    else if(index==(n_thr - 1)) then

        call polint(corr_hz(index), corr_dB(index), 2, freqHz, thr_corr, dy)
    else
        call polint(corr_hz(index), corr_dB(index), 3, freqHz, thr_corr, dy)
    endif
end
!---------------------------------------------
! This subroutine picks the numbers over which the interpolation is to be
! performed.  See Numerical Recipes p. 89-90.
subroutine locate(xx, n, x, j)
    implicit none
    integer ju, jl, jm, j, n
    real(kind = 8) xx(n), x

    jl = 0
    ju = n + 1
    10    if(ju - jl>1) then
        jm = (ju + jl) / 2
        if((xx(n)>xx(1)).eqv.(x>xx(jm))) then
            jl = jm
        else
            ju = jm
        endif
        goto 10
    endif
    j = jl
end
!-------------------------------------
! Subroutine polint performs polynomial interpolation.  See Numerical
! Recipes p. 80-83.
subroutine polint(xa, ya, n, x, y, dy)
    implicit none
    integer nmax, i, m, ns, n
    parameter (nmax = 10)
    real(kind = 8) x, y, dy, xa(n), ya(n), c(nmax), d(nmax)
    real(kind = 8) dif, dift, ho, hp, w, den

    ns = 1
    dif = abs(x - xa(1))
    do i = 1, n
        dift = abs(x - xa(i))
        if(dift<dif) then
            ns = i
            dif = dift
        endif
        c(i) = ya(i)
        d(i) = ya(i)
    end do
    y = ya(ns)
    ns = ns - 1
    do m = 1, n - 1
        do i = 1, n - m
            ho = xa(i) - x
            hp = xa(i + m) - x
            w = c(i + 1) - d(i)
            den = ho - hp
            !if(den==0.0d0) pause 'polint : den = 0.0 '
            write(6,*)'    polint : den = 0.0 ... Continuing'
            den = w / den
            d(i) = hp * den
            c(i) = ho * den
        end do
        if(2 * ns<n - m) then
            dy = c(ns + 1)
        else
            dy = d(ns)
            ns = ns - 1
        endif
        y = y + dy
    end do
end
!--------------------------------------------
!  Modname   simplex.f   from 2nd edition
!	Numerical Recipes: Press et al 2nd ed. p404

subroutine simplex(p, y, mp, np, ndim, ftol, funk, iter)
    implicit none
    integer iter, mp, ndim, np, NMAX, ITMAX
    real(kind = 8) ftol, p(mp, np), y(mp), funk, alpha, beta, gamma
    parameter (NMAX = 20, ITMAX = 5000, alpha = -1.0d0, beta = 0.5d0, &
            gamma = 2.0d0)
    external funk

    integer i, ihi, ilo, inhi, j, m, n
    real(kind = 8) rtol, sum, swap, ysave, ytry, psum(NMAX), amotry

    iter = 0
    1    do n = 1, ndim
        sum = 0.0d0
        do m = 1, ndim + 1
            sum = sum + p(m, n)
        end do
        psum(n) = sum
    end do

    2    ilo = 1
    if(y(1)>y(2)) then
        ihi = 1
        inhi = 2
    else
        ihi = 2
        inhi = 1
    endif
    do i = 1, ndim + 1
        if(y(i)<=y(ilo)) ilo = i
        if(y(i)>y(ihi)) then
            inhi = ihi
            ihi = i
        else if(y(i)>y(inhi)) then
            if(i/=ihi) inhi = i
        endif
    end do
    rtol = 2. * abs(y(ihi) - y(ilo)) / (abs(y(ihi)) + abs(y(ilo)))
    if (rtol<ftol) then
        swap = y(1)
        y(1) = y(ilo)
        y(ilo) = swap
        do n = 1, ndim
            swap = p(1, n)
            p(1, n) = p(ilo, n)
            p(ilo, n) = swap
        end do
        return
    endif
    if(iter>=ITMAX) print*, ' simplex exceeding maximum iterations'
    iter = iter + 2
    ytry = amotry(p, y, psum, mp, np, ndim, funk, ihi, alpha)
    if(ytry<=y(ilo)) then
        ytry = amotry(p, y, psum, mp, np, ndim, funk, ihi, gamma)
    else if(ytry>=y(inhi)) then
        ysave = y(ihi)
        ytry = amotry(p, y, psum, mp, np, ndim, funk, ihi, beta)
        if(ytry>=ysave) then
            do i = 1, ndim + 1
                if(i/=ilo)then
                    do j = 1, ndim
                        psum(j) = 0.5 * (p(i, j) + p(ilo, j))
                        p(i, j) = psum(j)
                    end do
                    y(i) = funk(psum)
                endif
            end do
            iter = iter + ndim
            goto 1
        endif
    else
        iter = iter - 1
    endif
    goto 2
end
!--------------------------------------
real(kind = 8) function amotry(p, y, psum, mp, np, ndim, funk, ihi, fac)
    implicit none
    integer ihi, mp, ndim, np, NMAX
    real(kind = 8) fac, p(mp, np), psum(np), y(mp), funk
    parameter(NMAX = 20)
    external funk
    integer j
    real(kind = 8) fac1, fac2, ytry, ptry(NMAX)

    fac1 = (1. - fac) / ndim
    fac2 = fac1 - fac
    do j = 1, ndim
        ptry(j) = psum(j) * fac1 - p(ihi, j) * fac2
    end do
    ytry = funk(ptry)
    if(ytry<y(ihi)) then
        y(ihi) = ytry
        do j = 1, ndim
            psum(j) = psum(j) - p(ihi, j) + ptry(j)
            p(ihi, j) = ptry(j)
        end do
    endif
    amotry = ytry
    return
end
!--------------------------------------------
! This is the subroutine brent which performs a one-dimensional minimization
! of the noise-to-signal ratio.  This is used for finding the optimum shift
! for each notch width. See Numerical Recipes, p. 283-286.
real(kind = 8) function brent(ax, bx, cx, f, tol, xmin)
    implicit none
    integer itmax, iter
    real(kind = 8) ax, bx, cx, f, tol, xmin
    real(kind = 8) cgold, zeps, a, b, v, x, e, fx, fu, fv, fw, xm
    real(kind = 8) tol1, tol2, w, q, p, r, etemp, d, u
    parameter (itmax = 100, cgold = 0.381966d0, zeps = 1.0d-10)
    external f

    a = min(ax, cx)
    b = max(ax, cx)
    v = bx
    w = v
    x = v
    e = 0.0d0
    fx = f(x)
    fv = fx
    fw = fx
    do iter = 1, itmax
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0d0 * tol1
        if(abs(x - xm)<=(tol2 - 0.5d0 * (b - a))) goto 3
        if(abs(e)>tol1) then
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0d0 * (q - r)
            if(q>0.0d0) p = -p
            q = abs(q)
            etemp = e
            e = d
            if(abs(p)>=abs(0.5d0 * q * etemp).or.p<=q * (a - x)&
                    .or.p>=q * (b - x)) goto 1
            d = p / q
            u = x + d
            if(u - a<tol2.or.b - u<tol2) d = sign(tol1, xm - x)
            goto 2
        endif
        1      if(x>=xm) then
            e = a - x
        else
            e = b - x
        endif
        d = cgold * e
        2      if(abs(d)>=tol1) then
            u = x + d
        else
            u = x + sign(tol1, d)
        endif
        fu = f(u)
        if(fu<=fx) then
            if(u>=x) then
                a = x
            else
                b = x
            endif
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
        else
            if(u<x) then
                a = u
            else
                b = u
            endif
            if(fu<=fw.or.w==x) then
                v = w
                fv = fw
                w = u
                fw = fu
            else if(fu<=fv.or.v==x.or.v==w) then
                v = u
                fv = fu
            endif
        endif
    end do
    !	print*,' brent exceeded maximum iterations'
    3    xmin = x
    brent = fx
end