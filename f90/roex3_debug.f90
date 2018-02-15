! roex3.for

!** Program for deriving auditory filter shapes from notched-noise ***
!********** data assuming filters with a roex(p,r) shape ************
! The language used is FORTRAN 77
! The program, subroutines and functions are as follows:
!       Main program
!	subroutine input(title,corr_file,No)
!	subroutine calc_abs(adjust,cf_corr,corr_file)
!	real*8 function sscalc(par)
!	real*8 function nscalc(ns_shift)
!	real*8 function p_erb(shift1)
!	real*8 function filter_wt(pg,r_lin)
!	real*8 function thresh(ns_shift)
!	subroutine qromb(a,b,ss,p,r_lin,hz1,hz2)
!	subroutine trapzd(a,b,s,n,p,r_lin,hz1,hz2)

! roexsub.for contains the routines:
!	subroutine openfiles()
!	character*(*) function noblanks(text)
!	real*8 function thr_corr(freqHz,n_thr,corr_hz,corr_dB)
!	subroutine locate(xx,n,x,j)
!	subroutine polint(xa,ya,n,x,y,dy)
!	subroutine simplex(p,y,mp,np,ndim,ftol,funk,iter)
!	real*8 function brent(ax,bx,cx,f,tol,xmin)

!---------------------------------------
!  The input file should have the following format:
!  1. Title
!  2. The number of data points (npts), how often the progress of the
!     fitting should be printed, e.g every ten iterations (n_print), and
!     the maximum allowed value of pl or pu (pmax)
!  3. The noise spectrum level in dB (No), the centre frequency in Hz
!     (cf), the noise bandwidth as a proportion of cf (width), and the
!      name of the correction file (corr_file)
!  4. Start values for the parameters & maximum shift: pl,pu,r,max_shift
!  5. The threshold values in dB
!  6. The lower & upper edges of the notch for each data point, in the
!     same order as the threshold values
!  See below for an example data file
!--------------------------------------
! main program
implicit none
include 'roex3.h'
integer k, j, iter, mp, np
real(kind = 8) ftol, st_par(ndim)
parameter(np = ndim, mp = np + 1, ftol = 1.0d-4)
real(kind = 8) params(mp, np), ssq(mp), sscalc, temp(np)
real(kind = 8) No, adjust, sscheck, cf_corr
real(kind = 8) erb, st_diff, r_lin
character title*60, corr_file*40
external sscalc

call openfile()
10    call input(title, corr_file, No)
adjust = No + 10.0d0 * log10(cf)
cferb = c1 * (c2 * cf / 1000.0d0 + 1.0d0)
call calc_abs(adjust, cf_corr, corr_file)
sscheck = 5001.0d0
ssq(1) = 5000.0d0
st_diff = 1.1d0
j = 1

! Label 30 is the return point for second and subsequent entries to simplex.
30    sscount = 0
st_par(1) = pl
st_par(2) = pu
st_par(3) = r
! Decide whether to repeat simplex.
if(ssq(j)<sscheck - 0.1.and.ssq(j)>0.5d0) then
    sscheck = ssq(j)
    if(sscheck<5000.0d0.and.n_print>0)&
            write(6, 9000) ssq(j)
    9000    format(' restart ssq = ', f9.3)
    ! Set starting vertices of simplex.
    do j = 1, mp
        do k = 1, np
            params(j, k) = st_par(k)
            if(j==k + 1) params(j, k) = params(j, k) / st_diff
            temp(k) = params(j, k)
        end do
        ! Calculate the sum of the squared deviations of the data from the
        ! fitted values (ssq) for each vertex of the simplex.  See later for
        ! function sscalc.
        ssq(j) = sscalc(temp)
    end do
    ! Perform multidimensional minimization using the downhill simplex method.
    call simplex(params, ssq, mp, np, ndim, ftol, sscalc, iter)
    st_diff = st_diff * 2.0d0
    j = 1
    do k = 2, mp
        if(ssq(k)<ssq(j)) j = k
    end do
    pl = abs(params(j, 1))
    pu = abs(params(j, 2))
    r = -abs(params(j, 3))
    goto 30
else
    temp(1) = sscalc(st_par)
    write(8, 9010)ssq(j)
    9010    format('   sum of squares =', f8.1, /)
    r_lin = 10.**(r / 10.)
    erb = (1. - r_lin) / pl * (2. - (2. + low_erb_lim * pl)&
            * exp(-low_erb_lim * pl)) + low_erb_lim * r_lin
    erb = erb + (1. - r_lin) / pu * (2. - (2. + up_erb_lim * pu)&
            * exp(-up_erb_lim * pu)) + up_erb_lim * r_lin
    !commented out	  erb = 2./pl+2./pu

    write(8, 9015)No, cf
    9015    format(' Noise spectrum level =', f5.1, '  cf =', f8.1, ' Hz')
    write(8, 9020)pl, pu, r, xk, erb, erb * cf
    9020    format(' pl =', f5.1, ' pu =', f5.1, ' r =', f6.1, 'dB  k(dB) =', f5.1, &
            ' erb =', f6.3, '( =', f6.1, ' Hz)', /)

    write(8, 9025)erb, low_erb_lim, up_erb_lim
    9025    format(' The erb ', f6.3, ' is calc with integration limits ', &
            'low side:', f5.2, ' high side:', f5.2)
    write(8, 9040)
    9040    format(6x, 'lower     upper     data      calc      resid', &
            '     shift     gain(dB)')

    do notch = 1, npts
        ! Return data and fitted values to initial form (before correction).
        data(notch) = data(notch) + adjust + cf_corr
        calcdb(notch) = calcdb(notch) + adjust + cf_corr
        write(8, 9050)el(notch), eu(notch), data(notch), calcdb(notch), &
                diff(notch), shift(notch), gain(notch)
        9050      format(1x, 5f10.2, 2f10.2)
    end do
    write(8, 9060)
    9060    format(10x, '-------------------------', /)
    goto 10
endif
end
!-------------------------------------
! This subroutine reads the data from the input file.
subroutine input(title, corr_file, No)
    implicit none
    include 'roex3.h'
    real(kind = 8) No
    character title*(*)
    character corr_file*(*)
    character noblanks*40

    read(7, 9000)title
    9000    format(a)
    read(7, 9010)npts, n_print, pmax
    9010    format(2i6, f6.0)
    if(npts>nmax) stop 'end of data'
    if(pmax==0.0d0) pmax = 50.0d0
    read(7, 9020)No, cf, width, corr_file
    9020    format(3f10.0, a)
    corr_file = noblanks(corr_file)
    if(corr_file==' ') corr_file = 'none'
    read(7, 9030)pl, pu, r, max_shift, low_erb_lim, up_erb_lim
    9030    format(7f10.0)
    ! If no value is specified for max_shift, default is 0.15.  If zero shift is
    ! desired, put a very small value for max_shift (0.01 or less) in data file.
    if(max_shift==0.0d0) max_shift = 0.15d0
    ! default of 1.0 for the integration limits when calculating the erb
    if(low_erb_lim==0.0) low_erb_lim = 1.0
    if(up_erb_lim==0.0)  up_erb_lim = 1.0
    read(7, *)(data(notch), notch = 1, npts)
    read(7, *)(el(notch), eu(notch), notch = 1, npts)
    if(corr_file=='none') then
        numeric = .false.
        write(8, 9060)
        9060    format(' roex3 : ( analytical integration)')
        write(8, 9070)
        9070    format(' No threshold correction used')
    else
        numeric = .true.
        write(8, 9080)
        9080    format(' roex3 : ( numerical integration)')
        write(8, 9090)corr_file
        9090    format(' Correction file : ', a)
    endif
    write(8, 9100) max_shift
    9100    format(' Maximum allowed shift = ', f6.2, /)
    write(8, 9110)title
    9110    format(1x, a)
    if(n_print>0) then
        write(6, 9120)
        9120    format(' roex3 ')
        write(6, 9130)corr_file
        9130    format(' correction file : ', a, /)
        write(6, *)title
    endif
end
!------------------------------------
! This subroutine uses the values in the correction file to determine
! correction values at 10-Hz spacing from 0 to 15000 Hz and stores them
! in array lin_corr as linear intensities.
subroutine calc_abs(adjust, cf_corr, corr_file)
    implicit none
    include 'roex3.h'
    integer i, n_thr, ioval, max_corr
    parameter(max_corr = 100)
    real(kind = 8) adjust, cf_corr, thr_corr, hz, corr_hz(max_corr)
    real(kind = 8) corr_dB(max_corr)
    character corr_file*(*)

    if(.not.numeric) then
        cf_corr = 0.0d0
        do notch = 1, npts
            data(notch) = data(notch) - adjust
        end do
        return
    else
        open (10, file = corr_file, err = 50, status = 'old', iostat = ioval)
        read (10, 8990) n_thr
        8990      format(i10)
        if(n_thr>max_corr) then
            write(6, 9000)corr_file
            9000        format(' too many points in file ', a)
            stop
        endif
        read (10, *)(corr_hz(i), i = 1, n_thr)
        read (10, *)(corr_dB(i), i = 1, n_thr)
        close (10)
        ! Calculate correction for the specified centre frequency.
        cf_corr = thr_corr(cf, n_thr, corr_hz, corr_dB)
        do notch = 1, npts
            data(notch) = data(notch) - adjust - cf_corr
        end do
        ! Convert index values to frequencies in Hz, calculate threshold correction
        ! at each frequency in linear intensity, and place results in array lin_corr.
        do i = 1, max_corr_data
            hz = real(i - 1) * 10.0d0
            lin_corr(i) = &
                    10.0d0 ** (-thr_corr(hz, n_thr, corr_hz, corr_dB) / 10.0d0)
        end do
    endif
    return
    50    stop ' error opening threshold correction file '
end
!--------------------------------------
! This function returns the total sum of squares difference between the
! threshold data and the calculated values using parameters held in par (as
! determined by the starting values, or, subsequently, by the simplex subroutine).
real(kind = 8) function sscalc(par)
    implicit none
    include 'roex3.h'
    real(kind = 8) par(ndim), ns, nscalc, sum, thresh
    real(kind = 8) tol1, neg_el, pos_eu, brent
    parameter(tol1 = 1.0d-3)
    external nscalc

    sscount = sscount + 1
    pl = abs(par(1))
    pu = abs(par(2))
    r = -abs(par(3))
    sum = 0.0d0
    do notch = 1, npts
        write(*,*)'notch=',notch,'*******************************'
        if(el(notch)==0.0d0.and.eu(notch)==0.0d0&
                .or.max_shift<0.01) then
            ns = thresh(0.0d0)
            shift(notch) = 0.0d0
        else
            ! Find the lowest noise-to-signal ratio (ns) and the shift at which it
            ! occurs for each notch width.
            neg_el = -el(notch)
            if(neg_el < -max_shift) neg_el = -max_shift
            pos_eu = eu(notch)
            if(pos_eu > max_shift)  pos_eu = max_shift
            write(*,*)'sscalc calling brent......'
            ns = brent(neg_el, 0.0d0, pos_eu, nscalc, tol1, shift(notch))
        endif
        write(*,*)'notch=',notch,' ns=',ns
        calcdb(notch) = 10. * log10(ns)
        gain(notch) = 10.0d0 * log10(thresh(0.0d0)) - calcdb(notch)
        sum = sum - calcdb(notch) + data(notch)
    end do
    xk = sum / dble(float(npts))
    sscalc = 0.0d0
    do notch = 1, npts
        calcdb(notch) = calcdb(notch) + xk
        diff(notch) = data(notch) - calcdb(notch)
        sscalc = sscalc + diff(notch)**2
    end do
    ! This is a safe way to set limits on the parameters.
    if(pl>pmax.or.pu>pmax)   sscalc = sscalc + 1000.0d0
    if(pl<0.1d0.or.pu<0.1d0) sscalc = sscalc + 1000.0d0
    if(n_print>0) then
        if(mod(sscount, n_print)==0.or.sscount==1)&
                write(6, 9000)sscount, pl, pu, r, sscalc
    endif
    9000    format(1x, i4, ' pl =', f7.3, '  pu =', f7.3, '  r = ', f9.2, &
            ' ssq =', f9.3)
end
!----------------------------------
! This function calculates the noise-to-signal ratio for a filter shifted by 
! ns_shift from the on-frequency filter.
real(kind = 8) function nscalc(ns_shift)
    implicit none
    include 'roex3.h'
    real(kind = 8) r_lin, sig_int, push, plsh, filterwt
    real(kind = 8) ns_shift, thresh, off_shift, pconst, p_const
    real(kind = 8) pu_shifted, pl_shifted

    write(*,*)'in nscalc ns_shift=',ns_shift
    r_lin = 10.**(r / 10.)

    ! The shift passed in is relative to the on-frequency filter but needs to be
    ! measured from the shifted (off-frequency) filter.
    ! sig_int is the weighting applied by the shifted filter at the signal frequency.

    pconst = p_const(ns_shift)
    pu_shifted = pu * pconst
    pl_shifted = pl * pconst
    off_shift = ns_shift / (1.0d0 + ns_shift)
    if(off_shift==0.0d0) then
        sig_int = 1.0d0
    else if(off_shift<0.0d0) then
        push = pu_shifted * (-off_shift)
        if(push>maxexp) push = maxexp
        sig_int = filterwt(push, r_lin)
    else if(off_shift>0.0d0) then
        plsh = pl_shifted * off_shift
        if(plsh>maxexp) plsh = maxexp
        sig_int = filterwt(plsh, r_lin)
    endif
    ! Calulate the noise-to-signal ratio at the output of the shifted filter.
    nscalc = thresh(ns_shift) / sig_int
    if(nscalc<=0.0d0) then
        write(6, 9000)
        9000        format(' nscalc .le.  zero !! ')
        nscalc = exp(-maxexp)
    endif
end
!--------------------------------------
! Equation of the filter weighting function: pg = p * g
real(kind = 8) function filterwt(pg, r_lin)
    implicit none
    real(kind = 8) pg, r_lin
    filterwt = (1.0d0 - r_lin) * (1.0d0 + pg) * exp(-pg) + r_lin
end
!---------------------------------------
! This function returns the amount of noise passing through the filter using
! either numerical integration (if a threshold correction file was specified) 
! or analytic integration (if no correction file was specified).
real(kind = 8) function thresh(ns_shift)
    implicit none
    include 'roex3.h'
    real(kind = 8) r_lin, cu, cl, hzl1, hzl2, hzu1, hzu2
    real(kind = 8) ns_shift, gl, gu, sum_low, sum_hih
    real(kind = 8) pgush, pglsh, pucu, plcl, thresh_l, thresh_u, thrint_l, thrint_u
    real(kind = 8) pl_shifted, pu_shifted, shifted_cf, pconst, p_const

    r_lin = 10.0d0**(r / 10.0d0)
    hzl1 = (1.0d0 - el(notch) - width) * cf
    if(hzl1<0.0d0) hzl1 = 0.0d0
    hzl2 = (1.0d0 - el(notch)) * cf
    hzu1 = (1.0d0 + eu(notch)) * cf
    hzu2 = hzu1 + width * cf

    shifted_cf = cf * (1.0d0 + ns_shift)
    gu = abs(hzu1 - shifted_cf) / shifted_cf
    cu = abs(hzu2 - shifted_cf) / shifted_cf
    gl = abs(shifted_cf - hzl2) / shifted_cf
    cl = abs(shifted_cf - hzl1) / shifted_cf

    pconst = p_const(ns_shift)
    pu_shifted = pu * pconst
    pl_shifted = pl * pconst

    if(numeric) then
        !  Integrate over equally spaced 'g' values.
        call qromb(gl, cl, sum_low, pl_shifted, r_lin, hzl2, hzl1)
        call qromb(gu, cu, sum_hih, pu_shifted, r_lin, hzu1, hzu2)
        thresh = sum_hih + sum_low
        ! correct for the cf of the shifted filter being used
        thresh = thresh * shifted_cf / cf
    else
        !  Analytic integration.
        pgush = pu_shifted * gu
        pglsh = pl_shifted * gl
        if(pgush>maxexp) pgush = maxexp
        if(pglsh>maxexp) pglsh = maxexp
        pucu = pu_shifted * cu
        plcl = pl_shifted * cl
        if(plcl>maxexp) plcl = maxexp
        if(pucu>maxexp) pucu = maxexp
        thrint_l = (r_lin - 1.0d0) * (2.0d0 + plcl) * exp(-plcl) / pl_shifted&
                + r_lin * cl
        thresh_l = (r_lin - 1.0d0) * (2.0d0 + pglsh) * exp(-pglsh) / pl_shifted&
                + r_lin * gl
        thrint_l = thrint_l - thresh_l
        thrint_u = (r_lin - 1.0d0) * (2.0d0 + pucu) * exp(-pucu) / pu_shifted&
                + r_lin * cu
        thresh_u = (r_lin - 1.0d0) * (2.0d0 + pgush) * exp(-pgush) / pu_shifted&
                + r_lin * gu
        thrint_u = thrint_u - thresh_u
        thresh = thrint_l + thrint_u
        ! correct for the cf of the shifted filter being used
        thresh = thresh * shifted_cf / cf
    endif
    if(thresh<=0.0d0) then
        write(6, 9000)
        9000        format(' thresh le 0.0d0')
        thresh = exp(-maxexp)
    endif
    write(*,*)'   thresh called with ns_shift=',ns_shift, '  thresh=', thresh
end
!--------------------------------------
! Subroutine for performing numerical integration by Romberg's method.
! See Numerical Recipes, p. 114-115.
! This routine has variables p,r_lin,hz1 and hz2 passed through it to trapzd.
subroutine qromb(a, b, ss, p, r_lin, hz1, hz2)
    implicit none
    integer jmax, j, jmaxp, k, km
    real(kind = 8) a, b, ss, dss, p, r_lin, hz1, hz2, eps
    ! eps = 1.0d-2 is almost certainly good enough and faster
    parameter (eps = 1.0d-5, jmax = 20, jmaxp = jmax + 1, k = 5, km = k - 1)
    real(kind = 8) s(jmaxp), h(jmaxp)

    h(1) = 1.0d0
    do j = 1, jmax
        call trapzd(a, b, s(j), j, p, r_lin, hz1, hz2)
        if(j>=k) then
            call polint(h(j - km), s(j - km), k, 0.0d0, ss, dss)
            if(abs(dss)<eps * abs(ss)) return
        endif
        s(j + 1) = s(j)
        h(j + 1) = 0.25 * h(j)
    end do
    ! pause ' qromb: Too many steps'
    write(6, *)' qromb: Too many steps... Continuing'
end
!----------------------------------------
! This routine decreases the fineness of spacing of points in the numerical
! integration until further changes make no significant difference.
! See Numerical Recipes, p. 111.
! This version of the routine calls function filterwt(pg,wlin) by name.
! It also calculates the threshold correction to be applied to each
! point before doing the integration.
subroutine trapzd(a, b, s, n, p, r_lin, hz1, hz2)
    implicit none
    include 'roex3.h'
    integer n, it, j, index1, index2
    real(kind = 8) filterwt, a, b, s, sum, x, tnm, del, p, pg, pgb, r_lin
    real(kind = 8) delhz, xhz, hz1, hz2

    if(n==1) then
        pg = p * a
        pgb = p * b
        xhz = hz1
        delhz = hz2
        index1 = ((xhz + 5.0d0) / 10.0d0) + 1
        index2 = ((delhz + 5.0d0) / 10.0d0) + 1
        if(index1>max_corr_data) index1 = max_corr_data
        if(index2>max_corr_data) index2 = max_corr_data
        s = 0.5 * (b - a) * (filterwt(pg, r_lin) * lin_corr(index1)&
                + filterwt(pgb, r_lin) * lin_corr(index2))
    else
        it = 2**(n - 2)
        tnm = it
        del = (b - a) / tnm
        delhz = (hz2 - hz1) / tnm
        x = a + 0.5 * del
        xhz = hz1 + 0.5 * delhz
        sum = 0
        do j = 1, it
            index1 = ((xhz + 5.0d0) / 10.0d0) + 1
            if(index1>max_corr_data) index1 = max_corr_data
            pg = p * x
            sum = sum + filterwt(pg, r_lin) * lin_corr(index1)
            x = x + del
            xhz = xhz + delhz
        end do
        s = 0.5 * (s + (b - a) * sum / tnm)
    endif
end
!---------------------------------------------
real(kind = 8) function p_const(ns_shift)
    ! This function returns a constant used to multiply the centre
    ! frequency pl or pl to get the shifted frequency pl or pu.

    implicit none
    include 'roex3.h'
    real(kind = 8) ns_shift, sh_freq

    sh_freq = 1.0d0 + ns_shift
    p_const = (sh_freq * cferb) / (c1 * (c2 * sh_freq * cf / 1000.0d0 + 1.0d0))
    ! write(*,*)'          p_const=',p_const
end