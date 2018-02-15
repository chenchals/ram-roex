! The include file roex2.h, used in compilation, is :-
	integer nmax,ndim,notch,npts,n_print,sscount
	integer max_corr_data
	parameter (ndim = 2,nmax = 100,max_corr_data=1501)
	real(kind=8) el(nmax),eu(nmax),data(nmax),calcdb(nmax)
	real(kind=8) gain(nmax),diff(nmax)
	real(kind=8) lin_corr(max_corr_data),xk,width
	real(kind=8) shift(nmax),low_erb_lim,up_erb_lim
	real(kind=8) pl,pu,r,pmax,max_shift
	real(kind=8) cf,cferb,pr_erb
	real(kind=8) c1,c2,maxexp
	parameter (c1 = 24.673, c2 = 4.368, maxexp = 700.0d0)
! Maxexp should be set to something near the largest number a negative
! exponential can have on the system being used; e.g. exp(-700.0)
! is close to 10**(-304) and exp(-182) is close to 10**(-78).
	logical numeric
	common el,eu,data,calcdb,gain,diff,xk,width,notch,npts,n_print,sscount
	common /tc/ lin_corr,numeric
	common /sh/ shift,max_shift
	common /param/ pl,pu,r,pmax
	common /erbscale/ cf,cferb,pr_erb
	common /erblimit/ low_erb_lim,up_erb_lim
