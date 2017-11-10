












      subroutine com_set
c
	use com
c
      !#deepcustom#	implicit none
      logical :: iex
c
        if (taskid.eq.0) write (6,*) 'enter com_set'
      allocate (luinit(nc+3))
c
      allocate (kx(nxhpad),ky(nypad),kz(nzpad))
c
      allocate (bk1(nxhpad,2))
      allocate (bk3(nzpad,2))
      allocate (bkk3(nzpad,2))
c
      allocate (sx(nxhpad,rkmethod),sy(nypad,rkmethod),sz(nzpad,rkmethod))
c
      allocate (cb(3,nc),svnorm(nc))
c

c#ifdef 1
      allocate (suo(nrproc),svo(nrproc))
c#endif
c
c
	allocate (fnwu(3+nc),luwu(3+nc)) !,fninit(3+nc))


	! to stop execution if tsp_max is larger than that specified in tsp_max
 	inquire(file='tsp_max',exist=iex)
	if (iex) then
	  open (999,file='tsp_max')
	  read (999,*) tsp_max
	  close (999)
	else
	  tsp_max=1.e6
	endif
c	
        if (taskid.eq.0) write (6,*) ' exit com_set'
      return
      end
