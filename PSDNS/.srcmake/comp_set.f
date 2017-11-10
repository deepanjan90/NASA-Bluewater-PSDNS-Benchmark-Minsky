












      subroutine comp_set
 
      use comp
      use force_mod
 
      !#deepcustom#	implicit none

      integer :: szux,szuy,szuz,padxm,padym,padzm



c      allocate (mask(ny,zjsz*xisz))
 
      allocate (xfac(nxh),yfac(ny),zfac(nz))
      allocate (xdif(nxh,nc),ydif(ny,nc),zdif(nz,nc))
 
!	print *, 'comp_set: xisz,zisz,zjsz,yjsz,ux,uy,uz = ',xisz,zisz,zjsz,yjsz,
!    >	nxh*zisz*yjsz, xisz*zjsz*ny ,xisz*nz*yjsz
! 	if (zisz .ne. 2*xisz) then
! 	  print *, 'comp_set: zisz .ne. 2*xisz. PROGRAM WILL GIVE WRONG RESULTS!'
!	endif  

!! Needs work !

! for uneven distribution of processors arrays are not the same size
! size of COMPLEX arrays used later in proc2a, proc2b, sptr, etc.
        szux=nxhpad*zisz*yjsz
	szuy=xisz*zjsz*nypad
	szuz=xisz*nzpad*yjsz

	szmax=max(szux,szuy,szuz)         ! max size
	padx=(szmax-szux)/(nxhpad*zisz)      ! added space required for ux
	pady=(szmax-szuy)/(nypad*zjsz)     ! added space required for uy
	padz=(szmax-szuz)/(xisz*nzpad)       ! added space required for uz

	if ( mod(szmax-szux,nxhpad*zisz)  .ne. 0 ) padx=padx+1
	if ( mod(szmax-szuy,nypad*zjsz) .ne. 0 ) pady=pady+1
	if ( mod(szmax-szuz,xisz*nzpad)   .ne. 0 ) padz=padz+1

	call MPI_REDUCE(padx,padxm,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE(pady,padym,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE(padz,padzm,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) print *,'comp_set: max pads=',padxm,padym,padzm

        szux=nxhpad*zisz*(yjsz+padx)
	szuy=nypad*zjsz*(xisz+pady)
	szuz=xisz*nzpad*(yjsz+padz)
!	print *, 'comp_set: szmax,szux,szuy,szuz=',szmax,szux,szuy,szuz
	szf=nxhpad*zisz*(yjsz+padx)
	
      allocate (un(nxpad,zisz,yjsz+padx,3+nc),stat=ierr)
	if (ierr.ne.0) then
	write (6,*) 'error allocating un, taskid=',taskid,nxpad,zisz,yjsz+padx,3+nc
	end if
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	if (ierr.ne.0) call abrt('abort: error allocating un in comp_set')

      allocate (u(nxpad,zisz,yjsz+padx,nu),stat=ierr)
	if (ierr.ne.0) call abrt('error allocating u in comp_set')
	un(:,:,:,:)=0.
	u(:,:,:,:)=0.
c
	if (rkmethod.eq.4) then
      allocate (u1(nxpad,zisz,yjsz+padx,nu),stat=ierr)
        if (ierr.ne.0) call abrt('error allocating u1 in comp_set')
        u1(:,:,:,:)=0.
	end if

 
	if  (kforce.gt.0.) allocate(for(kfor,k2fo,k2fo,3))
	if  (kforce.gt.0.) allocate(velinc(2,kfor,k2fo,k2fo))
      allocate (efki(nxh,3),efkz(nxh,3))
      allocate (efk(nxh))
 
!	allocate (bk2i(nxh))
 

	if (ivstep.gt.1) allocate (cvx(nx,zisz,yjsz,3,2,rkmethod))
 
      return
      end
