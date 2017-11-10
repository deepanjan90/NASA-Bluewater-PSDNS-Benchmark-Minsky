












c     
c     Line made redundant upon introduction of new inpen.f 
c     routine, PK Yeung, Dec 2007 have been removed
c     
c     routine to specify initial conditions
c     
      subroutine incond
      use comsp
      USE IO
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none
     
      character*2 str
      character*5 numer
      character*16 fn
      logical pic(ncd)
      integer i,izr,lstart
      real xnorm,ynorm,znorm,tmpvmaxz,vel
      integer m
      character(len=32) :: str1,str2

      m = 2


	if (taskid.eq.0) write (6,*) 'enter incond: b11(m)=',b11(m)
c      write (6,*) 'incond, entry, taskid=',taskid, 'irz=',irz
c     
c     
c----------------specified initial fields (un)  ----------------------
c     
      time0=0.
      istep0=0
      kstep=2
      if( istart .eq. 0 .or. any(kinit.eq.-1).or.any(kinit.eq.-3) ) then
c     
         un(:,:,:,1:3+nc)=0.
c     
c     
         go to 11
 9       write (6,*) 'all scalars with kinit=-2 must be placed first'
         write (6,*) 'edit ebdata.f and try again'
         stop 'aborts in sub. incond'
 11      continue
c     
c     initialize in Fourier space or read in previous data
c     
         call field(un)
c     
      end if
c     
c     
      if (taskid.eq.0) write (6,*) 'incond: irz=',irz
      if (any(kinit.gt.0)) then
         if (all(kinit.gt.0)) then
            if (hdf5_input_flag) then
               CALL RESTART_HDF5_DRIVER(TRIM(ADJUSTL(indir_fn)),hdf5_init,un)
            else
               if (fninit(1)(4:4).eq.'r') then
                  call incomm (un)
               else
                  call inpen (un)
               end if
            end if
         else
            if (hdf5_input_flag) then
               CALL RESTART_HDF5_DRIVER(TRIM(ADJUSTL(indir_fn)),hdf5_init,un)
            else
               if (fninit(1)(4:4).eq.'r') then
                  call incomm (un)
               else
                  call inpen (un)
               end if
            end if
         end if
         call time_stamp ('incond: after inpen')
      endif
c
	call time_stamp ('incond: before do 110')
      do 110 i=1,3
         u(:,:,:,i)=un(:,:,:,i)
 110  continue
	call time_stamp ('incond:  after do 110')

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	call time_stamp ('incond: at barrier')
      
!     option for differential diffusion: set one scalar to be
!     equal to the other at t=0
      
      do i=4,3+nc
         if (kinit(i).le.-10) then
            iscref=abs(kinit(i))-10
            u(:,:,:,i)=u(:,:,:,3+iscref)
            un(:,:,:,i)=un(:,:,:,3+iscref)
         end if
      end do
      
c     call inpvf to specify velocity field in physical space
c     if kinit(1)=-2 (in this case kshift must be equal to 4,
c     i.e., zero-shifts)
c     
      if (kinit(1).eq.-2) then
c     
         if (kshift.ne.4.and.nsteps.ge.1) then
            write (6,*) 'initial velocity field to be in physical space'
            write (6,*) 'but zero shifts must be selected,'
            write (6,*) 'i.e., use kshift=4 in ebdata'
            stop 'program aborts in sub. incond'
         end if
c     

c
         call inpvf(un,u)
c
      end if
c     
c     
c     
c----------------restart from checkpointed data ----------------------
c     
	call time_stamp ('incond: before store')
	
c     if (nx.ge.64) write (6,*) ' before call store'
                  if( istart .ne. 0 ) then
                     lstart = iabs( istart )
                     call store( 3 ,lstart )
c                     if( istart .gt. 0 ) then
                    if( istart .gt. 0.and.all(kinit.ge.0) ) then
                        
                        if (taskid.eq.0) rewind (luran1)
                        call ranseq( -1 , luran1 )
                     endif
                  endif
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
                   if (taskid.eq.0) write (6,*) 'initial data fields in place'    
                  if (istart.eq.0.and.tfiuo.lt.0.) then
                     call store (4,1)
                  end if
c     
c     initialize the time
c     
                  time=time0
                  istep=istep0
c     
c     produce spectra
c     
	if (taskid.eq.0) write (6,*) 'incond: b11(2)=',b11(2)
	if (b11(2).ne.beta1**2) then
	call init
	call waveno
	un(:,:,:,1)=un(:,:,:,1)/beta1
	un(:,:,:,2)=un(:,:,:,2)/beta2
	un(:,:,:,3)=un(:,:,:,3)/beta3
	if (taskid.eq.0) write (6,*) 'incond: b11(2)=',b11(2)
	end if


! spectrum calculated by sptvar routine is at present not correct
! except for the case of a 2pi^3 grid. Also must have nx >= ny & nz in
! order for the arrays in sptvar to not overflow.
      if (beta1.eq.1.0.and.beta2.eq.1.0.and.beta3.eq.1.0
     1    .and.nx.ge.ny.and.nx.ge.nz) then
	      call sptvar (un)
      end if
      str2 = 'AXI_SPTR'
	if (ioaxi(1).gt.0) then
      str1 = 'axi_kx'
      CALL AXI_SPTR(un,1,str1,str2)
	end if
	if (ioaxi(2).gt.0) then
      str1 = 'axi_ky'
      CALL AXI_SPTR(un,2,str1,str2)
	end if
	if (ioaxi(3).gt.0) then
      str1 = 'axi_kz'
      CALL AXI_SPTR(un,3,str1,str2)
	end if

        call sptr(un,1,1)
	call check_consym2 (un,3,'incond',3)
      return
      end
