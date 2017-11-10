












!################################################
!#      module timestep_module
!#
!#  contains time_stepping routines, in particular, timestep
!################################################
	module timestep_module	
	use rungekutta_module
	use wavespace_module
	use comsp
	use timers_rkstep
	use timers_tran
	!#deepcustom#	implicit none
	
	public  :: timestep

!################################################
        contains
!################################################	
!#	
!#     Subroutine TIMESTEP
!#      Executes one full RK step
!#      Number of substeps depends on RKscheme
!#
!################################################
	subroutine timestep()
          integer i,iost

	scgflag=.false.

	flagc=.false.
          do i=1, rkscheme%rkstages
             if(taskid .eq. 0) then
		write(6,*) 'RK stage: ', i, 'flagc=',flagc
	     endif
             call rksubstep(i) !execute one substep
             !rkstuff

          end do
	
c
	!done with entire timestep
	end subroutine timestep

!################################################
!#   Subroutine KSTEP
!#
!# Executes one Runge Kutta substep, currently for Rk2 only
!# Replacing predic.f and correc.f in original 
!# homogeneous, isotropic code
!#
!# Requires input for which part of the step is executing
!# i.e. predictor or corrector
c
c PKY: attempted extension  to RK4
!#
!################################################
        subroutine rksubstep(rkstep)
          !#deepcustom#	implicit none
          !#deepcustom# 	implicit none
          real strof
          integer i,ii,rkstep,iost,test,j,k,m,a,ymax,almax,jj
          logical flagt
          character*20 auxiliar         
	  real xnorm,ynorm,znorm,tmpvmaxz,vel,vm,velmaxi(3)
	  integer ymaxi(3),almaxi(3)
	  real(b8), allocatable :: ut1(:,:,:)

          data iost/0/
          data strof/1.e-6/

	integer ikstep
	real(8) rtime0,rtime1,rtime2,rtime3,rtime4
	real(8) tt(10),t_sptr
	save tt
c
	logical vorflag
      character(len=32) :: str1,str2
c
	if (rkstep.eq.1) then
	tt(:)=0.
	end if
c
	rtime0=MPI_WTIME()
	rtime1=MPI_WTIME()
	
	  m = 2
          !Need to know what substep it is.
          kstep=rkstep

	vorflag=.false.
 

c obtain velocity field on shifted grid used at subsequent R-K steps
	if (ivstep.gt.1) then
	if (kstep.eq.1.and.mod(jstep-1,ivstep).eq.0) then
	do ikstep=2,rkmethod
c	call phshift (u,cvx(1,1,1,1,1,ikstep),3,ikstep)
	call phshift (u,cvx(1,1,1,1,1,ikstep),3,3,ikstep)

        do i=1,3
	  call kxtran (cvx(1,1,1,i,1,ikstep),cvx(1,1,1,i,1,ikstep),
     1                 cvx(1,1,1,i,1,ikstep),1)
        end do

	end do
	end if
	end if
 
c
c
	if (vorflag) then 
	allocate (uom(nxpad,zisz,yjsz+padx,3),stat=ierr)
        if (ierr.ne.0) call abrt('error allocating uom in timestep_module.f')
	call phase1_om (u,uom,2)
	else
          call phase1 (u,2)
          allocate (ut(nxhpad,zisz,nut),stat=ierr)
          ut(:,:,:)=0.
	end if



          if(rkstep.eq.1)then
! flag for scalar histograms
             schflag=.false.
             if (nc.gt.0) then
          if (jstep.eq.1.and.(kinit(3+nc).ne.0.or.istep.ne.1)) schflag=.true.
                if (ioflag.eq.1.and.mod(iocount,scfreq).eq.0) schflag=.true.
             end if
          endif

	tt(1)=tt(1)+MPI_WTIME()-rtime1
	rtime1=MPI_WTIME()

	if(nzpad*xisz .ne. nxpad*zisz/2) then
	   allocate(ut1(nxpad,zisz,nut))
	   call itransform (u,u,u,ut,ut1,2)
	   deallocate(ut1)
	else


        if (taskid.eq.1) write (6,*) 'before itransform, vorflag=',vorflag
	if (nc.eq.0.and.(rkstep.ne.1.or.ioflag.ne.1)) then
        if (taskid.eq.1.and.istep.le.2) write (6,*) 'call itransform_vel'
   	   call itransform_vel (u,u,u)
	else
	if (vorflag) then
        if (taskid.eq.1.and.istep.le.2) write (6,*) 'call itransform_om'
  	   call itransform_om (u,u,u,uom,uom,uom)
	else
        if (taskid.eq.1.and.istep.le.2) write (6,*) 'call itransform'
	   call itransform (u,u,u,ut,ut,2)
	end if
	end if
	endif
c
          deallocate (ut,stat=ierr)
!
	tt(2)=tt(2)+MPI_WTIME()-rtime1
	rtime1=MPI_WTIME()
!
!#        All Velocity Components and Scalar Fluctuations in Phys Space
c
!
c

	if (nc.eq.0) then
           call realspace_vel (u,2,rkstep)
	else
          call realspace(u,2,rkstep)
	end if

	tt(3)=tt(3)+MPI_WTIME()-rtime1
	rtime1=MPI_WTIME()
!
	if (nc.eq.0) then
c          call transform_vel (u,u,u,2)
          call transform (u,u,u,2)
	else
          call transform (u,u,u,2)
	end if

	tt(4)=tt(4)+MPI_WTIME()-rtime1
	rtime1=MPI_WTIME()
	rtime2=MPI_WTIME()
	rtime3=MPI_WTIME()

         !Now all fields back to fourier 
          if(rkstep.eq.1)then

c
             call step

             call force
          else if(rkstep.eq.rkmethod) then
             if (mod(istep-istep0,iostep).eq.0.and.ioflag.ne.-1) ioflag=1
          endif
      
          !This subroutine performs the last of the 
          !Transforms (y) and then handles the rest
          !of the nonlinear terms (for homogeneous)
ccccc
	tt(7)=tt(7)+MPI_WTIME()-rtime3
	rtime3=MPI_WTIME()

	if (nc.eq.0) then
           call wavespace_vel (u,2)
	else
          call wavespace(u,2)
	end if
	tt(8)=tt(8)+MPI_WTIME()-rtime3
	rtime3=MPI_WTIME()
ccccc

 
 
	if (rkmethod.eq.4) then
c
	if (rkstep.eq.1) then
	call advanc_1 (u,un,u1)
	else if (rkstep.eq.2) then
	call advanc_2 (u,un,u1)
	else if (rkstep.eq.3) then
	call advanc_3 (u,un,u1)
	else if (rkstep.eq.4) then
	call advanc_4 (u,un,u1)
	end if
c
	else
c	

	if (rkstep.eq.1) then
	call proc4a (u,un,rkstep)
c	   call sptr(un,1,1)
	else if (rkstep.eq.2) then
	call proc4b (u,un,rkstep)
	end if
c
	end if


	tt(9)=tt(9)+MPI_WTIME()-rtime3
	rtime3=MPI_WTIME()
c
!  determine if this is the last step
!
!  in the case of shear flows, we should always stop at st=stend,
!  even if it means more than nsteps time-steps
 
	t_sptr=0.
	rtime4=MPI_WTIME()
c
!1st step
          if(rkstep.eq.1)then
             if (shear.eq.0.) then
                if (entime.gt.0..and.
     1          (time.ge.entime*(1.-strof).or.jstep.ge.nsteps)) istop=1
                if (entime.le.0..and.jstep.eq.nsteps) istop=1
             else if (shear.gt.0.) then
                if (stend.gt.0..and.shear*time.ge.stend*(1.-strof)) istop=1
                if (stend.le.0..and.jstep.eq.nsteps) istop=1
             end if
             ! make the last step an output step 
             if (istop.eq.1) then
                iost=iostep
                ioflag=1
             end if
!last substep
          else if(rkstep.eq.rkmethod)then
c#ifndef TIMERS	
 	if (ioflag.eq.1) then
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
	   if(taskid .eq. 0) then
	      print *,'Calling sptr, ioflag=',ioflag
	   endif

	if (taskid.eq.0) write (6,*) 'istep,iocount=',istep,iocount

      call sptr(un,2,1)
      str2 = 'AXI_SPTR'
	if (ioaxi(1).gt.0.and.mod(iocount,ioaxi(1)).eq.0) then
      str1 = 'axi_kx'
      CALL AXI_SPTR(un,1,str1,str2)
	end if
	if (ioaxi(2).gt.0.and.mod(iocount,ioaxi(2)).eq.0) then
      str1 = 'axi_ky'
      CALL AXI_SPTR(un,2,str1,str2)
	end if
	if (ioaxi(3).gt.0.and.mod(iocount,ioaxi(3)).eq.0) then
      str1 = 'axi_kz'
      CALL AXI_SPTR(un,3,str1,str2)
	end if

	   call check_consym2 (un,3,'rk',3)
c

	endif
c#endif	
          endif
c
	t_sptr=t_sptr+MPI_WTIME()-rtime4
c	if (taskid.eq.0) write (6,*) 'istep,t_sptr=',istep,t_sptr

	rtime2=MPI_WTIME()

	tt(5)=tt(5)+MPI_WTIME()-rtime1
	tt(6)=tt(6)+MPI_WTIME()-rtime0
	tt(10)=tt(10)+MPI_WTIME()-rtime2
c

c	if (jstep.ne.1.and.rkstep.eq.rkmethod.and.ioflag.ne.1) then
	if (jstep.ne.1.and.rkstep.eq.rkmethod) then
	if (nc.gt.0.and.(schflag.or.scgflag)) go to 81
c next line changed on 2/2/14
c	if (ioflag.eq.0) then
	if (ioflag.le.0) then
	jj=1
	ncpusteps=ncpusteps+1
	t_rks(:,1)=t_rks(:,1)+tt(:)
	if (taskid.eq.0) write (6,"('jstep, ncpusteps=',3x,2i4,1p,2e11.3)") jstep,ncpusteps,tt(6),t_rks(6,1)
	else
	jj=3
	ncpusteps_io=ncpusteps_io+1
	t_rks(:,2)=t_rks(:,2)+tt(:)
	if (taskid.eq.0) write (6,"('jstep, ncpusteps_io=',2i4,1p,2e11.3)") jstep,ncpusteps_io,tt(6),t_rks(6,2)
	end if

	t_itrans(:,jj)=t_itrans(:,jj)+t_itrans(:,2)
	t_trans(:,jj)=t_trans(:,jj)+t_trans(:,2)
	t_kxcomm1(:,jj)=t_kxcomm1(:,jj)+t_kxcomm1(:,2)
	if (nc.eq.0) then
	t_kxcomm2(:,jj)=t_kxcomm2(:,jj)+t_kxcomm2(:,2)
	else
	t_kxcomm2t(:,jj)=t_kxcomm2t(:,jj)+t_kxcomm2t(:,2)
	end if
	t_xkcomm1(:,jj)=t_xkcomm1(:,jj)+t_xkcomm1(:,2)
	t_xkcomm2(:,jj)=t_xkcomm2(:,jj)+t_xkcomm2(:,2)
 81	continue
	end if

          return
        end subroutine rksubstep

	end module timestep_module
