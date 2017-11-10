












!################################################
!#   Subroutine ITRANSFORM
!#
!#   Performs an inverse transform on X and Z 
!#        (return to physical space)
!#
!#
!################################################
	subroutine itransform(ux,uy,uz,utx,utz,m)
          use comsp
	use timers_rkstep
	use timers_tran
          !#deepcustom#	implicit none
          !#deepcustom# 	implicit none
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
          real(b8) :: ux(nxpad,zisz,yjsz,nu)
          complex(b8) :: uz(nzpad,xisz,yjsz,nu)
          complex(b8) :: uy(nypad,zjsz*xisz,nu)
          complex(b8) bk3i
          complex(b8) :: utz(nzpad,xisz,nut)
          complex(b8), allocatable :: buf(:,:),buf1(:,:,:)
          complex(b8), allocatable :: bk1i(:)
          integer :: m,i,j,zc,iz,int,ifc,ilc,idflag,x2,ix,dnz,iy,l
          real(b8) term1,term2,term3
          real(b8) factor
          real(b8) :: utx(nxpad,zisz,nut)          
          real(b8), allocatable :: scgy(:,:,:)
	  real xnorm,ynorm,znorm,tmpvmaxz,vel
c
	integer nddug
	real(b8), allocatable :: gradhist(:,:,:,:)
	real(b8), allocatable :: gradmom(:,:,:,:)

	real*8 rtime1,rtime2
c
	integer lstep
	real(b8) frac1,frac2
	real(b8), allocatable :: cvxint(:,:,:)

	integer i1
        integer iy1,iy2,ithr,omp_get_thread_num,ix1
c


        if (kstep.eq.1) then 
	t_itrans(:,2)=0.
	t_kxcomm1(:,2)=0.
	t_kxcomm2(:,2)=0.
	t_kxcomm2t(:,2)=0.
	end if

	if (nc.gt.0) then
          allocate (scgy(ncgd,3,0:num_thr-1))
		nddug=2.*gpdflim/0.2+1
		allocate (gradhist(nddug,nc,2,0:num_thr-1))
		allocate (gradmom(6,nc,2,0:num_thr-1))
	end if
          
	rtime1=MPI_WTIME()
c
        rtime2=MPI_WTIME()

        call divide_thr (xisz,ixp1,ixpi,'itransform: xisz')
        call divide_thr (zisz,izp1,izpi,'itransform: zisz')
c
       
          idflag=0
          int=nc+ncps
          if  (nc.eq.0.and.kstep.eq.1.and.ioflag.eq.1) then
             idflag=1
             int=2
          end if
	if (idflag.eq.1) then 
	call vgmomt (utx,0,0,0)
	end if

          allocate (bk1i(xisz))

          do xp=1,xisz
             x=xist+xp-1
             bk1i(xp)=imagi*bk1(x,m)
          end do

 	if (kstep.eq.1.and.(mod(jstep-1,iostep).eq.0.or.ioflag.eq.1)) then
	if (nc.gt.0) then
          scgy(:,:,:)=0.
	end if
	end if

        rtime2=MPI_WTIME()

! Transform in Y, transpose to get data in z-pencils, add extra elements if needed
	  if(jproc .gt. 1) then
	  

! Watch for array overlap issues ! 

	i1=4
	     if (mod(jstep-1,ivstep).eq.0.or.
     1        (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) i1=1
c
c     do i=i1,3+2*nc
		call kxcomm1_trans_cyl (uy(1,1,1),uz(1,1,1,1),3+2*nc)
c		call kxcomm1_trans_cyl (uy(1,1,i),uz(1,1,1,i),1)
c	     end do

	  else
	     
	     call i_reorder_yz_trans_y(uy,uz,3+2*nc)
	     
	  endif
	  

        rtime2=MPI_WTIME()-rtime2
        t_itrans(1,2)=t_itrans(1,2)+rtime2


          ! fft preliminaries

	allocate (cvxint(nx,zisz,3))
c
	  do 10 yp=1,yjsz

        ithr=0
!$OMP PARALLEL private (ithr,ix1,bk3i,i1)


             


             !npm converting loops so I can read this block
             if (nc.gt.0) then
!$OMP MASTER
	rtime2=MPI_WTIME()
!$OMP END MASTER

                do i=1,nc
                   j=i+3
                   do xp=ixp1(ithr),ixp1(ithr)+ixpi(ithr)-1
                      do z=1,nzpad
                         if (z.ne.nzhppad) then
			    bk3i=imagi*bk3(z,m)
			    utz(z,xp,i)=bk3i*uz(z,xp,yp,j)
			    utz(z,xp,i+nc)=uz(z,xp,yp,j)
			    uz(z,xp,yp,j)=bk1i(xp)*uz(z,xp,yp,j)
                         endif
                      enddo
                   enddo
                enddo
!$OMP MASTER
	t_itrans(4,2)=t_itrans(4,2)+MPI_WTIME()-rtime2
!$OMP END MASTER
             else !(nc.lt.0)
                if (idflag.eq.1) then
                   do xp=ixp1(ithr),ixp1(ithr)+ixpi(ithr)-1
                      do z=1,nzpad
                         if (z.ne.nzhppad) then
                            bk3i=imagi*bk3(z,m)
                            utz(z,xp,1)=bk1i(xp)*uz(z,xp,yp,1)
                            utz(z,xp,2)=bk3i*uz(z,xp,yp,3)
                         endif
                      enddo
                   enddo
                end if !(idflag.eq.1)
             end if

             !c
             !c set highest z-wavenumber modes to zero
             !c then transform to physical space (in z)
             !c for both uz and utz
             !c

	i1=4
	     if (mod(jstep-1,ivstep).eq.0.or.
     1        (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) i1=1
c
        ix1=ixp1(ithr)
c
!$OMP MASTER
	rtime2=MPI_WTIME()
!$OMP END MASTER
	     do i=i1,3+2*nc
       call fftw_execute(plan2_p2a(ithr),uz(1,ix1,yp,i),uz(1,ix1,yp,i))
	     enddo
!$OMP MASTER
	t_itrans(2,2)=t_itrans(2,2)+MPI_WTIME()-rtime2
!$OMP END MASTER

!$OMP END PARALLEL


	rtime2=MPI_WTIME()

	     if(iproc .eq. 1) then
                
! In case of 1D decomposition do only FFT in Z without data exchange,
! interchange indices, pad arrays
		
		do i=1,3+2*nc
		   call i_reorder_xz_trans_x(uz(1,1,yp,i),ux(1,1,yp,i))
                enddo

	     else
! if iproc > 1: Transpose X <-> Z and transform C2R in X
! Do a multivar transform in this case no matter if the directive 
! is defined, since this in this case only one Y plane is transformed
! so memory usage is not that high, and also in this way we avoid 
! problems with uz and ux array overlap in case of uneven division



	     if (mod(jstep-1,ivstep).eq.0.or.
     1        (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) then
               do i=1,3
                 call kxcomm2t_trans (uz(1,1,1,i),ux(1,1,1,i),yjsz,yp,1)
               end do
             endif
             do i=1,2*nc
      	       call kxcomm2t_trans (uz(1,1,1,i+3),ux(1,1,1,i+3),yjsz,yp,1)
             end do

	     endif

	t_itrans(3,2)=t_itrans(3,2)+MPI_WTIME()-rtime2
	rtime2=MPI_WTIME()

c       call MPI_FINALIZE (MPI_COMM_WORLD,mpierr)
c       stop 'force stop in itransform'

        ithr=0
!$OMP PARALLEL private (ithr,ix1)
c

        ix1=ixp1(ithr)

	     do i=1,int

!c In case of 1D decomposition do only FFT in Z without data exchange

	do ix=ixp1(ithr),ixp1(ithr)+ixpi(ithr)-1
	utz(nzhppad,ix,i)=cmplx(0.,0.)
	end do
		
		call fftw_execute(plan2_p2a(ithr),utz(1,ix1,i),utz(1,ix1,i))
	     end do

!$OMP END PARALLEL

	t_itrans(2,2)=t_itrans(2,2)+MPI_WTIME()-rtime2
	rtime2=MPI_WTIME()

	  if(iproc .eq. 1) then

	     do i=1,int
		call i_reorder_xz_trans_x(utz(1,1,i),utx(1,1,i))
	     enddo
	  else

        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	     call kxcomm2t_trans (utz,utx,1,1,int)
	  endif

	t_itrans(3,2)=t_itrans(3,2)+MPI_WTIME()-rtime2
	rtime2=MPI_WTIME()

!c at this point, all velocity components and scalar
!c fluctuations are in physical space
!c
        ithr=0

!$OMP PARALLEL private (ithr,term1,term2,term3,lstep,frac1,frac2,ifc,ilc,scgflag)


	if (idflag.eq.1) call vgmomt (utx,yp,1,ithr)
!c
	if (ivstep.gt.1) then
c
c (see sub. rksubstep for "cvx" array as well)
c
c
	if (kstep.eq.1) then

		if (mod(jstep-1,ivstep).eq.0) then
		do i=1,3
       		do zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
		do x=1,nxpad
		cvx(x,zp,yp,i,1,kstep)=ux(x,zp,yp,i)
		end do
		end do
		end do
		end if

		if (mod(jstep-2,ivstep).eq.0) then
		do i=1,3
       		do zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
		do x=1,nxpad
		cvx(x,zp,yp,i,2,kstep)=ux(x,zp,yp,i)
		end do
		end do
		end do
		end if
c
	else

		if (mod(jstep-1,ivstep).eq.0) then
		do i=1,3
       		do zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
		do x=1,nxpad
		cvx(x,zp,yp,i,2,kstep)=ux(x,zp,yp,i)
		end do
		end do
		end do
		end if
c
	end if

	end if


c some lines below (affecting output files called scgrad, scgpll*, scgpp* and 
c scdlpdf*) were fixed by P.K Yeung, 11/27/09

 	if (kstep.eq.1.and.(mod(jstep-1,iostep).eq.0.or.ioflag.eq.1)) then
	scgflag=.true.
	if (nc.gt.0) call scgrad (utx,ux,yp,scgy,ithr)
	end if
c
	do i=1,nc
      if (schflag.and.kstep.eq.1) then
	ifc=0 ; ilc=0
	if (yp.eq.1) ifc=1
	if (yp.eq.yjsz) ilc=1
        if (schflag) then
 	call scgpdf (ifc,ilc,i,yp,ux(1,1,yp,i+3),ux(1,1,yp,i+3+nc), utx(1,1,i),
     1              ithr,nddug,gradhist(1,1,1,ithr),gradmom(1,1,1,ithr),
     1              gradmom(1,1,2,ithr))
	end if
	end if
	end do
c
!$OMP BARRIER


	if (ivstep.eq.1) then
c
	do 70 i=1,nc
        do 71 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
	do 72 x=1,nxpad
	term1=ux(x,zp,yp,1)*(cb(1,i)+ux(x,zp,yp,i+3))
	term2=ux(x,zp,yp,2)*(cb(2,i)+ux(x,zp,yp,i+3+nc))
	term3=ux(x,zp,yp,3)*(cb(3,i)+utx(x,zp,i))
	ux(x,zp,yp,i+3)=term1+term2+term3
	ux(x,zp,yp,i+3+nc)=utx(x,zp,i+nc)
 72	continue
 71	continue
 70	continue
c
	else
c
	lstep=jstep-(jstep/ivstep)*ivstep-1
	if (lstep.lt.0) lstep=lstep+ivstep
	if (kstep.eq.2) lstep=lstep+1
	frac2=float(lstep)/ivstep
	frac1=1.-frac2
c
	do i=1,3
        do zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
	do x=1,nxpad
	cvxint(x,zp,i)=frac1*cvx(x,zp,yp,i,1,kstep)
     1	              +frac2*cvx(x,zp,yp,i,2,kstep)
	end do
	end do
	end do
	
c	cvxint(:,:,:)=frac1*cvx(:,:,yp,:,1,kstep)+frac2*cvx(:,:,yp,:,2,kstep)
c         write (140+taskid,911) istep,kstep,yp,cvxint(4,4,2),ux(4,4,yp,3+nc), lstep
 911	format ('istep,kstep,yp,cvx=',3i3,1p,2e12.4,'  lstep=',i2)
	
	do i=1,nc
        do zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
	do x=1,nxpad
	term1=cvxint(x,zp,1)*(cb(1,i)+ux(x,zp,yp,i+3))
	term2=cvxint(x,zp,2)*(cb(2,i)+ux(x,zp,yp,i+3+nc))
	term3=cvxint(x,zp,3)*(cb(3,i)+utx(x,zp,i))
	ux(x,zp,yp,i+3)=term1+term2+term3
	ux(x,zp,yp,i+3+nc)=utx(x,zp,i+nc)
	end do
	end do
	end do
c
	end if


!c


!$OMP END PARALLEL

	t_itrans(4,2)=t_itrans(4,2)+MPI_WTIME()-rtime2
	rtime2=MPI_WTIME()

 10	continue
c
	deallocate (cvxint)

              
	if (idflag.eq.1) call vgmomt (utx,0,2,0)

      if (kstep.eq.1.and.(mod(jstep-1,iostep).eq.0.or.ioflag.eq.1)) then
         if (nc.gt.0) call scgrad (utx,ux,0,scgy,0)
         if (schflag) call scgpdf_lc (nddug,gradhist,gradmom)
      end if
	deallocate (bk1i)

	if (nc.gt.0) deallocate (scgy,gradhist,gradmom)
          
c
	rtime1=MPI_WTIME()-rtime1

	return
        end subroutine itransform

!-------------------------------------------------------	
	subroutine i_reorder_yz_trans_y(uy,uz,int)

	use comp
	!#deepcustom#	implicit none
	!#deepcustom# 	implicit none

	integer i,dnz,iz,iy,int
	complex(b8) uz(nzpad,xisz,yjsz,int)
	complex(b8) uy(ny,zjsz*xisz,int)



	return
	end subroutine

!-------------------------------------------------------	
	subroutine i_reorder_xz_trans_x(uz,ux)

	use comp
	!#deepcustom#	implicit none
	!#deepcustom# 	implicit none

	complex(b8) uz(nzpad,xisz)
	real(b8) ux(nxpad,zisz)
	integer iz,ix,x2
	complex(b8), allocatable :: buf(:,:)
c


	allocate(buf(nxhppad,zisz))
	buf = 0

! Transpose for buf x <-> Z and pad in X, zero highest mode

	do z=1,nzpad,NB1_Z
	   z2 = min(z+NB1_Z-1,real(nzpad))
	   
	   do x=1,nxhp,NB1_X
	      x2 = min(x+NB1_X-1,real(xisz))
	      
	      do iz=z,z2
		 do ix=x,x2
		    buf(ix,iz) = uz(iz,ix)
		 enddo
	      enddo
	   enddo
	enddo
	   
	do z=1,nzpad
	   do x=nxhp,nxhppad
	      buf(x,z) = 0.
	   enddo
	enddo

!       transform C2R in X

	call fftw_execute_c2r(plan3_p2a,buf,ux)

	deallocate(buf)


	return
	end subroutine
