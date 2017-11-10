












!######################################################################
!#       Transform Routines
!# Takes as input the velocity fields
!# Transforms to fourier space
!#
!# Only transforming / transposing two components of velocity (x,z)
!#
!######################################################################

       subroutine transform(ux,uy,uz,m) 
	use comp
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
	real(b8) factor
	complex(b8), allocatable :: buf(:,:),buf1(:,:,:)
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)
 	complex(b8) :: uz(nzpad,xisz,yjsz,nu)
 	complex(b8) :: uy(nypad,zjsz*xisz,nu)
	integer :: m,i, is1,is2,novt,rkstep,i2f,x2,iz,ix,iz2,dnz,iy,l
        real(b8) s1,s2,s3
        complex(b8), allocatable :: bk1i(:,:),bk3i(:,:)

c 8/21/2015: fix for do 75 in ROTD mode which caused incorrect
c results on non-2pi^3 domains
c
	integer ithr,omp_get_thread_num,yp1,yp2

	integer i1,ii
	real(8) rtime1

        if (taskid.eq.0.and.istep.eq.1) write (6,*) 'call transform, kstep=',kstep
        if (kstep.eq.1) then
	t_trans(:,2)=0.
        t_xkcomm1(:,2)=0.
        t_xkcomm2(:,2)=0.
	end if

        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
        is1=4+nc
        is2=5+nc
!factor=1./nx
	factor=1./real(nxpad,b8)

! transform to x-wavenumber space, and take transpose 
!
	novt=5+nc

	rtime1=MPI_WTIME()
c
! Transform R2C in X, transpose x <-> Z, truncate in X if needed
        if (iproc.gt.1) then

! Watch for array overlap in case of uneven division

	i1=1
	if (mod(jstep-1,ivstep).ne.0) i1=4
           do i=i1,novt
              call xkcomm1_trans(ux(1,1,1,i),uz(1,1,1,i),1)
           enddo
        else
           call f_reorder_xz_trans_xz(ux,uz,novt)
        endif

        t_trans(1,2)=t_trans(1,2)+(MPI_WTIME()-rtime1)

	if (iproc .gt. 1) then
!
! transform in z
!
	i1=1
	if (mod(jstep-1,ivstep).ne.0) i1=4
c
	call divide_thr (yjsz,iyp1,iypi,'transform: yjsz')
c
	rtime1=MPI_WTIME()

	ithr=0

!$OMP PARALLEL private (ithr,yp1,yp2)
c
c
           do 60 i=i1,novt
c
        yp1=iyp1(ithr)

        call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,i),uz(1,1,yp1,i))
 60           continue

!$OMP END PARALLEL

        t_trans(2,2)=t_trans(2,2)+(MPI_WTIME()-rtime1)
        rtime1=MPI_WTIME()
c
 	endif


!#################################################
! operations for convective terms
!#################################################

	if (mod(jstep-1,ivstep).ne.0) go to 77
c
        allocate (bk1i(xist:xien,2))
        allocate (bk3i(nz,2))
c
	ithr=0
!$OMP PARALLEL private (ithr,y,x,bk1i,bk3i,yp1,yp2)
c
c
!,kstep these new arrays introduced by PKY, 4/9/2012


        do x=xist,xien
        bk1i(x,1)=imagi*kx(x)
        bk1i(x,2)=imagi*bk1(x,m)
        end do
c
        do z=1,nz
        bk3i(z,1)=imagi*bk3(z,m)
        bk3i(z,2)=imagi*kz(z)
        end do

	yp1=ithr*yjsz/num_thr+1
	yp2=yp1+yjsz/num_thr-1
c
c       do 70 yp=yp1,yp2
        do 70 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
           y=yjst+yp-1
           
           do 71 xp=1,xisz
              x=xist+xp-1

              do 76 z=1,nzpad
                 if (z.eq.nzhppad) go to 76
                 uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,is2)
                 uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,is2)+bk3i(z,2)*uz(z,xp,yp,3)
 76           continue

 71        continue

 70   continue


!$OMP END PARALLEL
	deallocate (bk1i,bk3i)

        t_trans(3,2)=t_trans(3,2)+(MPI_WTIME()-rtime1)
        rtime1=MPI_WTIME()

 77	continue 
!###################################################
        
      i2f=4+nc
      

! Truncate in Z if needed, Transpose ZXY -> YZX (Z- to Y-pencils), transform in Y 
        rtime1=MPI_WTIME()

      if(jproc .gt. 1) then
         if(num_al_i(0) .gt. 0) then
            do i=1,i2f
               call xkcomm2_trans_cyl (uz(1,1,1,i),uy(1,1,i),1)
            end do
         endif

      else
! 1D with jproc=1

         call f_reorder_yz_trans_y(uz,uy,i2f)

      endif

        t_trans(4,2)=t_trans(4,2)+(MPI_WTIME()-rtime1)


      return
      end subroutine transform

!-----------------------------------------------------------

      subroutine f_reorder_xz_trans_xz(ux,uz,int)

      use comp
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none

      complex(b8) uz(nzpad,xisz,yjsz,int)
      real(b8) ux(nxpad,zisz,yjsz,int)
      integer i,iz,ix,x2,int
      complex(b8), allocatable :: buf(:,:)
      real(b8) factor

      factor=1./real(nxpad,b8)

      allocate(buf(nxhppad,zisz))
      buf = 0

      do i=1,int
         do yp=1,yjsz
                 
            call sfftw_execute_dft_r2c(plan1_p2b,ux(1,1,yp,i),buf)


            do z=1,nzpad,NB1_Z
               z2 = min(z+NB1_Z-1,real(nzpad))
               
               do x=1,nxh,NB1_X
                  x2 = min(x+NB1_X-1,real(nxh))
                  
                  do iz=z,z2
                     do ix=x,x2
                        uz(iz,ix,yp,i) = buf(ix,iz) * factor
                     enddo
                  enddo
               enddo
            enddo
            
            call fftw_execute(plan2_p2b,uz(1,1,yp,i),uz(1,1,yp,i))
         enddo

      enddo

      deallocate(buf)

      return
      end subroutine

      
      subroutine f_reorder_yz_trans_y(uz,uy,int)

      use comp
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none

      integer i,iz,iy,dnz,int
      complex(b8) uz(nzpad,xisz,yjsz,int)

      complex(b8) uy(ny,zjsz*xisz,int)


      
      
      return
      end subroutine
