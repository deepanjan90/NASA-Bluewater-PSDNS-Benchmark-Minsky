












!  3D FFT inverse transform with 2D domain decomposition
!
!  This version uses MPI_Alltoallv to exchange data between processors
!  In the second step, y planes are sent separately
!  The order of array elements in memory is the same in all stages: (x,z,y) 

! Input: XZYg - comlpex array, with y dimension contained entirely,
!               while x and z are block-distributed among processors in 2D grid
! XZgY - complex array (auxiliary), z dimension contained entirely
!               while x and y are block-distributed among processors in 2D grid
! Output: XgZY - an array of real, x dimension is contained entirely within processors memory  
!               while z and y are block-distributed among processors in 2D grid

! !!! CAUTION: In this version: all arrays occupy the same memory space
!

      subroutine kxtran (XgZY,XZgY,XZYg,nv)
      use com
      !#deepcustom#	implicit none
	!#deepcustom# 	implicit none
!
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

      complex(b8) XgZY(nxpad/2,zisz,yjsz,nv)
      complex(b8) XZgY(nz,xisz,yjsz,nv)
      complex(b8) XZYg(ny,xisz,zjsz,nv)

      real(b8) factor
	integer :: i,nv
c
        integer iy1,iy2,ithr,omp_get_thread_num

c
! Transform in y dimension for all x and z
!
	XZYg(nyhp,:,:,1:nv)=0.

c transform in y direction and switch into z-lines
c
         call kxcomm1_trans_cyl(XZYg,XZgY,nv)

	call divide_thr (yjsz,iyp1,iypi,'kxtran: yjsz')

! Transform in z dimension for all x, one y-plane at a time
c using 1 or all y-planes in one call (ESSL)
!
c
	XZgY(nzhp,:,:,1:nv)=cmplx(0.,0.)
c
        ithr=0
!$OMP PARALLEL private(ithr,iy1)

c
	do 10 i=1,nv

        iy1=iyp1(ithr)
c
c
            call fftw_execute(plan2_kx(ithr),XZgY(1,1,iy1,i),XZgY(1,1,iy1,i))
c
 10	continue


!$OMP END PARALLEL
c
c
c switch into x-lines and take complex-real transform in x
c
         call kxcomm2_trans (XZgY(1,1,1,1),XgZY(1,1,1,1),nv)
c
c
      return
      end
