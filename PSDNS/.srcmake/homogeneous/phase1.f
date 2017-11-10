












          subroutine phase1 (uy,m)
c
            use comp
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

            complex(b8) :: uy(nypad,zjsz*xisz,nu)
            integer :: i,j,m
            complex(b8), allocatable :: shift(:)
            complex(b8) :: syz,sxz
            integer iy,a,xz(2)
            real(b8) factor,rk1,rk1sq

        integer ithr,omp_get_thread_num,ia

      allocate (shift(nypad))
c
	if (mod(jstep-1,ivstep).eq.0.or.
     1      (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) then
c	call phshift (uy(1,1,1),uy(1,1,1),nu,3,kstep)
 	call phshift_inplace (uy(1,1,1),3,kstep)
	end if

!$OMP PARALLEL private (ithr,sxz,bk2i,x,xp,z)

        ithr=0


      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 10 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c

         x = xp + xist-1
         
         sxz=sz(z,kstep)*sx(x,kstep)

! shift scalars and form y gradients
!
         do 30 j=4,3+nc
            do 30 y=1,nypad
               if(mask(y,a)) then
               shift(y) = sxz *sy(y,kstep)
                  bk2i=imagi*b22(m)*ky(y)
                  
                  uy(y,a,j)    = shift(y) * uy(y,a,j)
                  uy(y,a,j+nc) =  bk2i * uy(y,a,j)
               else
                  uy(y,a,j) = 0.0
                  uy(y,a,j+nc) = 0.0
               endif
 30         continue

            call next_xz(xp,z)
 10	continue

!$OMP END PARALLEL
!
c Note by PKY: in cylindrical truncation version the actual
c y-transform from wavenumber space to physical space is
c now absorbed into kxcomm1_trans_cyl (called from itransform.f)
c
! transform in y-direction to physical space
!
	deallocate(shift)
c
        return 
      end subroutine phase1
