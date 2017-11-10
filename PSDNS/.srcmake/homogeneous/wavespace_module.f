












! Version with Cylindrical representation 
!
!################################################
!#      module wavespace_module
!#
!#  contains routines related to operations in wavespace
!#  
!#  NOTE: also contains routines for homogeneous transforms
!#  I.E. the Y inverse and forward transform
!# 
!################################################
      module wavespace_module
      use comp
      !#deepcustom#	implicit none
  

      public  :: wavespace_sub, initialize_wavespace, yitransform_sub
!      private :: 

!############################################
        contains
!###########################################
!#     yitransform_module
!#
!#   transform the y-component 
!#   to physical space
!#
!#   Built on proc1.f
!############################################
          subroutine yitransform_sub(uy,m)
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

      allocate (shift(nypad))

c


      xp = mystart_x
      z = mystart_z
      do a=1,num_al

         x = xp + xist-1
         
         sxz=sz(z,kstep)*sx(x,kstep)

	if (mod(jstep-1,ivstep).eq.0.or.
     1      (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) then
         do 15 y=1,nypad
c      if (y.eq.nyhppad) go to 15
            if(mask(y,a)) then
               shift(y) = sxz *sy(y,kstep)

               uy(y,a,1)=shift(y)*uy(y,a,1)
               uy(y,a,2)=shift(y)*uy(y,a,2)
               uy(y,a,3)=shift(y)*uy(y,a,3)

            else
               uy(y,a,1:3) = 0.0
            endif
 15      continue
	end if
!
! shift scalars and form y gradients
!
         do 30 j=4,3+nc
            do 30 y=1,nypad
c              if (y.eq.nyhppad) go to 30
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
      enddo

!
c Note by PKY: in cylindrical truncation version the actual
c y-transform from wavenumber space to physical space is
c now absorbed into kxcomm1_trans_cyl (called from itransform.f)
c
! transform in y-direction to physical space
!
!      do 60 i=1,3+2*nc
!
!      do 20 xp=1,xisz
!      do 20 zp=1,zjsz
!      z=zjst+zp-1
!c      if (z.eq.nzhp) go to 20
!	
!      uy(nyhp,a,i)=0.
!
! 20   continue
!
!#ifdef 1
!      call fftw_execute(plan2_p1,uy(1,1,1,i),uy(1,1,1,i))
!#else
!      call escfft (uy(1,1,1,i),1,xisz*zjsz,ny,xisz*zjsz,1)
!#endif
!
! 60   continue

	deallocate(shift)
c
        return 
      end subroutine yitransform_sub

!############################################
!#           Wavespace
!#
!#     Program that operates when in wavespace
!#     Therefore, for homogeneous turbulence
!#     Program will transform y to finish forming 
!#     Nonlinear terms              
!#
!#     Built on proc3.f
!############################################
        subroutine wavespace_sub(uy,m)
          !#deepcustom#	implicit none
!          use comp
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
          integer :: i,j,m, is1, i2f,a,xz(2)
          complex(b8) :: uy(ny,zjsz*xisz,nu)
          complex(b8), allocatable :: shift(:)
          complex(b8) :: utmp1
          complex(b8) :: syz,sxz,sxyz
          real(b8) norm,xnorm,ynorm,znorm,vm,vel

        integer ithr,omp_get_thread_num
	integer ii

        
	i2f=4+nc

 
c	do i=1,i2f
c	   call xkcomm2_trans (uy(1,1,1,i),uy(1,1,1,i))
c	end do


!#############################################################
!# transpose to y-pencils and take inverse transform 
 

c	allocate (shift(ny))
 
      norm=dt2/ny/nz
 
      is1=4+nc

        ithr=0
!$OMP PARALLEL private (ithr,x,xp,z,sxz,sxyz,utmp1)


      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1


         x = xp + xist-1

         
 
            sxz=sz(z,kstep)*sx(x,kstep)


! form shift, shift, and complete formation of convective terms
 
 	if (mod(jstep-1,ivstep).ne.0) go to 47
            do 42 y=1,ny
               if( mask(y,a)) then

                  sxyz=-norm*conjg(sxz*sy(y,kstep))
c help get predictor estimate at 'ivstep's forward
	sxyz=sxyz*ivstep
 
                  utmp1=imagi*b22(m)*ky(y)*sxyz
c
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
 47	continue
 


! shift the (convective terms of) scalars
 
            do 50 j=4,3+nc
               do 50 y=1,ny
c ?            if (y.eq.nyhp) go to 50
                  if( mask(y,a)) then
                     sxyz=-norm*conjg(sxz*sy(y,kstep))
                     
                     uy(y,a,j)=sxyz*uy(y,a,j)
                  endif
 50            continue


      call next_xz(xp,z)

 100	continue

!$OMP END PARALLEL




      return
      end subroutine wavespace_sub
!############################################
!#           Initialize Wavespace
!#
!#     Program that sets up all variables needed
!#     By the wavespace module
!#
!############################################
        subroutine initialize_wavespace
          !#deepcustom#	implicit none
      

        end subroutine initialize_wavespace
!############################################
        end module wavespace_module
