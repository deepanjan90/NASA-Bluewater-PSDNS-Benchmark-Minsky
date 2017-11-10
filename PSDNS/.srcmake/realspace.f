












!#######################################################
!#     Nonlinear
!#
!#     Program to handle all nonlinear terms
!#
!#######################################################
	subroutine realspace(ux,m,rkstep)
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
!	real(b8) factor
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)

        real(b8), allocatable :: bv2(:,:)
	integer :: m,i, is1,is2,novt,rkstep,icall
        real(b8) s1,s2,s3,xnorm,ynorm,znorm,vel,k3z,bk3z
	real tmpu,tmpv,tmpw
c
        integer ithr,OMP_GET_THREAD_NUM,yp1,yp2
	real(b8), allocatable :: vmax(:)
	real sum,sum0

	integer ii
c
	kstep=rkstep

          !1st Step
          if(rkstep.eq.1)then
             if (taskid.eq.0.and.schflag)
     1       write (6,*) 'scalar PDFs generated at step=',istep
!
! extract PDFs
!
             do i=1,nc
                if (schflag) then
                   call phymom (u,i+3+nc,i+3+nc)
                   call scnpdf (u(1,1,1,i+3+nc),i)
                end if
             end do 

        if (mod(jstep-1,ivstep).eq.0.or.
     1    (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) go to 100
        return
  100  continue


          endif !done with first step only execution

      is1=4+nc
      is2=5+nc
!
      s1=sqrt(b11(m))
      s2=sqrt(b22(m))
      s3=sqrt(b33(m))
!
      xnorm=b11(m)*nxpad
      ynorm=b22(m)*nypad
      znorm=b33(m)*nzpad
!
      allocate (bv2(nxpad,0:num_thr-1))
!
	if (rkstep.eq.1) allocate (vmax(0:num_thr-1))

c
	call divide_thr (yjsz,iyp1,iypi,'yjsz: realspace')
c
! Courant number and convective terms in physical space
!
      if (rkstep.eq.1) then

	ithr=0
!$OMP PARALLEL private (ithr,vel,velmax,tmpu,tmpv,tmpw)
	velmax=0.
c
c     do 120 yp=1,yjsz
      do 120 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 120 zp=1,zisz

      do 125 x=1,nxpad
!
!
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)

 125   continue
 120   continue
	vmax(ithr)=velmax
!$OMP END PARALLEL
c
	velmax=0.
	do ithr=0,num_thr-1
	velmax=max(velmax,vmax(ithr))
	end do
c
	deallocate (vmax)
c  	print *,taskid,': realspace velmax=',velmax
c
      else if (kstep.ge.2) then
!
c
	ithr=0
!$OMP PARALLEL private (ithr)
c
c     do 130 yp=1,yjsz
      do 130 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 130 zp=1,zisz
!
      do 135 x=1,nxpad

!
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
!
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)
!
 135   continue
 130   continue
!$OMP END PARALLEL
!
      end if

!################
!#                   Done with realspace, now fft forward
!#################
      deallocate (bv2)
c
	return
	end
