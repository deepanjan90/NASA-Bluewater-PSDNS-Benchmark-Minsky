












! Input: X-pencils (real)
! Output: Z-pencils (complex)
!
! This routine performs real-to-complex FFT,
! truncates the X dimension from nxpad to nx, and
! transposes the data to arrange Z-pencils while 
! interchanging the order of X and Z indices (using loop blocking)
!
! Multivariable version
!

      subroutine xkcomm1_trans(source,dest,nv)

      use com
      use timers_comm
      use timers_comp
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

      real(b8) source(nxpad,zisz,yjsz,nv)
      complex(b8) dest(nzpad,xisz,yjsz,nv)
      complex(b8), allocatable :: buf(:,:,:) 

      integer position,pos0,pos1,pos2
      integer, allocatable :: pos(:)
      integer i,j,k,n,x2,ix,iz,ii,nv
      integer :: iarg
      real(b8) factor
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm
c
	real*8 rtime0,rtime1,rtime2
c
      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

c
	integer ithr,omp_get_thread_num,iy1,yp1,yp2,i1,i2,iz1,iz2
	integer, allocatable :: iip1(:),iipi(:)
c
c next line moved to sub. transform, 12/24/2012
c	if (kstep.eq.1) t_xkcomm1(:,2)=0.

      rtime0=MPI_WTIME()

      if(IfCntUneven) then
         allocate (pos(0:iproc-1))
      endif
      iarg = IfCntMax/(b8*2)
      allocate (sendbuf(iarg,nv,0:iproc-1))
      allocate (recvbuf(iarg,nv,0:iproc-1))

      allocate(buf(nxhppad,zisz,yjsz))
      factor=1./real(nxpad,b8)

      tp1_comm = tp1_comm - MPI_Wtime() 
c
	call divide_thr (yjsz,iyp1,iypi,'xkcomm1: yjsz')
	if (iproc.ge.num_thr) then
	allocate (iip1(0:num_thr-1),iipi(0:num_thr-1))
	call divide_thr (iproc,iip1,iipi,'xkcomm1: iproc')
	iip1(:)=iip1(:)-1
	end if
c
        t_xkcomm1(3,2)=t_xkcomm1(3,2)+(MPI_WTIME()-rtime0)

	ithr=0
!$OMP PARALLEL private (ithr,iy1,ii,pos0,z1,pos1,x2,pos2,position,n,i,j,k,iz,yp1,yp2,i1,i2,iz1,iz2,z2)
	

      do 100 j=1,nv

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER

! Transform R2C

	iy1=iyp1(ithr)
c
         call fftw_execute_r2c(plan1_p2b(ithr),source(1,1,iy1,j),buf(1,1,iy1))
c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_xkcomm1(2,2)=t_xkcomm1(2,2)+(rtime2-rtime1)
!$OMP END MASTER

c

c Pack the send buffer for exchanging z and x (within a given y plane ) into sendbuf
         
!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
c

! If using MPI_Alltoall

!$OMP DO
      do yp=1,yjsz

          do i=0,iproc-1
             ii = mymap(i)
             position = (yp-1)*iisz(ii)*zisz+1
             do z=1,zisz
! Pack only those components with kx < nxh
                do x=iist(ii),iien(ii)
                   sendbuf(position,j,i) = buf(x,z,yp) * factor
                position = position +1
                enddo
             enddo
             if(IfCntUneven) then
                pos(i) = position
             endif
          enddo

	end do
!$OMP END DO

c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_other=tcpu_other+(rtime2-rtime1)
        t_xkcomm1(3,2)=t_xkcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER
c
 100	continue
c
!$OMP BARRIER 
! Not needed if sendbuf is preallocated
       if(IfCntUneven) then
!$OMP DO
          do i=0,iproc-1
             ii = mymap(i)
             position = pos(i)
             n = IfCntMax/(b8*2) - iisz(ii)*zisz*yjsz
             do j=1,nv
             do k=1,n
                sendbuf(position,j,i) = 0.0
                position = position +1
             enddo
          enddo
          enddo
!$OMP END DO
       endif

!$OMP BARRIER

!$OMP MASTER


      tp1_comm = tp1_comm + MPI_Wtime() 
      ip_comm = ip_comm + 1

c Exchange the z-x buffers


c
      t_alltoall = t_alltoall - MPI_Wtime()
      t1_comm = t1_comm - MPI_Wtime() 

	rtime1=MPI_WTIME()

      iarg = IfCntMax*nv
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_row)

	rtime2=MPI_WTIME()

      t1_comm = t1_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
      t_xkcomm1(1,2) = t_xkcomm1(1,2) + rtime2-rtime1


!$OMP END MASTER


!Next line (the barrier) appears necessary to enure correct results!
!$OMP BARRIER

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
c

! Unpack the receive buffer, interchanging indices, with loop blocking NB1

      do 200 j=1,nv
!$OMP DO
      do yp=1,yjsz

         do i =0,iproc-1
            ii = mymap(i)
            pos0 =(yp-1)*kisz(ii)*xisz

            do z=kist(ii),kien(ii),NB1_Z
               z2 = min(z+NB1_Z-1,real(kien(ii))) !#deepcustom#
               pos1 = pos0 + (z-kist(ii))*xisz                

               do x=1,xisz,NB1_X
                  x2 = min(x+NB1_X-1,real(xisz)) !#deepcustom#
                  
                  pos2 = pos1 +x
                  do iz=z,z2
                     position = pos2
                     do ix=x,x2
                        dest(iz,ix,yp,j) = recvbuf(position,j,i)
                        position = position +1
                     enddo
                     pos2 = pos2 + xisz
                  enddo
               enddo
            enddo
         enddo
      enddo

!$OMP END DO
 200	continue
c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_other=tcpu_other+(rtime2-rtime1)
      t_xkcomm1(3,2) = t_xkcomm1(3,2) + rtime2-rtime1
!$OMP END MASTER
c
!$OMP END PARALLEL

	rtime2=MPI_WTIME()

      deallocate(buf)
c
      deallocate(recvbuf)
      deallocate (sendbuf)

      if(IfCntUneven) deallocate (pos)

	if (allocated(iipi)) deallocate (iip1,iipi)
c
      t_xkcomm1(3,2) = t_xkcomm1(3,2) + (MPI_WTIME()-rtime2)
c
      i1_comm = i1_comm + 1
c

      return
      end
