












! Update of kxcomm2, by PKY Jan 14, 2012, DP Jan 13, 2012
! all y-planes, but 1 alltoallv call per variable  being transposed

!
! This routine reorders X,Z indices while packing the send buffer
! (using loop blocking),
! exchanges data to arrange X pencils, expandsX dimension to nxhppad
! and performs inverse complex-to-real FFT
! Input: Z pencils, complex
! Output: X-pencils, real 
!
! Multivariable version
!

      subroutine kxcomm2_trans(source,dest,nv)

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

c        complex(b8) dest(nxpad/2,zisz,yjsz,nv)
        real(b8) dest(nxpad,zisz,yjsz,nv)
c
      complex(b8) source(nzpad,xisz,yjsz,nv)

      complex(b8), allocatable :: recvbuf(:,:,:),buf12(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)


      integer :: iarg
      integer position,pos0,pos1,pos2
      integer i,j,k,n,ix,iz,x2,ii,nv

        integer ostep
        data ostep/0/
	real*8 rtime1,rtime2,rtime3,rtime0
c
      real(b8) norm
c
        integer ithr,omp_get_thread_num,iz1,iz2
!        integer, allocatable :: jzp1(:,:),jzpi(:,:)
c
	rtime0=MPI_WTIME()
        rtime1=MPI_WTIME()
c
c already allocated in main program
c       allocate (iyp1(0:num_thr-1))
c        allocate (iypi(0:num_thr-1))

	call divide_thr (yjsz,iyp1,iypi,'kxcomm2: yjsz')
c
      norm = 1.

        if (jstep.gt.ostep.and.kstep.eq.1) then
	ostep=jstep
c        t_kxcomm2(:,2)=0.
        end if
c


      iarg = KrCntMax*yjsz/(b8*2)
      allocate (sendbuf(iarg,nv,0:iproc-1))
      allocate (recvbuf(iarg,nv,0:iproc-1))


      allocate(buf12(nxhppad,zisz,0:num_thr-1))

        rtime2=MPI_WTIME()
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)

        ithr=0

!$OMP PARALLEL private (ithr,ii,z2,pos0,x2,pos1,position,z,x,iz,n,pos2,iz1,j,ix,iz2,i)



! Pack send buffer
! Interchange indices, using loop blocking NB1

! Use MPI_Alltoall
c
!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

      do 100 j=1,nv
c
        do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1

         do i=0,iproc-1
            ii = mymap(i)

                pos0 = xisz * kisz(ii) * (yp-1)

               do z=kist(ii),kien(ii),NB1_Z
                  z2 = min(z+NB1_Z-1,real(kien(ii))) !#deepcustom#
                  pos1 = pos0 + (z-kist(ii))*xisz

                  do x=1,xisz,NB1_X
                     x2 = min(x+NB1_X-1,real(xisz)) !#deepcustom#
                     pos2 = pos1 + x

                     do iz=z,z2
                        position = pos2
                        do ix=x,x2
                           sendbuf(position,j,i) = source(iz,ix,yp,j)
                           position = position +1
                        enddo
                        pos2 = pos2 + xisz
                     enddo
                  enddo
               enddo
            enddo
         enddo
 100  continue

!$OMP BARRIER
!$OMP MASTER



! Exchange x-z buffers

	rtime3=MPI_WTIME()


      iarg = KrCntMax*nv*yjsz
      call mpi_alltoall(sendbuf,iarg,mpi_byte,
     &                  recvbuf,iarg,mpi_byte,mpi_comm_row,ierr)
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_row)

	rtime3=MPI_WTIME()-rtime3
	t_kxcomm2(1,2)=t_kxcomm2(1,2)+rtime3

!$OMP END MASTER
!$OMP BARRIER


! Unpack receive buffers 



      do 300 j=1,nv
            do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1

!$OMP MASTER
            rtime1=MPI_WTIME()
!$OMP END MASTER

            do i=0,iproc-1
               ii = mymap(i)
               position = 1
               position = position + (yp-1)*iisz(ii)*zisz 
               do z=1,zisz
                  do x=iist(ii),iien(ii)
                     buf12(x,z,ithr) = recvbuf(position,j,i)
                     position = position+1
                  enddo
               enddo
            enddo
         
         
! Add and zero extra elements in X
!
            do z=1,zisz
               do x=nxhp,nxhppad
                  buf12(x,z,ithr) = 0.0
               enddo
            enddo
!$OMP MASTER
            rtime2=MPI_WTIME()
            tcpu_other=tcpu_other+(rtime2-rtime1)
            t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)
! C2R Transform    
            rtime1=MPI_WTIME()
!$OMP END MASTER
             call fftw_execute_c2r(plan3_p2c(ithr),buf12(1,1,ithr),dest(1,1,yp,j))
!$OMP MASTER
            rtime2=MPI_WTIME()
            tcpu_fft=tcpu_fft+(rtime2-rtime1)
            t_kxcomm2(2,2)=t_kxcomm2(2,2)+(rtime2-rtime1)
!$OMP END MASTER
         enddo
 300  continue
         
!$OMP END PARALLEL

        rtime1=MPI_WTIME()

      deallocate(sendbuf)
      deallocate (recvbuf)
      deallocate (buf12)

        rtime2=MPI_WTIME()
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)

      i4_comm = i4_comm + 1

	t_kxcomm2(4,2)=t_kxcomm2(4,2)+(MPI_WTIME()-rtime0)
c
      return
      end
