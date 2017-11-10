













! Version with Cylindrical representation 
!
! Updates by PKY, 1/16/2012: assume jproc.gt.num_thr
!
! Starting with Y-cylinders, wavenumber-space in all three dimensions
! - Inverse-Transform in Y, 
! - transpose to get data in z-pencils with stride-1 in Z
!
! Multivariable version
!
!     subroutine kxcomm1_trans_cyl (source,dest,nv)
      subroutine kxcomm1_clean (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
	use timers_rkstep

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

      complex(b8) source(ny,zjsz*xisz,nv)
      complex(b8) dest(nz,xisz,yjsz,nv)

      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)
      integer position,pos0,pos1,n,a
      integer, allocatable :: pos(:,:)
      integer i,j,iz,dnz,iy,ia,a1,a2,c,nv
      real vmax
	integer :: iarg
c
	real*8 rtime0,rtime1,rtime2
c
        integer ithr,omp_get_thread_num,iy1,iy2,num_fft,ia1,jthr,NB2Y
        integer, allocatable :: ja_st(:)

       integer, allocatable :: jpp1(:),jppi(:)

        if (kstep.eq.1) t_kxcomm1(:,2)=0.
c

	rtime1=MPI_WTIME()

      if(num_al_i(0) .gt. 0) then
c
        allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
        allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
c
c ASSUME jproc.ge.num_thr
c
        allocate (jpp1(0:num_thr-1),jppi(0:num_thr-1))
        call divide_thr (jproc,jpp1,jppi,'xkcomm2: jproc')
        jpp1(:)=jpp1(:)-1

! This loop does the following:
! - for each x, it executes zjsz transforms in Y, to reduce memory references
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

        allocate(pos(0:jproc-1,0:num_thr-1))
	
        allocate (ja_st(0:num_thr-1))
        ia=1
        do jthr=0,num_thr-1
        ja_st(jthr)=ia
        ia=ia+num_fft0(jthr)
        end do
c
      end if

      rtime2=MPI_WTIME()
      tcpu_other=tcpu_other+(rtime2-rtime1)
      t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)

      call divide_thr (yjsz,iyp1,iypi,'kxcomm1: yjsz')

      ithr=0
c
!$OMP PARALLEL private (ithr,y,a2,pos1,ia,ia1,position,jthr,i,n,c,a1,iy1,iy2,x,z,a)


      if(num_al_i(0) .gt. 0) then

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
        do j=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
        sendbuf (:,:,j)=0.
        end do
!$OMP MASTER
        rtime2=MPI_WTIME()
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

        if(num_al .gt. 0) then
c
c
         do 100 j=1,nv

!$OMP MASTER
          rtime1=MPI_WTIME()
!$OMP END MASTER

          ia1=ja_st(ithr)
            call fftw_execute(plan4_p1(ithr),source(1,ia1,j),source(1,ia1,j))

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_kxcomm1(2,2)=t_kxcomm1(2,2)+(rtime2-rtime1)
!$OMP END MASTER


!$OMP BARRIER

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

            do i=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1

                do y=jjst(i),jjen(i),NB2_Y


                  y2 = min(y+NB2_Y-1,real(jjen(i))) !#deepcustom#
               
                  do a=1,num_al,NB2_Z
                     a2 = min(a+NB2_Z-1,num_al)
                  
                     pos1 = a + (y-jjst(i))*num_al
               
                     do iy=y,y2
                        position = pos1
                        do ia=a,a2
                           sendbuf(position,j,i) = source(iy,ia,j)
                           position = position +1
                        enddo
                        pos1 = pos1 + num_al
                     enddo
                  enddo
               enddo
            
            enddo

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

 100	continue
      
c
      endif

!$OMP BARRIER
!$OMP MASTER

      t_alltoall = t_alltoall - MPI_Wtime()
      t3_comm = t3_comm - MPI_Wtime() 

	iarg=jjsz(0)*num_al_i(0)*b8*nv*2
c
	rtime1=MPI_WTIME()

!      call mpi_alltoall(sendbuf,iarg,mpi_byte,
!     &      		recvbuf,iarg,mpi_byte,mpi_comm_col,ierr)

      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_col)

	rtime2=MPI_WTIME()

      t3_comm = t3_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
        t_kxcomm1(1,2)=t_kxcomm1(1,2)+(rtime2-rtime1)
!$OMP END MASTER
!$OMP BARRIER


!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
	do j=1,nv
	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	do x=max_al_x+1,xisz
	do z=1,nz
	dest(z,x,y,j)=0.
	end do
	end do
	end do
	end do
c

	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	
	do i=0,jproc-1
	pos(i,ithr)=(y-1)*num_al_i(i)+1
	end do

            i = 0
            n=0
            do x=1,max_al_x
               c = cut_z(x)
               a1 = nzhp - c -1
               a2 = nzhp + c +1
               do z=1,a1
                  dest(z,x,y,1:nv) = recvbuf(pos(i,ithr),1:nv,i) 
                  pos(i,ithr) = pos(i,ithr)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
	do z=a1+1,a2-1
	dest(z,x,y,1:nv)=0.
	end do
	
               
               do z=a2,nz
                  dest(z,x,y,1:nv) = recvbuf(pos(i,ithr),1:nv,i) 

                  pos(i,ithr) = pos(i,ithr)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
            enddo
         enddo
c
         
!$OMP MASTER
      i3_comm = i3_comm + 1
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

      else
! if num_al_i(0) = 0
c
!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
         do j=1,nv
	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
               do x=1,xisz
                  do z=1,nz
                     dest(z,x,y,j) = 0.0
                  enddo
               enddo
            enddo
         enddo

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER
c
      endif


!$OMP END PARALLEL

c

      if(num_al_i(0) .gt. 0) then
c
      deallocate (pos)
c
	deallocate (ja_st)
      deallocate (recvbuf)
      deallocate (jpp1,jppi)
c
	end if

!RAF Must allocate coarrays on all images, so all ranks must deallocate them
      if (allocated(sendbuf)) deallocate(sendbuf)
c
      return
      end
