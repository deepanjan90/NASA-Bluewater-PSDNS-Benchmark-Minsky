












! Version with Cylindrical representation 
!
! Updates by PKY, 1/16/2012: assume jproc.gt.num_thr
!
! Starting with data in Z-pencils, X and Z in wavenumber space, 
! Y in real space:
! - transpose data into Y-cylinders, changing indices ordering as well
!   to arrange the data in stride-1 in Y
! - Then do a forward transform in Y
!
c This version has improvements made by Dmitry, Dec 22, 2011
c
! Multivariable version
!
      subroutine xkcomm2_trans_cyl(source,dest,nv)

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

      integer i,j,iz,dnz,iy,ia,a1,a2,c,nv
      complex(b8) source(nzpad,xisz,yjsz,nv)
      complex(b8) dest(nypad,xisz*zjsz,nv)

      integer position,pos0,pos1,pos2,n,a
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

      integer, allocatable :: pos(:,:)


      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

	integer :: iarg
	integer :: ii,jj
c
	real*8 rtime0,rtime1,rtime2,rtime3,rtime4
c
	integer lu
c
        integer ithr,omp_get_thread_num,iy1,iy2,ia1
 	integer, allocatable :: ja_st(:)

	integer ip,index,jthr,jzthr,NB2Y
c
	integer, allocatable :: jpp1(:),jppi(:)
c
c	if (kstep.eq.1) t_xkcomm2(:,2)=0.

	rtime1=MPI_WTIME()
c
      if(num_al_i(0) .gt. 0) then
        allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
      endif
      if(num_al_i(0) .gt. 0) then
        allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
        sendbuf = 0.
c

	rtime1=MPI_WTIME()-rtime1
	tcpu_other=tcpu_other+rtime1
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime1

        rtime1=MPI_WTIME()

      allocate(pos(0:jproc-1,0:num_thr-1))
c
c ASSUME jproc.ge.num_thr
c
	allocate (jpp1(0:num_thr-1),jppi(0:num_thr-1))
	call divide_thr (jproc,jpp1,jppi,'xkcomm2: jproc')
	jpp1(:)=jpp1(:)-1

	ithr=0
c
!$OMP PARALLEL private(ithr,i,j,x,y,z,n,c,a1,a2,iy1,iy2,ip,index,jthr,pos1,position,ia,ia1,y1,y2)
c


      do i=0,jproc-1
         pos(i,ithr) = 1+num_al_i(i)*(iyp1(ithr)-1)
      enddo
	
	iy1=iyp1(ithr)
	iy2=iyp1(ithr)+iypi(ithr)-1

      do y=iy1,iy2
         i = 0
         n=0
         do x=1,max_al_x
            c = cut_z(x) 
	    a1 = nzhp - c -1
            a2 = nzhp + c +1
            do z=1,a1

	sendbuf(pos(i,ithr),1:nv,i) = source(z,x,y,1:nv)
               pos(i,ithr) = pos(i,ithr)+1
               n = n+1
               if(n .eq. num_al_i(i)) then
                  n = 0
                  i = i+1
               endif
            enddo
            
            do z=a2,nz
	sendbuf(pos(i,ithr),1:nv,i) = source(z,x,y,1:nv)
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
!$OMP BARRIER
!$OMP MASTER
c

        rtime1=MPI_WTIME()-rtime1
        tcpu_other=tcpu_other+rtime1
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime1

      t_alltoall = t_alltoall - MPI_Wtime()
      t2_comm = t2_comm - MPI_Wtime() 


	iarg=jjsz(0)*num_al_i(0)*b8*nv*2
c
	rtime1=MPI_WTIME()

      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_col)
        rtime1=MPI_WTIME()-rtime1
        t_xkcomm2(1,2)=t_xkcomm2(1,2)+rtime1
c
      t2_comm = t2_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()

      deallocate (pos)
c
!$OMP END MASTER
!$OMP BARRIER

c
      if(num_al .gt. 0) then


!$OMP MASTER
	allocate (ja_st(0:num_thr-1))
	ia=1
	do jthr=0,num_thr-1
	ja_st(jthr)=ia
	ia=ia+num_fft0(jthr)
	end do

!$OMP END MASTER
!$OMP BARRIER

c
           do j=1,nv

!$OMP MASTER
        rtime2=MPI_WTIME()
!$OMP END MASTER

         do i=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
c
            do y=jjst(i),jjen(i),NB2_Y
               y2 = min(y+NB2_Y-1,real(jjen(i))) !#deepcustom#
               do a=1,num_al,NB2_Z
                  a2 = min(a+NB2_Z-1,num_al)
                  pos1 = a + (y-jjst(i))*num_al
                  do iy=y,y2
                     position = pos1
                     do ia=a,a2
                         dest(iy,ia,j) = recvbuf(position,j,i)
                        position = position +1
                     enddo
                     pos1 = pos1 + num_al
                  enddo
               enddo
            enddo
            
	end do
c	end if

!$OMP BARRIER

!$OMP MASTER
        rtime2=MPI_WTIME()-rtime2
        tcpu_other=tcpu_other+rtime2
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime2

        rtime3=MPI_WTIME()
!$OMP END MASTER

	ia1=ja_st(ithr)
            call fftw_execute(plan4_p3(ithr),dest(1,ia1,j),dest(1,ia1,j))
        

!$OMP BARRIER
!$OMP MASTER
        rtime3=MPI_WTIME()-rtime3
        tcpu_fft=tcpu_fft+rtime3
        t_xkcomm2(2,2)=t_xkcomm2(2,2)+rtime3
!$OMP END MASTER

      enddo
c
!$OMP MASTER
	deallocate (ja_st)
!$OMP END MASTER

      endif
         
!$OMP END PARALLEL

      deallocate (recvbuf)

      i2_comm = i2_comm + 1

	if (allocated(jppi)) deallocate (jpp1,jppi)
c
      endif
! num_al_i(0) > 0

      if (allocated(sendbuf)) deallocate (sendbuf)
c
      return
      end

!-----------------------------------------------------------------
! Square pencils version (used in dwrtu2)
!
! Starting with data in Z-pencils (X-dimension leading), X and Z in 
! wavenumber space, 
! - transpose data into Y-pencils, changing indices ordering as well
!   to arrange the data in stride-1 in Y
! (only transpose - no transform)
!-----------------------------------------------------------------

      subroutine xkcomm2_pen(source,dest,nv)

      use com
      use timers_comm
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none

      integer i,j,iz,dnz,iy,ia,a,a1,a2,c,nv
      complex(b8) source(xisz,nz,yjsz,nv)
      complex(b8) dest(xisz,ny,zjsz,nv)

      integer position,pos0,pos1,pos2,n,cnt
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

      complex(b8), allocatable :: sendbuf(:,:,:),recvbuf(:,:,:)
      integer, allocatable :: pos(:)

      t2_comm = t2_comm - MPI_Wtime() 

      allocate (sendbuf(xisz*yjsz*nz,nv,1))
      allocate (recvbuf(xisz*zjsz*ny,nv,1))

      t2_comm = t2_comm - MPI_Wtime() 

      do j=1,nv
         
         position = 1
         do i=0,jproc-1
            do y=1,yjsz
               do z=kjst(i),kjen(i)
                  do x=1,xisz
                     sendbuf(position,j,1) = source(x,z,y,j)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo
         
      enddo
      
! Now send the buffers (exchange data)

        cnt = xisz * yjsz * zjsz * b8 * 2


      call compi_alltoall(sendbuf,recvbuf,cnt,mpi_comm_col)


      do j=1,nv
         
         position = 1
         do y=1,ny
            do z=1,zjsz
               do x=1,xisz
                  dest(x,y,z,j) = recvbuf(position,j,1)
                  position = position +1
               enddo
            enddo
         enddo

      enddo
         
      deallocate (sendbuf)
      deallocate (recvbuf)
      
      t2_comm = MPI_Wtime() + t2_comm
      i2_comm = i2_comm + 1
c
      
      return
      end


!-----------------------------------------------------------------
! Version with Cylindrical representation (used in inpen_read)
!
! Starting with data in Z-pencils, X and Z in wavenumber space, 
! - transpose data into Y-cylinders, changing indices ordering as well
!   to arrange the data in stride-1 in Y
! (only transpose - no transform)
!-----------------------------------------------------------------

      subroutine xkcomm2_sq2cyl(source,dest,nv)

      use com
      use timers_comm
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none

      integer position,pos0,pos1,pos2,n,nv
      complex(b8) source(nz,xisz,yjsz,nv)
      complex(b8) dest(ny,xisz*zjsz,nv)

      integer i,j,iz,dnz,iy,ia,a,a1,a2,c
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

      integer, allocatable :: pos(:)
      complex(b8), allocatable :: sendbuf(:,:,:),recvbuf(:,:,:)

      t2_comm = t2_comm - MPI_Wtime() 

      allocate(pos(0:jproc-1))
      allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
      allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))

      sendbuf = 0
      recvbuf = 0

      do j=1,nv

      do i=0,jproc-1
         pos(i) = 1
      enddo

      do y=1,yjsz
         i = 0
         n=0
         do x=1,max_al_x
            c = cut_z(x)
            a1 = nzhp - c -1
            a2 = nzhp + c +1
            do z=1,a1
               sendbuf(pos(i),j,i) = source(z,x,y,j)
               pos(i) = pos(i)+1
               n = n+1
               if(n .eq. num_al_i(i)) then
                  n = 0
                  i = i+1
               endif
            enddo
            
            do z=a2,nz
               sendbuf(pos(i),j,i) = source(z,x,y,j)
               pos(i) = pos(i)+1
               n = n+1
               if(n .eq. num_al_i(i)) then
                  n = 0
                  i = i+1
               endif
            enddo
            
         enddo
      enddo
      enddo
      
      deallocate (pos)


      call compi_alltoall(sendbuf,recvbuf,jjsz(0)*num_al_i(0)*b8*2*nv,
     &  mpi_comm_col)


      if(num_al .gt. 0) then

      do j=1,nv
         do i=0,jproc-1
            
            do y=jjst(i),jjen(i),NB2_Y
               y2 = min(y+NB2_Y-1,real(jjen(i))) !#deepcustom#
               
               do a=1,num_al,NB2_Z
                  a2 = min(a+NB2_Z-1,num_al)
                  
                  pos1 = a + (y-jjst(i))*num_al
                  
                  do iy=y,y2
                     position = pos1
                     do ia=a,a2
                        dest(iy,ia,j) = recvbuf(position,j,i)
                        position = position +1
                     enddo
                     pos1 = pos1 + num_al
                  enddo
               enddo
            enddo
            
         enddo

      enddo
      endif

         
      deallocate (sendbuf)
      deallocate (recvbuf)
      
      t2_comm = MPI_Wtime() + t2_comm
      i2_comm = i2_comm + 1
      
      return
      end

!-----------------------------------------------------------------
! Square pencils version (used in inhomogeneous version)
!
! Starting with data in Z-pencils, X and Z in wavenumber space,
! Y in real space and padded:
! - remove extra elements from the center in Z,
! - transpose data into Y-pencils, changing indices ordering as well
!   to arrange the data in stride-1 in Y
! - Then do a forward transform in Y
!-----------------------------------------------------------------

      subroutine xkcomm2_trans_pen(source,dest)

      use com
      use timers_comm
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

      complex(b8) source(nzpad,xisz,yjsz)
      complex(b8) dest(nypad,zjsz,xisz)

      integer position,pos0,pos1,pos2
      integer i,j,n,iz,dnz,iy
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

! If using MPI_Alltoall
      complex(b8), allocatable :: sendbuf(:,:),recvbuf(:,:)

      allocate (sendbuf(KfCntMax/(b8*2),0:jproc-1))
      allocate (recvbuf(KfCntMax/(b8*2),0:jproc-1))

      t2_comm = t2_comm - MPI_Wtime() 


! Pack the sendbuf, omitting the center nzpad-nz elements in Z dimension
 
      dnz = nzpad-nz
      do i=0,jproc-1
         position = 1

! If clearly in the first half of nz
         if(kjen(i) .le. nzh) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),kjen(i)
                     sendbuf(position,i) = source(z,x,y)
                     position = position+1
                  enddo
               enddo
            enddo
! If clearly in the second half of nz
         else if (kjst(i) .ge. nzh+1) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i)+dnz,kjen(i)+dnz
                     sendbuf(position,i) = source(z,x,y)
                     position = position +1
                  enddo
               enddo
            enddo         

! If spanning the first and second half of nz (e.g. iproc is odd)  
         else
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),nzh
                     sendbuf(position,i) = source(z,x,y)
                     position = position +1
                  enddo
                  do z=nzpad-nzh+1,kjen(i)+dnz
                     sendbuf(position,i) = source(z,x,y)
                     position = position +1
                  enddo
               enddo
            enddo
         endif
         
! Not needed if preallocating send buffer
         n = KfCntMax/(b8*2) - kjsz(i)*xisz*yjsz
         do j=1,n
            sendbuf(position,i) = 0.0
            position = position +1
         enddo
!
      enddo

! Now send the buffers (exchange data)


      call compi_alltoall(sendbuf,recvbuf,KfCntMax,
     &  mpi_comm_col)


! Unpack data from recvbuf using loop blocking, change the order of array indices 
! and do the Y transpose in the same loop to save on memory references

      do x=1,xisz

         do i=0,jproc-1
            
            pos0 = (x-1)*zjsz 
            do y=jjst(i),jjen(i),NB2_Y
               y2 = min(y+NB2_Y-1,real(jjen(i))) !#deepcustom#

               do z=1,zjsz,NB2_Z
                  z2 = min(z+NB2_Z-1,real(zjsz)) !#deepcustom#
                  pos1 = pos0 +z
                  
                  do iy=y,y2
                     position = pos1 
                     do iz=z,z2
                        dest(iy,iz,x) = recvbuf(position,i)
                        position = position +1
                     enddo
                     pos1 = pos1 + zjsz * xisz
                  enddo
               enddo
               pos0 = pos0 + xisz*zjsz*NB2_Y
            enddo
         enddo

      call fftw_execute(plan2_p3,dest(1,1,x),dest(1,1,x))

      enddo


      deallocate (sendbuf)
      deallocate (recvbuf)
      
      t2_comm = MPI_Wtime() + t2_comm
      i2_comm = i2_comm + 1
      
      return
      end

