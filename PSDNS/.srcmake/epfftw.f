












	subroutine epfftw(ux,uxc,uy,uz)
	use comp
	!#deepcustom#	implicit none

 	real(b8) :: ux(nxpad,zisz,yjsz,nu)
 	complex(b8) :: uxc(nxhpad,zisz,yjsz,nu)
	complex(b8) :: uy(nypad,zjsz*xisz,nu)
 	complex(b8) :: uz(nzpad,xisz,yjsz,nu)

	complex(b8), allocatable :: buf(:,:,:)
	complex(b8), allocatable :: buf2(:,:,:)
c

	integer num_fft,num_fft1,num,l
	integer ithr

	real*8 rtime1,rtime2


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

	if (fftw_flag.eq.FFTW_ESTIMATE) then
	if (taskid.eq.0) print *,'enter epfftw: ESTIMATE'
	elseif (fftw_flag.eq.FFTW_MEASURE) then
	if (taskid.eq.0) print *,'enter epfftw: MEASURE'
	endif

!$OMP MASTER
	rtime1=MPI_WTIME(ierr)
!$OMP END MASTER
c
	call epfftw_wisdom(0,taskid,2) ! try to read wisdom

      num_thr = 1

!      allocate (num_fft0(0:num_thr-1))!moved to main.f Shine 20150313
	allocate (plan4_p1(0:num_thr-1))
	allocate (plan4_p3(0:num_thr-1))
	allocate (plan2_p2a(0:num_thr-1))
        allocate (plan3_p2a(0:num_thr-1))
        allocate (plan3_p2c(0:num_thr-1))
c
	num=num_al/num_thr
	num_fft0(:)=num
	l=mod(num_al,num_thr)
	if (l.ne.0) then
	do ithr=0,l-1
	num_fft0(ithr)=num+1
	end do
	end if
c	write (6,"('epfftw: taskid,num_al,num_fft0=',4i5)") 
c     1       taskid,num_al,num_fft0

! proc1	
	do ithr=0,num_thr-1
      call fftw_plan(plan4_p1(ithr),1,nypad,num_fft0(ithr),uy,
     &    NULL,1,nypad, uy,NULL,1,nypad,
     &    FFTW_BACKWARD,fftw_flag)
	end do

! proc2a	
	call divide_thr (xisz, ixp1,ixpi,'epfftw:plan2_p2a')
	do ithr=0,num_thr-1
	call fftw_plan(plan2_p2a(ithr),1,nzpad,ixpi(ithr),
     &    uz,NULL,1,nzpad, uz,NULL,1,nzpad,
     &    FFTW_BACKWARD,fftw_flag)
	end do
c
c       call divide (zisz, num_thr, num_fft)
c       allocate (buf(nxhppad,zisz,yjsz))
c       call fftw_plan_c2r(plan3_p2a,1,nxpad,num_fft,
c     &     buf,NULL,1,nxhppad, uxc,NULL,1,nxpad,fftw_flag)
c       deallocate (buf)
c
        call divide_thr (zisz, izp1,izpi,'epfftw: zisz')
        allocate (buf(nxhppad,zisz,yjsz))
        do ithr=0,num_thr-1
        call fftw_plan_c2r(plan3_p2a(ithr),1,nxpad,izpi(ithr),
     &     buf,NULL,1,nxhppad, uxc,NULL,1,nxpad,fftw_flag)
        end do
        deallocate (buf)


        allocate (buf(nxhppad,zisz,1))
        do ithr=0,num_thr-1
        call fftw_plan_c2r(plan3_p2c(ithr),1,nxpad,zisz,
     &     buf,NULL,1,nxhppad, ux,NULL,1,nxpad,fftw_flag)
        end do
        deallocate (buf)

!proc2b
c

	allocate(buf2(nxhppad,zisz,yjsz))
c
        allocate (plan1_p2b(0:num_thr-1))
c
        call divide_thr (yjsz,iyp1,iypi,'epfftw: yjsz')
        do ithr=0,num_thr-1
        call fftw_plan_r2c(plan1_p2b(ithr),1,nxpad,iypi(ithr)*zisz,
     &     ux(1,1,1,1),nxpad,1,nxpad, buf2,nxhppad,1,nxhppad,fftw_flag)
        end do
        if (taskid.eq.0) then
        write (6,"('iyp1:',6i8)") iyp1
        write (6,"('iypi:',6i8)") iypi
        end if

	call divide (xisz, num_thr, num_fft)
	call fftw_plan(plan2_p2b,1,nzpad,num_fft,uz(1,1,1,1),
     &    NULL,1,nzpad, uz(1,1,1,1),NULL,1,nzpad,
     &    FFTW_FORWARD,fftw_flag)
c
	deallocate(buf2)

! proc3     
c
c
	do ithr=0,num_thr-1
      call fftw_plan(plan4_p3(ithr),1,nypad,num_fft0(ithr),uy,
     &    NULL,1,nypad, uy,NULL,1,nypad,
     &    FFTW_FORWARD,fftw_flag)
	end do

! kxtran     
	if (jproc.gt.1) then
c
        allocate (plan2_kx(0:num_thr-1))
       call divide_thr (yjsz,iyp1,iypi,'epfftw: yjsz')
        do ithr=0,num_thr-1
        call fftw_plan(plan2_kx(ithr),1,nzpad,iypi(ithr)*xisz,uz,
     &     NULL,1,nzpad, uz,NULL,1,nzpad,FFTW_BACKWARD,fftw_flag)
        end do

      else
	call divide (nz, num_thr, num_fft)
         call fftw_plan(plan1_kx,1,nypad,num_fft,
     &      uy,NULL,1,nypad, uz,NULL,1,nypad,
     &      FFTW_BACKWARD,fftw_flag)
	endif

      allocate(buf(nxhppad,zist:zien,yjst:yjen))
c

	call divide (zisz*yjsz, num_thr, num_fft)
      call fftw_plan_c2r(plan3_kx,1,nxpad,num_fft,
     &  buf,NULL,1,nxhp, ux,NULL,1,nxpad,fftw_flag)
c
      deallocate (buf)

! xktran     
c
	allocate(buf2(nxhppad,zisz,yjsz),stat=ierr)
c
c	call divide (xisz*yjsz, num_thr, num_fft)
c	call fftw_plan(plan2_xk,1,nzpad,num_fft,uz,
c    & 	   NULL,1,nzpad, uz,NULL,1,nzpad,FFTW_FORWARD,fftw_flag)
c
        allocate (plan2_xk(0:num_thr-1))
c
       call divide_thr (yjsz,iyp1,iypi,'epfftw: yjsz')
        do ithr=0,num_thr-1
        call fftw_plan(plan2_xk(ithr),1,nzpad,iypi(ithr)*xisz,uz,
     &     NULL,1,nzpad, uz,NULL,1,nzpad,FFTW_FORWARD,fftw_flag)
        end do

c
	if (jproc.eq.1) then
	   allocate(buf(nypad,nzpad,xisz))
	   buf = 0.
	call divide (num_al, num_thr, num_fft)
	   call fftw_plan(plan3_xk,1,nypad,num_fft,
     &	      buf,NULL,1,nypad, uy,NULL,1,nypad, FFTW_FORWARD,fftw_flag)
	   deallocate(buf)
	endif
c
	deallocate (buf2)

!$OMP MASTER
	rtime2=MPI_WTIME(ierr)
!$OMP END MASTER

	if (taskid.eq.0) 
     1 write (6,*) 'exit epfftw: taskid, time(secs) = ',taskid,rtime2-rtime1

	call epfftw_wisdom(1,taskid,2) ! save wisdom
	return
	end
c
	subroutine divide (n,m,num)
	integer l,n,m,num
	integer ithr,omp_get_thread_num
	num = n/m
	l = mod(n,m)
	if (l.ne.0) num=num+1
	return
	end
c
	subroutine divide2 (ithr,n,m,num)
	integer l,n,m,num
	num = n/m
	l = mod(n,m)
	if (l.ne.0.and.ithr.lt.m-1) num=num+1
	return
	end
