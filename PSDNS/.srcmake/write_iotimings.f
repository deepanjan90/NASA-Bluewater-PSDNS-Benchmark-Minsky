












	subroutine write_iotimings 
c
       	use comp
	use timers_comm
	use timers_comp
	use timers_io
	use timers_rkstep
	use timers_tran

        !#deepcustom#	implicit none

	integer itask,i,i2t
	real(8) sum,comm2
        real(8), allocatable :: t_rks_all(:,:),t_itrans_all(:,:)
        integer, allocatable :: numal_all(:)
        real(8), allocatable :: t_xkcomm1_all(:,:)
        real(8), allocatable :: t_xkcomm2_all(:,:)
        real(8), allocatable :: t_kxcomm1_all(:,:)
        real(8), allocatable :: t_kxcomm2_all(:,:)
        real(8), allocatable :: t_kxcomm2t_all(:,:)
        real(8), allocatable :: t_trans_all(:,:)
c
	integer itask1,itask2
	character*20 string
	character*6 filepos
c
c detailed instrumentation timings, added by PKY, 2/3/2012
c
	if (ncpusteps_io.eq.0) return
c
	allocate (t_rks_all(10,0:numtasks-1))
	allocate (t_itrans_all(4,0:numtasks-1))
	allocate (t_trans_all(4,0:numtasks-1))
	allocate (t_kxcomm1_all(4,0:numtasks-1))
	allocate (t_xkcomm1_all(4,0:numtasks-1))
	allocate (t_xkcomm2_all(4,0:numtasks-1))

	t_rks(:,2)=t_rks(:,2)/ncpusteps_io
        call  MPI_ALLGATHER (t_rks(1,2),10,MPI_DOUBLE_PRECISION,
     1                       t_rks_all,10,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	t_itrans(:,3)=t_itrans(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_itrans(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_itrans_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)

	t_trans(:,3)=t_trans(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_trans(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_trans_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)

	t_kxcomm1(:,3)=t_kxcomm1(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_kxcomm1(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm1_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
c#ifdef 1
	allocate (t_kxcomm2_all(4,0:numtasks-1))
	t_kxcomm2(:,3)=t_kxcomm2(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_kxcomm2(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm2_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c#else
	allocate (t_kxcomm2t_all(4,0:numtasks-1))
	t_kxcomm2t(:,3)=t_kxcomm2t(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_kxcomm2t(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm2t_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c#endif
c
	t_xkcomm1(:,3)=t_xkcomm1(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_xkcomm1(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_xkcomm1_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
	t_xkcomm2(:,3)=t_xkcomm2(:,3)/ncpusteps_io
        call  MPI_ALLGATHER (t_xkcomm2(1,3),4,MPI_DOUBLE_PRECISION,
     1                       t_xkcomm2_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
	t_rks(:,2)=t_rks(:,2)*ncpusteps_io
	t_itrans(:,3)=t_itrans(:,3)*ncpusteps_io
	t_trans(:,3)=t_trans(:,3)*ncpusteps_io
	t_kxcomm1(:,3)=t_kxcomm1(:,3)*ncpusteps_io
	t_kxcomm2(:,3)=t_kxcomm2(:,3)*ncpusteps_io
	t_kxcomm2t(:,3)=t_kxcomm2t(:,3)*ncpusteps_io
	t_xkcomm1(:,3)=t_xkcomm1(:,3)*ncpusteps_io
	t_xkcomm2(:,3)=t_xkcomm2(:,3)*ncpusteps_io
c
	if (taskid.eq.1) then
c
        write (string,"('USE_EVEN, kxcomm2')")
c
	itask1=2*iproc-1
	itask2=numtasks-2*iproc

	filepos='rewind'
	if (ichkpt.gt.1) filepos='append'

 	open (77,file='iostep_timings/cpu_rk_all',position=filepos)
      write (77,711) nx,ny,nz,nc,iproc, jproc
 711  format ('Pencils code, nx,ny,nz,nc=',3i6,i3,'  iproc, jproc=',i4,
     1        ' x ',i6)
	num_thr=0
        write (77,"('num_thr=',i3,'  double prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	num_thr=1

	write (77,"(a20)") string
	write (77,"('Detailed breakdown of CPU costs per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid           itransform realspace',
     1             ' transform           overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	sum=0.
	do i=1,5
	sum=sum+t_rks_all(i,itask)
	end do
        write (77,630) itask,ipid_all(itask),jpid_all(itask),
     1             (t_rks_all(i,itask),i=1,6)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
 630    format (i6,2i5,1x,1p,6e10.3)
	close (77)
c
	open (77,file='iostep_timings/aftrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, after sub. transform, per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid    step     wavespace   advanc     barrier    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_rks_all(i,itask),i=7,10),t_rks_all(5,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
	end if
	end if
	end do
 631    format (i6,2i5,1p,5e11.3)
	close (77)
c
	open (77,file='iostep_timings/t_itrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. itransform_vel, per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
        if (nc.gt.0) then
        write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    scalars    total')")
        else
        write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    velgrad    total')")
        end if
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_itrans_all(i,itask),i=1,4),
     1             t_itrans_all(1,itask)+t_itrans_all(2,itask)
     1             +t_itrans_all(3,itask)+t_itrans_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='iostep_timings/t_trans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. transform_vel, per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid   xkcomm1    fft(z)     convec     xkcomm2    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_trans_all(i,itask),i=1,4),
     1             t_trans_all(1,itask)+t_trans_all(2,itask)
     1             +t_trans_all(3,itask)+t_trans_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='iostep_timings/t_kxcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm1_all(i,itask),i=1,3),
     1             t_kxcomm1_all(1,itask)+t_kxcomm1_all(2,itask)
     1             +t_kxcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	if (nc.eq.0) then
	open (77,file='iostep_timings/t_kxcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2 per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total      overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm2_all(i,itask),i=1,3),
     1             t_kxcomm2_all(1,itask)+t_kxcomm2_all(2,itask)
     1             +t_kxcomm2_all(3,itask),t_kxcomm2_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	else
	open (77,file='iostep_timings/t_kxcomm2t_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2t per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total      overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm2t_all(i,itask),i=1,3),
     1             t_kxcomm2t_all(1,itask)+t_kxcomm2t_all(2,itask)
     1             +t_kxcomm2t_all(3,itask),t_kxcomm2t_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	end if
	close (77)
c
	open (77,file='iostep_timings/t_xkcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_xkcomm1_all(i,itask),i=1,3),
     1             t_xkcomm1_all(1,itask)+t_xkcomm1_all(2,itask)
     1             +t_xkcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='iostep_timings/t_xkcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm2_clean per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_xkcomm2_all(i,itask),i=1,3),
     1             t_xkcomm2_all(1,itask)+t_xkcomm2_all(2,itask)
     1             +t_xkcomm2_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
	close (77)
c
	i2t=2
	if (nc.eq.0) i2t=1

	open (77,file='iostep_timings/t_comm_all',position=filepos)
	write (77,"('Detailed breakdown of alltoall communication costs, per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	if (i2t.eq.1) then
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")
	else
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2t  xkcomm1   xkcomm2    total     %Code')")
	end if

        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	if (i2t.eq.1) then
	comm2=t_kxcomm2_all(1,itask)
	else
	comm2=t_kxcomm2t_all(1,itask)
	end if
        sum=t_kxcomm1_all(1,itask)+comm2+
     1             t_xkcomm1_all(1,itask)+t_xkcomm2_all(1,itask)
        write (77,632) itask,ipid_all(itask),jpid_all(itask),
     1             t_kxcomm1_all(1,itask),comm2,
     1             t_xkcomm1_all(1,itask),t_xkcomm2_all(1,itask),
     1		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
 632    format (i6,2i5,2x,1p,5e10.3,0p,f7.1,'%')
c
	open (77,file='iostep_timings/t_loctsp_all',position=filepos)
	write (77,"('Detailed breakdown of local-transpose costs, per time step')")
	write (77,"('averaged over', i4,'  output steps')") ncpusteps_io
	if (i2t.eq.1) then
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")
	else
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2t  xkcomm1   xkcomm2    total     %Code')")
	end if
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	if (i2t.eq.1) then
	comm2=t_kxcomm2_all(3,itask)
	else
	comm2=t_kxcomm2t_all(3,itask)
	end if
        sum=t_kxcomm1_all(3,itask)+comm2+
     1             t_xkcomm1_all(3,itask)+t_xkcomm2_all(3,itask)
        write (77,632) itask,ipid_all(itask),jpid_all(itask),
     1             t_kxcomm1_all(3,itask),comm2,
     1             t_xkcomm1_all(3,itask),t_xkcomm2_all(3,itask),
     1		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	end if

	deallocate (t_rks_all)
	deallocate (t_itrans_all)
	deallocate (t_trans_all)
	deallocate (t_kxcomm1_all)
	deallocate (t_xkcomm1_all)
	deallocate (t_xkcomm2_all)
	deallocate (t_kxcomm2_all)
	deallocate (t_kxcomm2t_all)

	return
	end
