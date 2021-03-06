      subroutine new_write (k,j1,j2)
c
c alternative to pop2_write4, 2/1/09
c 
c output action corresponding to wrt2pop4
c
c let every processor handle some of the output, 7/4/08
c (desirable if no. of particles and no. of processors are both very large)
c
c output of particle properties
c
c this routine is used if there is only 1 "group" of particles
c (otherwise, use wrtpop2 or wrtpop3)
c
c if iopt=0, simply close the files without writing new data
c
#ifdef LAG
#ifdef LAGSC2
c
	use compart
	implicit none
c
      character*15 fn
      character*6 name,numer
      character*11 unform
      logical secvel
      data unform/'unformatted'/
      character*110 string
	character*30 caux
c
	real(p8), allocatable :: pptmp(:,:),rbuf(:)
	integer, allocatable :: ppindtmp(:,:),rindbuf(:)
c
 	integer j,j1,j2,k,mp,ic,indx,ic2
c
	integer ipfst,iplst,nl,i,i1,i2,il,mm,npr,lu,npset
c
	real*8 rtime1,rtime2,rtime0
	real cpumpi,cpuio,cpu
c
cccc 
	integer iraddr,itask,jtask,irtask1,irtask2,dest,msgids,imt1
        integer, allocatable :: msgidr(:),msgtype(:,:)
        integer, allocatable :: mpistr(:,:),mpists(:)
	integer, allocatable :: ipout(:),ipin(:),ipall(:),ipall2(:),idisp(:)
c
	integer isubset,nid,nvar,iset,jj,itemp
c

	if (taskid.eq.0.and.jstep.eq.1) write (6,*) 'enter new_write, 
     1                    istep,j1,j2=',istep,j1,j2

	rtime0 = MPI_WTIME()

	cpumpi=0.
	cpuio=0.

	npset = ndpart2/nsubset 

! index of particle in pp2 is now used from pp2ind array
!	nvar = j2-j1+2
	nvar = j2-j1+1

	isubset=(taskid+1)/(numtasks/nsubset)+1
c
! npr = no. of particles per record
c
      npr=min0(8192,npset)
      nl=npset/npr
      if (mod(npset,npr).ne.0) nl=nl+1
c

	sync memory
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

	allocate (pptmp(nump2*nvar+1,nsubset))
	allocate (ppindtmp(nump2+1,nsubset))
	allocate (rbuf(npset*nvar + numtasks))
	allocate (rindbuf(npset + numtasks))
	allocate (ipout(nsubset),ipall(numtasks))
	allocate (ipall2(numtasks),idisp(numtasks))

	pptmp = 0.

	ipout(:) = 0

	do i=1,nump2
!	itemp = int(pp2(i,0))
!	iset = (itemp-1)/npset + 1
	iset = (pp2ind(i) - 1)/npset + 1
!	pptmp(nvar*ipout(iset) + 1,iset) = pp2(i,0)
	ppindtmp(ipout(iset) + 1,iset) = pp2ind(i)
	do j=j1,j2
	pptmp(nvar*ipout(iset) + j-j1+1,iset) = pp2(i,j)
	enddo
	ipout(iset) = ipout(iset) + 1 
	enddo

	ipall(:)=0


	rtime1=MPI_WTIME()

	do jj=1,nsubset

	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHER (ipout(jj), 1, MPI_INTEGER, ipall, 1, MPI_INTEGER,
     1                   nid, MPI_COMM_WORLD,mpierr)

	enddo

	if(sum(ipall).ne.0) then
	write(6,*) 'in new write: nump2,sum ipall,taskid=',nump2,sum(ipall),taskid
	endif


! send the variables from j1 to j2
	idisp(1)=0
	ipall2(1) = ipall(1)*nvar + 1
	do i=2,numtasks
	ipall2(i) = ipall(i)*nvar + 1
	idisp(i) = idisp(i-1) + ipall2(i-1)
	enddo
	do jj=1,nsubset
	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHERV (pptmp(1,jj), ipout(jj)*nvar+1, pptype, rbuf,
     1                   ipall2, idisp, pptype, nid, MPI_COMM_WORLD,mpierr)
	enddo

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

! send the particle indices
	idisp(1)=0
	ipall2(1) = ipall(1) + 1
	do i=2,numtasks
	ipall2(i) = ipall(i) + 1
	idisp(i) = idisp(i-1) + ipall2(i-1)
	enddo
	do jj=1,nsubset
	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHERV (ppindtmp(1,jj), ipout(jj)+1, MPI_INTEGER, rindbuf,
     1                   ipall2, idisp, MPI_INTEGER, nid, MPI_COMM_WORLD,mpierr)
	enddo

	deallocate (pptmp)
	deallocate (ppindtmp)

	rtime2=MPI_WTIME()
	cpumpi=cpumpi+rtime2-rtime1


!	allocate (pptmp(npset,0:nvar-1))
	allocate (pptmp(npset,nvar))


        if (mod(taskid,numtasks/nsubset).eq.0) then

	ic = 1
	ic2 = 1

	do i=1,numtasks

	if(ipall(i).ne.0) then
	do j=1,ipall(i)
!	indx = mod(int(rbuf(ic))-1,npset) + 1
!	pptmp(indx,0:nvar-1) = rbuf(ic:ic+nvar-1)
	indx = mod(rindbuf(ic2) - 1, npset) + 1
	pptmp(indx,1:nvar) = rbuf(ic:ic+nvar-1)
	ic = ic + nvar 
	ic2 = ic2 + 1
	enddo
	endif

	ic = ic + 1
	ic2 = ic2 + 1

	enddo


#ifdef DEBUG
	write(7000+taskid,*) 'j1,j2,nvar=',j1,j2,nvar
	do i=1,npset
	write(7000+taskid,*) i, pptmp(i,0)
	enddo
#endif



	do 100 j=j1,j2

	rtime1=MPI_WTIME()
      lu=lupop+k-1+1000
c
      if (j.eq.j1.and.jstep.eq.1) then 
	write (6,903) taskid,lu,istep-1,time,j1,j2
 903  format ('new_write: taskid,lu,istep-1,time,j1,j2=', 3i6,1p,e12.4,2i3)
	end if
c
	if (j.eq.j1) write (lu) istep-1,time

      do 15 il=1,nl
      i1=(il-1)*npr+1
      i2=min0(i1+npr-1,npset)
      write (lu) (pptmp(i,j-j1+1),i=i1,i2)
 15   continue
c
	rtime2=MPI_WTIME()
	cpuio=cpuio+rtime2-rtime1
c
 100   continue

	endif    ! if (mod(taskid,numtasks/nsubset).eq.0) then
c

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

	cpu = MPI_WTIME() -rtime0
	if (taskid.eq.0) then
	if (jstep.eq.1) open (71,file='lag_timings/timings_new_write')
	write (71,"('istep,cpumpi,cpuio,cpu=',i6,2i3,1p,3e11.3)") istep,j1,j2,cpumpi,cpuio,cpu
	end if 

	deallocate (pptmp,rbuf,rindbuf)


	return
c
#endif
#endif
      end
