        subroutine ranseq(k,lu)
c
c  purpose:
c       this routine permits different random numbers sequences within
c       the same program.  each sequence is designated a k-number.
c       ( k is the first argument to ranseq.  1 .le. k .le. nk=20 .)
c       random numbers from the sequence can be generated by ranu2,
c       or by any routine calling ranu2, or by any routine using the
c       seed  iseed  in common block  rancom.
c  use:
c       call ranseq(k) causes random numbers to be generated from sequen
c       this sequence will be used until the next call.
c
c       the argument lu is relevant only if k is -1 or -2 (see below)
c
c  initialization:
c       by default, all k seeds are initially set to the same value.
c       ranseq must be called to set initial k value before any other
c       random number routines are called.
c       if ranseq is called with the argument k= -1, then the seeds are
c       read (formatted) from logical unit lu.
c       alternatively, seeds can be specified by calling seqput.
c
c  checkpointing:
c       if ranseq is called with the argument k= -2, then the seeds are
c       written (formatted) on logical unit lu.
c       alternatively, seeds can be obtained by calling seqget.
c
c-----------------------------------------------------------------------
c
#ifdef MOL
	use compart, only: ngm,nom
#endif
	use mpicom
	  use com, only: seed_input
c
#ifdef RANDOM_SEED
	use ranseed
#endif
	implicit none
c
	integer iseq,nk,klast
        parameter ( nk=20 )
        common/seqcom/iseq(nk),klast
        integer iseqdf(nk)
	integer seed,ifst,k,lu,i,iset
c
	integer no_err,no_end
c
        data ifst/0/
	data iseqdf/ nk* 0 /
c
#ifdef MOL
	integer, allocatable :: alliseq(:),jseed(:),allseeds(:)
        integer seedsize
#endif

c  first call only
c
c	if (taskid.eq.0) write(6,*)'enter ranseq: ifst,k=',ifst,k

        if( ifst .eq. 0 ) then
           ifst = 1

	 if (taskid.eq.0) then
	 write (6,*) 'taskid,1st call of ranseq: k,lu=',taskid,k,lu
 	 write (6,*) 'ranseq: taskid,seed_input=',taskid,seed_input
	 end if

!	if (iseqdf(1).eq.0.and.k.ne.-1) then
!	write (6,*) 'specify default random no. seed'
!	read (5,*) seed

	if (iseqdf(1).eq.0.and.k.ne.-1) then
	seed=seed_input
	do 5 i=1,nk
 5 	iseqdf(i)=seed
        if (taskid.eq.0) write (6,*) 'ranseq: taskid,seed from stdin:',taskid,seed
	end if
c
c  initial value of seeds
c
           if( k .eq. -1 ) then
c
		no_err=0
	       if (taskid.eq.0) then
		read (lu,100,err=15,end=15) iseq
		end if
		go to 16
 15		no_err=1
 16		continue
	       call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
       if (no_err.gt.0) then
        if (taskid.eq.0) write (6,*)
     1          'Code stops because of error in reading ''ranin'' file'
        stop
	end if
		
	       call MPI_BCAST (iseq,nk,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c
           else
c
c  have they been set in seqput ?  if not, use default values.
c
              call seqset(1,iset)
              if( iset .eq. 0 ) then
                 do 10 i=1,nk
10               iseq(i) = iseqdf(i)

#ifdef PARTIN_AR
		iseq(12) = iseq(12) + taskid
#endif

#ifdef MOL
		do i=13,12+ngm
		iseq(i) = iseq(i) + taskid
		enddo
#endif
              endif

           endif
c
           klast = 1
           iseed = iseq(1)
           if( k .eq. -1 ) return
c
        endif
c
c  general call
c
        if( k .gt. nk ) then
           write(6,*)' k=',k,' must be less than nk=',nk
           write(6,*)' parameter nk can be increased in ranseq '
           write(6,*)' stopped in ranseq '
           stop
c
        else if( k .gt. 0 ) then
           iseq(klast) = iseed
           klast = k
           iseed = iseq(klast)
           continue
c
        else if( k .eq. -1 ) then
           read(lu,100)iseq

        else if( k .eq. -4 ) then
#ifdef MOLOLD
	if (nom.gt.0) then
	allocate (alliseq(numtasks*ngm))


        no_err=0
        no_end=0
	if(taskid.eq.0) then
        read(lu,100,err=91,end=92) alliseq
        go to 30
 91     no_err=1
        go to 30
 92     no_end=1
 30     continue
        end if
        call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_BCAST (no_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if (no_err.gt.0) then
        if (taskid.eq.0) write (6,*)
     1          'Code stops because of error in reading ''mranin'' file'
        stop
        end if
        if (no_end.gt.0) then
        if (taskid.eq.0) write (6,*)
     1   'Code stops because of end-of-file in reading ''mranin'' file'
        stop
        end if

	call MPI_SCATTER (alliseq,ngm,MPI_INTEGER,iseq(13:12+ngm),ngm,
     1                 MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	
	deallocate (alliseq)
	end if
#endif
c
        else if( k .eq. -2 ) then
           iseq(klast) = iseed
!          rewind lu
!       call fclose1 (lu)
!       call fopen1 (lu,'ranout ','formatted  ')
!           write(lu,100)iseq
!       call fclose1 (lu)
c
	if (taskid.eq.0) write(lu,100)iseq

#ifdef MOLOLD
	if (nom.gt.0) then
	allocate (alliseq(numtasks*ngm))
	call MPI_GATHER (iseq(13:12+ngm),ngm,MPI_INTEGER,alliseq,ngm,
     1                 MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if(taskid.eq.0) write(lu+1,100) alliseq
	deallocate (alliseq)
	end if

#elif defined MOL
        if(nom.gt.0) then
        
        call random_seed (size=seedsize)
        allocate(jseed(seedsize),allseeds(seedsize*numtasks))
        call random_seed (get=jseed)
        call MPI_GATHER (jseed,seedsize,MPI_INTEGER,allseeds,seedsize,
     1                 MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

        if(taskid.eq.0) then
        open(unit=2301,file='molseedout')      
        write(2301,*) seedsize, numtasks
        write(2301,100) allseeds
	close(2301)
        endif

        deallocate(jseed,allseeds)
        endif


#endif
c
        else
           write(6,*)' invalid argument: k=',k
           write(6,*)' stopped in ranseq '
           stop
        endif
c 
        return
c
100     format( (5i15) )
c
        end