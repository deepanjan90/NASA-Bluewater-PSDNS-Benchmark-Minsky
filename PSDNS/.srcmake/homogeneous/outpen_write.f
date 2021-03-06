












      subroutine outpen_write (buf,iwgid)
c
c Revised by PK Yeung, 3/24/09, for batches/relay mode
c called by sub. outpen
c
c Clean-up of previous directives:
c --- Implemented JUl1308 and TIMERS_IO
c --- Removed BGW, WRTALL
c Also removed old choice of iwfld (always operate with iwfld=2)
c
c
	use comp
	!#deepcustom#	implicit none
	!#deepcustom# 	implicit none
c
	complex(b8) :: buf(xisz,ny,zjsz,3+nc)
c
c  routine to write all u-fields on the appropriate logical units
c
 	integer maxsav,lu,i,icall
      parameter (maxsav=100)
      character*6 numer
      character*6 head
      character*16 fn
c
        character*40 longname
        character*4 cpid
        character*2 char_icall
        integer nchar_pid,nchar_icall
c
	character*8 fname
	integer nchar,nf
	integer iwgid
c
      save icall
      data icall/0/
c
        write (numer,"(i6)") taskid
        call blanks (numer,nchar)

        write (cpid,"(i3,'/')") iwgid
        call blanks (cpid,nchar_pid)
c
c With fnmu(i)=outpen for i=1 to 3+nc, 
c files to be called outpen.pXXXX if expecting only 1 checkpoint,
c but outpenYY.pXXXX if more than 1 checkpoint
c
      if (isave.eq.0.and.dtsave.eq.0.) then
	fname=fnwu(1)
        else
	icall=icall+1
	head=fnwu(1)
	write (char_icall,"(i2)") icall
	call blanks (char_icall,nchar_icall)
	fname=head//char_icall(1:nchar_icall)
	end if
c
	call blanks (fname,nf)
      longname='outpen/'//cpid(1:nchar_pid)//fname(1:nf)//'.p'//numer(1:nchar)
	call blanks (longname,nchar)
c
      lu=luwu(1)
        open (lu,file=longname(1:nchar),form='unformatted')
c
	call time_stamp ('outpen_write: before do 100')
c
      do 100 i = 1 , 3+nc
      write( lu ) i , nx , ny , nz
      write( lu ) xist,xisz,zjst,zjsz,1,ny
c
      do 120 zp=1,zjsz
      z=zjst+zp-1
      if( z .eq. nzhp ) go to 120
      write (lu) ((buf(xp,y,zp,i),xp=1,xisz),y=1,nyh)
      write (lu) ((buf(xp,y,zp,i),xp=1,xisz),y=nyhp+1,ny)
120   continue
100   continue
c
	call time_stamp ('outpen_write:  after do 100')
c
	close (luwu(1))
c
c
      return
      end
