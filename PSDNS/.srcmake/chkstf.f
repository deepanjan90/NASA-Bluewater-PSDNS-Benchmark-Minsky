












      subroutine chkstf (icall)
c
c used in MPL code, based on xlf 3.2
c 
c close and re-open statistics files with append option
c
	use com
      character*11 form,unform
      data form,unform/'formatted','unformatted'/
c
      character*1 tail
	character*8 caux
c
        write (6,*) 'enter chkstf, icall=',icall
      go to (1,2,3), icall
c
 1    continue
c
c files updated at the end of output step
c
      call clreop (7,'log     ',form)
      call clreop (8,'run     ',form)
      if (shear.gt.0.) call clreop (12,'e12sptr',form)
      if (io1d.eq.1) call clreop (13,'sptr1d ',form)
      if (io1d.eq.1) call clreop (14,'corr1d ',form)
      call clreop (15,'eulstat',form)
      CALL clreop(996, 'rstensor', form)
      if (nc.gt.0) call clreop (16,'escout ',form)
      if (nc.gt.0) then
      do i=1,nc
      if(grad(1,i).ne.0..or.grad(2,i).ne.0..or.grad(3,i).ne.0.)go to 21
      end do
      go to 22
 21   call clreop (17,'scflxsp',form)
      call clreop (94,'scgmom',form)
 22   if (ictf.eq.1) call clreop (17,'scflxsp',form)
      end if
      call clreop (30,'corrmat',form)
      if (icsptr.eq.1) call clreop (32,'radsptr',form)
      if (icsptr.eq.1) call clreop (33,'dijsptr',form)
      if (iovor.eq.1) call clreop (35,'vorstat',form)
      if (kforce.gt.0.) call clreop (50,'forcinp',form)
c
      if (nc.gt.0) call clreop (38,'scgrad ',form)
      if (nc.eq.2) call clreop (24,'difpdf1',form)
      if (nc.eq.2) call clreop (27,'cosptr1',form)
      if (nc.eq.3) then
      call clreop (27,'cosptr1',form)
      call clreop (28,'cosptr2',form)
      call clreop (29,'cosptr3',form)
      end if


	call clreop (751,'dns_pdf_epsenst',form)
	call clreop (752,'dns_max_epsenst',form)
c
c
      return
c
 2    continue
c
c files updated at the beginning of next step after an output step
c (e.g. physical space output)
c
      if (nc.eq.0) call clreop (18,'vgskew ',form)
      call clreop (31,'trfsptr',form)
      if (icvc.eq.1) call clreop (34,'tijsptr',form)
      if (nc.gt.0.and.ictf.eq.1) 
     1    call clreop (37,'sctrfsp',form)
      if (nc.gt.0) call clreop (38,'scgrad ',form)
c
      return
c
 3    continue
c
c scalar histograms
c
	write (6,*) 'chkstf, at 3'
      do i=1,nc
        write (tail,"(i1)") i
c     lu=lusc+(i-1)*5
c     call clreop (lu,'scpdf'//tail,form)
      lu=luscn+(i-1)*6
	caux='scnpdf'//tail
      call clreop (lu,caux,form)
      lu=luscd+(i-1)*6
	caux='scdlpdf'//tail
      call clreop (lu,caux,form)
      end do
c
c
c
 90     return
        end
