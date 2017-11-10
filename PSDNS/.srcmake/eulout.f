













	subroutine eulout
	
	use comsp
	!#deepcustom#	implicit none       
	integer lu,i,j
	real uiui,sumjeps
	integer :: k,ios


	integer nshell
	real(b8) shell

	lu=15

	if (taskid.ne.0) return
c
        uiui=corr(1)+corr(4+nc)+corr(6+2*nc)

      write (lu,210) istep,time,sqrt(b11(2)),sqrt(b22(2)),sqrt(b33(2))
c
      write (lu,211) 'u',(laak(1,j),j=1,3),ekl
      write (lu,211) 'v',(laak(2,j),j=1,3)
      write (lu,211) 'w',(laak(3,j),j=1,3)

	write (lu,212) (taak(i,i),i=1,3)
      write (lu,213) (raak(i),i=1,3)
      write (lu,214) (tmre(i),i=1,3),(tmre(1)+tmre(2)+tmre(3))/3.
      write (lu,215) (rms(i),i=1,3),corr(2)

	write (lu,216) tke
 	write (lu,217) epslon, kmax*klen
      write (lu,218) skew
      write (lu,219) klen,kvel,ktime
      write (lu,231) (ett(i),i=1,3)

      shell=min(beta1,beta2,beta3)
c     nshell=nxh
c     nshell=mxyz
        nshell=max(nxh,nyh,nzh)
c
	write (lu,2101) shell
c
	write (lu,220) (ek(k),k=1,nshell)
	write (lu,2102)
	write (lu,220) (dk(k),k=1,nshell)

      close (lu)
      open (lu,file='eulstat',form='formatted',iostat=ios,
     1            position='append')
      ! Write out the Reynolds stresses.
      WRITE(996,997) corr(1), corr(2), corr(3),
     1               corr(4+nc), corr(5+nc),
     2               corr(6+2*nc)
 997  FORMAT(1P,6E13.5)
      CLOSE(UNIT=996)
      OPEN(UNIT=996,FILE='rstensor',FORM='FORMATTED',POSITION='APPEND')

      ! Write out grid metrics.
      WRITE(999,*) istep, time, beta1, beta2, beta3
      CLOSE(UNIT=999)
      OPEN(UNIT=999,FILE='metrics',FORM='FORMATTED',POSITION='APPEND')

      ! Write out the 3D spectra, and their trace.
      WRITE(997,240) istep,time,3+nc
      WRITE(997,242) '11', (eijk(i,1),i=1,10)
      WRITE(997,243)       (eijk(i,1),i=11,mxyz)
      WRITE(997,242) '12', (eijk(i,2),i=1,10)
      WRITE(997,243)       (eijk(i,2),i=11,mxyz)
      WRITE(997,242) '13', (eijk(i,3),i=1,10)
      WRITE(997,243)       (eijk(i,3),i=11,mxyz)
      WRITE(997,242) '22', (eijk(i,4+nc),i=1,10)
      WRITE(997,243)       (eijk(i,4+nc),i=11,mxyz)
      WRITE(997,242) '23', (eijk(i,5+nc),i=1,10)
      WRITE(997,243)       (eijk(i,5+nc),i=11,mxyz)
      WRITE(997,242) '33', (eijk(i,6+2*nc),i=1,10)
      WRITE(997,243)       (eijk(i,6+2*nc),i=11,mxyz)
      WRITE(997,242) 'ii', ((eijk(i,1)+eijk(i,4+nc)+eijk(i,6+2*nc)),i=1,10)
      WRITE(997,243)       ((eijk(i,1)+eijk(i,4+nc)+eijk(i,6+2*nc)),i=11,mxyz)
      CLOSE(UNIT=997)
      OPEN(UNIT=997,FILE='sptr3d',FORM='FORMATTED',POSITION='APPEND')
 240  FORMAT (/,1X,' istep,time,3+nc=',I5,1P,E14.5,2X,I4)
 241  FORMAT (1X,A)
 242  FORMAT (1X,A,1P,10E13.5)
 243  FORMAT (3X,1P,10E13.5)

      WRITE(991,240) istep,time,3+nc
      write (991,214) (tmre(i),i=1,3),(tmre(1)+tmre(2)+tmre(3))/3.
	write (991,216) tke
 	write (991,217) epslon, kmax*klen
      write (991,219) klen,kvel,ktime
	write (991,2101) shell
        do k=1,nshell
        write (991,*) k,ek(k)
        end do
      CLOSE(UNIT=991)
      OPEN(UNIT=991,FILE='spectrum.free',FORM='FORMATTED',POSITION='APPEND')



 210  format (/'istep=',i6,2x,'time=',1p,e13.5,'  betas=',1p,3e13.5)
 211  format (a1,'integral length scales:',1p,4e12.4)
 212  format (' taylor microscales    :',1p,3e12.4)
 213  format (' int. scale reynolds no:',1p,3e12.4)
 214  format (' taylor reynolds no.   :',1p,4e12.4)
 215  format (' rms velocities & uv   :',1p,4e12.4)
 216  format (' turb. kinetic energy  :',1p,e12.4)
 217  format (' dissipation rate      :',1p,e12.4,2x,'kmax.eta=',1p,e12.4)
 218  format (' dissipation skewness  :',1p,e12.4)
 219  format (' kol. scales (l,v,t)   :',1p,3e12.4)
 220  format ((3x,1p,10e13.5))
 231  format (' eddy turnover times   :',1p,3e12.4)
 2101 format ('3-d energy spectrum, shell thickness=',1p,e12.4)
 2102 format ('3-d dissipation spectrum')

	return 
	end subroutine eulout
