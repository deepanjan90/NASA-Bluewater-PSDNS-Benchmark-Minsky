












!############################################
!#           Wavespace
!#
!#     Program that operates when in wavespace
!#     Therefore, for homogeneous turbulence
!#     Program will transform y to finish forming 
!#     Nonlinear terms              
!#
!#     Built on proc3.f
!############################################
        subroutine wavespace_vel(uy,m)
           use comp
          !#deepcustom#	implicit none
          !#deepcustom# 	implicit none
 
          integer :: i,j,m, is1, i2f,a,xz(2)

          complex(b8) :: uy(ny,zjsz*xisz,nu)
          complex(b8), allocatable :: shift(:)

          complex(b8) :: utmp1
          complex(b8) :: syz,sxz,sxyz
	real(b8) norm

        integer ithr,omp_get_thread_num
	real(b8) upy_factor
	integer ii
	
	
        real(b8) rk1,rk1sq
        real s1,s2,s3


        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
c


	i2f=4+nc
	upy_factor=.5



!#############################################################
!# transpose to y-pencils and take inverse transform 
 

 
      norm=dt2/ny/nz
 
      is1=4+nc

        ithr=0
!$OMP PARALLEL private (ithr,x,xp,z,sxz,sxyz,utmp1)


      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1


         x = xp + xist-1

         rk1=s1*kx(x)
         rk1sq = rk1**2

            sxz=sz(z,kstep)*sx(x,kstep)


! form shift, shift, and complete formation of convective terms

            do 42 y=1,ny
               if( mask(y,a)) then
                  sxyz=-norm*conjg(sxz*sy(y,kstep))
                  utmp1=imagi*b22(m)*ky(y)*sxyz
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
 

 

      call next_xz(xp,z)

 100	continue

!$OMP END PARALLEL


c

      return
	end subroutine wavespace_vel
