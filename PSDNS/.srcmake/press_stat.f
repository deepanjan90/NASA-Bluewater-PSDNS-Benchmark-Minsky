












      subroutine press_stat

! Last deallocate statement in this routine added on 10/3/2012
c

      use comsp
      !#deepcustom#	implicit none
      !#deepcustom# 	implicit none

	real, allocatable :: eky(:),eks(:)
c
	real term,rk2,sum
	integer a,ik
c
	integer ithr,omp_get_thread_num,iy1,iy2,ii
c
	integer, allocatable :: ixp(:),iz(:)
        real, allocatable :: eky_thr(:,:)
	integer icall
	save icall
	data icall/0/
c
        real(b8) beta_min

        beta_min=min(beta1,beta2,beta3)

	return
	end
