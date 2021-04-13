#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SETCOVERED(
     &           vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,m
     &           ,mask
     &           ,imasklo0,imasklo1,imasklo2
     &           ,imaskhi0,imaskhi1,imaskhi2
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2)
      REAL_T m
      integer imasklo0,imasklo1,imasklo2
      integer imaskhi0,imaskhi1,imaskhi2
      BYTE_T mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           imasklo2:imaskhi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      if ( m .eq. 0) then
       
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         vel(i,j,k)= vel(i,j,k)*mask(i,j,k)
       
      enddo
      enddo
      enddo
       ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxhi2- iboxlo2+1)
      else
       
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
 
         vel(i,j,k)= vel(i,j,k)*mask(i,j,k)
     &                     + (1-mask(i,j,k))*m
       
      enddo
      enddo
      enddo
       ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxhi2- iboxlo2+1)*4
      endif
      return
      end
