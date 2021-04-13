#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine DIVERGESTENF(
     &           idcalclo0,idcalclo1,idcalclo2
     &           ,idcalchi0,idcalchi1,idcalchi2
     &           ,divf
     &           ,idivflo0,idivflo1,idivflo2
     &           ,idivfhi0,idivfhi1,idivfhi2
     &           ,ndivfcomp
     &           ,flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2
     &           ,ifluxhi0,ifluxhi1,ifluxhi2
     &           ,nfluxcomp
     &           ,facedir
     &           ,nconserved
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer ndivfcomp
      integer idivflo0,idivflo1,idivflo2
      integer idivfhi0,idivfhi1,idivfhi2
      REAL_T divf(
     &           idivflo0:idivfhi0,
     &           idivflo1:idivfhi1,
     &           idivflo2:idivfhi2,
     &           0:ndivfcomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL_T dx
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = CH_SPACEDIM
      do iv = 0,nconserved - 1
         
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0

         divf(i,j,k,iv) = divf(i,j,k,iv) +
     &        (flux(i+ioff,j+joff,k+koff,iv)
     &        -flux(i     ,j     ,k     ,iv))/dx
         
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*(idcalchi1- idcalclo1+1)*(idcalchi2- idcalclo2+1)*3*nconserved
      return
      end
