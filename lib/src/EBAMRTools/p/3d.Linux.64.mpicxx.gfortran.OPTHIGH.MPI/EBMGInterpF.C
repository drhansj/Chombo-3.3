#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

        subroutine REGPROLONG(
     &           phithislevel
     &           ,iphithislevello0,iphithislevello1,iphithislevello2
     &           ,iphithislevelhi0,iphithislevelhi1,iphithislevelhi2
     &           ,correctcoarse
     &           ,icorrectcoarselo0,icorrectcoarselo1,icorrectcoarselo2
     &           ,icorrectcoarsehi0,icorrectcoarsehi1,icorrectcoarsehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,irefboxlo0,irefboxlo1,irefboxlo2
     &           ,irefboxhi0,irefboxhi1,irefboxhi2
     &           ,reftocoar
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphithislevello0,iphithislevello1,iphithislevello2
      integer iphithislevelhi0,iphithislevelhi1,iphithislevelhi2
      REAL_T phithislevel(
     &           iphithislevello0:iphithislevelhi0,
     &           iphithislevello1:iphithislevelhi1,
     &           iphithislevello2:iphithislevelhi2)
      integer icorrectcoarselo0,icorrectcoarselo1,icorrectcoarselo2
      integer icorrectcoarsehi0,icorrectcoarsehi1,icorrectcoarsehi2
      REAL_T correctcoarse(
     &           icorrectcoarselo0:icorrectcoarsehi0,
     &           icorrectcoarselo1:icorrectcoarsehi1,
     &           icorrectcoarselo2:icorrectcoarsehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer irefboxlo0,irefboxlo1,irefboxlo2
      integer irefboxhi0,irefboxhi1,irefboxhi2
      integer reftocoar
        integer iic,jjc,kkc
        integer iie,jje,kke
        integer iif,jjf,kkf
        real_t fineval, coarval
        
      do kkc = icoarboxlo2,icoarboxhi2
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0

        
      do kke = irefboxlo2,irefboxhi2
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0

        
        iif  =  reftocoar*iic  + iie
        jjf  =  reftocoar*jjc  + jje
        kkf  =  reftocoar*kkc  + kke
        fineval = phithislevel(iif,jjf,kkf)
        coarval = correctcoarse(iic,jjc,kkc)
        phithislevel(iif,jjf,kkf) = fineval + coarval
        
      enddo
      enddo
      enddo
        
      enddo
      enddo
      enddo
        ch_flops=ch_flops+(icoarboxhi0- icoarboxlo0+1)*(icoarboxhi1- icoarboxlo1+1)*(icoarboxhi2- icoarboxlo2+1)*(irefboxhi0- irefboxlo0+1)*(irefboxhi1- irefboxlo1+1)*(irefboxhi2- irefboxlo2+1)
        return
        end
        subroutine PROLONGADDSLOPE(
     &           phithislevel
     &           ,iphithislevello0,iphithislevello1,iphithislevello2
     &           ,iphithislevelhi0,iphithislevelhi1,iphithislevelhi2
     &           ,correctcoarse
     &           ,icorrectcoarselo0,icorrectcoarselo1,icorrectcoarselo2
     &           ,icorrectcoarsehi0,icorrectcoarsehi1,icorrectcoarsehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,irefboxlo0,irefboxlo1,irefboxlo2
     &           ,irefboxhi0,irefboxhi1,irefboxhi2
     &           ,idir
     &           ,dxf
     &           ,dxc
     &           ,reftocoar
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iphithislevello0,iphithislevello1,iphithislevello2
      integer iphithislevelhi0,iphithislevelhi1,iphithislevelhi2
      REAL_T phithislevel(
     &           iphithislevello0:iphithislevelhi0,
     &           iphithislevello1:iphithislevelhi1,
     &           iphithislevello2:iphithislevelhi2)
      integer icorrectcoarselo0,icorrectcoarselo1,icorrectcoarselo2
      integer icorrectcoarsehi0,icorrectcoarsehi1,icorrectcoarsehi2
      REAL_T correctcoarse(
     &           icorrectcoarselo0:icorrectcoarsehi0,
     &           icorrectcoarselo1:icorrectcoarsehi1,
     &           icorrectcoarselo2:icorrectcoarsehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer irefboxlo0,irefboxlo1,irefboxlo2
      integer irefboxhi0,irefboxhi1,irefboxhi2
      integer idir
      REAL_T dxf
      REAL_T dxc
      integer reftocoar
        integer iic,jjc,kkc
        integer iie,jje,kke
        integer iif,jjf,kkf
        integer ioff,joff,koff
        real_t coarslope, hival, loval, phifold, phifnew
        real_t dxcoar, dxfine, coarloc, fineloc, midval
        integer iindexc, iindexf,cbox
        dxcoar = dxc
        dxfine = dxf
        ioff = chf_id(0,idir)
                  joff = chf_id(1,idir)
                  koff = chf_id(2,idir)
        cbox=(icoarboxhi0- icoarboxlo0+1)*(icoarboxhi1- icoarboxlo1+1)*(icoarboxhi2- icoarboxlo2+1)
        
      do kkc = icoarboxlo2,icoarboxhi2
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0

        iindexc = ioff*iic + joff*jjc +koff*kkc
        coarloc = dxcoar*(iindexc + half)
        hival   = correctcoarse(iic + ioff,jjc + joff,kkc + koff)
        loval   = correctcoarse(iic - ioff,jjc - joff,kkc - koff)
        midval  = correctcoarse(iic,jjc,kkc)
        coarslope = (hival -loval)/(two*dxcoar)
        
      do kke = irefboxlo2,irefboxhi2
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0

        
        iif  =  reftocoar*iic  + iie
        jjf  =  reftocoar*jjc  + jje
        kkf  =  reftocoar*kkc  + kke
        iindexf = ioff*iif + joff*jjf +koff*kkf
        fineloc = dxfine*(iindexf + half)
        phifold = phithislevel(iif,jjf,kkf)
        phifnew = phifold + coarslope*(fineloc - coarloc)
        phithislevel(iif,jjf,kkf) = phifnew
        
      enddo
      enddo
      enddo
        
      enddo
      enddo
      enddo
        ch_flops=ch_flops+cbox*(irefboxhi0- irefboxlo0+1)*(irefboxhi1- irefboxlo1+1)*(irefboxhi2- irefboxlo2+1)*3+cbox*5
        return
        end
