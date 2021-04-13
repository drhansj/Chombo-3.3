        subroutine COPYCFFAB(
     & finelevel
     & ,ifinelevello0,ifinelevello1,ifinelevello2
     & ,ifinelevelhi0,ifinelevelhi1,ifinelevelhi2
     & ,coarlevel
     & ,icoarlevello0,icoarlevello1,icoarlevello2
     & ,icoarlevelhi0,icoarlevelhi1,icoarlevelhi2
     & ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     & ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     & ,irefboxlo0,irefboxlo1,irefboxlo2
     & ,irefboxhi0,irefboxhi1,irefboxhi2
     & ,reftocoar
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelevello0,ifinelevello1,ifinelevello2
      integer ifinelevelhi0,ifinelevelhi1,ifinelevelhi2
      REAL*8 finelevel(
     & ifinelevello0:ifinelevelhi0,
     & ifinelevello1:ifinelevelhi1,
     & ifinelevello2:ifinelevelhi2)
      integer icoarlevello0,icoarlevello1,icoarlevello2
      integer icoarlevelhi0,icoarlevelhi1,icoarlevelhi2
      REAL*8 coarlevel(
     & icoarlevello0:icoarlevelhi0,
     & icoarlevello1:icoarlevelhi1,
     & icoarlevello2:icoarlevelhi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer irefboxlo0,irefboxlo1,irefboxlo2
      integer irefboxhi0,irefboxhi1,irefboxhi2
      integer reftocoar
        integer iic,jjc,kkc
        integer iie,jje,kke
        integer iif,jjf,kkf
        REAL*8 coarval
      do kkc = icoarboxlo2,icoarboxhi2
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0
          coarval = coarlevel(iic,jjc,kkc)
      do kke = irefboxlo2,irefboxhi2
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0
            iif = reftocoar*iic + iie
            jjf = reftocoar*jjc + jje
            kkf = reftocoar*kkc + kke
            finelevel(iif,jjf,kkf) = coarval
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
        ch_flops=ch_flops+(irefboxhi0- irefboxlo0+1)*(irefboxhi1- irefbo
     &xlo1+1)*(irefboxhi2- irefboxlo2+1)*(icoarboxhi0- icoarboxlo0+1)*(i
     &coarboxhi1- icoarboxlo1+1)*(icoarboxhi2- icoarboxlo2+1)
        return
        end
