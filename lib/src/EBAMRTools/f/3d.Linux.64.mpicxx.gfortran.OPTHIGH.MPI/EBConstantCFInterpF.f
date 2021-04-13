        subroutine REGCONSTANTINTERP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,ivalidboxlo0,ivalidboxlo1,ivalidboxlo2
     & ,ivalidboxhi0,ivalidboxhi1,ivalidboxhi2
     & ,dir
     & ,nghost
     & ,side
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer ivalidboxlo0,ivalidboxlo1,ivalidboxlo2
      integer ivalidboxhi0,ivalidboxhi1,ivalidboxhi2
      integer dir
      integer nghost
      integer side
        integer ioff,joff,koff
        integer ighost, jghost, kghost
        integer ivalid, jvalid, kvalid
        integer whichghost
        if((side.ne.1).and.(side.ne.-1)) then
           call MAYDAY_ERROR()
        endif
        ioff = chf_id(0,dir)*side
        joff = chf_id(1,dir)*side
        koff = chf_id(2,dir)*side
      do kvalid = ivalidboxlo2,ivalidboxhi2
      do jvalid = ivalidboxlo1,ivalidboxhi1
      do ivalid = ivalidboxlo0,ivalidboxhi0
        do whichghost = 1, nghost
           ighost = ivalid + whichghost*ioff
           jghost = jvalid + whichghost*joff
           kghost = kvalid + whichghost*koff
           phi(ighost,jghost,kghost) = phi(ivalid,jvalid,kvalid)
        enddo
      enddo
      enddo
      enddo
        return
        end
