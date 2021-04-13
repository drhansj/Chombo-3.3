      subroutine EBAVECELLTOFACE(
     & facedata
     & ,ifacedatalo0,ifacedatalo1,ifacedatalo2
     & ,ifacedatahi0,ifacedatahi1,ifacedatahi2
     & ,celldata
     & ,icelldatalo0,icelldatalo1,icelldatalo2
     & ,icelldatahi0,icelldatahi1,icelldatahi2
     & ,facedir
     & ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     & ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacedatalo0,ifacedatalo1,ifacedatalo2
      integer ifacedatahi0,ifacedatahi1,ifacedatahi2
      REAL*8 facedata(
     & ifacedatalo0:ifacedatahi0,
     & ifacedatalo1:ifacedatahi1,
     & ifacedatalo2:ifacedatahi2)
      integer icelldatalo0,icelldatalo1,icelldatalo2
      integer icelldatahi0,icelldatahi1,icelldatahi2
      REAL*8 celldata(
     & icelldatalo0:icelldatahi0,
     & icelldatalo1:icelldatahi1,
     & icelldatalo2:icelldatahi2)
      integer facedir
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k
      integer ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      facedata(i,j,k) =
     & ( celldata(i ,j ,k )
     & + celldata(i-ioff,j-joff,k-koff)
     & )*(0.500d0)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*(ifaceboxhi1- iface
     &boxlo1+1)*(ifaceboxhi2- ifaceboxlo2+1)*2
      return
      end
