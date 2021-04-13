#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
        subroutine SETDIRICHLETFACEFLUX(
     &           faceFlux
     &           ,ifaceFluxlo0,ifaceFluxlo1,ifaceFluxlo2
     &           ,ifaceFluxhi0,ifaceFluxhi1,ifaceFluxhi2
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,value
     &           ,dx
     &           ,idir
     &           ,iside
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifaceFluxlo0,ifaceFluxlo1,ifaceFluxlo2
      integer ifaceFluxhi0,ifaceFluxhi1,ifaceFluxhi2
      REAL_T faceFlux(
     &           ifaceFluxlo0:ifaceFluxhi0,
     &           ifaceFluxlo1:ifaceFluxhi1,
     &           ifaceFluxlo2:ifaceFluxhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      REAL_T value
      REAL_T dx(0:2)
      integer idir
      integer iside
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k
        real_t ihdx
        ihdx = two/dx(idir)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          faceFlux(i,j,k) = iside * ihdx * (phi(i,j,k) - value)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine SETHODIRICHLETFACEFLUX(
     &           faceFlux
     &           ,ifaceFluxlo0,ifaceFluxlo1,ifaceFluxlo2
     &           ,ifaceFluxhi0,ifaceFluxhi1,ifaceFluxhi2
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,value
     &           ,dx
     &           ,idir
     &           ,iside
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifaceFluxlo0,ifaceFluxlo1,ifaceFluxlo2
      integer ifaceFluxhi0,ifaceFluxhi1,ifaceFluxhi2
      REAL_T faceFlux(
     &           ifaceFluxlo0:ifaceFluxhi0,
     &           ifaceFluxlo1:ifaceFluxhi1,
     &           ifaceFluxlo2:ifaceFluxhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      REAL_T value
      REAL_T dx(0:2)
      integer idir
      integer iside
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k, ioff,joff,koff
        real_t dxinv
        ioff = iside*chf_id(0,idir)
                  joff = iside*chf_id(1,idir)
                  koff = iside*chf_id(2,idir)
        dxinv = one/dx(idir)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          faceFlux(i,j,k) =
     &       iside * dxinv * (three*phi(i,j,k)
     &                  - one/three*phi(i+ioff,j+joff,k+koff)
     &                  - eight/three*value)
        
      enddo
      enddo
      enddo
        return
        end
