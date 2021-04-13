#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine EBAVERAGE(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,refrat
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2)
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer refrat
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      real_t refscale, coarsesum, frefine
      refscale = one
      frefine = refrat
      do ii = 1, CH_SPACEDIM
         refscale = refscale/frefine
      enddo
      
      do kcc = icoarboxlo2,icoarboxhi2
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
      
      do kk = ibreflo2,ibrefhi2
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      kff = kcc*refrat + kk 
      coarsesum = coarsesum + fine(iff,jff,kff)
      
      enddo
      enddo
      enddo
      coarse(icc,jcc,kcc) = coarsesum * refscale
      
      enddo
      enddo
      enddo
      return
      end
      subroutine EBAVERARZ(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,refrat
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           ,dxcoar
     &           ,dxfine
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2)
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer refrat
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      REAL_T dxcoar
      REAL_T dxfine
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      real_t coarsesum, coarsevol, finevol, sumfinevol
      
      do kcc = icoarboxlo2,icoarboxhi2
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
      sumfinevol = zero
      
      do kk = ibreflo2,ibrefhi2
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      kff = kcc*refrat + kk 
      finevol = (iff + half)*dxfine*dxfine*dxfine
      sumfinevol = sumfinevol + finevol
      coarsesum = coarsesum + finevol*fine(iff,jff,kff)
      
      enddo
      enddo
      enddo
      coarsevol = (icc + half)*dxcoar*dxcoar*dxcoar
      coarse(icc,jcc,kcc) = coarsesum/coarsevol
      
      enddo
      enddo
      enddo
      return
      end
      subroutine EBAVERAGEFACE(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,refrat
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           ,idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2)
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer refrat
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer idir
      integer idoloop, jdoloop, kdoloop
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      real_t refscale, coarsesum, frefine
      
      idoloop = 1-CHF_ID(idir,0)
      jdoloop = 1-CHF_ID(idir,1)
      kdoloop = 1-CHF_ID(idir,2)
      refscale = one
      frefine = refrat
      do ii = 1, (CH_SPACEDIM - 1)
         refscale = refscale/frefine
      enddo
      
      do kcc = icoarboxlo2,icoarboxhi2
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      coarsesum = zero
#if (CH_SPACEDIM == 3)
      do kk = 0,(refrat-1)*kdoloop
#endif
         do jj = 0,(refrat-1)*jdoloop
            do ii = 0,(refrat-1)*idoloop
               
               iff = icc*refrat + ii
               jff = jcc*refrat + jj
               kff = kcc*refrat + kk 
               coarsesum = coarsesum + fine(iff,jff,kff)
            enddo
         enddo
#if (CH_SPACEDIM == 3)
      enddo
#endif
      coarse(icc,jcc,kcc) = coarsesum * refscale
      
      enddo
      enddo
      enddo
      return
      end
      subroutine EBCOARSEN(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,h2laplfine
     &           ,ih2laplfinelo0,ih2laplfinelo1,ih2laplfinelo2
     &           ,ih2laplfinehi0,ih2laplfinehi1,ih2laplfinehi2
     &           ,icoarboxlo0,icoarboxlo1,icoarboxlo2
     &           ,icoarboxhi0,icoarboxhi1,icoarboxhi2
     &           ,refrat
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2)
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2)
      integer ih2laplfinelo0,ih2laplfinelo1,ih2laplfinelo2
      integer ih2laplfinehi0,ih2laplfinehi1,ih2laplfinehi2
      REAL_T h2laplfine(
     &           ih2laplfinelo0:ih2laplfinehi0,
     &           ih2laplfinelo1:ih2laplfinehi1,
     &           ih2laplfinelo2:ih2laplfinehi2)
      integer icoarboxlo0,icoarboxlo1,icoarboxlo2
      integer icoarboxhi0,icoarboxhi1,icoarboxhi2
      integer refrat
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      real_t refscale,coarsesum,frefine,coarseave,laplsum, laplave
      refscale = one
      frefine = refrat
      do ii = 1, CH_SPACEDIM
         refscale = refscale/2
      enddo
      
      do kcc = icoarboxlo2,icoarboxhi2
      do jcc = icoarboxlo1,icoarboxhi1
      do icc = icoarboxlo0,icoarboxhi0

      laplsum   = zero
      coarsesum = zero
      
      do kk = ibreflo2,ibrefhi2
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0

      
      iff = icc*refrat + ii
      jff = jcc*refrat + jj
      kff = kcc*refrat + kk 
      coarsesum = coarsesum +       fine(iff,jff,kff)
      laplsum   = laplsum   + h2laplfine(iff,jff,kff)
      
      enddo
      enddo
      enddo
      laplave   =   laplsum * refscale
      coarseave = coarsesum * refscale
      coarse(icc,jcc,kcc) = coarseave - laplave/eight
      
      enddo
      enddo
      enddo
      return
      end
      subroutine H2LAPL1DADDITIVE(
     &           opphidir
     &           ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     &           ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,idir
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,haslo
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,hashi
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iopphidirlo0,iopphidirlo1,iopphidirlo2
      integer iopphidirhi0,iopphidirhi1,iopphidirhi2
      REAL_T opphidir(
     &           iopphidirlo0:iopphidirhi0,
     &           iopphidirlo1:iopphidirhi1,
     &           iopphidirlo2:iopphidirhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i,j,k
      integer ioff,joff,koff
      
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

      opphidir(i,j,k) = opphidir(i,j,k) +
     $     (    phi(i+ioff,j+joff,k+koff)
     &     -two*phi(i     ,j     ,k     )
     &     +    phi(i-ioff,j-joff,k-koff))
      
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
         
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         opphidir(i,j,k) = opphidir(i,j,k) +
     $        (    phi(i       ,j       ,k       )
     &        -two*phi(i+  ioff,j+  joff,k+  koff)
     &        +    phi(i+2*ioff,j+2*joff,k+2*koff))
         
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
         
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         opphidir(i,j,k) = opphidir(i,j,k) +
     $        (    phi(i       ,j       ,k       )
     &        -two*phi(i-  ioff,j-  joff,k-  koff)
     &        +    phi(i-2*ioff,j-2*joff,k-2*koff))
         
      enddo
      enddo
      enddo
      endif
      return
      end
