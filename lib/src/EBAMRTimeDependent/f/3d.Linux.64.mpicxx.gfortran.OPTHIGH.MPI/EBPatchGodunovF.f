      subroutine MINFLAT(
     & flattening
     & ,iflatteninglo0,iflatteninglo1,iflatteninglo2
     & ,iflatteninghi0,iflatteninghi1,iflatteninghi2
     & ,zetadir
     & ,izetadirlo0,izetadirlo1,izetadirlo2
     & ,izetadirhi0,izetadirhi1,izetadirhi2
     & ,nzetadircomp
     & ,du
     & ,idulo0,idulo1,idulo2
     & ,iduhi0,iduhi1,iduhi2
     & ,nducomp
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflatteninglo0,iflatteninglo1,iflatteninglo2
      integer iflatteninghi0,iflatteninghi1,iflatteninghi2
      REAL*8 flattening(
     & iflatteninglo0:iflatteninghi0,
     & iflatteninglo1:iflatteninghi1,
     & iflatteninglo2:iflatteninghi2)
      integer nzetadircomp
      integer izetadirlo0,izetadirlo1,izetadirlo2
      integer izetadirhi0,izetadirhi1,izetadirhi2
      REAL*8 zetadir(
     & izetadirlo0:izetadirhi0,
     & izetadirlo1:izetadirhi1,
     & izetadirlo2:izetadirhi2,
     & 0:nzetadircomp-1)
      integer nducomp
      integer idulo0,idulo1,idulo2
      integer iduhi0,iduhi1,iduhi2
      REAL*8 du(
     & idulo0:iduhi0,
     & idulo1:iduhi1,
     & idulo2:iduhi2,
     & 0:nducomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      integer iv
      REAL*8 sumdu,minflattot,minzetadir
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sumdu = (0.0d0)
      do iv = 0,nducomp - 1
         sumdu = sumdu + du(i,j,k,iv)
      enddo
      if (sumdu .lt. (0.0d0)) then
         minflattot = zetadir(i,j,k,0)
         do iv = 1,nducomp - 1
            minzetadir = zetadir(i,j,k,iv)
            minflattot = min(minflattot,minzetadir)
         enddo
         flattening(i,j,k) = minflattot
      else
         flattening(i,j,k) = (1.0d0)
      endif
      enddo
      enddo
      enddo
      return
      end
      subroutine GETDPTWO(
     & delta2p
     & ,idelta2plo0,idelta2plo1,idelta2plo2
     & ,idelta2phi0,idelta2phi1,idelta2phi2
     & ,delta1p
     & ,idelta1plo0,idelta1plo1,idelta1plo2
     & ,idelta1phi0,idelta1phi1,idelta1phi2
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idelta2plo0,idelta2plo1,idelta2plo2
      integer idelta2phi0,idelta2phi1,idelta2phi2
      REAL*8 delta2p(
     & idelta2plo0:idelta2phi0,
     & idelta2plo1:idelta2phi1,
     & idelta2plo2:idelta2phi2)
      integer idelta1plo0,idelta1plo1,idelta1plo2
      integer idelta1phi0,idelta1phi1,idelta1phi2
      REAL*8 delta1p(
     & idelta1plo0:idelta1phi0,
     & idelta1plo1:idelta1phi1,
     & idelta1plo2:idelta1phi2)
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k
      integer ioff,joff,koff
      REAL*8 dp1hi, dp1lo
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      dp1hi = delta1p(i+ioff,j+joff,k+koff)
      dp1lo = delta1p(i-ioff,j-joff,k-koff)
      delta2p(i,j,k) = dp1hi + dp1lo
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         dp1hi = delta1p(i+ioff,j+joff,k+koff)
         dp1lo = delta1p(i ,j ,k )
         delta2p(i,j,k) = dp1hi + dp1lo
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         dp1hi = delta1p(i ,j ,k )
         dp1lo = delta1p(i-ioff,j-joff,k-koff)
         delta2p(i,j,k) = dp1hi + dp1lo
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine GETFLAT(
     & zetatwiddle
     & ,izetatwiddlelo0,izetatwiddlelo1,izetatwiddlelo2
     & ,izetatwiddlehi0,izetatwiddlehi1,izetatwiddlehi2
     & ,delta1p
     & ,idelta1plo0,idelta1plo1,idelta1plo2
     & ,idelta1phi0,idelta1phi1,idelta1phi2
     & ,delta2p
     & ,idelta2plo0,idelta2plo1,idelta2plo2
     & ,idelta2phi0,idelta2phi1,idelta2phi2
     & ,bulkmin
     & ,ibulkminlo0,ibulkminlo1,ibulkminlo2
     & ,ibulkminhi0,ibulkminhi1,ibulkminhi2
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer izetatwiddlelo0,izetatwiddlelo1,izetatwiddlelo2
      integer izetatwiddlehi0,izetatwiddlehi1,izetatwiddlehi2
      REAL*8 zetatwiddle(
     & izetatwiddlelo0:izetatwiddlehi0,
     & izetatwiddlelo1:izetatwiddlehi1,
     & izetatwiddlelo2:izetatwiddlehi2)
      integer idelta1plo0,idelta1plo1,idelta1plo2
      integer idelta1phi0,idelta1phi1,idelta1phi2
      REAL*8 delta1p(
     & idelta1plo0:idelta1phi0,
     & idelta1plo1:idelta1phi1,
     & idelta1plo2:idelta1phi2)
      integer idelta2plo0,idelta2plo1,idelta2plo2
      integer idelta2phi0,idelta2phi1,idelta2phi2
      REAL*8 delta2p(
     & idelta2plo0:idelta2phi0,
     & idelta2plo1:idelta2phi1,
     & idelta2plo2:idelta2phi2)
      integer ibulkminlo0,ibulkminlo1,ibulkminlo2
      integer ibulkminhi0,ibulkminhi1,ibulkminhi2
      REAL*8 bulkmin(
     & ibulkminlo0:ibulkminhi0,
     & ibulkminlo1:ibulkminhi1,
     & ibulkminlo2:ibulkminhi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      REAL*8 d,r0,r1,ratio,strength
      REAL*8 d1pvof, d2pvof, smallp
      data d /0.33d0/
      data r0 /0.75d0/
      data r1 /0.85d0/
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      strength = abs(delta1p(i,j,k)/bulkmin(i,j,k))
      smallp = 1.0e-7
      if (strength .ge. d) then
         d1pvof = abs(delta1p(i,j,k))
         d2pvof = max(abs(delta2p(i,j,k)),smallp)
         ratio = d1pvof/d2pvof
         if (ratio .le. r0) then
            zetatwiddle(i,j,k) = (1.0d0)
         else if (ratio .ge. r1) then
            zetatwiddle(i,j,k) = (0.0d0)
         else
            zetatwiddle(i,j,k) = (1.0d0) - (ratio - r0)/(r1 - r0)
         endif
      else
         zetatwiddle(i,j,k) = (1.0d0)
      endif
      enddo
      enddo
      enddo
      return
      end
      subroutine GETGRAD(
     & du
     & ,idulo0,idulo1,idulo2
     & ,iduhi0,iduhi1,iduhi2
     & ,u
     & ,iulo0,iulo1,iulo2
     & ,iuhi0,iuhi1,iuhi2
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idulo0,idulo1,idulo2
      integer iduhi0,iduhi1,iduhi2
      REAL*8 du(
     & idulo0:iduhi0,
     & idulo1:iduhi1,
     & idulo2:iduhi2)
      integer iulo0,iulo1,iulo2
      integer iuhi0,iuhi1,iuhi2
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & iulo2:iuhi2)
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k
      integer ioff,joff,koff
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      du(i,j,k) = (0.500d0)*
     & ( u(i+ioff,j+joff,k+koff)
     & - u(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         du(i,j,k) =
     & ( u(i+ioff,j+joff,k+koff)
     & - u(i ,j ,k ))
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         du(i,j,k) =
     & ( u(i ,j ,k )
     & - u(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine GETRELATIVEGRAD(
     & du
     & ,idulo0,idulo1,idulo2
     & ,iduhi0,iduhi1,iduhi2
     & ,u
     & ,iulo0,iulo1,iulo2
     & ,iuhi0,iuhi1,iuhi2
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idulo0,idulo1,idulo2
      integer iduhi0,iduhi1,iduhi2
      REAL*8 du(
     & idulo0:iduhi0,
     & idulo1:iduhi1,
     & idulo2:iduhi2)
      integer iulo0,iulo1,iulo2
      integer iuhi0,iuhi1,iuhi2
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & iulo2:iuhi2)
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k
      integer ioff,joff,koff
      REAL*8 diff, ave
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      diff = (0.500d0)*
     & ( u(i+ioff,j+joff,k+koff)
     & - u(i-ioff,j-joff,k-koff))
      ave = (0.500d0)*
     & ( u(i+ioff,j+joff,k+koff)
     & + u(i-ioff,j-joff,k-koff))
      du(i,j,k) = diff/ave
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         diff =
     & ( u(i+ioff,j+joff,k+koff)
     & - u(i ,j ,k ))
         ave =(0.500d0)*
     & ( u(i+ioff,j+joff,k+koff)
     & + u(i ,j ,k ))
         du(i,j,k) = diff/ave
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         diff =
     & ( u(i ,j ,k )
     & - u(i-ioff,j-joff,k-koff))
         ave = (0.500d0)*
     & ( u(i ,j ,k )
     & + u(i-ioff,j-joff,k-koff))
         du(i,j,k) = diff/ave
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine MAGNITUDE(
     & magdata
     & ,imagdatalo0,imagdatalo1,imagdatalo2
     & ,imagdatahi0,imagdatahi1,imagdatahi2
     & ,data
     & ,idatalo0,idatalo1,idatalo2
     & ,idatahi0,idatahi1,idatahi2
     & ,ndatacomp
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imagdatalo0,imagdatalo1,imagdatalo2
      integer imagdatahi0,imagdatahi1,imagdatahi2
      REAL*8 magdata(
     & imagdatalo0:imagdatahi0,
     & imagdatalo1:imagdatahi1,
     & imagdatalo2:imagdatahi2)
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2
      integer idatahi0,idatahi1,idatahi2
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & idatalo2:idatahi2,
     & 0:ndatacomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      integer iv
      REAL*8 cur,sum
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sum = (0.0d0)
      do iv = 0,ndatacomp-1
         cur = data(i,j,k,iv)
         sum = sum + cur*cur
      enddo
      magdata(i,j,k) = sqrt(sum)
      enddo
      enddo
      enddo
      return
      end
      subroutine MIN3PTS(
     & mindata
     & ,imindatalo0,imindatalo1,imindatalo2
     & ,imindatahi0,imindatahi1,imindatahi2
     & ,data
     & ,idatalo0,idatalo1,idatalo2
     & ,idatahi0,idatahi1,idatahi2
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer imindatalo0,imindatalo1,imindatalo2
      integer imindatahi0,imindatahi1,imindatahi2
      REAL*8 mindata(
     & imindatalo0:imindatahi0,
     & imindatalo1:imindatahi1,
     & imindatalo2:imindatahi2)
      integer idatalo0,idatalo1,idatalo2
      integer idatahi0,idatahi1,idatahi2
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & idatalo2:idatahi2)
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k
      integer ioff,joff,koff
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      mindata(i,j,k) = min(
     & data(i ,j ,k ),
     & data(i+ioff,j+joff,k+koff),
     & data(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      if (haslo .ne. 0) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         mindata(i,j,k) = min(
     & data(i ,j ,k ),
     & data(i+ioff,j+joff,k+koff))
      enddo
      enddo
      enddo
      endif
      if (hashi .ne. 0) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         mindata(i,j,k) = min(
     & data(i ,j ,k ),
     & data(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine SECONDSLOPEDIFFS(
     & deltawc
     & ,ideltawclo0,ideltawclo1,ideltawclo2
     & ,ideltawchi0,ideltawchi1,ideltawchi2
     & ,ndeltawccomp
     & ,deltawl
     & ,ideltawllo0,ideltawllo1,ideltawllo2
     & ,ideltawlhi0,ideltawlhi1,ideltawlhi2
     & ,ndeltawlcomp
     & ,deltawr
     & ,ideltawrlo0,ideltawrlo1,ideltawrlo2
     & ,ideltawrhi0,ideltawrhi1,ideltawrhi2
     & ,ndeltawrcomp
     & ,w
     & ,iwlo0,iwlo1,iwlo2
     & ,iwhi0,iwhi1,iwhi2
     & ,nwcomp
     & ,numslopes
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndeltawccomp
      integer ideltawclo0,ideltawclo1,ideltawclo2
      integer ideltawchi0,ideltawchi1,ideltawchi2
      REAL*8 deltawc(
     & ideltawclo0:ideltawchi0,
     & ideltawclo1:ideltawchi1,
     & ideltawclo2:ideltawchi2,
     & 0:ndeltawccomp-1)
      integer ndeltawlcomp
      integer ideltawllo0,ideltawllo1,ideltawllo2
      integer ideltawlhi0,ideltawlhi1,ideltawlhi2
      REAL*8 deltawl(
     & ideltawllo0:ideltawlhi0,
     & ideltawllo1:ideltawlhi1,
     & ideltawllo2:ideltawlhi2,
     & 0:ndeltawlcomp-1)
      integer ndeltawrcomp
      integer ideltawrlo0,ideltawrlo1,ideltawrlo2
      integer ideltawrhi0,ideltawrhi1,ideltawrhi2
      REAL*8 deltawr(
     & ideltawrlo0:ideltawrhi0,
     & ideltawrlo1:ideltawrhi1,
     & ideltawrlo2:ideltawrhi2,
     & 0:ndeltawrcomp-1)
      integer nwcomp
      integer iwlo0,iwlo1,iwlo2
      integer iwhi0,iwhi1,iwhi2
      REAL*8 w(
     & iwlo0:iwhi0,
     & iwlo1:iwhi1,
     & iwlo2:iwhi2,
     & 0:nwcomp-1)
      integer numslopes
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k ,lvar
      integer ioff,joff,koff
      REAL*8 dwr,dwl
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do lvar = 0, numslopes - 1
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         dwr = w(i+ioff,j+joff,k+koff,lvar)
     & - w(i ,j ,k ,lvar)
         dwl = w(i ,j ,k ,lvar)
     & - w(i-ioff,j-joff,k-koff,lvar)
         deltawr(i,j,k,lvar) = dwr
         deltawl(i,j,k,lvar) = dwl
         deltawc(i,j,k,lvar) = (0.500d0)*(dwr + dwl)
      enddo
      enddo
      enddo
         if (haslo .ne. 0) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            dwr = w(i+ioff,j+joff,k+koff,lvar)
     & - w(i ,j ,k ,lvar)
            deltawc(i,j,k,lvar) = dwr
            deltawl(i,j,k,lvar) = dwr
            deltawr(i,j,k,lvar) = dwr
      enddo
      enddo
      enddo
         endif
         if (hashi .ne. 0) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            dwl = w(i ,j ,k ,lvar)
     & - w(i-ioff,j-joff,k-koff,lvar)
            deltawc(i,j,k,lvar) = dwl
            deltawl(i,j,k,lvar) = dwl
            deltawr(i,j,k,lvar) = dwl
      enddo
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine FORTHSLOPEDIFFS(
     & delta4wc
     & ,idelta4wclo0,idelta4wclo1,idelta4wclo2
     & ,idelta4wchi0,idelta4wchi1,idelta4wchi2
     & ,ndelta4wccomp
     & ,w
     & ,iwlo0,iwlo1,iwlo2
     & ,iwhi0,iwhi1,iwhi2
     & ,nwcomp
     & ,delta2w
     & ,idelta2wlo0,idelta2wlo1,idelta2wlo2
     & ,idelta2whi0,idelta2whi1,idelta2whi2
     & ,ndelta2wcomp
     & ,numslopes
     & ,idir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndelta4wccomp
      integer idelta4wclo0,idelta4wclo1,idelta4wclo2
      integer idelta4wchi0,idelta4wchi1,idelta4wchi2
      REAL*8 delta4wc(
     & idelta4wclo0:idelta4wchi0,
     & idelta4wclo1:idelta4wchi1,
     & idelta4wclo2:idelta4wchi2,
     & 0:ndelta4wccomp-1)
      integer nwcomp
      integer iwlo0,iwlo1,iwlo2
      integer iwhi0,iwhi1,iwhi2
      REAL*8 w(
     & iwlo0:iwhi0,
     & iwlo1:iwhi1,
     & iwlo2:iwhi2,
     & 0:nwcomp-1)
      integer ndelta2wcomp
      integer idelta2wlo0,idelta2wlo1,idelta2wlo2
      integer idelta2whi0,idelta2whi1,idelta2whi2
      REAL*8 delta2w(
     & idelta2wlo0:idelta2whi0,
     & idelta2wlo1:idelta2whi1,
     & idelta2wlo2:idelta2whi2,
     & 0:ndelta2wcomp-1)
      integer numslopes
      integer idir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i ,j ,k ,lvar
      integer ioff,joff,koff
      REAL*8 dwr,dwl, vall, slol, valr, slor
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do lvar = 0, numslopes - 1
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         valr = w(i+ioff,j+joff,k+koff,lvar)
         slor = delta2w(i+ioff,j+joff,k+koff,lvar)
         vall = w(i-ioff,j-joff,k-koff,lvar)
         slol = delta2w(i-ioff,j-joff,k-koff,lvar)
         dwl = vall + (0.250d0)*slol
         dwr = valr - (0.250d0)*slor
         delta4wc(i,j,k,lvar) = (2.000d0 / 3.000d0)*(dwr - dwl)
      enddo
      enddo
      enddo
         if (haslo .ne. 0) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            delta4wc(i,j,k,lvar) = delta2w(i,j,k,lvar)
      enddo
      enddo
      enddo
         endif
         if (hashi .ne. 0) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            delta4wc(i,j,k,lvar) = delta2w(i,j,k,lvar)
      enddo
      enddo
      enddo
         endif
      enddo
      return
      end
      subroutine APPLYFLAT(
     & dw
     & ,idwlo0,idwlo1,idwlo2
     & ,idwhi0,idwhi1,idwhi2
     & ,ndwcomp
     & ,flattening
     & ,iflatteninglo0,iflatteninglo1,iflatteninglo2
     & ,iflatteninghi0,iflatteninghi1,iflatteninghi2
     & ,numslopes
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndwcomp
      integer idwlo0,idwlo1,idwlo2
      integer idwhi0,idwhi1,idwhi2
      REAL*8 dw(
     & idwlo0:idwhi0,
     & idwlo1:idwhi1,
     & idwlo2:idwhi2,
     & 0:ndwcomp-1)
      integer iflatteninglo0,iflatteninglo1,iflatteninglo2
      integer iflatteninghi0,iflatteninghi1,iflatteninghi2
      REAL*8 flattening(
     & iflatteninglo0:iflatteninghi0,
     & iflatteninglo1:iflatteninghi1,
     & iflatteninglo2:iflatteninghi2)
      integer numslopes
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k,lvar
      do lvar = 0, numslopes - 1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         dw(i,j,k,lvar) = flattening(i,j,k)
     & * dw(i,j,k,lvar)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INCSOURCE(
     & prim
     & ,iprimlo0,iprimlo1,iprimlo2
     & ,iprimhi0,iprimhi1,iprimhi2
     & ,nprimcomp
     & ,source
     & ,isourcelo0,isourcelo1,isourcelo2
     & ,isourcehi0,isourcehi1,isourcehi2
     & ,nsourcecomp
     & ,scale
     & ,idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nprimcomp
      integer iprimlo0,iprimlo1,iprimlo2
      integer iprimhi0,iprimhi1,iprimhi2
      REAL*8 prim(
     & iprimlo0:iprimhi0,
     & iprimlo1:iprimhi1,
     & iprimlo2:iprimhi2,
     & 0:nprimcomp-1)
      integer nsourcecomp
      integer isourcelo0,isourcelo1,isourcelo2
      integer isourcehi0,isourcehi1,isourcehi2
      REAL*8 source(
     & isourcelo0:isourcehi0,
     & isourcelo1:isourcehi1,
     & isourcelo2:isourcehi2,
     & 0:nsourcecomp-1)
      REAL*8 scale
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer i,j,k, iv
      REAL*8 increment
      do iv = 0,nprimcomp - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         increment = scale*source(i,j,k, iv)
         prim(i,j,k, iv) = prim(i,j,k, iv) + increment
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine VLLIMITER(
     & slopeprim
     & ,islopeprimlo0,islopeprimlo1,islopeprimlo2
     & ,islopeprimhi0,islopeprimhi1,islopeprimhi2
     & ,nslopeprimcomp
     & ,slopeleft
     & ,islopeleftlo0,islopeleftlo1,islopeleftlo2
     & ,islopelefthi0,islopelefthi1,islopelefthi2
     & ,nslopeleftcomp
     & ,sloperigh
     & ,isloperighlo0,isloperighlo1,isloperighlo2
     & ,isloperighhi0,isloperighhi1,isloperighhi2
     & ,nsloperighcomp
     & ,idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nslopeprimcomp
      integer islopeprimlo0,islopeprimlo1,islopeprimlo2
      integer islopeprimhi0,islopeprimhi1,islopeprimhi2
      REAL*8 slopeprim(
     & islopeprimlo0:islopeprimhi0,
     & islopeprimlo1:islopeprimhi1,
     & islopeprimlo2:islopeprimhi2,
     & 0:nslopeprimcomp-1)
      integer nslopeleftcomp
      integer islopeleftlo0,islopeleftlo1,islopeleftlo2
      integer islopelefthi0,islopelefthi1,islopelefthi2
      REAL*8 slopeleft(
     & islopeleftlo0:islopelefthi0,
     & islopeleftlo1:islopelefthi1,
     & islopeleftlo2:islopelefthi2,
     & 0:nslopeleftcomp-1)
      integer nsloperighcomp
      integer isloperighlo0,isloperighlo1,isloperighlo2
      integer isloperighhi0,isloperighhi1,isloperighhi2
      REAL*8 sloperigh(
     & isloperighlo0:isloperighhi0,
     & isloperighlo1:isloperighhi1,
     & isloperighlo2:isloperighhi2,
     & 0:nsloperighcomp-1)
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer i,j,k, iv
      REAL*8 dql, dqr, dqlim
      do iv = 0,nslopeprimcomp - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         dql = slopeleft(i,j,k, iv)
         dqr = sloperigh(i,j,k, iv)
         dqlim = slopeprim(i,j,k, iv)
         call pointvllimiter(dqlim, dql, dqr)
         slopeprim(i,j,k,iv) = dqlim
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine POINTVLLIMITER(
     & dqlim
     & ,dql
     & ,dqr
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dqlim
      REAL*8 dql
      REAL*8 dqr
      REAL*8 dqc
      dqc = dqlim
      dqlim = min((2.0d0)*abs(dql),(2.0d0)*abs(dqr))
      dqlim = min(dqlim, abs(dqc))
      if (dql*dqr .lt. (0.0d0)) then
         dqlim = (0.0d0)
      else
         dqlim = dqlim*sign((1.0d0), dql)
         ch_flops=ch_flops+2
      endif
      ch_flops=ch_flops+8
      return
      end
        subroutine DIVUEDGE(
     & divu
     & ,idivulo0,idivulo1,idivulo2
     & ,idivuhi0,idivuhi1,idivuhi2
     & ,facedir
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1,idivulo2
      integer idivuhi0,idivuhi1,idivuhi2
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1,
     & idivulo2:idivuhi2)
      integer facedir
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer i,j,k,ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      if (haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         divu(i,j,k) = divu(i+ioff,j+joff,k+koff)
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         divu(i,j,k) = divu(i-ioff,j-joff,k-koff)
      enddo
      enddo
      enddo
      endif
      return
      end
        subroutine AVEFLUXTOFACE(
     & faceflux
     & ,ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
     & ,ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
     & ,ccflux
     & ,iccfluxlo0,iccfluxlo1,iccfluxlo2
     & ,iccfluxhi0,iccfluxhi1,iccfluxhi2
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
      integer ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
      integer ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
      REAL*8 faceflux(
     & ifacefluxlo0:ifacefluxhi0,
     & ifacefluxlo1:ifacefluxhi1,
     & ifacefluxlo2:ifacefluxhi2)
      integer iccfluxlo0,iccfluxlo1,iccfluxlo2
      integer iccfluxhi0,iccfluxhi1,iccfluxhi2
      REAL*8 ccflux(
     & iccfluxlo0:iccfluxhi0,
     & iccfluxlo1:iccfluxhi1,
     & iccfluxlo2:iccfluxhi2)
      integer facedir
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k,ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      faceflux(i,j,k) = (0.500d0)*(
     & ccflux(i ,j ,k ) +
     & ccflux(i - ioff,j - joff,k - koff))
      enddo
      enddo
      enddo
      return
      end
        subroutine DIVUONED(
     & divu
     & ,idivulo0,idivulo1,idivulo2
     & ,idivuhi0,idivuhi1,idivuhi2
     & ,velnorm
     & ,ivelnormlo0,ivelnormlo1,ivelnormlo2
     & ,ivelnormhi0,ivelnormhi1,ivelnormhi2
     & ,facedir
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1,idivulo2
      integer idivuhi0,idivuhi1,idivuhi2
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1,
     & idivulo2:idivuhi2)
      integer ivelnormlo0,ivelnormlo1,ivelnormlo2
      integer ivelnormhi0,ivelnormhi1,ivelnormhi2
      REAL*8 velnorm(
     & ivelnormlo0:ivelnormhi0,
     & ivelnormlo1:ivelnormhi1,
     & ivelnormlo2:ivelnormhi2)
      integer facedir
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i,j,k,ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      divu(i,j,k) =
     & velnorm(i ,j ,k ) -
     & velnorm(i - ioff,j - joff,k - koff)
      enddo
      enddo
      enddo
      return
      end
        subroutine DIVUTRAN(
     & divu
     & ,idivulo0,idivulo1,idivulo2
     & ,idivuhi0,idivuhi1,idivuhi2
     & ,slopevel
     & ,islopevello0,islopevello1,islopevello2
     & ,islopevelhi0,islopevelhi1,islopevelhi2
     & ,facedir
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivulo0,idivulo1,idivulo2
      integer idivuhi0,idivuhi1,idivuhi2
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1,
     & idivulo2:idivuhi2)
      integer islopevello0,islopevello1,islopevello2
      integer islopevelhi0,islopevelhi1,islopevelhi2
      REAL*8 slopevel(
     & islopevello0:islopevelhi0,
     & islopevello1:islopevelhi1,
     & islopevello2:islopevelhi2)
      integer facedir
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer i,j,k,ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      divu(i,j,k) = divu(i,j,k)+ (0.500d0)*(
     & slopevel(i ,j ,k ) +
     & slopevel(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      return
      end
        subroutine ARTVISC(
     & f
     & ,iflo0,iflo1,iflo2
     & ,ifhi0,ifhi1,ifhi2
     & ,nfcomp
     & ,u
     & ,iulo0,iulo1,iulo2
     & ,iuhi0,iuhi1,iuhi2
     & ,nucomp
     & ,divu
     & ,idivulo0,idivulo1,idivulo2
     & ,idivuhi0,idivuhi1,idivuhi2
     & ,coeff
     & ,idir
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,numcons
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1,
     & iflo2:ifhi2,
     & 0:nfcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2
      integer iuhi0,iuhi1,iuhi2
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & iulo2:iuhi2,
     & 0:nucomp-1)
      integer idivulo0,idivulo1,idivulo2
      integer idivuhi0,idivuhi1,idivuhi2
      REAL*8 divu(
     & idivulo0:idivuhi0,
     & idivulo1:idivuhi1,
     & idivulo2:idivuhi2)
      REAL*8 coeff
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer numcons
      REAL*8 dx
        integer i , j , k
        integer ioff, joff, koff
        integer iv
        REAL*8 fc,dv,s1,s2
        ioff = chf_id(0,idir)
        joff = chf_id(1,idir)
        koff = chf_id(2,idir)
        do iv = 0,numcons - 1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
            fc = f (i ,j ,k ,iv)
            dv = divu(i ,j ,k )
            s1 = u (i ,j ,k ,iv)
            s2 = u (i-ioff,j-joff,k-koff,iv)
            f(i,j,k,iv) = fc - coeff*max(-dv, 0.d0)*(s1-s2)
      enddo
      enddo
      enddo
        enddo
        return
        end
      subroutine UPDATE(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,state
     & ,istatelo0,istatelo1,istatelo2
     & ,istatehi0,istatehi1,istatehi2
     & ,nstatecomp
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,facedir
     & ,nconserved
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer nstatecomp
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & istatelo2:istatehi2,
     & 0:nstatecomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dtbydx
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         state(i,j,k,iv) = state(i,j,k,iv) -
     & dtbydx *
     & ( flux(i+ioff,j+joff,k+koff,iv)
     & - flux(i ,j ,k ,iv))
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine DIVERGEF(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,divf
     & ,idivflo0,idivflo1,idivflo2
     & ,idivfhi0,idivfhi1,idivfhi2
     & ,ndivfcomp
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,facedir
     & ,nconserved
     & ,dx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer ndivfcomp
      integer idivflo0,idivflo1,idivflo2
      integer idivfhi0,idivfhi1,idivfhi2
      REAL*8 divf(
     & idivflo0:idivfhi0,
     & idivflo1:idivfhi1,
     & idivflo2:idivfhi2,
     & 0:ndivfcomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dx
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         divf(i,j,k,iv) = divf(i,j,k,iv) +
     & (flux(i+ioff,j+joff,k+koff,iv)
     & -flux(i ,j ,k ,iv))/dx
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine REGUPDATE(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,consstate
     & ,iconsstatelo0,iconsstatelo1,iconsstatelo2
     & ,iconsstatehi0,iconsstatehi1,iconsstatehi2
     & ,nconsstatecomp
     & ,divf
     & ,idivflo0,idivflo1,idivflo2
     & ,idivfhi0,idivfhi1,idivfhi2
     & ,ndivfcomp
     & ,nconserved
     & ,dt
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer nconsstatecomp
      integer iconsstatelo0,iconsstatelo1,iconsstatelo2
      integer iconsstatehi0,iconsstatehi1,iconsstatehi2
      REAL*8 consstate(
     & iconsstatelo0:iconsstatehi0,
     & iconsstatelo1:iconsstatehi1,
     & iconsstatelo2:iconsstatehi2,
     & 0:nconsstatecomp-1)
      integer ndivfcomp
      integer idivflo0,idivflo1,idivflo2
      integer idivfhi0,idivfhi1,idivfhi2
      REAL*8 divf(
     & idivflo0:idivfhi0,
     & idivflo1:idivfhi1,
     & idivflo2:idivfhi2,
     & 0:ndivfcomp-1)
      integer nconserved
      REAL*8 dt
      integer i, j, k
      integer iv
      do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         consstate(i,j,k,iv) = consstate(i,j,k,iv)
     & - dt*divf(i,j,k,iv)
      enddo
      enddo
      enddo
      enddo
      return
      end
