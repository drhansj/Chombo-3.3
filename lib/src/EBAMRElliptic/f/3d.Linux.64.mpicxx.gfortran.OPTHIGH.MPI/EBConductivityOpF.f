      subroutine CONDUCTIVITYGSRB(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,relcoef
     & ,irelcoeflo0,irelcoeflo1,irelcoeflo2
     & ,irelcoefhi0,irelcoefhi1,irelcoefhi2
     & ,acoef
     & ,iacoeflo0,iacoeflo1,iacoeflo2
     & ,iacoefhi0,iacoefhi1,iacoefhi2
     & ,b0
     & ,ib0lo0,ib0lo1,ib0lo2
     & ,ib0hi0,ib0hi1,ib0hi2
     & ,b1
     & ,ib1lo0,ib1lo1,ib1lo2
     & ,ib1hi0,ib1hi1,ib1hi2
     & ,b2
     & ,ib2lo0,ib2lo1,ib2lo2
     & ,ib2hi0,ib2hi1,ib2hi2
     & ,alpha
     & ,beta
     & ,dx
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,redBlack
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2)
      integer irelcoeflo0,irelcoeflo1,irelcoeflo2
      integer irelcoefhi0,irelcoefhi1,irelcoefhi2
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1,
     & irelcoeflo2:irelcoefhi2)
      integer iacoeflo0,iacoeflo1,iacoeflo2
      integer iacoefhi0,iacoefhi1,iacoefhi2
      REAL*8 acoef(
     & iacoeflo0:iacoefhi0,
     & iacoeflo1:iacoefhi1,
     & iacoeflo2:iacoefhi2)
      integer ib0lo0,ib0lo1,ib0lo2
      integer ib0hi0,ib0hi1,ib0hi2
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1,
     & ib0lo2:ib0hi2)
      integer ib1lo0,ib1lo1,ib1lo2
      integer ib1hi0,ib1hi1,ib1hi2
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1,
     & ib1lo2:ib1hi2)
      integer ib2lo0,ib2lo1,ib2lo2
      integer ib2hi0,ib2hi1,ib2hi2
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1,
     & ib2lo2:ib2hi2)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer redBlack
      integer i,j,k
      REAL*8 laplphi, dx0
      integer imin,imax,indtot
      dx0 = beta/(dx * dx)
      do k=iregionlo2, iregionhi2
         do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               laplphi =
     & (b0(i+1,j ,k )*(phi(i+1,j ,k ) - phi(i ,j ,k ))
     & - b0(i ,j ,k )*(phi(i ,j ,k ) - phi(i-1,j ,k )))*dx0
     & +(b1(i ,j+1,k )*(phi(i ,j+1,k ) - phi(i ,j ,k ))
     & - b1(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j-1,k )))*dx0
     & +(b2(i ,j ,k+1)*(phi(i ,j ,k+1) - phi(i ,j ,k ))
     & - b2(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j ,k-1)))*dx0
               laplphi = laplphi + alpha * acoef(i,j,k) * phi(i,j,k)
               phi(i,j,k) = phi(i,j,k) + relcoef(i,j,k)*(rhs(i,j,k) - la
     &plphi)
            enddo
         enddo
      enddo
      return
      end
        subroutine EBCOREGAPPLYDOMAINFLUX(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,faceflux
     & ,ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
     & ,ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
     & ,bc
     & ,ibclo0,ibclo1,ibclo2
     & ,ibchi0,ibchi1,ibchi2
     & ,dx
     & ,side
     & ,idir
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
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
      integer ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
      integer ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
      REAL*8 faceflux(
     & ifacefluxlo0:ifacefluxhi0,
     & ifacefluxlo1:ifacefluxhi1,
     & ifacefluxlo2:ifacefluxhi2)
      integer ibclo0,ibclo1,ibclo2
      integer ibchi0,ibchi1,ibchi2
      REAL*8 bc(
     & ibclo0:ibchi0,
     & ibclo1:ibchi1,
     & ibclo2:ibchi2)
      REAL*8 dx
      integer side
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k, ioff,joff,koff
        REAL*8 scaledflux, tol
        tol = 1.0d-15
        ioff = chf_id(0,idir)
        joff = chf_id(1,idir)
        koff = chf_id(2,idir)
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        if (side.eq.1) then
           if (ABS(bc(i-ioff,j-joff,k-koff)) .GT. tol) then
              scaledflux = faceflux(i-ioff,j-joff,k-koff)/bc(i-ioff,j-jo
     &ff,k-koff)
              phi(i,j,k) = phi(i-ioff,j-joff,k-koff) + scaledflux*dx
           else
              phi(i,j,k) = phi(i-ioff,j-joff,k-koff)
           endif
        else
           if(ABS(bc(i+ioff,j+joff,k+koff)) .GT. tol) then
              scaledflux = faceflux(i+ioff,j+joff,k+koff)/bc(i+ioff,j+jo
     &ff,k+koff)
              phi(i,j,k) = phi(i+ioff,j+joff,k+koff) - scaledflux*dx
           else
              phi(i,j,k) = phi(i+ioff,j+joff,k+koff)
           endif
        endif
      enddo
      enddo
      enddo
        return
        end
      subroutine APPLYOPEBCONDNOBCS(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,aco
     & ,iacolo0,iacolo1,iacolo2
     & ,iacohi0,iacohi1,iacohi2
     & ,b0
     & ,ib0lo0,ib0lo1,ib0lo2
     & ,ib0hi0,ib0hi1,ib0hi2
     & ,b1
     & ,ib1lo0,ib1lo1,ib1lo2
     & ,ib1hi0,ib1hi1,ib1hi2
     & ,b2
     & ,ib2lo0,ib2lo1,ib2lo2
     & ,ib2hi0,ib2hi1,ib2hi2
     & ,dx
     & ,alpha
     & ,beta
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iopphidirlo0,iopphidirlo1,iopphidirlo2
      integer iopphidirhi0,iopphidirhi1,iopphidirhi2
      REAL*8 opphidir(
     & iopphidirlo0:iopphidirhi0,
     & iopphidirlo1:iopphidirhi1,
     & iopphidirlo2:iopphidirhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer iacolo0,iacolo1,iacolo2
      integer iacohi0,iacohi1,iacohi2
      REAL*8 aco(
     & iacolo0:iacohi0,
     & iacolo1:iacohi1,
     & iacolo2:iacohi2)
      integer ib0lo0,ib0lo1,ib0lo2
      integer ib0hi0,ib0hi1,ib0hi2
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1,
     & ib0lo2:ib0hi2)
      integer ib1lo0,ib1lo1,ib1lo2
      integer ib1hi0,ib1hi1,ib1hi2
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1,
     & ib1lo2:ib1hi2)
      integer ib2lo0,ib2lo1,ib2lo2
      integer ib2hi0,ib2hi1,ib2hi2
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1,
     & ib2lo2:ib2hi2)
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      REAL*8 laplphi, dxinvsq
      dxinvsq = (1.0d0)/(dx * dx)
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & (b0(i+1,j ,k )*(phi(i+1,j ,k ) - phi(i ,j ,k ))
     & - b0(i ,j ,k )*(phi(i ,j ,k ) - phi(i-1,j ,k )))*dxinvsq
     & +(b1(i ,j+1,k )*(phi(i ,j+1,k ) - phi(i ,j ,k ))
     & - b1(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j-1,k )))*dxinvsq
     & +(b2(i ,j ,k+1)*(phi(i ,j ,k+1) - phi(i ,j ,k ))
     & - b2(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j ,k-1)))*dxinvsq
      opphidir(i,j,k) = alpha*aco(i,j,k)*phi(i,j,k) + beta*laplphi
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*(4 + 6*3)
      return
      end
      subroutine GSCOLOREBCONDOP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,relco
     & ,irelcolo0,irelcolo1,irelcolo2
     & ,irelcohi0,irelcohi1,irelcohi2
     & ,aco
     & ,iacolo0,iacolo1,iacolo2
     & ,iacohi0,iacohi1,iacohi2
     & ,b0
     & ,ib0lo0,ib0lo1,ib0lo2
     & ,ib0hi0,ib0hi1,ib0hi2
     & ,b1
     & ,ib1lo0,ib1lo1,ib1lo2
     & ,ib1hi0,ib1hi1,ib1hi2
     & ,b2
     & ,ib2lo0,ib2lo1,ib2lo2
     & ,ib2hi0,ib2hi1,ib2hi2
     & ,dx
     & ,alpha
     & ,beta
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2)
      integer irelcolo0,irelcolo1,irelcolo2
      integer irelcohi0,irelcohi1,irelcohi2
      REAL*8 relco(
     & irelcolo0:irelcohi0,
     & irelcolo1:irelcohi1,
     & irelcolo2:irelcohi2)
      integer iacolo0,iacolo1,iacolo2
      integer iacohi0,iacohi1,iacohi2
      REAL*8 aco(
     & iacolo0:iacohi0,
     & iacolo1:iacohi1,
     & iacolo2:iacohi2)
      integer ib0lo0,ib0lo1,ib0lo2
      integer ib0hi0,ib0hi1,ib0hi2
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1,
     & ib0lo2:ib0hi2)
      integer ib1lo0,ib1lo1,ib1lo2
      integer ib1hi0,ib1hi1,ib1hi2
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1,
     & ib1lo2:ib1hi2)
      integer ib2lo0,ib2lo1,ib2lo2
      integer ib2hi0,ib2hi1,ib2hi2
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1,
     & ib2lo2:ib2hi2)
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k,ncolors
      REAL*8 laplphi, dxinvsq, opphipt, rhspt, relcopt, phipt
      dxinvsq = (1.0d0)/(dx * dx)
      do k = iboxlo2,iboxhi2,2
      do j = iboxlo1,iboxhi1,2
      do i = iboxlo0,iboxhi0,2
      laplphi =
     & (b0(i+1,j ,k )*(phi(i+1,j ,k ) - phi(i ,j ,k ))
     & - b0(i ,j ,k )*(phi(i ,j ,k ) - phi(i-1,j ,k )))*dxinvsq
     & +(b1(i ,j+1,k )*(phi(i ,j+1,k ) - phi(i ,j ,k ))
     & - b1(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j-1,k )))*dxinvsq
     & +(b2(i ,j ,k+1)*(phi(i ,j ,k+1) - phi(i ,j ,k ))
     & - b2(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j ,k-1)))*dxinvsq
      opphipt = alpha*aco(i,j,k)*phi(i,j,k) + beta*laplphi
      phipt = phi(i,j,k)
      rhspt = rhs(i,j,k)
      relcopt = relco(i,j,k)
      phi(i,j,k) = phipt + relcopt*(rhspt - opphipt)
      enddo
      enddo
      enddo
      ncolors = 2 *2 *2
      ch_flops=ch_flops+((iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(ibox
     &hi2- iboxlo2+1)*(4 + 6*3))/ncolors
      return
      end
      subroutine CONDUCTIVITYINPLACE(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,b0
     & ,ib0lo0,ib0lo1,ib0lo2
     & ,ib0hi0,ib0hi1,ib0hi2
     & ,b1
     & ,ib1lo0,ib1lo1,ib1lo2
     & ,ib1hi0,ib1hi1,ib1hi2
     & ,b2
     & ,ib2lo0,ib2lo1,ib2lo2
     & ,ib2hi0,ib2hi1,ib2hi2
     & ,beta
     & ,dx
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iopphidirlo0,iopphidirlo1,iopphidirlo2
      integer iopphidirhi0,iopphidirhi1,iopphidirhi2
      REAL*8 opphidir(
     & iopphidirlo0:iopphidirhi0,
     & iopphidirlo1:iopphidirhi1,
     & iopphidirlo2:iopphidirhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer ib0lo0,ib0lo1,ib0lo2
      integer ib0hi0,ib0hi1,ib0hi2
      REAL*8 b0(
     & ib0lo0:ib0hi0,
     & ib0lo1:ib0hi1,
     & ib0lo2:ib0hi2)
      integer ib1lo0,ib1lo1,ib1lo2
      integer ib1hi0,ib1hi1,ib1hi2
      REAL*8 b1(
     & ib1lo0:ib1hi0,
     & ib1lo1:ib1hi1,
     & ib1lo2:ib1hi2)
      integer ib2lo0,ib2lo1,ib2lo2
      integer ib2hi0,ib2hi1,ib2hi2
      REAL*8 b2(
     & ib2lo0:ib2hi0,
     & ib2lo1:ib2hi1,
     & ib2lo2:ib2hi2)
      REAL*8 beta
      REAL*8 dx
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      REAL*8 laplphi, dx0
      dx0 = beta/(dx * dx)
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & (b0(i+1,j ,k )*(phi(i+1,j ,k ) - phi(i ,j ,k ))
     & - b0(i ,j ,k )*(phi(i ,j ,k ) - phi(i-1,j ,k )))*dx0
     & +(b1(i ,j+1,k )*(phi(i ,j+1,k ) - phi(i ,j ,k ))
     & - b1(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j-1,k )))*dx0
     & +(b2(i ,j ,k+1)*(phi(i ,j ,k+1) - phi(i ,j ,k ))
     & - b2(i ,j ,k )*(phi(i ,j ,k ) - phi(i ,j ,k-1)))*dx0
      opphidir(i,j,k) = opphidir(i,j,k) + laplphi
      enddo
      enddo
      enddo
      return
      end
      subroutine INCRAPPLYEBCO(
     & lhs
     & ,ilhslo0,ilhslo1,ilhslo2
     & ,ilhshi0,ilhshi1,ilhshi2
     & ,interiorflux
     & ,iinteriorfluxlo0,iinteriorfluxlo1,iinteriorfluxlo2
     & ,iinteriorfluxhi0,iinteriorfluxhi1,iinteriorfluxhi2
     & ,domainfluxlo
     & ,idomainfluxlolo0,idomainfluxlolo1,idomainfluxlolo2
     & ,idomainfluxlohi0,idomainfluxlohi1,idomainfluxlohi2
     & ,domainfluxhi
     & ,idomainfluxhilo0,idomainfluxhilo1,idomainfluxhilo2
     & ,idomainfluxhihi0,idomainfluxhihi1,idomainfluxhihi2
     & ,beta
     & ,dx
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & ,haslo
     & ,hashi
     & ,facedir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1,
     & ilhslo2:ilhshi2)
      integer iinteriorfluxlo0,iinteriorfluxlo1,iinteriorfluxlo2
      integer iinteriorfluxhi0,iinteriorfluxhi1,iinteriorfluxhi2
      REAL*8 interiorflux(
     & iinteriorfluxlo0:iinteriorfluxhi0,
     & iinteriorfluxlo1:iinteriorfluxhi1,
     & iinteriorfluxlo2:iinteriorfluxhi2)
      integer idomainfluxlolo0,idomainfluxlolo1,idomainfluxlolo2
      integer idomainfluxlohi0,idomainfluxlohi1,idomainfluxlohi2
      REAL*8 domainfluxlo(
     & idomainfluxlolo0:idomainfluxlohi0,
     & idomainfluxlolo1:idomainfluxlohi1,
     & idomainfluxlolo2:idomainfluxlohi2)
      integer idomainfluxhilo0,idomainfluxhilo1,idomainfluxhilo2
      integer idomainfluxhihi0,idomainfluxhihi1,idomainfluxhihi2
      REAL*8 domainfluxhi(
     & idomainfluxhilo0:idomainfluxhihi0,
     & idomainfluxhilo1:idomainfluxhihi1,
     & idomainfluxhilo2:idomainfluxhihi2)
      REAL*8 beta
      REAL*8 dx
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer haslo
      integer hashi
      integer facedir
      integer ii,i,jj,j,kk,k
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      lhs(i,j,k) = lhs(i,j,k)
     $ +beta*
     $ (interiorflux(i+ii,j+jj,k+kk)
     $ -interiorflux(i ,j ,k ))/dx
      enddo
      enddo
      enddo
      if(haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         lhs(i,j,k) = lhs(i,j,k)
     $ + beta*
     $ (interiorflux(i+ii,j+jj,k+kk)
     $ -domainfluxlo(i ,j ,k ))/dx
      enddo
      enddo
      enddo
      endif
      if(hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         lhs(i,j,k) = lhs(i,j,k)
     $ + beta*
     $ (domainfluxhi(i+ii,j+jj,k+kk)
     $ -interiorflux(i ,j ,k ))/dx
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine DECRINVRELCOEFEBCO(
     & relcoef
     & ,irelcoeflo0,irelcoeflo1,irelcoeflo2
     & ,irelcoefhi0,irelcoefhi1,irelcoefhi2
     & ,bcoef
     & ,ibcoeflo0,ibcoeflo1,ibcoeflo2
     & ,ibcoefhi0,ibcoefhi1,ibcoefhi2
     & ,beta
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,dx
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer irelcoeflo0,irelcoeflo1,irelcoeflo2
      integer irelcoefhi0,irelcoefhi1,irelcoefhi2
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1,
     & irelcoeflo2:irelcoefhi2)
      integer ibcoeflo0,ibcoeflo1,ibcoeflo2
      integer ibcoefhi0,ibcoefhi1,ibcoefhi2
      REAL*8 bcoef(
     & ibcoeflo0:ibcoefhi0,
     & ibcoeflo1:ibcoefhi1,
     & ibcoeflo2:ibcoefhi2)
      REAL*8 beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      REAL*8 dx
      integer idir
      integer i,j,k
      integer ii,jj,kk
      ii = chf_id(idir, 0)
      jj = chf_id(idir, 1)
      kk = chf_id(idir, 2)
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      relcoef(i,j,k) = relcoef(i,j,k)
     $ - beta*(
     $ bcoef(i+ii,j+jj,k+kk) +
     $ bcoef(i ,j ,k ))/(dx*dx)
      enddo
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAEBCO(
     & lambda
     & ,ilambdalo0,ilambdalo1,ilambdalo2
     & ,ilambdahi0,ilambdahi1,ilambdahi2
     & ,safety
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & ilambdalo2:ilambdahi2)
      REAL*8 safety
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         lambda(i,j,k) = safety/(lambda(i,j,k))
      enddo
      enddo
      enddo
      return
      end
      subroutine GETFLUXEBCO(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,bcoef
     & ,ibcoeflo0,ibcoeflo1,ibcoeflo2
     & ,ibcoefhi0,ibcoefhi1,ibcoefhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
     & ,iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
     & ,dx
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2)
      integer ibcoeflo0,ibcoeflo1,ibcoeflo2
      integer ibcoefhi0,ibcoefhi1,ibcoefhi2
      REAL*8 bcoef(
     & ibcoeflo0:ibcoefhi0,
     & ibcoeflo1:ibcoefhi1,
     & ibcoeflo2:ibcoefhi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
      integer iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
      REAL*8 dx
      integer idir
      integer i,j,k
      integer ioff,joff,koff
      REAL*8 oneoverdx
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = iopphiboxlo2,iopphiboxhi2
      do j = iopphiboxlo1,iopphiboxhi1
      do i = iopphiboxlo0,iopphiboxhi0
      oneoverdx = bcoef(i,j,k)/dx
      flux(i,j,k) =
     $ oneoverdx*(
     $ phi(i ,j ,k ) -
     $ phi(i-ioff,j-joff,k-koff) )
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBEBCO(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,lph
     & ,ilphlo0,ilphlo1,ilphlo2
     & ,ilphhi0,ilphhi1,ilphhi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,lam
     & ,ilamlo0,ilamlo1,ilamlo2
     & ,ilamhi0,ilamhi1,ilamhi2
     & ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     & ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer ilphlo0,ilphlo1,ilphlo2
      integer ilphhi0,ilphhi1,ilphhi2
      REAL*8 lph(
     & ilphlo0:ilphhi0,
     & ilphlo1:ilphhi1,
     & ilphlo2:ilphhi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2)
      integer ilamlo0,ilamlo1,ilamlo2
      integer ilamhi0,ilamhi1,ilamhi2
      REAL*8 lam(
     & ilamlo0:ilamhi0,
     & ilamlo1:ilamhi1,
     & ilamlo2:ilamhi2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer i,j,k
      REAL*8 phio, lamo, rhso, lpho
      do k = icoloredboxlo2,icoloredboxhi2,2
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2
         phio = phi(i,j,k)
         lamo = lam(i,j,k)
         rhso = rhs(i,j,k)
         lpho = lph(i,j,k)
         phi(i,j,k) =
     $ phi(i,j,k) +
     $ lam(i,j,k)*(
     $ rhs(i,j,k) -
     $ lph(i,j,k))
      enddo
      enddo
      enddo
      return
      end
