      subroutine EBAMRPDOTPROD(
     & sum
     & ,aone
     & ,iaonelo0,iaonelo1,iaonelo2
     & ,iaonehi0,iaonehi1,iaonehi2
     & ,atwo
     & ,iatwolo0,iatwolo1,iatwolo2
     & ,iatwohi0,iatwohi1,iatwohi2
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 sum
      integer iaonelo0,iaonelo1,iaonelo2
      integer iaonehi0,iaonehi1,iaonehi2
      REAL*8 aone(
     & iaonelo0:iaonehi0,
     & iaonelo1:iaonehi1,
     & iaonelo2:iaonehi2)
      integer iatwolo0,iatwolo1,iatwolo2
      integer iatwohi0,iatwohi1,iatwohi2
      REAL*8 atwo(
     & iatwolo0:iatwohi0,
     & iatwolo1:iatwohi1,
     & iatwolo2:iatwohi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sum = sum + aone(i,j,k)*atwo(i,j,k)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*2
      return
      end
      subroutine GETINVDIAGRHS(
     & lhs
     & ,ilhslo0,ilhslo1,ilhslo2
     & ,ilhshi0,ilhshi1,ilhshi2
     & ,nlhscomp
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,nrhscomp
     & ,scale
     & ,ncomp
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlhscomp
      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1,
     & ilhslo2:ilhshi2,
     & 0:nlhscomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & 0:nrhscomp-1)
      REAL*8 scale
      integer ncomp
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k, ivar
      do ivar = 0, ncomp-1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         lhs(i,j,k, ivar) = scale*rhs(i,j,k, ivar)
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*ncomp
      return
      end
      subroutine MAXNORM(
     & m
     & ,vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 m
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        m=max(m,abs(vel(i,j,k)))
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*2
      return
      end
      subroutine MAXNORMMASK(
     & m
     & ,vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,mask
     & ,imasklo0,imasklo1,imasklo2
     & ,imaskhi0,imaskhi1,imaskhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 m
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer imasklo0,imasklo1,imasklo2
      integer imaskhi0,imaskhi1,imaskhi2
      integer*1 mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1,
     & imasklo2:imaskhi2)
      integer i,j,k
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         m=max(m,abs(vel(i,j,k)*mask(i,j,k)))
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*2
      return
      end
      subroutine AMRPZEROSUB(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,ioverlayboxlo0,ioverlayboxlo1,ioverlayboxlo2
     & ,ioverlayboxhi0,ioverlayboxhi1,ioverlayboxhi2
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer ioverlayboxlo0,ioverlayboxlo1,ioverlayboxlo2
      integer ioverlayboxhi0,ioverlayboxhi1,ioverlayboxhi2
      integer ncomp
      integer i,j,k
      integer ivar
      do ivar = 0, ncomp-1
      do k = ioverlayBoxlo2,ioverlayBoxhi2
      do j = ioverlayBoxlo1,ioverlayBoxhi1
      do i = ioverlayBoxlo0,ioverlayBoxhi0
         phi(i,j,k, ivar) = (0.0d0)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine DOALLREGULARMULTICOLOR(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,weight
     & ,alpha
     & ,beta
     & ,dx
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
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2)
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer i,j,k
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2
        laplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1)
     $ -(2.0d0)*phi(i ,j ,k ))*dx2
        laplphi = laplphi + alpha * phi(i,j,k)
        phi(i,j,k) = phi(i,j,k) +
     & weight*(rhs(i,j,k) - laplphi)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(icoloredBoxhi0- icoloredBoxlo0+1)*(icoloredBoxh
     &i1- icoloredBoxlo1+1)*(icoloredBoxhi2- icoloredBoxlo2+1)/2*(4*3 + 
     &5)
      return
      end
      subroutine DOALLREGULARUPDATE(
     & phinew
     & ,iphinewlo0,iphinewlo1,iphinewlo2
     & ,iphinewhi0,iphinewhi1,iphinewhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,weight
     & ,alpha
     & ,beta
     & ,dx
     & ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     & ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphinewlo0,iphinewlo1,iphinewlo2
      integer iphinewhi0,iphinewhi1,iphinewhi2
      REAL*8 phinew(
     & iphinewlo0:iphinewhi0,
     & iphinewlo1:iphinewhi1,
     & iphinewlo2:iphinewhi2)
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
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer i,j,k
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2
        laplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1)
     $ -(2.0d0)*phi(i ,j ,k ))*dx2
        laplphi = laplphi + alpha * phi(i,j,k)
        phinew(i,j,k) = phi(i,j,k) +
     & weight*(rhs(i,j,k) - laplphi)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(icoloredBoxhi0- icoloredBoxlo0+1)*(icoloredBoxh
     &i1- icoloredBoxlo1+1)*(icoloredBoxhi2- icoloredBoxlo2+1)/2*(4*3 + 
     &5)
      return
      end
      subroutine DOALLREGULARGSRB(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,weight
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
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer redBlack
      integer i,j,k
      integer imin,imax,indtot
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do k=iregionlo2, iregionhi2
         do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               laplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1)
     $ -(2.0d0)*phi(i ,j ,k ))*dx2
               laplphi = laplphi + alpha * phi(i,j,k)
               phi(i,j,k) = phi(i,j,k) +
     & weight*(rhs(i,j,k) - laplphi)
            enddo
         enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*(iregionhi2- iregionlo2+1)/2*(4*3 + 5)
      return
      end
      subroutine SLOWGSRBEBAMRPO(
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
      subroutine DOALLREGULARJACOBI(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,weight
     & ,alpha
     & ,beta
     & ,dx
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
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        laplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1)
     $ -(2.0d0)*phi(i ,j ,k ))*dx2
        laplphi = laplphi + alpha * phi(i,j,k)
        phi(i,j,k) = phi(i,j,k) +
     & weight*(rhs(i,j,k) - laplphi)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*(4*3 + 5)
      return
      end
      subroutine UNDOREGULARGS(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,weight
     & ,alpha
     & ,beta
     & ,dx
     & ,iv
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
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:2)
      integer iv(0:2)
      integer i,j,k
      REAL*8 sublaplphi, dx0,dx1,dx2
      REAL*8 bigk, sumtwooverh2, numer
      i = iv(0)
      j = iv(1)
      k = iv(2)
      dx0 = (1.0d0)/(dx(0) * dx(0))
      dx1 = (1.0d0)/(dx(1) * dx(1))
      dx2 = (1.0d0)/(dx(2) * dx(2))
      sumtwooverh2 = (2.0d0)*dx0 + (2.0d0)*dx1 + (2.0d0)*dx2
      bigk = (1.0d0) + weight*alpha - beta*weight*sumtwooverh2
      sublaplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1))*dx2
      numer = phi(i,j,k)
     $ + weight*rhs(i,j,k) - weight*beta*sublaplphi
      phi(i,j,k) = numer/bigk
      return
      end
        subroutine REGAPPLYDOMAINFLUX_INPLACE(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,faceflux
     & ,ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
     & ,ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
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
      REAL*8 dx
      integer side
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k, ioff,joff,koff
        ioff = chf_id(0,idir)
                  joff = chf_id(1,idir)
                  koff = chf_id(2,idir)
        if (side.eq.1) then
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          phi(i,j,k) =
     & phi( i-ioff,j-joff,k-koff) +
     & faceflux(i-ioff,j-joff,k-koff)*dx
      enddo
      enddo
      enddo
        else
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          phi(i,j,k) =
     & phi( i+ioff,j+joff,k+koff) -
     & faceflux(i+ioff,j+joff,k+koff)*dx
      enddo
      enddo
      enddo
        endif
        ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(ibo
     &xhi2- iboxlo2+1)*2
        return
        end
      subroutine MVOPERATORLAP(
     & lph
     & ,ilphlo0,ilphlo1,ilphlo2
     & ,ilphhi0,ilphhi1,ilphhi2
     & ,nlphcomp
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,beta
     & ,ncomps
     & ,dx
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlphcomp
      integer ilphlo0,ilphlo1,ilphlo2
      integer ilphhi0,ilphhi1,ilphhi2
      REAL*8 lph(
     & ilphlo0:ilphhi0,
     & ilphlo1:ilphhi1,
     & ilphlo2:ilphhi2,
     & 0:nlphcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      REAL*8 beta
      integer ncomps
      REAL*8 dx(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k, icomp
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do icomp = 0, ncomps-1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & ( phi(i+1,j ,k , icomp)
     & + phi(i-1,j ,k , icomp)
     $ -(2.0d0)*phi(i ,j ,k , icomp))*dx0
     $ +( phi(i ,j+1,k , icomp)
     & + phi(i ,j-1,k , icomp)
     $ -(2.0d0)*phi(i ,j ,k , icomp))*dx1
     $ +( phi(i ,j ,k+1, icomp)
     & + phi(i ,j ,k-1, icomp)
     $ -(2.0d0)*phi(i ,j ,k , icomp))*dx2
      lph(i,j,k,icomp) = laplphi
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+ncomps*(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)
     &*(iboxhi2- iboxlo2+1)*(4*3 + 1)
      return
      end
      subroutine REGINCRLAPLACIAN_INPLACE(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
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
      REAL*8 beta
      REAL*8 dx(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      REAL*8 laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      laplphi =
     & ( phi(i+1,j ,k )
     & + phi(i-1,j ,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx0
     $ +( phi(i ,j+1,k )
     & + phi(i ,j-1,k )
     $ -(2.0d0)*phi(i ,j ,k ))*dx1
     $ +( phi(i ,j ,k+1)
     & + phi(i ,j ,k-1)
     $ -(2.0d0)*phi(i ,j ,k ))*dx2
      opphidir(i,j,k) = opphidir(i,j,k) + laplphi
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iboxhi0- iboxlo0+1)*(iboxhi1- iboxlo1+1)*(iboxh
     &i2- iboxlo2+1)*(4*3 + 1)
      return
      end
      subroutine REGGET1DLAPLACIAN(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,dx
     & ,beta
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
      REAL*8 dx
      REAL*8 beta
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
      REAL*8 bdx2
      bdx2 = beta/dx/dx
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      opphidir(i,j,k) =
     & bdx2 *
     $ ( phi(i+ioff,j+joff,k+koff)
     & -(2.0d0)*phi(i ,j ,k )
     & + phi(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
         opphidir(i,j,k) =
     & bdx2 *
     $ ( phi(i+ioff,j+joff,k+koff)
     & - phi(i ,j ,k ))
      enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
         opphidir(i,j,k) =
     & bdx2 *
     $ ( phi(i ,j ,k )
     & - phi(i-ioff,j-joff,k-koff))
      enddo
      enddo
      enddo
      endif
      return
      end
        subroutine REGAPPLYDOMAINFLUX(
     & opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,faceflux
     & ,ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
     & ,ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
     & ,dx
     & ,beta
     & ,idir
     & ,side
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
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
      integer ifacefluxlo0,ifacefluxlo1,ifacefluxlo2
      integer ifacefluxhi0,ifacefluxhi1,ifacefluxhi2
      REAL*8 faceflux(
     & ifacefluxlo0:ifacefluxhi0,
     & ifacefluxlo1:ifacefluxhi1,
     & ifacefluxlo2:ifacefluxhi2)
      REAL*8 dx
      REAL*8 beta
      integer idir
      integer side
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
        integer i,j,k
        REAL*8 idx
        idx = (1.0d0)/dx
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
        opphidir(i,j,k) = idx * side *
     & (beta*faceflux(i,j,k) - opphidir(i,j,k))
      enddo
      enddo
      enddo
        return
        end
      subroutine REGSUMLAPLACIAN(
     & opphi
     & ,iopphilo0,iopphilo1,iopphilo2
     & ,iopphihi0,iopphihi1,iopphihi2
     & ,opphidir
     & ,iopphidirlo0,iopphidirlo1,iopphidirlo2
     & ,iopphidirhi0,iopphidirhi1,iopphidirhi2
     & ,iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
     & ,iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iopphilo0,iopphilo1,iopphilo2
      integer iopphihi0,iopphihi1,iopphihi2
      REAL*8 opphi(
     & iopphilo0:iopphihi0,
     & iopphilo1:iopphihi1,
     & iopphilo2:iopphihi2)
      integer iopphidirlo0,iopphidirlo1,iopphidirlo2
      integer iopphidirhi0,iopphidirhi1,iopphidirhi2
      REAL*8 opphidir(
     & iopphidirlo0:iopphidirhi0,
     & iopphidirlo1:iopphidirhi1,
     & iopphidirlo2:iopphidirhi2)
      integer iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
      integer iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
        integer i,j,k
      do k = iopphiboxlo2,iopphiboxhi2
      do j = iopphiboxlo1,iopphiboxhi1
      do i = iopphiboxlo0,iopphiboxhi0
        opphi(i,j,k) =
     $ opphi( i,j,k) +
     $ opphidir(i,j,k)
      enddo
      enddo
      enddo
        return
        end
        subroutine REGMULTICOLORGS(
     & newphi
     & ,inewphilo0,inewphilo1,inewphilo2
     & ,inewphihi0,inewphihi1,inewphihi2
     & ,weight
     & ,resid
     & ,iresidlo0,iresidlo1,iresidlo2
     & ,iresidhi0,iresidhi1,iresidhi2
     & ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     & ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer inewphilo0,inewphilo1,inewphilo2
      integer inewphihi0,inewphihi1,inewphihi2
      REAL*8 newphi(
     & inewphilo0:inewphihi0,
     & inewphilo1:inewphihi1,
     & inewphilo2:inewphihi2)
      REAL*8 weight
      integer iresidlo0,iresidlo1,iresidlo2
      integer iresidhi0,iresidhi1,iresidhi2
      REAL*8 resid(
     & iresidlo0:iresidhi0,
     & iresidlo1:iresidhi1,
     & iresidlo2:iresidhi2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
        integer i,j,k
      do k = icoloredboxlo2,icoloredboxhi2,2
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2
          newphi(i,j,k) = newphi(i,j,k)
     & + weight * resid(i,j,k)
      enddo
      enddo
      enddo
        return
        end
      subroutine REGGSRB(
     & newphi
     & ,inewphilo0,inewphilo1,inewphilo2
     & ,inewphihi0,inewphihi1,inewphihi2
     & ,resid
     & ,iresidlo0,iresidlo1,iresidlo2
     & ,iresidhi0,iresidhi1,iresidhi2
     & ,weight
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,color
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer inewphilo0,inewphilo1,inewphilo2
      integer inewphihi0,inewphihi1,inewphihi2
      REAL*8 newphi(
     & inewphilo0:inewphihi0,
     & inewphilo1:inewphihi1,
     & inewphilo2:inewphihi2)
      integer iresidlo0,iresidlo1,iresidlo2
      integer iresidhi0,iresidhi1,iresidhi2
      REAL*8 resid(
     & iresidlo0:iresidhi0,
     & iresidlo1:iresidhi1,
     & iresidlo2:iresidhi2)
      REAL*8 weight
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer color
      integer i,j,k
      integer imin,imax,indtot
      do k=iboxlo2, iboxhi2
         do j=iboxlo1, iboxhi1
            imin = iboxlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot + color, 2))
            imax = iboxhi0
            do i = imin, imax, 2
               newphi(i,j,k) = newphi(i,j,k)
     & + weight*resid(i,j,k)
            enddo
         enddo
      enddo
      return
      end
      subroutine REGGETFLUX(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
     & ,iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
     & ,beta
     & ,dx
     & ,idir
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer iopphiboxlo0,iopphiboxlo1,iopphiboxlo2
      integer iopphiboxhi0,iopphiboxhi1,iopphiboxhi2
      REAL*8 beta
      REAL*8 dx(0:2)
      integer idir
      integer ncomp
      integer i,j,k
      integer ioff,joff,koff
      integer ivar
      REAL*8 oneoverdx
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      oneoverdx = beta/dx(idir)
      do ivar = 0, ncomp-1
      do k = iopphiboxlo2,iopphiboxhi2
      do j = iopphiboxlo1,iopphiboxhi1
      do i = iopphiboxlo0,iopphiboxhi0
         flux(i,j,k, ivar) =
     $ oneoverdx*(
     $ phi(i ,j ,k , ivar) -
     $ phi(i-ioff,j-joff,k-koff, ivar) )
      enddo
      enddo
      enddo
      enddo
      return
      end
