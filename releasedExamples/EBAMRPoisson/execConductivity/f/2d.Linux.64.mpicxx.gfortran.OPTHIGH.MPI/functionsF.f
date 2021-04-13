        subroutine GETETARESIST(
     & eta
     & ,ietalo0,ietalo1
     & ,ietahi0,ietahi1
     & ,freq
     & ,dx
     & ,problo
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & ,etaval
     & ,eps
     & ,whicheta
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ietalo0,ietalo1
      integer ietahi0,ietahi1
      REAL*8 eta(
     & ietalo0:ietahi0,
     & ietalo1:ietahi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      REAL*8 etaval
      REAL*8 eps
      integer whicheta
        integer i,j,jdir
        REAL*8 x(0:2 -1)
        integer iv(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
           iv(0) = i
           iv(1) = j
           do jdir = 0, 2 -1
              if(idir .eq. jdir) then
                 x(jdir) = iv(jdir)*dx(jdir) + problo(jdir)
              else
                 x(jdir) = (iv(jdir)+(0.500d0))*dx(jdir) + problo(jdir)
              endif
           enddo
          call getetapointresist(eta(i,j),freq,x,
     $ etaval, eps, whicheta)
      enddo
      enddo
        return
        end
        subroutine GETETAPOINTRESIST(
     & eta
     & ,freq
     & ,xval
     & ,etaval
     & ,eps
     & ,whicheta
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 eta
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      REAL*8 etaval
      REAL*8 eps
      integer whicheta
        REAL*8 x, y
        if(whicheta.eq. 1) then
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           eta = etaval*((1.0d0) + eps*(sin(x) + sin(y)))
        elseif(whicheta.eq.0) then
           eta = (1.0d0)
        elseif(whicheta.eq.3) then
           eta = (0.500d0)
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBRESIST(
     & klb
     & ,iklblo0,iklblo1
     & ,iklbhi0,iklbhi1
     & ,freq
     & ,dx
     & ,problo
     & ,alpha
     & ,beta
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,icomp
     & ,eps
     & ,whichmag
     & ,whicheta
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL*8 klb(
     & iklblo0:iklbhi0,
     & iklblo1:iklbhi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL*8 eps
      integer whichmag
      integer whicheta
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getklbpointresist(klb(i,j),freq,x,
     $ alpha, beta, icomp, eps, whichmag, whicheta)
      enddo
      enddo
        return
        end
        subroutine GETMAGRESIST(
     & mag
     & ,imaglo0,imaglo1
     & ,imaghi0,imaghi1
     & ,freq
     & ,dx
     & ,problo
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,icomp
     & ,whichmag
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imaglo0,imaglo1
      integer imaghi0,imaghi1
      REAL*8 mag(
     & imaglo0:imaghi0,
     & imaglo1:imaghi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      integer whichmag
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getmagpointresist(mag(i,j),freq,x,
     $ icomp, whichmag)
      enddo
      enddo
        return
        end
        subroutine GETMAGPOINTRESIST(
     & mag
     & ,freq
     & ,xval
     & ,icomp
     & ,whichmag
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 mag
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      integer icomp
      integer whichmag
        REAL*8 x,y, time
        integer i,j
        i = icomp
        j = max(1-icomp, 0)
        if(whichmag.eq. 2) then
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
           mag = sin(y)
        else if (whichmag .eq.3) then
           time = (0.0d0)
           call getphipoint(mag, freq, xval, time)
        elseif(whichmag.eq. 1) then
           x = freq(icomp)*xval(icomp)
           mag = sin(x)
        elseif(whichmag.eq.0) then
           x = xval(icomp)
           mag = x*x
        elseif(whichmag.eq.4) then
           x = xval(icomp)
           mag = x
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETDVDXPOINTRESIST(
     & mag
     & ,freq
     & ,xval
     & ,icomp
     & ,ideriv
     & ,whichmag
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 mag
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      integer icomp
      integer ideriv
      integer whichmag
        REAL*8 x,y, time, gradphi(0:2 -1)
        integer i,j
        i = icomp
        j = max(1-icomp, 0)
        if(whichmag.eq. 2) then
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
           if(ideriv.eq.j) then
              mag = freq(j)*cos(y)
           else
              mag = (0.0d0)
           endif
        else if (whichmag .eq.3) then
           time = (0.0d0)
           call getgradphipoint(gradphi, freq, xval, time)
           mag = gradphi(ideriv)
        elseif(whichmag.eq. 1) then
           if(ideriv.eq.icomp) then
              x = freq(icomp)*xval(icomp)
              mag = freq(icomp)*cos(x)
           else
              mag = (0.0d0)
           endif
        elseif(whichmag.eq. 4) then
           if(ideriv.eq.icomp) then
              mag = (1.0d0)
           else
              mag = (0.0d0)
           endif
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBPOINTRESIST(
     & klb
     & ,freq
     & ,xvec
     & ,alpha
     & ,beta
     & ,icomp
     & ,eps
     & ,whichmag
     & ,whicheta
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 klb
      REAL*8 freq(0:1)
      REAL*8 xvec(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer icomp
      REAL*8 eps
      integer whichmag
      integer whicheta
        REAL*8 fx,fy, time
        REAL*8 x,y, termone, etaval
        REAL*8 freqx,freqy, mag, divf, eta
        integer i,j
        i = icomp
        j = max(1-icomp, 0)
        etaval = (1.0d0)
        call getetapointresist(eta,freq,xvec, etaval, eps, whicheta)
        call getmagpointresist(mag,freq,xvec, icomp, whichmag)
        freqx = freq(i)
        freqy = freq(j)
        x = xvec(i)
        y = xvec(j)
        if((whichmag.eq. 2).and.(whicheta.eq.0)) then
           divf = -(freqy*freqy*sin(freqy*y))
        elseif((whichmag.eq. 3).and.(whicheta.eq.0)) then
           time = 0.0
           call getlofphipoint(klb, freq, xvec, alpha, beta, time)
           goto 123
        elseif((whichmag.eq. 2).and.(whicheta.eq.1)) then
           fx = freqx*x
           fy = freqy*y
           divf = (freqy*cos(fy) - freqx*cos(fx))*(eps*freqy*cos(fy)) - 
     &freqy*freqy*eta*sin(fy)
        elseif((whichmag.eq. 1).and.(whicheta.eq.1)) then
           termone =
     $ freqx*cos(freqx*x) +
     $ freqy*cos(freqy*y)
          divf = eps*freqx*cos(freqx*x)*termone
     $ -freqx*freqx*sin(freqx*x)*eta
        elseif((whichmag.eq.0).and.(whicheta.eq.0)) then
           divf = (2.0d0)
        elseif((whichmag.eq.4).and.(whicheta.eq.0)) then
           divf = (0.0d0)
        elseif((whichmag.eq.1).and.(whicheta.eq.0)) then
           divf = -freqx*freqx*sin(freqx*x)
        else
           call MayDay_Error()
        endif
        klb = alpha*mag + beta*divf
  123 continue
        return
        end
        subroutine GETKLVVISCOUS(
     & klb
     & ,iklblo0,iklblo1
     & ,iklbhi0,iklbhi1
     & ,freq
     & ,dx
     & ,problo
     & ,alpha
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,icomp
     & ,eps
     & ,whichvel
     & ,whicheta
     & ,whichlambda
     & ,beta
     & ,lambdafactor
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL*8 klb(
     & iklblo0:iklbhi0,
     & iklblo1:iklbhi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 alpha
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL*8 eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL*8 beta
      REAL*8 lambdafactor
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getklvpointviscous(klb(i,j),
     $ freq,x, alpha, icomp, eps,
     $ whichvel, whicheta, whichlambda, beta, lambdafactor)
      enddo
      enddo
        return
        end
        subroutine GETKLVPOINTVISCOUS(
     & klv
     & ,freq
     & ,xvec
     & ,alpha
     & ,icomp
     & ,eps
     & ,whichvel
     & ,whicheta
     & ,whichlambda
     & ,beta
     & ,lambdafactor
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 klv
      REAL*8 freq(0:1)
      REAL*8 xvec(0:1)
      REAL*8 alpha
      integer icomp
      REAL*8 eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL*8 beta
      REAL*8 lambdafactor
        REAL*8 x,y, lambda
        REAL*8 freqx,freqy, vel, divf, eta
        REAL*8 fx,fy, etaval
        integer i,j
        i = icomp
        j = max(1-icomp, 0)
        etaval = (1.0d0)
        call getetapointresist( eta,freq,xvec, etaval, eps, whicheta)
        if(whichlambda .eq. 2) then
           lambda = -lambdafactor*eta
        else
           call getetapointresist(lambda,freq, xvec, etaval, eps, whichl
     &ambda)
        endif
        call getmagpointresist( vel,freq,xvec, icomp, whichvel)
        freqx = freq(i)
        freqy = freq(j)
        x = xvec(i)
        y = xvec(j)
        fx = freqx*x
        fy = freqy*y
        if((whichvel.eq.2).and.(whicheta.eq.0)) then
           divf = -freqy*freqy*sin(fy)
        else if((whichvel.eq.1).and.(whicheta.eq.3).and.(whichlambda.eq.
     &3)) then
           divf = -(3.0d0)*(0.500d0)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.
     &0)) then
           divf = -(3.0d0)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.
     &2)) then
           divf = -((2.0d0) - lambdafactor)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.
     &1)) then
           divf = -(3.0d0)*eta*freqx*freqx*sin(fx)
           divf = divf + (3.0d0)*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf + eps*freqx*freqy*cos(fx)*cos(fy)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.
     &2)) then
           divf = -((2.0d0) - lambdafactor)*eta*freqx*freqx*sin(fx)
           divf = divf + ((2.0d0) - lambdafactor)*eps*freqx*freqx*cos(fx
     &)*cos(fx)
           divf = divf - (lambdafactor)*eps*freqx*freqy*cos(fx)*cos(fy)
        else if((whichvel.eq.2).and.(whicheta.eq.1)) then
           divf = - eta*freqy*freqy*sin(fy)
           divf = divf + eps*freqy*cos(fy)*(freqy*cos(fy) + freqx*cos(fx
     &))
        else if(whichvel.eq.3) then
           divf = 0
           call getlofphipoint(klv,freq,xvec,alpha,beta,divf)
        else if(whichvel.eq.4) then
           divf = (0.0d0)
        else
           call MayDay_Error()
        endif
        klv = alpha*vel + beta*divf
        return
        end
      REAL*8 function getphirzfunc(radius)
      implicit none
      REAL*8 radius
      getphirzfunc = radius*radius
      return
      end
      REAL*8 function getgradphirzfunc(radius)
      implicit none
      REAL*8 radius
      getgradphirzfunc = (2.0d0)*radius
      return
      end
      REAL*8 function getlaplphirzfunc(radius)
      implicit none
      REAL*8 radius
      getlaplphirzfunc = (4.0d0)
      return
      end
        subroutine GETPHI(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,freq
     & ,dx
     & ,time
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 time
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getphipoint(phi(i,j),freq,x,time)
      enddo
      enddo
        return
        end
        subroutine GETPHIPOINT(
     & phi
     & ,freq
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 phi
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        phi = sin(freq(0)*x(0))
     & * sin(freq(1)*x(1))
        phi = phi*cos(time)
        return
        end
        subroutine GETSHPHI(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,lmp
     & ,dx
     & ,time
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      REAL*8 lmp(0:1)
      REAL*8 dx(0:1)
      REAL*8 time
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getshphipoint(phi(i,j),lmp,x,time)
      enddo
      enddo
        return
        end
        subroutine GETSHPHIPOINT(
     & phi
     & ,lmp
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 phi
      REAL*8 lmp(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        external normalization
        REAL*8 normalization
        external algndr
        REAL*8 algndr
        external bessel
        REAL*8 bessel
        external trigphi
        REAL*8 trigphi
        REAL*8 r(0:2 -1)
        integer l,m
        REAL*8 z,phase
        call convert_spherical(x,r)
        l = int(lmp(0))
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     & (r(0)**z)*
     & trigphi(l,r(1),phase)
        return
        end
        subroutine GETLOFPHIRZPOLY(
     & lofphi
     & ,x
     & ,alpha
     & ,beta
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 lofphi
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 beta
        REAL*8 phi, laplphi
        REAL*8 dist
        external getlaplphirzfunc
        REAL*8 getlaplphirzfunc
        external getphirzfunc
        REAL*8 getphirzfunc
        dist = abs(x(0))
        phi = getphirzfunc(dist)
        laplphi = getlaplphirzfunc(dist)
        lofphi = alpha*phi + beta*laplphi
        return
        end
        subroutine GETPHIRZPOLY(
     & phi
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 phi
      REAL*8 x(0:1)
        REAL*8 dist
        external getphirzfunc
        REAL*8 getphirzfunc
        dist =abs(x(0))
        phi = getphirzfunc(dist)
        return
        end
      subroutine GETGRADPHIRZPOLY(
     & gradphi
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 x(0:1)
        REAL*8 dist
        external getgradphirzfunc
        REAL*8 getgradphirzfunc
        dist = abs(x(0))
        gradphi(0) = getgradphirzfunc(dist)
        gradphi(1) = (0.0d0)
        return
        end
        subroutine GETGRADPHIPOINT(
     & gradphi
     & ,freq
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        integer idir
        gradphi(0) = freq(0) * cos(freq(0)*x(0)) * sin(freq(1)*x(1))
        gradphi(1) = freq(1) * sin(freq(0)*x(0)) * cos(freq(1)*x(1))
        do idir = 0, 2 -1
           gradphi(idir) = gradphi(idir)*cos(time)
        enddo
        return
        end
        subroutine GETGRADSHPHIPOINT(
     & gradphi
     & ,lmp
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 lmp(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        external normalization
        REAL*8 normalization
        external algndr
        REAL*8 algndr
        external bessel
        REAL*8 bessel
        external d_algndr
        REAL*8 d_algndr
        external d_bessel
        REAL*8 d_bessel
        external trigphi
        REAL*8 trigphi
        external d_trigphi
        REAL*8 d_trigphi
        REAL*8 gr,gtheta
        REAL*8 rad,st
        REAL*8 ir,ct
        REAL*8 r(0:2 -1)
        integer l,m,idir
        REAL*8 z,phase
        call convert_spherical(x,r)
                  rad = r(0)
                  st = sin(r(1))
        if (rad.eq.0.) call MayDay_Error()
                  ir = 1./rad
                  ct = cos(r(1))
        l = int(lmp(0))
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        gr = normalization(l,m)*
     & z*(rad**(z-1.))*
     & trigphi(l,r(1),phase)
        gtheta = normalization(l,m)*
     & (rad**z)*
     & ir*d_trigphi(l,r(1),phase)
        gradphi(0) = ct*gr - st*gtheta
        gradphi(1) = st*gr + ct*gtheta
        do idir = 0, 2 -1
           gradphi(idir) = gradphi(idir)*cos(time)
        enddo
        return
        end
        subroutine GETMARSHAGRADPHIPOINT(
     & gradphi
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 x(0:1)
        gradphi(0) = cos(x(0))*exp(x(1))
        gradphi(1) = sin(x(0))*exp(x(1))
        return
        end
        subroutine GETLOFPHI(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,freq
     & ,dx
     & ,problo
     & ,probhi
     & ,alpha
     & ,beta
     & ,time
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getlofphipoint(lofphi(i,j),freq,x,alpha,beta,time)
      enddo
      enddo
        return
        end
        subroutine GETLOFSHPHI(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,lmp
     & ,dx
     & ,problo
     & ,probhi
     & ,alpha
     & ,beta
     & ,time
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1)
      REAL*8 lmp(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getlofshphipoint(lofphi(i,j),lmp,x,alpha,beta,time)
      enddo
      enddo
        return
        end
        subroutine GETMARSHALOFPHI(
     & lofphi
     & ,ilofphilo0,ilofphilo1
     & ,ilofphihi0,ilofphihi1
     & ,dx
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getmarshalofphipoint(lofphi(i,j),x)
      enddo
      enddo
        return
        end
        subroutine GETMARSHALOFPHIPOINT(
     & lofphi
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 lofphi
      REAL*8 x(0:1)
        lofphi = (0.0d0)
        return
        end
        subroutine GETMARSHAPHI(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,dx
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getmarshaphipoint(phi(i,j),x)
      enddo
      enddo
        return
        end
        subroutine GETMARSHAPHIPOINT(
     & phi
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 phi
      REAL*8 x(0:1)
        phi = sin(x(0))
     $ *exp(x(1))
        return
        end
        subroutine GETLOFPHIPOINT(
     & lofphi
     & ,freq
     & ,x
     & ,alpha
     & ,beta
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 lofphi
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 time
        REAL*8 fac,phi
        fac = -(freq(0)**2
     & + freq(1)**2)
        call getphipoint(phi, freq, x, time)
        lofphi = fac*phi
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETLOFSHPHIPOINT(
     & lofphi
     & ,lmp
     & ,x
     & ,alpha
     & ,beta
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 lofphi
      REAL*8 lmp(0:1)
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 beta
      REAL*8 time
        external normalization
        REAL*8 normalization
        external algndr
        REAL*8 algndr
        external bessel
        REAL*8 bessel
        external trigphi
        REAL*8 trigphi
        REAL*8 r(0:2 -1)
        integer l,m
        REAL*8 z,phase,phi
        call convert_spherical(x,r)
        l = int(lmp(0))
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     & (r(0)**z)*
     & trigphi(l,r(1),phase)
        phi = phi*cos(time)
        lofphi = -phi
        lofphi = lofphi/(r(0)**2)
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETDBGPHI(
     & dbgphi
     & ,idbgphilo0,idbgphilo1
     & ,idbgphihi0,idbgphihi1
     & ,beta
     & ,ibetalo0,ibetalo1
     & ,ibetahi0,ibetahi1
     & ,freq
     & ,dx
     & ,problo
     & ,probhi
     & ,alpha
     & ,time
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idbgphilo0,idbgphilo1
      integer idbgphihi0,idbgphihi1
      REAL*8 dbgphi(
     & idbgphilo0:idbgphihi0,
     & idbgphilo1:idbgphihi1)
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL*8 beta(
     & ibetalo0:ibetahi0,
     & ibetalo1:ibetahi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      REAL*8 alpha
      REAL*8 time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getdbgphipoint(dbgphi(i,j),
     & beta(i,j),freq,x,alpha,time)
      enddo
      enddo
        return
        end
        subroutine GETDBGPHIPOINT(
     & dbgphi
     & ,beta
     & ,freq
     & ,x
     & ,alpha
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dbgphi
      REAL*8 beta
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 time
        REAL*8 gradphi(0:2 -1),gradbeta(0:2 -1)
        REAL*8 alphaphiplusbetalapphi,gradbetadotgradphi
        call getbetapoint(beta,freq,x,time)
        call getlofphipoint(alphaphiplusbetalapphi,freq,x,alpha,beta,tim
     &e)
        call getgradbetapoint(gradbeta,freq,x,time)
        call getgradphipoint(gradphi,freq,x,time)
        gradbetadotgradphi = gradbeta(0)*gradphi(0)
     & + gradbeta(1)*gradphi(1)
        dbgphi = alphaphiplusbetalapphi
        dbgphi = dbgphi + gradbetadotgradphi
        return
        end
        subroutine GETBETAPOINT(
     & beta
     & ,freq
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 beta
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        beta = 2 +2
     & +sin((2.0d0)*(3.14159265358979323846264338327950288D0)*x(0))
     & + sin((2.0d0)*(3.14159265358979323846264338327950288D0)*x(1))
        return
        end
        subroutine GETGRADBETAPOINT(
     & gradbeta
     & ,freq
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradbeta(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        integer idir
        do idir = 0, 2 -1
            gradbeta(idir) = (2.0d0)*(3.14159265358979323846264338327950
     &288D0)*cos((2.0d0)*(3.14159265358979323846264338327950288D0)*x(idi
     &r))
        enddo
        return
        end
        subroutine GETBETAGRADPHIPOINT(
     & gradphi
     & ,freq
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        integer idir
        REAL*8 beta
        call getbetapoint(beta,freq,x,time)
        call getgradphipoint(gradphi,freq,x,time)
        do idir = 0, 2 -1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETBETAGRADSHPHIPOINT(
     & gradphi
     & ,lmp
     & ,x
     & ,time
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 gradphi(0:1)
      REAL*8 lmp(0:1)
      REAL*8 x(0:1)
      REAL*8 time
        integer idir
        REAL*8 beta
        call getbetapoint(beta,lmp,x,time)
        call getgradshphipoint(gradphi,lmp,x,time)
        do idir = 0, 2 -1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETSRC(
     & src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,freq
     & ,dx
     & ,time
     & ,diffconst
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 time
      REAL*8 diffconst
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getsrcpoint(src(i,j),freq,x,time,diffconst)
      enddo
      enddo
        return
        end
        subroutine GETSRCPOINT(
     & src
     & ,freq
     & ,x
     & ,time
     & ,diffconst
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 src
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 time
      REAL*8 diffconst
        REAL*8 fac,phi
        fac = -(freq(0)**2
     & + freq(1)**2)
        phi = (sin(freq(0)*x(0))
     & * sin(freq(1)*x(1)))
        src = (-fac*diffconst*cos(time) - sin(time))*phi
        return
        end
        subroutine GETSHSRC(
     & src
     & ,isrclo0,isrclo1
     & ,isrchi0,isrchi1
     & ,lmp
     & ,dx
     & ,time
     & ,diffconst
     & ,problo
     & ,probhi
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     & isrclo0:isrchi0,
     & isrclo1:isrchi1)
      REAL*8 lmp(0:1)
      REAL*8 dx(0:1)
      REAL*8 time
      REAL*8 diffconst
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getshsrcpoint(src(i,j),lmp,x,time,diffconst)
      enddo
      enddo
        return
        end
        subroutine GETSHSRCPOINT(
     & src
     & ,lmp
     & ,x
     & ,time
     & ,diffconst
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 src
      REAL*8 lmp(0:1)
      REAL*8 x(0:1)
      REAL*8 time
      REAL*8 diffconst
        external normalization
        REAL*8 normalization
        external algndr
        REAL*8 algndr
        external bessel
        REAL*8 bessel
        external trigphi
        REAL*8 trigphi
        REAL*8 rad,z
        REAL*8 r(0:2 -1)
        integer l,m
        call convert_spherical(x,r)
        rad = r(0)
        l = int(lmp(0))
        m = 0
        z = cos(r(1))
        src = -diffconst*normalization(l,m)*
     & bessel(l,rad)*
     & algndr(l,m,z)
        return
        end
        REAL*8 function normalization(l,m)
        implicit none
        integer l,m
        external factorial
        integer factorial
        normalization = (2*l+1.)*factorial(l)*0.5/factorial(l)/
     & ((2.0d0)*(3.14159265358979323846264338327950288D0))
        normalization = sqrt(normalization)
        return
        end
        REAL*8 function algndr(l,m,z)
        implicit none
        integer l,m
        REAL*8 z
        integer i,ll
        REAL*8 fac,pll,pmm,pmmp1,somz2
        if (l.lt.0 .or. m.lt.0 .or. m.gt.l) call MayDay_Error()
        pmm=1.
        if (m.gt.0) then
           somz2 = sqrt((1.-z)*(1.+z))
           fac = 1.
           do i=1,m
              pmm = -pmm*fac*somz2
              fac = fac+2.
           enddo
        endif
        if (l.eq.m) then
           algndr = pmm
        else
           pmmp1 = z*(2*m+1)*pmm
           if (l.eq.m+1) then
              algndr = pmmp1
           else
              do ll=m+2,l
                 pll = (z*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                 pmm = pmmp1
                 pmmp1 = pll
              enddo
              algndr = pll
           endif
        endif
        return
        end
        REAL*8 function d_algndr(l,m,z)
        implicit none
        external algndr
        REAL*8 algndr
        integer l,m
        REAL*8 z
        if (z.eq.1) then
           if (l.eq.0) then
              d_algndr = 0.
           elseif (l.eq.1) then
              if (m.eq.0) then
                 d_algndr = 0.
              elseif (m.eq.1) then
                 d_algndr = -1.
              else
                 call MayDay_Error()
              endif
           elseif (l.eq.2) then
              if (m.eq.0) then
                 d_algndr = 0.
              elseif (m.eq.1) then
                 d_algndr = -3.
              elseif (m.eq.2) then
                 d_algndr = 0.
              else
                 call MayDay_Error()
              endif
           elseif (l.eq.3) then
              if (m.eq.0) then
                 d_algndr = 0.
              elseif (m.eq.1) then
                 d_algndr = -6.
              elseif (m.eq.2) then
                 d_algndr = 0.
              elseif (m.eq.3) then
                 d_algndr = 0.
              else
                 call MayDay_Error()
              endif
           else
              call MayDay_Error()
           endif
        else
           d_algndr = l*z*algndr(l,m,z)
           if (l.gt.0 .and. l.gt.m) then
              d_algndr = d_algndr - (l+m)*algndr(l-1,m,z)
           endif
           d_algndr = d_algndr/sqrt(1.-z**2)
        endif
        return
        end
        subroutine CONVERT_SPHERICAL(
     & x
     & ,r
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 x(0:1)
      REAL*8 r(0:1)
        r(0) = sqrt(x(0)**2 + x(1)**2)
        r(1) = 0.
        if (x(0).ne.0. .or. x(1).ne.0.) then
           r(1) = atan2(x(1),x(0))
           if (r(1).lt.0.) r(1) = r(1) + (2.0d0)*(3.14159265358979323846
     &264338327950288D0)
        endif
        return
        end
        integer function factorial(n)
        implicit none
        integer n
        integer i
        factorial=1
        do i=n,2,-1
           factorial = factorial*i
        enddo
        return
        end
        REAL*8 function trigphi(m, azphi, phase)
        implicit none
        integer m
        REAL*8 azphi, phase
        trigphi = cos(m*azphi - phase)/sqrt((3.1415926535897932384626433
     &8327950288D0))
        return
        end
        REAL*8 function d_trigphi(m, azphi, phase)
        implicit none
        integer m
        REAL*8 azphi, phase
        d_trigphi = -m*sin(m*azphi - phase)/sqrt((3.14159265358979323846
     &264338327950288D0))
        return
        end
        REAL*8 function bessel(l, r)
        implicit none
        integer l
        REAL*8 r
        REAL*8 ir,s,c
        if (r.eq.0.) then
           if (l.eq.0) then
              bessel = 1.
           elseif (l.ge.1 .and. l.le.3) then
              bessel = 0.
           else
              call MayDay_Error()
           endif
        else
           ir = 1./r
           s = sin(r)
           if (l.eq.0) then
              bessel = s*ir
           else
              c = cos(r)
              if (l.eq.1) then
                 bessel = s*(ir**2) - c*ir
              elseif (l.eq.2) then
                 bessel = (3.*ir**2-1.)*s*ir - 3.*c*(ir**2)
              elseif (l.eq.3) then
                 bessel = (15.*ir**2-6.)*s*(ir**2) -
     & (15.*ir**2-1.)*c*ir
              else
                 call MayDay_Error()
              endif
           endif
        endif
        return
        end
        REAL*8 function d_bessel(l, r)
        implicit none
        integer l
        REAL*8 r
        REAL*8 ir,irsq,s,c
        if (r.eq.0.) then
           if (l.eq.0) then
              d_bessel = 0.
           elseif (l.eq.1) then
              d_bessel = 1./3.
           elseif (l.eq.2) then
              d_bessel = 0.
           elseif (l.eq.3) then
              d_bessel = 2.
           else
              call MayDay_Error()
           endif
        else
           ir = 1./r
           irsq = ir**2
           s = sin(r)
           c = cos(r)
           if (l.eq.0) then
              d_bessel = -s*irsq + c*ir
           elseif (l.eq.1) then
              d_bessel = (-2.*irsq+1.)*s*ir + 2.*c*irsq
           elseif (l.eq.2) then
              d_bessel = (-9.*irsq+4.)*s*irsq +
     & ( 9.*irsq-1.)*c*ir
           elseif (l.eq.3) then
              d_bessel = (-60.*irsq**2+15.*irsq-1.)*s*ir +
     & (60.*irsq-7.)*c*irsq
           else
              call MayDay_Error()
           endif
        endif
        return
        end
