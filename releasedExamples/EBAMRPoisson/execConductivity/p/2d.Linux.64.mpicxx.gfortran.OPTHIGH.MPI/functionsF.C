#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
        subroutine GETETARESIST(
     &           eta
     &           ,ietalo0,ietalo1
     &           ,ietahi0,ietahi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,idir
     &           ,etaval
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ietalo0,ietalo1
      integer ietahi0,ietahi1
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      REAL_T etaval
      REAL_T eps
      integer whicheta
        integer i,j,jdir
        real_t x(0:CH_SPACEDIM-1)
        integer iv(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

           
           iv(0) = i
           iv(1) = j
           do jdir = 0, CH_SPACEDIM-1
              if(idir .eq. jdir) then
                 x(jdir) = iv(jdir)*dx(jdir) + problo(jdir)
              else
                 x(jdir) = (iv(jdir)+half)*dx(jdir) + problo(jdir)
              endif
           enddo
          call getetapointresist(eta(i,j),freq,x,
     $          etaval, eps, whicheta)
        
      enddo
      enddo
        return
        end
        subroutine GETETAPOINTRESIST(
     &           eta
     &           ,freq
     &           ,xval
     &           ,etaval
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T eta
      REAL_T freq(0:1)
      REAL_T xval(0:1)
      REAL_T etaval
      REAL_T eps
      integer whicheta
        REAL_T x, y
        if(whicheta.eq. 1) then
           
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           eta  = etaval*(one + eps*(sin(x) + sin(y)))
        elseif(whicheta.eq.0) then
           eta = one
        elseif(whicheta.eq.3) then
           eta = half
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBRESIST(
     &           klb
     &           ,iklblo0,iklblo1
     &           ,iklbhi0,iklbhi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL_T klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getklbpointresist(klb(i,j),freq,x,
     $         alpha, beta, icomp, eps, whichmag, whicheta)
        
      enddo
      enddo
        return
        end
        subroutine GETMAGRESIST(
     &           mag
     &           ,imaglo0,imaglo1
     &           ,imaghi0,imaghi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imaglo0,imaglo1
      integer imaghi0,imaghi1
      REAL_T mag(
     &           imaglo0:imaghi0,
     &           imaglo1:imaghi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      integer whichmag
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getmagpointresist(mag(i,j),freq,x,
     $         icomp, whichmag)
        
      enddo
      enddo
        return
        end
        subroutine GETMAGPOINTRESIST(
     &           mag
     &           ,freq
     &           ,xval
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T mag
      REAL_T freq(0:1)
      REAL_T xval(0:1)
      integer icomp
      integer whichmag
        REAL_T  x,y, time
        integer i,j
        
        i = icomp
        j = max(1-icomp, 0)
        if(whichmag.eq. 2) then
           
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
#if CH_SPACEDIM==2
           mag = sin(y)
#elif CH_SPACEDIM==3
           mag = sin(y) + sin(z)
#else
           bogus_spacedim()
#endif
        else if (whichmag .eq.3) then
           time = zero
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
     &           mag
     &           ,freq
     &           ,xval
     &           ,icomp
     &           ,ideriv
     &           ,whichmag
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T mag
      REAL_T freq(0:1)
      REAL_T xval(0:1)
      integer icomp
      integer ideriv
      integer whichmag
        REAL_T  x,y, time, gradphi(0:CH_SPACEDIM-1)
        integer i,j
        
        i = icomp
        j = max(1-icomp, 0)
        if(whichmag.eq. 2) then
           
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
#if CH_SPACEDIM==2
           if(ideriv.eq.j) then
              mag = freq(j)*cos(y)
           else
              mag = zero
           endif
#elif CH_SPACEDIM==3
           if(ideriv.eq.j) then
              mag = freq(j)*cos(y)
           else if (ideriv.eq.k) then
              mag = freq(k)*cos(z)
           else
              mag = zero
           endif
#else
           bogus_spacedim()
#endif
        else if (whichmag .eq.3) then
           time = zero
           call getgradphipoint(gradphi, freq, xval, time)
           mag = gradphi(ideriv)
        elseif(whichmag.eq. 1) then
           if(ideriv.eq.icomp) then
              x = freq(icomp)*xval(icomp)
              mag = freq(icomp)*cos(x)
           else
              mag = zero
           endif
        elseif(whichmag.eq. 4) then
           if(ideriv.eq.icomp) then
              mag = one
           else
              mag = zero
           endif
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBPOINTRESIST(
     &           klb
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,beta
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T klb
      REAL_T freq(0:1)
      REAL_T xvec(0:1)
      REAL_T alpha
      REAL_T beta
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        REAL_T  fx,fy, time
        REAL_T  x,y,  termone, etaval
        REAL_T  freqx,freqy, mag, divf, eta
        integer i,j
        
        i = icomp
        j = max(1-icomp, 0)
        etaval = one
        call getetapointresist(eta,freq,xvec,  etaval, eps, whicheta)
        call getmagpointresist(mag,freq,xvec,  icomp, whichmag)
        
        freqx = freq(i)
        freqy = freq(j)
        
        x = xvec(i)
        y = xvec(j)
        if((whichmag.eq. 2).and.(whicheta.eq.0)) then
           divf = -(freqy*freqy*sin(freqy*y))
#if CH_SPACEDIM==3
           divf = divf - (freqz*freqz*sin(freqz*z))
#endif
        elseif((whichmag.eq. 3).and.(whicheta.eq.0)) then
           time = 0.0
           call getlofphipoint(klb, freq, xvec, alpha, beta, time)
           goto  123
        elseif((whichmag.eq. 2).and.(whicheta.eq.1)) then
           
           fx = freqx*x
           fy = freqy*y
           divf = (freqy*cos(fy) - freqx*cos(fx))*(eps*freqy*cos(fy)) - freqy*freqy*eta*sin(fy)
#if CH_SPACEDIM==3
           divf = eps*freqy*cos(fy)*(freqy*cos(fy) - freqx*cos(fx))  - eta*freqy*freqy*sin(fy)
     $          + eps*freqz*cos(fz)*(freqz*cos(fz) - freqx*cos(fx))  - eta*freqz*freqz*sin(fz)
#endif
        elseif((whichmag.eq. 1).and.(whicheta.eq.1)) then
           termone =  
     $          freqx*cos(freqx*x) +
     $          freqy*cos(freqy*y)
          divf = eps*freqx*cos(freqx*x)*termone
     $         -freqx*freqx*sin(freqx*x)*eta
        elseif((whichmag.eq.0).and.(whicheta.eq.0)) then
           divf = two
        elseif((whichmag.eq.4).and.(whicheta.eq.0)) then
           divf = zero
        elseif((whichmag.eq.1).and.(whicheta.eq.0)) then
           divf = -freqx*freqx*sin(freqx*x)
        else
           call MayDay_Error()
        endif
        klb = alpha*mag  + beta*divf
  123   continue
        return
        end
        subroutine GETKLVVISCOUS(
     &           klb
     &           ,iklblo0,iklblo1
     &           ,iklbhi0,iklbhi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,beta
     &           ,lambdafactor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL_T klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T alpha
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T beta
      REAL_T lambdafactor
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getklvpointviscous(klb(i,j),
     $         freq,x, alpha, icomp, eps,
     $         whichvel, whicheta, whichlambda, beta, lambdafactor)
        
      enddo
      enddo
        return
        end
        subroutine GETKLVPOINTVISCOUS(
     &           klv
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,beta
     &           ,lambdafactor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T klv
      REAL_T freq(0:1)
      REAL_T xvec(0:1)
      REAL_T alpha
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T beta
      REAL_T lambdafactor
        REAL_T  x,y,   lambda
        REAL_T  freqx,freqy, vel, divf, eta
        real_t fx,fy, etaval
        integer i,j
        
        i = icomp
        j = max(1-icomp, 0)
        etaval = one
        call getetapointresist(   eta,freq,xvec, etaval, eps, whicheta)
        if(whichlambda .eq. 2) then
           lambda = -lambdafactor*eta
        else
           call getetapointresist(lambda,freq, xvec, etaval, eps, whichlambda)
        endif
        call getmagpointresist(   vel,freq,xvec, icomp, whichvel)
        
        freqx = freq(i)
        freqy = freq(j)
        
        x = xvec(i)
        y = xvec(j)
        
        fx = freqx*x
        fy = freqy*y
        if((whichvel.eq.2).and.(whicheta.eq.0)) then
           divf = -freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf -freqz*freqz*sin(fz)
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.3).and.(whichlambda.eq.3)) then
           divf = -three*half*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.0)) then
           divf = -three*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.2)) then
           divf = -(two - lambdafactor)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.1)) then
           divf =       -three*eta*freqx*freqx*sin(fx)
           divf = divf + three*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf +       eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf +       eps*freqx*freqz*cos(fx)*cos(fz)
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.2)) then
           divf =       -(two - lambdafactor)*eta*freqx*freqx*sin(fx)
           divf = divf + (two - lambdafactor)*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf -       (lambdafactor)*eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf -       (lambdafactor)*eps*freqx*freqz*cos(fx)*cos(fz)
#endif
        else if((whichvel.eq.2).and.(whicheta.eq.1)) then
           divf =       - eta*freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf  - eta*freqz*freqz*sin(fz)
#endif
           divf = divf + eps*freqy*cos(fy)*(freqy*cos(fy) + freqx*cos(fx))
#if CH_SPACEDIM==3
           divf = divf + eps*freqz*cos(fz)*(freqz*cos(fz) + freqx*cos(fx))
#endif
        else if(whichvel.eq.3) then
           divf = 0
           call getlofphipoint(klv,freq,xvec,alpha,beta,divf)
        else if(whichvel.eq.4) then
           divf = zero
        else
           call MayDay_Error()
        endif
        klv = alpha*vel  + beta*divf
        return
        end
      real_t  function getphirzfunc(radius)
      implicit none
      real_t radius
      getphirzfunc = radius*radius
      return
      end
      real_t  function getgradphirzfunc(radius)
      implicit none
      real_t radius
      getgradphirzfunc = two*radius
      return
      end
      real_t  function getlaplphirzfunc(radius)
      implicit none
      real_t radius
      getlaplphirzfunc = four
      return
      end
        subroutine GETPHI(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,freq
     &           ,dx
     &           ,time
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T time
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getphipoint(phi(i,j),freq,x,time)
        
      enddo
      enddo
        return
        end
        subroutine GETPHIPOINT(
     &           phi
     &           ,freq
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
        phi = sin(freq(0)*x(0))
     &                * sin(freq(1)*x(1))
        phi = phi*cos(time)
        return
        end
        subroutine GETSHPHI(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,lmp
     &           ,dx
     &           ,time
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T lmp(0:1)
      REAL_T dx(0:1)
      REAL_T time
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getshphipoint(phi(i,j),lmp,x,time)
        
      enddo
      enddo
        return
        end
        subroutine GETSHPHIPOINT(
     &           phi
     &           ,lmp
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T lmp(0:1)
      REAL_T x(0:1)
      REAL_T time
        external normalization
        real_t normalization
        external algndr
        real_t algndr
        external bessel
        real_t bessel
        external trigphi
        real_t trigphi
        real_t r(0:CH_SPACEDIM-1)
        integer l,m
        real_t z,phase
        call convert_spherical(x,r)
        l = int(lmp(0))
#if CH_SPACEDIM == 1
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     &        trigphi(l,r(0),phase)
#elif CH_SPACEDIM == 2
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     &        (r(0)**z)*
     &        trigphi(l,r(1),phase)
#elif CH_SPACEDIM == 3
        m = int(lmp(1))
        phase = PI*lmp(2)
        z = cos(r(1))
        phi = normalization(l,m)*
     &        bessel(l,r(0))*
     &        algndr(l,m,z)*
     &        trigphi(m,r(2),phase)
#else
        bogus_spacedim()
#endif
        return
        end
        subroutine GETLOFPHIRZPOLY(
     &           lofphi
     &           ,x
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T x(0:1)
      REAL_T alpha
      REAL_T beta
        real_t phi, laplphi
        real_t dist
        external getlaplphirzfunc
        real_t getlaplphirzfunc
        external getphirzfunc
        real_t getphirzfunc
        dist = abs(x(0))
        phi = getphirzfunc(dist)
        laplphi = getlaplphirzfunc(dist)
        lofphi = alpha*phi + beta*laplphi
        return
        end
        subroutine GETPHIRZPOLY(
     &           phi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T x(0:1)
        real_t dist
        external getphirzfunc
        real_t getphirzfunc
        dist =abs(x(0))
        phi = getphirzfunc(dist)
        return
        end
      subroutine GETGRADPHIRZPOLY(
     &           gradphi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T x(0:1)
        real_t dist
        external getgradphirzfunc
        real_t getgradphirzfunc
        dist = abs(x(0))
        
        gradphi(0) = getgradphirzfunc(dist)
        gradphi(1) = zero
        return
        end
        subroutine GETGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
        integer idir
        
        gradphi(0) = freq(0) * cos(freq(0)*x(0)) * sin(freq(1)*x(1))
        gradphi(1) = freq(1) * sin(freq(0)*x(0)) * cos(freq(1)*x(1))                   
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*cos(time)
        enddo
        return
        end
        subroutine GETGRADSHPHIPOINT(
     &           gradphi
     &           ,lmp
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T lmp(0:1)
      REAL_T x(0:1)
      REAL_T time
        external normalization
        real_t normalization
        external algndr
        real_t algndr
        external bessel
        real_t bessel
        external d_algndr
        real_t d_algndr
        external d_bessel
        real_t d_bessel
        external trigphi
        real_t trigphi
        external d_trigphi
        real_t d_trigphi
        real_t gr,gtheta
        real_t rad,st
        real_t ir,ct
        real_t r(0:CH_SPACEDIM-1)
        integer l,m,idir
        real_t z,phase
        call convert_spherical(x,r)
        
                  rad = r(0)
                  st = sin(r(1))
        if (rad.eq.0.) call MayDay_Error()
        
                  ir = 1./rad
                  ct = cos(r(1))
        l = int(lmp(0))
#if CH_SPACEDIM == 1
        m = 0
        phase = 0.
        gr = normalization(l,m)*
     &       d_trigphi(l,rad,phase)
#elif CH_SPACEDIM == 2
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        gr = normalization(l,m)*
     &       z*(rad**(z-1.))*
     &       trigphi(l,r(1),phase)
        gtheta = normalization(l,m)*
     &           (rad**z)*
     &           ir*d_trigphi(l,r(1),phase)
#elif CH_SPACEDIM == 3
        m = int(lmp(1))
        phase = PI*lmp(2)
        gr = normalization(l,m)*
     &       d_bessel(l,rad)*
     &       algndr(l,m,ct)*
     &       trigphi(m,r(2),phase)
        gtheta = normalization(l,m)*
     &           bessel(l,rad)*
     &           ir*d_algndr(l,m,ct)*
     &           trigphi(m,r(2),phase)
        if (st.ne.0.) then
           gphi = normalization(l,m)*
     &            bessel(l,rad)*
     &            algndr(l,m,ct)*
     &            (ir/st)*d_trigphi(m,r(2),phase)
        else
           gphi = 0.
        endif
#else
        bogus_spacedim();
#endif
        
        gradphi(0) = ct*gr - st*gtheta
        gradphi(1) = st*gr + ct*gtheta                
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*cos(time)
        enddo
        return
        end
        subroutine GETMARSHAGRADPHIPOINT(
     &           gradphi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T x(0:1)
#if CH_SPACEDIM == 1
        gradphi(0) = zero
#elif CH_SPACEDIM == 2
        gradphi(0) = cos(x(0))*exp(x(1))
        gradphi(1) = sin(x(0))*exp(x(1))
#elif CH_SPACEDIM == 3
        gradphi(0) = cos(x(0))*exp(x(1))
        gradphi(1) = sin(x(0))*exp(x(1))
        gradphi(2) = zero
#else
        bogus_spacedim();
#endif
        return
        end
        subroutine GETLOFPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,beta
     &           ,time
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      REAL_T alpha
      REAL_T beta
      REAL_T time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getlofphipoint(lofphi(i,j),freq,x,alpha,beta,time)
        
      enddo
      enddo
        return
        end
        subroutine GETLOFSHPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,lmp
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,beta
     &           ,time
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1)
      REAL_T lmp(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      REAL_T alpha
      REAL_T beta
      REAL_T time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getlofshphipoint(lofphi(i,j),lmp,x,alpha,beta,time)
        
      enddo
      enddo
        return
        end
        subroutine GETMARSHALOFPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getmarshalofphipoint(lofphi(i,j),x)
        
      enddo
      enddo
        return
        end
        subroutine GETMARSHALOFPHIPOINT(
     &           lofphi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T x(0:1)
        lofphi = zero
        return
        end
        subroutine GETMARSHAPHI(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getmarshaphipoint(phi(i,j),x)
        
      enddo
      enddo
        return
        end
        subroutine GETMARSHAPHIPOINT(
     &           phi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T x(0:1)
        
        phi = sin(x(0))
     $       *exp(x(1))
        return
        end
        subroutine GETLOFPHIPOINT(
     &           lofphi
     &           ,freq
     &           ,x
     &           ,alpha
     &           ,beta
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T alpha
      REAL_T beta
      REAL_T time
        real_t fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2)
        call getphipoint(phi, freq, x, time)
        lofphi = fac*phi
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETLOFSHPHIPOINT(
     &           lofphi
     &           ,lmp
     &           ,x
     &           ,alpha
     &           ,beta
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T lmp(0:1)
      REAL_T x(0:1)
      REAL_T alpha
      REAL_T beta
      REAL_T time
        external normalization
        real_t normalization
        external algndr
        real_t algndr
        external bessel
        real_t bessel
        external trigphi
        real_t trigphi
        real_t r(0:CH_SPACEDIM-1)
        integer l,m
        real_t z,phase,phi
        call convert_spherical(x,r)
        l = int(lmp(0))
#if CH_SPACEDIM == 1
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     &        (-l**2)*trigphi(l,r(0),phase)
#elif CH_SPACEDIM == 2
        z = -sqrt(l**2-1.)
        m = 0
        phase = 0.
        phi = normalization(l,m)*
     &        (r(0)**z)*
     &        trigphi(l,r(1),phase)
#elif CH_SPACEDIM == 3
        m = int(lmp(1))
        phase = PI*lmp(2)
        z = cos(r(1))
        phi = normalization(l,m)*
     &        bessel(l,r(0))*
     &        algndr(l,m,z)*
     &        trigphi(m,r(2),phase)
#else
        bogus_spacedim()
#endif
        phi = phi*cos(time)
        lofphi = -phi
#if CH_SPACEDIM == 2
        lofphi = lofphi/(r(0)**2)
#endif
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETDBGPHI(
     &           dbgphi
     &           ,idbgphilo0,idbgphilo1
     &           ,idbgphihi0,idbgphihi1
     &           ,beta
     &           ,ibetalo0,ibetalo1
     &           ,ibetahi0,ibetahi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,time
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idbgphilo0,idbgphilo1
      integer idbgphihi0,idbgphihi1
      REAL_T dbgphi(
     &           idbgphilo0:idbgphihi0,
     &           idbgphilo1:idbgphihi1)
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL_T beta(
     &           ibetalo0:ibetahi0,
     &           ibetalo1:ibetahi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      REAL_T alpha
      REAL_T time
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getdbgphipoint(dbgphi(i,j),
     &        beta(i,j),freq,x,alpha,time)
        
      enddo
      enddo
        return
        end
        subroutine GETDBGPHIPOINT(
     &           dbgphi
     &           ,beta
     &           ,freq
     &           ,x
     &           ,alpha
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dbgphi
      REAL_T beta
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T alpha
      REAL_T time
        real_t gradphi(0:CH_SPACEDIM-1),gradbeta(0:CH_SPACEDIM-1)
        real_t alphaphiplusbetalapphi,gradbetadotgradphi
        call getbetapoint(beta,freq,x,time)
        call getlofphipoint(alphaphiplusbetalapphi,freq,x,alpha,beta,time)
        call getgradbetapoint(gradbeta,freq,x,time)
        call getgradphipoint(gradphi,freq,x,time)
        gradbetadotgradphi = gradbeta(0)*gradphi(0)
     &                               + gradbeta(1)*gradphi(1)
        dbgphi = alphaphiplusbetalapphi
        dbgphi = dbgphi + gradbetadotgradphi
        return
        end
        subroutine GETBETAPOINT(
     &           beta
     &           ,freq
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T beta
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
        beta = CH_SPACEDIM+2
     &       +sin(two*PI*x(0))
     &                + sin(two*PI*x(1))
        return
        end
        subroutine GETGRADBETAPOINT(
     &           gradbeta
     &           ,freq
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradbeta(0:1)
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
        integer idir
        do idir = 0, CH_SPACEDIM-1
            gradbeta(idir) = two*PI*cos(two*PI*x(idir))
        enddo
        return
        end
        subroutine GETBETAGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
        integer idir
        real_t beta
        call getbetapoint(beta,freq,x,time)
        call getgradphipoint(gradphi,freq,x,time)
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETBETAGRADSHPHIPOINT(
     &           gradphi
     &           ,lmp
     &           ,x
     &           ,time
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:1)
      REAL_T lmp(0:1)
      REAL_T x(0:1)
      REAL_T time
        integer idir
        real_t beta
        call getbetapoint(beta,lmp,x,time)
        call getgradshphipoint(gradphi,lmp,x,time)
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETSRC(
     &           src
     &           ,isrclo0,isrclo1
     &           ,isrchi0,isrchi1
     &           ,freq
     &           ,dx
     &           ,time
     &           ,diffconst
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1)
      REAL_T freq(0:1)
      REAL_T dx(0:1)
      REAL_T time
      REAL_T diffconst
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getsrcpoint(src(i,j),freq,x,time,diffconst)
        
      enddo
      enddo
        return
        end
        subroutine GETSRCPOINT(
     &           src
     &           ,freq
     &           ,x
     &           ,time
     &           ,diffconst
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T src
      REAL_T freq(0:1)
      REAL_T x(0:1)
      REAL_T time
      REAL_T diffconst
        real_t fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2)
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1)))
        src = (-fac*diffconst*cos(time) - sin(time))*phi
        return
        end
        subroutine GETSHSRC(
     &           src
     &           ,isrclo0,isrclo1
     &           ,isrchi0,isrchi1
     &           ,lmp
     &           ,dx
     &           ,time
     &           ,diffconst
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1)
      REAL_T lmp(0:1)
      REAL_T dx(0:1)
      REAL_T time
      REAL_T diffconst
      REAL_T problo(0:1)
      REAL_T probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        real_t x(0:CH_SPACEDIM-1)
        
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          call getshsrcpoint(src(i,j),lmp,x,time,diffconst)
        
      enddo
      enddo
        return
        end
        subroutine GETSHSRCPOINT(
     &           src
     &           ,lmp
     &           ,x
     &           ,time
     &           ,diffconst
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T src
      REAL_T lmp(0:1)
      REAL_T x(0:1)
      REAL_T time
      REAL_T diffconst
        external normalization
        real_t normalization
        external algndr
        real_t algndr
        external bessel
        real_t bessel
        external trigphi
        real_t trigphi
        real_t rad,z
        real_t r(0:CH_SPACEDIM-1)
        integer l,m
        call convert_spherical(x,r)
        rad = r(0)
        l = int(lmp(0))
#if CH_SPACEDIM == 1
        m = 0
        src = -diffconst*normalization(l,m)*
     &                   bessel(l,rad)
#elif CH_SPACEDIM == 2
        m = 0
        z = cos(r(1))
        src = -diffconst*normalization(l,m)*
     &                   bessel(l,rad)*
     &                   algndr(l,m,z)
#elif CH_SPACEDIM == 3
        m = int(lmp(1))
        phase = PI*lmp(2)
        z = cos(r(1))
        src = -diffconst*normalization(l,m)*
     &                   bessel(l,rad)*
     &                   algndr(l,m,z)*
     &                   trigphi(m,r(2),phase)
#else
        bogus_spacedim()
#endif
        return
        end
        real_t  function normalization(l,m)
        implicit none
        integer l,m
        external factorial
        integer factorial
#if CH_SPACEDIM == 1
        normalization = 1./(four*PI)
#elif CH_SPACEDIM == 2
        normalization = (2*l+1.)*factorial(l)*0.5/factorial(l)/
     &                  (two*PI)
#elif CH_SPACEDIM == 3
        normalization = (2*l+1.)*factorial(l-m)*0.5/factorial(l+m)/
     &                  (two*PI)
        if (m.gt.0) then
           normalization = two*normalization
        endif
#else
        bogus_spacedim()
#endif
        normalization = sqrt(normalization)
        return
        end
        real_t  function algndr(l,m,z)
        implicit none
        integer l,m
        real_t z
        integer i,ll
        real_t fac,pll,pmm,pmmp1,somz2
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
        real_t  function d_algndr(l,m,z)
        implicit none
        external algndr
        real_t algndr
        integer l,m
        real_t z
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
     &           x
     &           ,r
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T x(0:1)
      REAL_T r(0:1)
#if CH_SPACEDIM == 1
        r(0) = abs(x(0))
#elif CH_SPACEDIM == 2
        r(0) = sqrt(x(0)**2 + x(1)**2)
        r(1) = 0.
        if (x(0).ne.0. .or. x(1).ne.0.) then
           r(1) = atan2(x(1),x(0))
           if (r(1).lt.0.) r(1) = r(1) + two*PI
        endif
#elif CH_SPACEDIM == 3
        r(0) = sqrt(x(0)**2 + x(1)**2 + x(2)**2)
        if (r(0)>0.) then
           r(1) = acos(x(2)/r(0))
        else
           r(1) = 0.
        endif
        r(2) = 0.
        if (x(0).ne.0. .or. x(1).ne.0.) then
           r(2) = atan2(x(1),x(0))
           if (r(2).lt.0.) r(2) = r(2) + two*PI
        endif
#else
#endif
        return
        end
        integer  function factorial(n)
        implicit none
        integer n
        integer i
        factorial=1
        do i=n,2,-1
           factorial = factorial*i
        enddo
        return
        end
        real_t  function trigphi(m, azphi, phase)
        implicit none
        integer m
        real_t azphi, phase
        trigphi = cos(m*azphi - phase)/sqrt(PI)
        return
        end
        real_t  function d_trigphi(m, azphi, phase)
        implicit none
        integer m
        real_t azphi, phase
        d_trigphi = -m*sin(m*azphi - phase)/sqrt(PI)
        return
        end
        real_t  function bessel(l, r)
        implicit none
        integer l
        real_t r
        real_t ir,s,c
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
     &                                 (15.*ir**2-1.)*c*ir
              else
                 call MayDay_Error()
              endif
           endif
        endif
        return
        end
        real_t  function d_bessel(l, r)
        implicit none
        integer l
        real_t r
        real_t ir,irsq,s,c
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
     &                                ( 9.*irsq-1.)*c*ir
           elseif (l.eq.3) then
              d_bessel = (-60.*irsq**2+15.*irsq-1.)*s*ir +
     &                                (60.*irsq-7.)*c*irsq
           else
              call MayDay_Error()
           endif
        endif
        return
        end
