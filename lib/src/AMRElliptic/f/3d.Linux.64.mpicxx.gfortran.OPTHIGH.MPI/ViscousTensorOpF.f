      subroutine EXTRAPTOGHOSTVEL(
     & vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,nvelcomp
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,facedir
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2,
     & 0:nvelcomp-1)
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      integer facedir
      integer ncomp
      integer ii,i,jj,j,kk,k, icomp
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      if(hashi.eq. 1) then
         do icomp = 0, ncomp-1
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            vel(i,j,k,icomp) =
     $ 2*vel(i- ii,j- jj,k- kk, icomp)
     $ - vel(i-2*ii,j-2*jj,k-2*kk, icomp)
      enddo
      enddo
      enddo
         enddo
      endif
      if(haslo .eq. 1) then
         do icomp = 0, ncomp-1
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            vel(i,j,k,icomp) =
     $ 2*vel(i+ ii,j+ jj,k+ kk, icomp)
     $ - vel(i+2*ii,j+2*jj,k+2*kk, icomp)
      enddo
      enddo
      enddo
         enddo
      endif
      return
      end
      subroutine APPLYOPVTOPNOBCS(
     & lphfab
     & ,ilphfablo0,ilphfablo1,ilphfablo2
     & ,ilphfabhi0,ilphfabhi1,ilphfabhi2
     & ,nlphfabcomp
     & ,phifab
     & ,iphifablo0,iphifablo1,iphifablo2
     & ,iphifabhi0,iphifabhi1,iphifabhi2
     & ,nphifabcomp
     & ,acofab
     & ,iacofablo0,iacofablo1,iacofablo2
     & ,iacofabhi0,iacofabhi1,iacofabhi2
     & ,eta0fab
     & ,ieta0fablo0,ieta0fablo1,ieta0fablo2
     & ,ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
     & ,eta1fab
     & ,ieta1fablo0,ieta1fablo1,ieta1fablo2
     & ,ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
     & ,eta2fab
     & ,ieta2fablo0,ieta2fablo1,ieta2fablo2
     & ,ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
     & ,lam0fab
     & ,ilam0fablo0,ilam0fablo1,ilam0fablo2
     & ,ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
     & ,lam1fab
     & ,ilam1fablo0,ilam1fablo1,ilam1fablo2
     & ,ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
     & ,lam2fab
     & ,ilam2fablo0,ilam2fablo1,ilam2fablo2
     & ,ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
     & ,dx
     & ,alpha
     & ,beta
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlphfabcomp
      integer ilphfablo0,ilphfablo1,ilphfablo2
      integer ilphfabhi0,ilphfabhi1,ilphfabhi2
      REAL*8 lphfab(
     & ilphfablo0:ilphfabhi0,
     & ilphfablo1:ilphfabhi1,
     & ilphfablo2:ilphfabhi2,
     & 0:nlphfabcomp-1)
      integer nphifabcomp
      integer iphifablo0,iphifablo1,iphifablo2
      integer iphifabhi0,iphifabhi1,iphifabhi2
      REAL*8 phifab(
     & iphifablo0:iphifabhi0,
     & iphifablo1:iphifabhi1,
     & iphifablo2:iphifabhi2,
     & 0:nphifabcomp-1)
      integer iacofablo0,iacofablo1,iacofablo2
      integer iacofabhi0,iacofabhi1,iacofabhi2
      REAL*8 acofab(
     & iacofablo0:iacofabhi0,
     & iacofablo1:iacofabhi1,
     & iacofablo2:iacofabhi2)
      integer ieta0fablo0,ieta0fablo1,ieta0fablo2
      integer ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
      REAL*8 eta0fab(
     & ieta0fablo0:ieta0fabhi0,
     & ieta0fablo1:ieta0fabhi1,
     & ieta0fablo2:ieta0fabhi2)
      integer ieta1fablo0,ieta1fablo1,ieta1fablo2
      integer ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
      REAL*8 eta1fab(
     & ieta1fablo0:ieta1fabhi0,
     & ieta1fablo1:ieta1fabhi1,
     & ieta1fablo2:ieta1fabhi2)
      integer ieta2fablo0,ieta2fablo1,ieta2fablo2
      integer ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
      REAL*8 eta2fab(
     & ieta2fablo0:ieta2fabhi0,
     & ieta2fablo1:ieta2fabhi1,
     & ieta2fablo2:ieta2fabhi2)
      integer ilam0fablo0,ilam0fablo1,ilam0fablo2
      integer ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
      REAL*8 lam0fab(
     & ilam0fablo0:ilam0fabhi0,
     & ilam0fablo1:ilam0fabhi1,
     & ilam0fablo2:ilam0fabhi2)
      integer ilam1fablo0,ilam1fablo1,ilam1fablo2
      integer ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
      REAL*8 lam1fab(
     & ilam1fablo0:ilam1fabhi0,
     & ilam1fablo1:ilam1fabhi1,
     & ilam1fablo2:ilam1fabhi2)
      integer ilam2fablo0,ilam2fablo1,ilam2fablo2
      integer ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
      REAL*8 lam2fab(
     & ilam2fablo0:ilam2fabhi0,
     & ilam2fablo1:ilam2fabhi1,
     & ilam2fablo2:ilam2fabhi2)
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 divf(0:3 -1)
      REAL*8 aphi(0:3 -1)
      REAL*8 etaL(0:3 -1)
      REAL*8 etaH(0:3 -1)
      REAL*8 lamL(0:3 -1)
      REAL*8 lamH(0:3 -1)
      REAL*8 fluxL(0:3 -1,0:3 -1)
      REAL*8 fluxH(0:3 -1,0:3 -1)
      REAL*8 gphiL(0:3 -1,0:3 -1,0:3 -1)
      REAL*8 gphiH(0:3 -1,0:3 -1,0:3 -1)
      REAL*8 divuL(0:3 -1)
      REAL*8 divuH(0:3 -1)
      integer i,j,k, facedir, derivdir, veldir
      integer iif,jjf,kkf
      integer iid,jjd,kkd
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      etaL(0) =eta0fab(i ,j ,k )
      etaL(1) =eta1fab(i ,j ,k )
      etaL(2) =eta2fab(i ,j ,k )
      etaH(0) =eta0fab(i+1 ,j ,k )
      etaH(1) =eta1fab(i ,j+1 ,k )
      etaH(2) =eta2fab(i ,j ,k+1 )
      lamL(0) =lam0fab(i ,j ,k )
      lamL(1) =lam1fab(i ,j ,k )
      lamL(2) =lam2fab(i ,j ,k )
      lamH(0) =lam0fab(i+1 ,j ,k )
      lamH(1) =lam1fab(i ,j+1 ,k )
      lamH(2) =lam2fab(i ,j ,k+1 )
      do facedir = 0, 3 -1
         iif = chf_id(facedir, 0)
         jjf = chf_id(facedir, 1)
         kkf = chf_id(facedir, 2)
         do derivdir = 0, 3 -1
            iid = chf_id(derivdir, 0)
            jjd = chf_id(derivdir, 1)
            kkd = chf_id(derivdir, 2)
            do veldir = 0, 3 -1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phifab(i+iid,j+jjd
     &,k+kkd,veldir) - phifab(i ,j ,k ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phifab(i ,j ,k ,ve
     &ldir) - phifab(i-iid,j-jjd,k-kkd,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = ((1.0d0)/((4.0d0)*d
     &x))*(
     $ phifab(i+iid+iif,j+jjd+jjf,k+kkd+kkf,veldir) - phifab(i-iid+iif,j
     &-jjd+jjf,k-kkd+kkf,veldir) +
     $ phifab(i+iid ,j+jjd ,k+kkd ,veldir) - phifab(i-iid ,j-jjd ,k-kkd 
     &,veldir) )
                  gphiL(veldir, derivdir ,facedir) = ((1.0d0)/((4.0d0)*d
     &x))*(
     $ phifab(i+iid-iif,j+jjd-jjf,k+kkd-kkf,veldir) - phifab(i-iid-iif,j
     &-jjd-jjf,k-kkd-kkf,veldir) +
     $ phifab(i+iid ,j+jjd ,k+kkd ,veldir) - phifab(i-iid ,j-jjd ,k-kkd 
     &,veldir) )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, 3 -1
         divuL(facedir) = (0.0d0)
         divuH(facedir) = (0.0d0)
         do veldir = 0, 3 -1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, faced
     &ir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, faced
     &ir)
         enddo
      enddo
      do facedir = 0, 3 -1
         do veldir = 0, 3 -1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldi
     &r, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldi
     &r, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(fa
     &cedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(fa
     &cedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, 3 -1
         divf(veldir) = (0.0d0)
         do facedir = 0, 3 -1
            divf(veldir) = divf(veldir)
     $ + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      aphi(0) = acofab(i,j,k)*phifab(i,j,k, 0)
      aphi(1) = acofab(i,j,k)*phifab(i,j,k, 1)
      aphi(2) = acofab(i,j,k)*phifab(i,j,k, 2)
      lphfab(i,j,k,0) = alpha*aphi(0) + beta*divf(0)
      lphfab(i,j,k,1) = alpha*aphi(1) + beta*divf(1)
      lphfab(i,j,k,2) = alpha*aphi(2) + beta*divf(2)
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionl
     &o1+1)*(iregionhi2- iregionlo2+1)*(4 + 6*3 + 12*(3*3))
      return
      end
      subroutine GETFLUXVTOPNOBCS(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,phifab
     & ,iphifablo0,iphifablo1,iphifablo2
     & ,iphifabhi0,iphifabhi1,iphifabhi2
     & ,nphifabcomp
     & ,etafab
     & ,ietafablo0,ietafablo1,ietafablo2
     & ,ietafabhi0,ietafabhi1,ietafabhi2
     & ,lamfab
     & ,ilamfablo0,ilamfablo1,ilamfablo2
     & ,ilamfabhi0,ilamfabhi1,ilamfabhi2
     & ,dx
     & ,facedir
     & ,beta
     & ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     & ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
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
      integer nphifabcomp
      integer iphifablo0,iphifablo1,iphifablo2
      integer iphifabhi0,iphifabhi1,iphifabhi2
      REAL*8 phifab(
     & iphifablo0:iphifabhi0,
     & iphifablo1:iphifabhi1,
     & iphifablo2:iphifabhi2,
     & 0:nphifabcomp-1)
      integer ietafablo0,ietafablo1,ietafablo2
      integer ietafabhi0,ietafabhi1,ietafabhi2
      REAL*8 etafab(
     & ietafablo0:ietafabhi0,
     & ietafablo1:ietafabhi1,
     & ietafablo2:ietafabhi2)
      integer ilamfablo0,ilamfablo1,ilamfablo2
      integer ilamfabhi0,ilamfabhi1,ilamfabhi2
      REAL*8 lamfab(
     & ilamfablo0:ilamfabhi0,
     & ilamfablo1:ilamfabhi1,
     & ilamfablo2:ilamfabhi2)
      REAL*8 dx
      integer facedir
      REAL*8 beta
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      REAL*8 gphi(0:3 -1,0:3 -1)
      REAL*8 divu, eta, lam, fluxpt
      integer i,j,k , derivdir, veldir
      integer iif,jjf,kkf
      integer iid,jjd,kkd
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      eta = etafab(i ,j ,k );
      lam = lamfab(i ,j ,k );
      iif = chf_id(facedir, 0)
      jjf = chf_id(facedir, 1)
      kkf = chf_id(facedir, 2)
      do derivdir = 0, 3 -1
         iid = chf_id(derivdir, 0)
         jjd = chf_id(derivdir, 1)
         kkd = chf_id(derivdir, 2)
         do veldir = 0, 3 -1
            if(facedir .eq. derivdir) then
               gphi(veldir, derivdir) = (phifab(i,j,k,veldir) - phifab(i
     &-iid,j-jjd,k-kkd,veldir))/dx
            else
               gphi(veldir, derivdir) = ((1.0d0)/((4.0d0)*dx))*(
     $ phifab(i+iid-iif,j+jjd-jjf,k+kkd-kkf,veldir) - phifab(i-iid-iif,j
     &-jjd-jjf,k-kkd-kkf,veldir) +
     $ phifab(i+iid ,j+jjd ,k+kkd ,veldir) - phifab(i-iid ,j-jjd ,k-kkd 
     &,veldir) )
            endif
         enddo
      enddo
      divu = (0.0d0)
      do veldir = 0, 3 -1
         divu = divu + gphi(veldir, veldir)
      enddo
      do veldir = 0, 3 -1
         fluxpt = eta*(gphi(facedir, veldir) + gphi(veldir, facedir))
         if(veldir .eq. facedir) then
            fluxpt = fluxpt + lam*divu
         endif
         flux(i,j,k, veldir) = beta*fluxpt
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONGVTOP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,ncoarsecomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,m
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
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2,
     & 0:ncoarsecomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer m
      integer ncomp, n
      integer i,j,k
      integer ii,jj,kk
      ncomp = nphicomp
      do n = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          ii = (i-abs(mod(i,m)))/m
          jj = (j-abs(mod(j,m)))/m
          kk = (k-abs(mod(k,m)))/m
          phi(i,j,k,n) = phi(i,j,k,n) +
     & coarse(ii,jj,kk,n)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESVTOP(
     & res
     & ,ireslo0,ireslo1,ireslo2
     & ,ireshi0,ireshi1,ireshi2
     & ,nrescomp
     & ,resfine
     & ,iresfinelo0,iresfinelo1,iresfinelo2
     & ,iresfinehi0,iresfinehi1,iresfinehi2
     & ,nresfinecomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0,ireslo1,ireslo2
      integer ireshi0,ireshi1,ireshi2
      REAL*8 res(
     & ireslo0:ireshi0,
     & ireslo1:ireshi1,
     & ireslo2:ireshi2,
     & 0:nrescomp-1)
      integer nresfinecomp
      integer iresfinelo0,iresfinelo1,iresfinelo2
      integer iresfinehi0,iresfinehi1,iresfinehi2
      REAL*8 resfine(
     & iresfinelo0:iresfinehi0,
     & iresfinelo1:iresfinehi1,
     & iresfinelo2:iresfinehi2,
     & 0:nresfinecomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer ncomp
      REAL*8 denom
      integer n
      integer i,j,k
      integer ii,jj,kk
      denom = (2.0d0) *(2.0d0) *(2.0d0)
      do n = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         ii = (i-abs(mod(i,2)))/2
         jj = (j-abs(mod(j,2)))/2
         kk = (k-abs(mod(k,2)))/2
         res(ii,jj,kk,n) = res(ii,jj,kk,n) + resfine(i,j,k,n)/denom
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine CELLGRADVTOP(
     & grad
     & ,igradlo0,igradlo1,igradlo2
     & ,igradhi0,igradhi1,igradhi2
     & ,vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,igridlo0,igridlo1,igridlo2
     & ,igridhi0,igridhi1,igridhi2
     & ,dx
     & ,divdir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1,
     & igradlo2:igradhi2)
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2)
      integer igridlo0,igridlo1,igridlo2
      integer igridhi0,igridhi1,igridhi2
      REAL*8 dx
      integer divdir
      integer ii,i,jj,j,kk,k
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      kk = chf_id(divdir, 2)
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
      grad(i,j,k) =
     $ ( vel(i+ii,j+jj,k+kk)
     $ - vel(i-ii,j-jj,k-kk) )/((2.0d0)*dx)
      enddo
      enddo
      enddo
      return
      end
      subroutine ADDGRADTOFLUXVTOP(
     & flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,eta
     & ,ietalo0,ietalo1,ietalo2
     & ,ietahi0,ietahi1,ietahi2
     & ,fluxcomp
     & ,grad
     & ,igradlo0,igradlo1,igradlo2
     & ,igradhi0,igradhi1,igradhi2
     & ,ngradcomp
     & ,gradcomp
     & ,gradtran
     & ,iregionfacelo0,iregionfacelo1,iregionfacelo2
     & ,iregionfacehi0,iregionfacehi1,iregionfacehi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer ietalo0,ietalo1,ietalo2
      integer ietahi0,ietahi1,ietahi2
      REAL*8 eta(
     & ietalo0:ietahi0,
     & ietalo1:ietahi1,
     & ietalo2:ietahi2)
      integer fluxcomp
      integer ngradcomp
      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1,
     & igradlo2:igradhi2,
     & 0:ngradcomp-1)
      integer gradcomp
      integer gradtran
      integer iregionfacelo0,iregionfacelo1,iregionfacelo2
      integer iregionfacehi0,iregionfacehi1,iregionfacehi2
      integer i,j,k
      REAL*8 gradcontrib, trancontrib, etafac
      do k = iregionfacelo2,iregionfacehi2
      do j = iregionfacelo1,iregionfacehi1
      do i = iregionfacelo0,iregionfacehi0
      etafac = eta(i,j,k)
      gradcontrib = grad(i,j,k, gradcomp)
      trancontrib = grad(i,j,k, gradtran)
      flux(i,j,k, fluxcomp) = flux(i,j,k, fluxcomp) +
     $ eta(i,j,k)*
     $ ( grad(i,j,k, gradcomp)
     $ + grad(i,j,k, gradtran) )
      enddo
      enddo
      enddo
      return
      end
      subroutine GETFACEGRADVTOP(
     & gradvelface
     & ,igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
     & ,igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
     & ,gradvelcell
     & ,igradvelcelllo0,igradvelcelllo1,igradvelcelllo2
     & ,igradvelcellhi0,igradvelcellhi1,igradvelcellhi2
     & ,velcomp
     & ,ivelcomplo0,ivelcomplo1,ivelcomplo2
     & ,ivelcomphi0,ivelcomphi1,ivelcomphi2
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,dx
     & ,facedir
     & ,divdir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
      integer igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
      REAL*8 gradvelface(
     & igradvelfacelo0:igradvelfacehi0,
     & igradvelfacelo1:igradvelfacehi1,
     & igradvelfacelo2:igradvelfacehi2)
      integer igradvelcelllo0,igradvelcelllo1,igradvelcelllo2
      integer igradvelcellhi0,igradvelcellhi1,igradvelcellhi2
      REAL*8 gradvelcell(
     & igradvelcelllo0:igradvelcellhi0,
     & igradvelcelllo1:igradvelcellhi1,
     & igradvelcelllo2:igradvelcellhi2)
      integer ivelcomplo0,ivelcomplo1,ivelcomplo2
      integer ivelcomphi0,ivelcomphi1,ivelcomphi2
      REAL*8 velcomp(
     & ivelcomplo0:ivelcomphi0,
     & ivelcomplo1:ivelcomphi1,
     & ivelcomplo2:ivelcomphi2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      REAL*8 dx
      integer facedir
      integer divdir
      integer ii,i,jj,j,kk,k
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      if (facedir .eq. divdir) then
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         gradvelface(i,j,k) =
     $ ( velcomp(i ,j ,k )
     $ - velcomp(i-ii,j-jj,k-kk) )/dx
      enddo
      enddo
      enddo
      else
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         gradvelface(i,j,k) =
     $ ( gradvelcell(i ,j ,k )
     $ + gradvelcell(i-ii,j-jj,k-kk) )/(2.0d0)
      enddo
      enddo
      enddo
         if(haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            gradvelface(i,j,k) =
     $ ((3.0d0)*gradvelcell(i ,j ,k )
     $ - gradvelcell(i+ii,j+jj,k+kk))/(2.0d0)
      enddo
      enddo
      enddo
         endif
         if(hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            gradvelface(i,j,k) =
     $ ((3.0d0)*gradvelcell(i- ii,j- jj,k- kk)
     $ - gradvelcell(i-2*ii,j-2*jj,k-2*kk))/(2.0d0)
      enddo
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine CELLDIVINCRVTOP(
     & divvel
     & ,idivvello0,idivvello1,idivvello2
     & ,idivvelhi0,idivvelhi1,idivvelhi2
     & ,vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,nvelcomp
     & ,dx
     & ,divdir
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivvello0,idivvello1,idivvello2
      integer idivvelhi0,idivvelhi1,idivvelhi2
      REAL*8 divvel(
     & idivvello0:idivvelhi0,
     & idivvello1:idivvelhi1,
     & idivvello2:idivvelhi2)
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2,
     & 0:nvelcomp-1)
      REAL*8 dx
      integer divdir
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer ii,i,jj,j,kk,k
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      kk = chf_id(divdir, 2)
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      divvel(i,j,k) = divvel(i,j,k) +
     $ ( vel(i+ii,j+jj,k+kk,divdir)
     $ - vel(i-ii,j-jj,k-kk,divdir) )/((2.0d0)*dx)
      enddo
      enddo
      enddo
      return
      end
      subroutine FACEDIVINCRVTOP(
     & divvel
     & ,idivvello0,idivvello1,idivvello2
     & ,idivvelhi0,idivvelhi1,idivvelhi2
     & ,vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,nvelcomp
     & ,gradvel
     & ,igradvello0,igradvello1,igradvello2
     & ,igradvelhi0,igradvelhi1,igradvelhi2
     & ,ngradvelcomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     & ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     & ,iloboxlo0,iloboxlo1,iloboxlo2
     & ,iloboxhi0,iloboxhi1,iloboxhi2
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     & ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     & ,hashi
     & ,dx
     & ,facedir
     & ,divdir
     & ,gradcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idivvello0,idivvello1,idivvello2
      integer idivvelhi0,idivvelhi1,idivvelhi2
      REAL*8 divvel(
     & idivvello0:idivvelhi0,
     & idivvello1:idivvelhi1,
     & idivvello2:idivvelhi2)
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2,
     & 0:nvelcomp-1)
      integer ngradvelcomp
      integer igradvello0,igradvello1,igradvello2
      integer igradvelhi0,igradvelhi1,igradvelhi2
      REAL*8 gradvel(
     & igradvello0:igradvelhi0,
     & igradvello1:igradvelhi1,
     & igradvello2:igradvelhi2,
     & 0:ngradvelcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      REAL*8 dx
      integer facedir
      integer divdir
      integer gradcomp
      integer ii,i,jj,j,kk,k
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      if (facedir .eq. divdir) then
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         divvel(i,j,k) = divvel(i,j,k) +
     $ ( vel(i ,j ,k ,facedir)
     $ - vel(i-ii,j-jj,k-kk,facedir) )/dx
      enddo
      enddo
      enddo
      else
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         divvel(i,j,k) = divvel(i,j,k) +
     $ ( gradvel(i ,j ,k , gradcomp)
     $ + gradvel(i-ii,j-jj,k-kk, gradcomp) )/(2.0d0)
      enddo
      enddo
      enddo
         if(haslo .eq. 1) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            divvel(i,j,k) = divvel(i,j,k) +
     $ ((3.0d0)*gradvel(i ,j ,k , gradcomp)
     $ - gradvel(i+ii,j+jj,k+kk, gradcomp))/(2.0d0)
      enddo
      enddo
      enddo
         endif
         if(hashi .eq. 1) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            divvel(i,j,k) = divvel(i,j,k) +
     $ ((3.0d0)*gradvel(i- ii,j- jj,k- kk, gradcomp)
     $ - gradvel(i-2*ii,j-2*jj,k-2*kk, gradcomp))/(2.0d0)
      enddo
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine DECRINVRELCOEFVTOP(
     & relcoef
     & ,irelcoeflo0,irelcoeflo1,irelcoeflo2
     & ,irelcoefhi0,irelcoefhi1,irelcoefhi2
     & ,nrelcoefcomp
     & ,eta
     & ,ietalo0,ietalo1,ietalo2
     & ,ietahi0,ietahi1,ietahi2
     & ,netacomp
     & ,lambda
     & ,ilambdalo0,ilambdalo1,ilambdalo2
     & ,ilambdahi0,ilambdahi1,ilambdahi2
     & ,nlambdacomp
     & ,beta
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
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
      integer nrelcoefcomp
      integer irelcoeflo0,irelcoeflo1,irelcoeflo2
      integer irelcoefhi0,irelcoefhi1,irelcoefhi2
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1,
     & irelcoeflo2:irelcoefhi2,
     & 0:nrelcoefcomp-1)
      integer netacomp
      integer ietalo0,ietalo1,ietalo2
      integer ietahi0,ietahi1,ietahi2
      REAL*8 eta(
     & ietalo0:ietahi0,
     & ietalo1:ietahi1,
     & ietalo2:ietahi2,
     & 0:netacomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & ilambdalo2:ilambdahi2,
     & 0:nlambdacomp-1)
      REAL*8 beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      REAL*8 dx
      integer idir
      integer ncomp
      integer ii,jj,kk
      integer i,j,k
      integer icomp
      REAL*8 lamh, laml,etah, etal, relcoold, relconew, incr
      ii = chf_id(idir, 0)
      jj = chf_id(idir, 1)
      kk = chf_id(idir, 2)
      do icomp = 0, ncomp-1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         relcoold = relcoef(i,j,k,icomp)
         etah = eta(i+ii,j+jj,k+kk,0)
         etal = eta(i ,j ,k ,0)
         lamh = lambda(i+ii,j+jj,k+kk,0)
         laml = lambda(i ,j ,k ,0)
         if(icomp .eq. idir) then
            incr =
     $ beta*(
     $ (2.0d0)*eta(i+ii,j+jj,k+kk,0) +
     $ lambda (i+ii,j+jj,k+kk,0) +
     $ (2.0d0)*eta(i ,j ,k ,0) +
     $ lambda (i ,j ,k ,0)
     $ )/(dx*dx)
         else
            incr =
     $ beta*(
     $ eta(i+ii,j+jj,k+kk,0) +
     $ eta(i ,j ,k ,0)
     $ )/(dx*dx)
         endif
         relcoef(i,j,k, icomp) = relcoef(i,j,k,icomp) - incr
         relconew = relcoef(i,j,k,icomp)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INITIALIZERELAXCOEF(
     & relcoef
     & ,irelcoeflo0,irelcoeflo1,irelcoeflo2
     & ,irelcoefhi0,irelcoefhi1,irelcoefhi2
     & ,nrelcoefcomp
     & ,acoef
     & ,iacoeflo0,iacoeflo1,iacoeflo2
     & ,iacoefhi0,iacoefhi1,iacoefhi2
     & ,alpha
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrelcoefcomp
      integer irelcoeflo0,irelcoeflo1,irelcoeflo2
      integer irelcoefhi0,irelcoefhi1,irelcoefhi2
      REAL*8 relcoef(
     & irelcoeflo0:irelcoefhi0,
     & irelcoeflo1:irelcoefhi1,
     & irelcoeflo2:irelcoefhi2,
     & 0:nrelcoefcomp-1)
      integer iacoeflo0,iacoeflo1,iacoeflo2
      integer iacoefhi0,iacoefhi1,iacoefhi2
      REAL*8 acoef(
     & iacoeflo0:iacoefhi0,
     & iacoeflo1:iacoefhi1,
     & iacoeflo2:iacoefhi2)
      REAL*8 alpha
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer ncomp
      integer i,j,k
      integer icomp
      do icomp = 0, ncomp-1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         relcoef(i,j,k, icomp) = alpha*acoef(i,j,k)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAVTOP(
     & lambda
     & ,ilambdalo0,ilambdalo1,ilambdalo2
     & ,ilambdahi0,ilambdahi1,ilambdahi2
     & ,nlambdacomp
     & ,safety
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & ilambdalo2:ilambdahi2,
     & 0:nlambdacomp-1)
      REAL*8 safety
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer ncomp
      integer i,j,k
      integer icomp,jcomp
      REAL*8 zeroval
      zeroval = 1.0e-20
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
         do icomp = 0, ncomp-1
            lambda(i,j,k, icomp) = safety/lambda(i,j,k, icomp)
         enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBVTOP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,lphi
     & ,ilphilo0,ilphilo1,ilphilo2
     & ,ilphihi0,ilphihi1,ilphihi2
     & ,nlphicomp
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,nrhscomp
     & ,lambda
     & ,ilambdalo0,ilambdalo1,ilambdalo2
     & ,ilambdahi0,ilambdahi1,ilambdahi2
     & ,nlambdacomp
     & ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     & ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
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
      integer nlphicomp
      integer ilphilo0,ilphilo1,ilphilo2
      integer ilphihi0,ilphihi1,ilphihi2
      REAL*8 lphi(
     & ilphilo0:ilphihi0,
     & ilphilo1:ilphihi1,
     & ilphilo2:ilphihi2,
     & 0:nlphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & 0:nrhscomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & ilambdalo2:ilambdahi2,
     & 0:nlambdacomp-1)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer ncomp
      integer i,j,k
      integer icomp
      REAL*8 phio, lamo, rhso, lphio
      do icomp = 0, ncomp-1
      do k = icoloredboxlo2,icoloredboxhi2,2
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2
         phio = phi( i,j,k,icomp)
         lamo = lambda(i,j,k,icomp)
         rhso = rhs( i,j,k,icomp)
         lphio = lphi( i,j,k,icomp)
         phi(i,j,k, icomp) =
     $ phi( i,j,k,icomp) +
     & lambda(i,j,k,icomp)*(
     $ rhs( i,j,k,icomp) -
     $ lphi( i,j,k,icomp))
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine ADDDIVFLUXDIRVTOP(
     & lhs
     & ,ilhslo0,ilhslo1,ilhslo2
     & ,ilhshi0,ilhshi1,ilhshi2
     & ,nlhscomp
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,dx
     & ,ncomp
     & ,facedir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlhscomp
      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1,
     & ilhslo2:ilhshi2,
     & 0:nlhscomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 dx
      integer ncomp
      integer facedir
      integer ii,i,jj,j,kk,k
      integer icomp
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      do icomp = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
         lhs(i,j,k, icomp) = lhs(i,j,k, icomp)
     $ +
     $ (flux(i+ii,j+jj,k+kk, icomp)
     $ -flux(i ,j ,k , icomp))/dx
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine LINEAREXTRAPVTOP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,ighostBoxlo0,ighostBoxlo1,ighostBoxlo2
     & ,ighostBoxhi0,ighostBoxhi1,ighostBoxhi2
     & ,dir
     & ,hiLo
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer ighostBoxlo0,ighostBoxlo1,ighostBoxlo2
      integer ighostBoxhi0,ighostBoxhi1,ighostBoxhi2
      integer dir
      integer hiLo
      integer n
      integer i0,i1,i2
      integer ii0,ii1,ii2
      ii0= hiLo*CHF_ID(0, dir)
      ii1= hiLo*CHF_ID(1, dir)
      ii2= hiLo*CHF_ID(2, dir)
      do n=0, nphicomp-1
      do i2 = ighostBoxlo2,ighostBoxhi2
      do i1 = ighostBoxlo1,ighostBoxhi1
      do i0 = ighostBoxlo0,ighostBoxhi0
           phi(i0,i1,i2,n) = (2.0d0)*phi(i0-ii0,i1-ii1,i2-ii2,n)
     & - phi(i0-2*ii0,i1-2*ii1,i2-2*ii2,n)
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SLOPESVTOP(
     & slopes
     & ,islopeslo0,islopeslo1,islopeslo2
     & ,islopeshi0,islopeshi1,islopeshi2
     & ,nslopescomp
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,ncoarsecomp
     & ,icBoxlo0,icBoxlo1,icBoxlo2
     & ,icBoxhi0,icBoxhi1,icBoxhi2
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nslopescomp
      integer islopeslo0,islopeslo1,islopeslo2
      integer islopeshi0,islopeshi1,islopeshi2
      REAL*8 slopes(
     & islopeslo0:islopeshi0,
     & islopeslo1:islopeshi1,
     & islopeslo2:islopeshi2,
     & 0:nslopescomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2,
     & 0:ncoarsecomp-1)
      integer icBoxlo0,icBoxlo1,icBoxlo2
      integer icBoxhi0,icBoxhi1,icBoxhi2
      integer dir
      integer n, i0,i1,i2, ii0,ii1,ii2
      ii0= 1*CHF_ID(0, dir)
      ii1= 1*CHF_ID(1, dir)
      ii2= 1*CHF_ID(2, dir)
      do n=0, ncoarsecomp-1
      do i2 = icBoxlo2,icBoxhi2
      do i1 = icBoxlo1,icBoxhi1
      do i0 = icBoxlo0,icBoxhi0
         slopes(i0,i1,i2,n) = 0.5*(coarse(i0+ii0,i1+ii1,i2+ii2,n)
     & -coarse(i0-ii0,i1-ii1,i2-ii2,n))
      enddo
      enddo
      enddo
      enddo
      return
      end
