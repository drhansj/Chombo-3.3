      subroutine ADVECTSLOPEDIFFS(
     & deltawc
     & ,ideltawclo0,ideltawclo1,ideltawclo2
     & ,ideltawchi0,ideltawchi1,ideltawchi2
     & ,ndeltawccomp
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
     & ,iuselimiting
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
      integer iuselimiting
      integer i ,j ,k ,lvar
      integer ioff,joff,koff
      REAL*8 dwr,dwl, dwc
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
         dwc = (0.500d0)*(dwr + dwl)
         if(iuselimiting .eq. 1) then
            call pointvllimiter(dwc, dwl, dwr)
         endif
         deltawc(i,j,k,lvar) = dwc
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
      enddo
      enddo
      enddo
         endif
      enddo
      ch_flops=ch_flops+((icenterboxhi0- icenterboxlo0+1)*(icenterboxhi1
     &- icenterboxlo1+1)*(icenterboxhi2- icenterboxlo2+1)*4+(iloboxhi0- 
     &iloboxlo0+1)*(iloboxhi1- iloboxlo1+1)*(iloboxhi2- iloboxlo2+1)+(ih
     &iboxhi0- ihiboxlo0+1)*(ihiboxhi1- ihiboxlo1+1)*(ihiboxhi2- ihiboxl
     &o2+1))*numslopes
      return
      end
      subroutine PREDADVECT(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,rho
     & ,irholo0,irholo1,irholo2
     & ,irhohi0,irhohi1,irhohi2
     & ,drho
     & ,idrholo0,idrholo1,idrholo2
     & ,idrhohi0,idrhohi1,idrhohi2
     & ,velcc
     & ,ivelcclo0,ivelcclo1,ivelcclo2
     & ,ivelcchi0,ivelcchi1,ivelcchi2
     & ,nvelcccomp
     & ,rholo
     & ,irhololo0,irhololo1,irhololo2
     & ,irholohi0,irholohi1,irholohi2
     & ,rhohi
     & ,irhohilo0,irhohilo1,irhohilo2
     & ,irhohihi0,irhohihi1,irhohihi2
     & ,normdir
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer irholo0,irholo1,irholo2
      integer irhohi0,irhohi1,irhohi2
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1,
     & irholo2:irhohi2)
      integer idrholo0,idrholo1,idrholo2
      integer idrhohi0,idrhohi1,idrhohi2
      REAL*8 drho(
     & idrholo0:idrhohi0,
     & idrholo1:idrhohi1,
     & idrholo2:idrhohi2)
      integer nvelcccomp
      integer ivelcclo0,ivelcclo1,ivelcclo2
      integer ivelcchi0,ivelcchi1,ivelcchi2
      REAL*8 velcc(
     & ivelcclo0:ivelcchi0,
     & ivelcclo1:ivelcchi1,
     & ivelcclo2:ivelcchi2,
     & 0:nvelcccomp-1)
      integer irhololo0,irhololo1,irhololo2
      integer irholohi0,irholohi1,irholohi2
      REAL*8 rholo(
     & irhololo0:irholohi0,
     & irhololo1:irholohi1,
     & irhololo2:irholohi2)
      integer irhohilo0,irhohilo1,irhohilo2
      integer irhohihi0,irhohihi1,irhohihi2
      REAL*8 rhohi(
     & irhohilo0:irhohihi0,
     & irhohilo1:irhohihi1,
     & irhohilo2:irhohihi2)
      integer normdir
      REAL*8 dtbydx
      REAL*8 veloc(0:3 -1)
      integer i,j,k, idir
      REAL*8 dense, denlo, denhi, denslope
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
      dense = rho(i,j,k)
      denslope = drho(i,j,k)
      do idir = 0, 3 -1
         veloc(idir) = velcc(i,j,k,idir)
      enddo
      call pointpredadvect(
     & dense, denlo, denhi, denslope, veloc,
     & normdir, dtbydx)
      rholo(i,j,k) = denlo
      rhohi(i,j,k) = denhi
      enddo
      enddo
      enddo
      return
      end
      subroutine POINTPREDADVECT(
     & dense
     & ,denlo
     & ,denhi
     & ,denslope
     & ,veloc
     & ,normdir
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dense
      REAL*8 denlo
      REAL*8 denhi
      REAL*8 denslope
      REAL*8 veloc(0:2)
      integer normdir
      REAL*8 dtbydx
      REAL*8 velpos, velneg
      REAL*8 tol
      denhi = dense + (0.500d0)*min(((1.0d0)-veloc(normdir)*dtbydx),(1.0
     &d0))*denslope
      denlo = dense - (0.500d0)*min(((1.0d0)+veloc(normdir)*dtbydx),(1.0
     &d0))*denslope
      ch_flops=ch_flops+12
      return
      end
      subroutine PREDADVECTTRANS(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,rho
     & ,irholo0,irholo1,irholo2
     & ,irhohi0,irhohi1,irhohi2
     & ,drho
     & ,idrholo0,idrholo1,idrholo2
     & ,idrhohi0,idrhohi1,idrhohi2
     & ,velcc
     & ,ivelcclo0,ivelcclo1,ivelcclo2
     & ,ivelcchi0,ivelcchi1,ivelcchi2
     & ,nvelcccomp
     & ,rholo
     & ,irhololo0,irhololo1,irhololo2
     & ,irholohi0,irholohi1,irholohi2
     & ,rhohi
     & ,irhohilo0,irhohilo1,irhohilo2
     & ,irhohihi0,irhohihi1,irhohihi2
     & ,tandir
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer irholo0,irholo1,irholo2
      integer irhohi0,irhohi1,irhohi2
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1,
     & irholo2:irhohi2)
      integer idrholo0,idrholo1,idrholo2
      integer idrhohi0,idrhohi1,idrhohi2
      REAL*8 drho(
     & idrholo0:idrhohi0,
     & idrholo1:idrhohi1,
     & idrholo2:idrhohi2)
      integer nvelcccomp
      integer ivelcclo0,ivelcclo1,ivelcclo2
      integer ivelcchi0,ivelcchi1,ivelcchi2
      REAL*8 velcc(
     & ivelcclo0:ivelcchi0,
     & ivelcclo1:ivelcchi1,
     & ivelcclo2:ivelcchi2,
     & 0:nvelcccomp-1)
      integer irhololo0,irhololo1,irhololo2
      integer irholohi0,irholohi1,irholohi2
      REAL*8 rholo(
     & irhololo0:irholohi0,
     & irhololo1:irholohi1,
     & irhololo2:irholohi2)
      integer irhohilo0,irhohilo1,irhohilo2
      integer irhohihi0,irhohihi1,irhohihi2
      REAL*8 rhohi(
     & irhohilo0:irhohihi0,
     & irhohilo1:irhohihi1,
     & irhohilo2:irhohihi2)
      integer tandir
      REAL*8 dtbydx
      REAL*8 veloc(0:3 -1)
      integer i,j,k, idir
      REAL*8 dense, denlo, denhi, denslope
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
      dense = rho(i,j,k)
      denslope = drho(i,j,k)
      denlo = rholo(i,j,k)
      denhi = rhohi(i,j,k)
      do idir = 0, 3 -1
         veloc(idir) = velcc(i,j,k,idir)
      enddo
      call pointpredadvecttrans(
     & dense, denlo, denhi, denslope, veloc,
     & tandir, dtbydx)
      rholo(i,j,k) = denlo
      rhohi(i,j,k) = denhi
      enddo
      enddo
      enddo
      return
      end
      subroutine POINTPREDADVECTTRANS(
     & dense
     & ,denlo
     & ,denhi
     & ,denslope
     & ,veloc
     & ,tandir
     & ,dtbydx
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 dense
      REAL*8 denlo
      REAL*8 denhi
      REAL*8 denslope
      REAL*8 veloc(0:2)
      integer tandir
      REAL*8 dtbydx
      REAL*8 velpos, velneg
      REAL*8 tol
      denhi = denhi - (0.500d0)*dtbydx*veloc(tandir)*denslope
      denlo = denlo - (0.500d0)*dtbydx*veloc(tandir)*denslope
      ch_flops=ch_flops+8
      return
      end
      subroutine GETDSDT(
     & dsdtplus
     & ,idsdtpluslo0,idsdtpluslo1,idsdtpluslo2
     & ,idsdtplushi0,idsdtplushi1,idsdtplushi2
     & ,ndsdtpluscomp
     & ,dsdtminu
     & ,idsdtminulo0,idsdtminulo1,idsdtminulo2
     & ,idsdtminuhi0,idsdtminuhi1,idsdtminuhi2
     & ,ndsdtminucomp
     & ,slopeprim
     & ,islopeprimlo0,islopeprimlo1,islopeprimlo2
     & ,islopeprimhi0,islopeprimhi1,islopeprimhi2
     & ,nslopeprimcomp
     & ,slopeupwi
     & ,islopeupwilo0,islopeupwilo1,islopeupwilo2
     & ,islopeupwihi0,islopeupwihi1,islopeupwihi2
     & ,nslopeupwicomp
     & ,normalvel
     & ,inormalvello0,inormalvello1,inormalvello2
     & ,inormalvelhi0,inormalvelhi1,inormalvelhi2
     & ,nnormalvelcomp
     & ,dx
     & ,dt
     & ,facedir
     & ,extrapdir
     & ,numslopes
     & ,ientireboxlo0,ientireboxlo1,ientireboxlo2
     & ,ientireboxhi0,ientireboxhi1,ientireboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndsdtpluscomp
      integer idsdtpluslo0,idsdtpluslo1,idsdtpluslo2
      integer idsdtplushi0,idsdtplushi1,idsdtplushi2
      REAL*8 dsdtplus(
     & idsdtpluslo0:idsdtplushi0,
     & idsdtpluslo1:idsdtplushi1,
     & idsdtpluslo2:idsdtplushi2,
     & 0:ndsdtpluscomp-1)
      integer ndsdtminucomp
      integer idsdtminulo0,idsdtminulo1,idsdtminulo2
      integer idsdtminuhi0,idsdtminuhi1,idsdtminuhi2
      REAL*8 dsdtminu(
     & idsdtminulo0:idsdtminuhi0,
     & idsdtminulo1:idsdtminuhi1,
     & idsdtminulo2:idsdtminuhi2,
     & 0:ndsdtminucomp-1)
      integer nslopeprimcomp
      integer islopeprimlo0,islopeprimlo1,islopeprimlo2
      integer islopeprimhi0,islopeprimhi1,islopeprimhi2
      REAL*8 slopeprim(
     & islopeprimlo0:islopeprimhi0,
     & islopeprimlo1:islopeprimhi1,
     & islopeprimlo2:islopeprimhi2,
     & 0:nslopeprimcomp-1)
      integer nslopeupwicomp
      integer islopeupwilo0,islopeupwilo1,islopeupwilo2
      integer islopeupwihi0,islopeupwihi1,islopeupwihi2
      REAL*8 slopeupwi(
     & islopeupwilo0:islopeupwihi0,
     & islopeupwilo1:islopeupwihi1,
     & islopeupwilo2:islopeupwihi2,
     & 0:nslopeupwicomp-1)
      integer nnormalvelcomp
      integer inormalvello0,inormalvello1,inormalvello2
      integer inormalvelhi0,inormalvelhi1,inormalvelhi2
      REAL*8 normalvel(
     & inormalvello0:inormalvelhi0,
     & inormalvello1:inormalvelhi1,
     & inormalvello2:inormalvelhi2,
     & 0:nnormalvelcomp-1)
      REAL*8 dx
      REAL*8 dt
      integer facedir
      integer extrapdir
      integer numslopes
      integer ientireboxlo0,ientireboxlo1,ientireboxlo2
      integer ientireboxhi0,ientireboxhi1,ientireboxhi2
      integer i ,j ,k ,lvar
      REAL*8 velpt, dw, velplus, velminu
      do lvar = 0, numslopes - 1
      do k = ientireboxlo2,ientireboxhi2
      do j = ientireboxlo1,ientireboxhi1
      do i = ientireboxlo0,ientireboxhi0
         velpt = normalvel(i,j,k,extrapdir)
         if(extrapdir .eq. facedir) then
            dw = slopeprim(i,j,k, lvar)
            velplus = max(velpt, (0.0d0))
            velminu = min(velpt, (0.0d0))
            ch_flops=ch_flops+2
         else
            dw = slopeupwi(i,j,k, lvar)
            velplus = velpt
            velminu = velpt
         endif
         dsdtplus(i,j,k,lvar) =
     $ dsdtplus(i,j,k,lvar) - velplus*dw/dx
         dsdtminu(i,j,k,lvar) =
     $ dsdtminu(i,j,k,lvar) - velminu*dw/dx
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(ientireboxhi0- ientireboxlo0+1)*(ientireboxhi1-
     & ientireboxlo1+1)*(ientireboxhi2- ientireboxlo2+1)*6*numslopes
      return
      end
      subroutine UPWINDDIFFS(
     & slopeupwi
     & ,islopeupwilo0,islopeupwilo1,islopeupwilo2
     & ,islopeupwihi0,islopeupwihi1,islopeupwihi2
     & ,nslopeupwicomp
     & ,primstate
     & ,iprimstatelo0,iprimstatelo1,iprimstatelo2
     & ,iprimstatehi0,iprimstatehi1,iprimstatehi2
     & ,nprimstatecomp
     & ,normalvel
     & ,inormalvello0,inormalvello1,inormalvello2
     & ,inormalvelhi0,inormalvelhi1,inormalvelhi2
     & ,nnormalvelcomp
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
      integer nslopeupwicomp
      integer islopeupwilo0,islopeupwilo1,islopeupwilo2
      integer islopeupwihi0,islopeupwihi1,islopeupwihi2
      REAL*8 slopeupwi(
     & islopeupwilo0:islopeupwihi0,
     & islopeupwilo1:islopeupwihi1,
     & islopeupwilo2:islopeupwihi2,
     & 0:nslopeupwicomp-1)
      integer nprimstatecomp
      integer iprimstatelo0,iprimstatelo1,iprimstatelo2
      integer iprimstatehi0,iprimstatehi1,iprimstatehi2
      REAL*8 primstate(
     & iprimstatelo0:iprimstatehi0,
     & iprimstatelo1:iprimstatehi1,
     & iprimstatelo2:iprimstatehi2,
     & 0:nprimstatecomp-1)
      integer nnormalvelcomp
      integer inormalvello0,inormalvello1,inormalvello2
      integer inormalvelhi0,inormalvelhi1,inormalvelhi2
      REAL*8 normalvel(
     & inormalvello0:inormalvelhi0,
     & inormalvello1:inormalvelhi1,
     & inormalvello2:inormalvelhi2,
     & 0:nnormalvelcomp-1)
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
      REAL*8 dwr,dwl,velpt, dw
      REAL*8 tol
      tol = 1.e-12
      ioff = chf_id(0,idir)
      joff = chf_id(1,idir)
      koff = chf_id(2,idir)
      do lvar = 0, numslopes - 1
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
         velpt = normalvel(i,j,k,idir)
         dwr = primstate(i+ioff,j+joff,k+koff,lvar)
     & - primstate(i ,j ,k ,lvar)
         dwl = primstate(i ,j ,k ,lvar)
     & - primstate(i-ioff,j-joff,k-koff,lvar)
         if(velpt .gt. tol) then
            dw = dwl
         else if(velpt .lt. -tol) then
            dw = dwr
         else
            dw = (0.500d0)*(dwl+dwr)
            ch_flops=ch_flops+2
         endif
         slopeupwi(i,j,k,lvar) = dw
      enddo
      enddo
      enddo
         if (haslo .ne. 0) then
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            velpt = normalvel(i,j,k,idir)
            dwl = 0.d0
            dwr = primstate(i+ioff,j+joff,k+koff,lvar)
     & - primstate(i ,j ,k ,lvar)
            dw = dwr
            slopeupwi(i,j,k,lvar) = dw
      enddo
      enddo
      enddo
         endif
         if (hashi .ne. 0) then
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            velpt = normalvel(i,j,k,idir)
            dwl = primstate(i ,j ,k ,lvar)
     & - primstate(i-ioff,j-joff,k-koff,lvar)
            dwr = 0.d0
            dw = dwl
            slopeupwi(i,j,k,lvar) = dw
      enddo
      enddo
      enddo
         endif
      enddo
      ch_flops=ch_flops+((icenterboxhi0- icenterboxlo0+1)*(icenterboxhi1
     &- icenterboxlo1+1)*(icenterboxhi2- icenterboxlo2+1)*2+(iloboxhi0- 
     &iloboxlo0+1)*(iloboxhi1- iloboxlo1+1)*(iloboxhi2- iloboxlo2+1)+(ih
     &iboxhi0- ihiboxlo0+1)*(ihiboxhi1- ihiboxlo1+1)*(ihiboxhi2- ihiboxl
     &o2+1))*numslopes
      return
      end
      subroutine ADVECTIVEF(
     & udelrho
     & ,iudelrholo0,iudelrholo1,iudelrholo2
     & ,iudelrhohi0,iudelrhohi1,iudelrhohi2
     & ,nudelrhocomp
     & ,facerho
     & ,ifacerholo0,ifacerholo1,ifacerholo2
     & ,ifacerhohi0,ifacerhohi1,ifacerhohi2
     & ,nfacerhocomp
     & ,facevel
     & ,ifacevello0,ifacevello1,ifacevello2
     & ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     & ,facedir
     & ,nconserved
     & ,dx
     & ,idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,doingvel
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nudelrhocomp
      integer iudelrholo0,iudelrholo1,iudelrholo2
      integer iudelrhohi0,iudelrhohi1,iudelrhohi2
      REAL*8 udelrho(
     & iudelrholo0:iudelrhohi0,
     & iudelrholo1:iudelrhohi1,
     & iudelrholo2:iudelrhohi2,
     & 0:nudelrhocomp-1)
      integer nfacerhocomp
      integer ifacerholo0,ifacerholo1,ifacerholo2
      integer ifacerhohi0,ifacerhohi1,ifacerhohi2
      REAL*8 facerho(
     & ifacerholo0:ifacerhohi0,
     & ifacerholo1:ifacerhohi1,
     & ifacerholo2:ifacerhohi2,
     & 0:nfacerhocomp-1)
      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1,
     & ifacevello2:ifacevelhi2)
      integer facedir
      integer nconserved
      REAL*8 dx
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer doingvel
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      REAL*8 uave, rhodiff, hival, loval
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      if(doingvel .eq. 1) then
         do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
            uave =
     & (facevel(i+ioff,j+joff,k+koff)
     & +facevel(i ,j ,k ))/(2.0d0)
            rhodiff =
     & (facerho(i+ioff,j+joff,k+koff,iv)
     & -facerho(i ,j ,k ,iv))/dx
            udelrho(i,j,k,iv) = udelrho(i,j,k,iv) + uave*rhodiff
      enddo
      enddo
      enddo
         enddo
         ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*(idcalchi1- idcalclo
     &1+1)*(idcalchi2- idcalclo2+1)*6*nconserved
      else
         do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
            hival = facevel(i+ioff,j+joff,k+koff)*facerho(i+ioff,j+joff,
     &k+koff,iv)
            loval = facevel(i ,j ,k )*facerho(i ,j ,k ,iv)
            udelrho(i,j,k,iv) = udelrho(i,j,k,iv) + (hival-loval)/dx
      enddo
      enddo
      enddo
         enddo
         ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*(idcalchi1- idcalclo
     &1+1)*(idcalchi2- idcalclo2+1)*5*nconserved
      endif
      return
      end
      subroutine EBAVEFACETOCELL(
     & cellvel
     & ,icellvello0,icellvello1,icellvello2
     & ,icellvelhi0,icellvelhi1,icellvelhi2
     & ,ncellvelcomp
     & ,facevel
     & ,ifacevello0,ifacevello1,ifacevello2
     & ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     & ,idir
     & ,idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ncellvelcomp
      integer icellvello0,icellvello1,icellvello2
      integer icellvelhi0,icellvelhi1,icellvelhi2
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & icellvello1:icellvelhi1,
     & icellvello2:icellvelhi2,
     & 0:ncellvelcomp-1)
      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1,
     & ifacevello2:ifacevelhi2)
      integer idir
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer i, j, k
      integer ii,jj,kk
      ii = chf_id(0,idir)
      jj = chf_id(1,idir)
      kk = chf_id(2,idir)
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
      cellvel(i,j,k,idir) = (0.500d0)*(
     & facevel(i ,j ,k ) +
     & facevel(i+ii,j+jj,k+kk))
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*(idcalchi1- idcalclo1+1
     &)*(idcalchi2- idcalclo2+1)*2
      return
      end
      subroutine ADVECTUPDATE(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,primminu
     & ,iprimminulo0,iprimminulo1,iprimminulo2
     & ,iprimminuhi0,iprimminuhi1,iprimminuhi2
     & ,nprimminucomp
     & ,primplus
     & ,iprimpluslo0,iprimpluslo1,iprimpluslo2
     & ,iprimplushi0,iprimplushi1,iprimplushi2
     & ,nprimpluscomp
     & ,primface
     & ,iprimfacelo0,iprimfacelo1,iprimfacelo2
     & ,iprimfacehi0,iprimfacehi1,iprimfacehi2
     & ,nprimfacecomp
     & ,normvel
     & ,inormvello0,inormvello1,inormvello2
     & ,inormvelhi0,inormvelhi1,inormvelhi2
     & ,nnormvelcomp
     & ,facedir
     & ,nprim
     & ,dtbydx
     & ,icellboxlo0,icellboxlo1,icellboxlo2
     & ,icellboxhi0,icellboxhi1,icellboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer nprimminucomp
      integer iprimminulo0,iprimminulo1,iprimminulo2
      integer iprimminuhi0,iprimminuhi1,iprimminuhi2
      REAL*8 primminu(
     & iprimminulo0:iprimminuhi0,
     & iprimminulo1:iprimminuhi1,
     & iprimminulo2:iprimminuhi2,
     & 0:nprimminucomp-1)
      integer nprimpluscomp
      integer iprimpluslo0,iprimpluslo1,iprimpluslo2
      integer iprimplushi0,iprimplushi1,iprimplushi2
      REAL*8 primplus(
     & iprimpluslo0:iprimplushi0,
     & iprimpluslo1:iprimplushi1,
     & iprimpluslo2:iprimplushi2,
     & 0:nprimpluscomp-1)
      integer nprimfacecomp
      integer iprimfacelo0,iprimfacelo1,iprimfacelo2
      integer iprimfacehi0,iprimfacehi1,iprimfacehi2
      REAL*8 primface(
     & iprimfacelo0:iprimfacehi0,
     & iprimfacelo1:iprimfacehi1,
     & iprimfacelo2:iprimfacehi2,
     & 0:nprimfacecomp-1)
      integer nnormvelcomp
      integer inormvello0,inormvello1,inormvello2
      integer inormvelhi0,inormvelhi1,inormvelhi2
      REAL*8 normvel(
     & inormvello0:inormvelhi0,
     & inormvello1:inormvelhi1,
     & inormvello2:inormvelhi2,
     & 0:nnormvelcomp-1)
      integer facedir
      integer nprim
      REAL*8 dtbydx
      integer icellboxlo0,icellboxlo1,icellboxlo2
      integer icellboxhi0,icellboxhi1,icellboxhi2
      REAL*8 unorm
      REAL*8 primdiff
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      do iv = 0,nprim - 1
      do k = icellboxlo2,icellboxhi2
      do j = icellboxlo1,icellboxhi1
      do i = icellboxlo0,icellboxhi0
         unorm = normvel(i,j,k, facedir)
         primdiff =
     $ ( primface(i+ioff,j+joff,k+koff,iv)
     $ - primface(i ,j ,k ,iv))
         primminu(i,j,k,iv) = primminu(i,j,k,iv) -
     & dtbydx * unorm * primdiff
         primplus(i,j,k,iv) = primplus(i,j,k,iv) -
     & dtbydx * unorm * primdiff
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(icellboxhi0- icellboxlo0+1)*(icellboxhi1- icell
     &boxlo1+1)*(icellboxhi2- icellboxlo2+1)*7*nprim
      return
      end
      subroutine ADVECTRIEMANN(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,primgdnv
     & ,iprimgdnvlo0,iprimgdnvlo1,iprimgdnvlo2
     & ,iprimgdnvhi0,iprimgdnvhi1,iprimgdnvhi2
     & ,nprimgdnvcomp
     & ,primleft
     & ,iprimleftlo0,iprimleftlo1,iprimleftlo2
     & ,iprimlefthi0,iprimlefthi1,iprimlefthi2
     & ,nprimleftcomp
     & ,primrigh
     & ,iprimrighlo0,iprimrighlo1,iprimrighlo2
     & ,iprimrighhi0,iprimrighhi1,iprimrighhi2
     & ,nprimrighcomp
     & ,advectvel
     & ,iadvectvello0,iadvectvello1,iadvectvello2
     & ,iadvectvelhi0,iadvectvelhi1,iadvectvelhi2
     & ,facedir
     & ,nprim
     & ,curcomp
     & ,doingvel
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer nprimgdnvcomp
      integer iprimgdnvlo0,iprimgdnvlo1,iprimgdnvlo2
      integer iprimgdnvhi0,iprimgdnvhi1,iprimgdnvhi2
      REAL*8 primgdnv(
     & iprimgdnvlo0:iprimgdnvhi0,
     & iprimgdnvlo1:iprimgdnvhi1,
     & iprimgdnvlo2:iprimgdnvhi2,
     & 0:nprimgdnvcomp-1)
      integer nprimleftcomp
      integer iprimleftlo0,iprimleftlo1,iprimleftlo2
      integer iprimlefthi0,iprimlefthi1,iprimlefthi2
      REAL*8 primleft(
     & iprimleftlo0:iprimlefthi0,
     & iprimleftlo1:iprimlefthi1,
     & iprimleftlo2:iprimlefthi2,
     & 0:nprimleftcomp-1)
      integer nprimrighcomp
      integer iprimrighlo0,iprimrighlo1,iprimrighlo2
      integer iprimrighhi0,iprimrighhi1,iprimrighhi2
      REAL*8 primrigh(
     & iprimrighlo0:iprimrighhi0,
     & iprimrighlo1:iprimrighhi1,
     & iprimrighlo2:iprimrighhi2,
     & 0:nprimrighcomp-1)
      integer iadvectvello0,iadvectvello1,iadvectvello2
      integer iadvectvelhi0,iadvectvelhi1,iadvectvelhi2
      REAL*8 advectvel(
     & iadvectvello0:iadvectvelhi0,
     & iadvectvello1:iadvectvelhi1,
     & iadvectvello2:iadvectvelhi2)
      integer facedir
      integer nprim
      integer curcomp
      integer doingvel
      REAL*8 velhi, vello, velface
      REAL*8 tol
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      tol = 1.e-12
      do iv = 0,nprim - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         velface = advectvel(i,j,k)
         if(velface.gt.tol) then
            primgdnv(i,j,k,iv) = primleft(i-ioff,j-joff,k-koff,iv)
         else if(velface.lt.-tol) then
            primgdnv(i,j,k,iv) = primrigh(i,j,k,iv)
         else
            if( (doingvel .eq. 1) .and. (curcomp.eq.facedir)) then
               primgdnv(i,j,k,iv) = (0.0d0)
            else
               primgdnv(i,j,k,iv) =
     $ (0.500d0)*(
     $ primrigh(i ,j ,k ,iv) +
     $ primleft(i-ioff,j-joff,k-koff,iv))
               ch_flops=ch_flops+2
            endif
         endif
      enddo
      enddo
      enddo
      enddo
      return
      end
