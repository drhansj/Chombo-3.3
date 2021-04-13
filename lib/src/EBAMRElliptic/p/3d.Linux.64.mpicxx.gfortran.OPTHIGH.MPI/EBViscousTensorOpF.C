#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NWOEBVTOPOINTLPH(
     &           lphi
     &           ,iv
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,acofab
     &           ,iacofablo0,iacofablo1,iacofablo2
     &           ,iacofabhi0,iacofabhi1,iacofabhi2
     &           ,eta0fab
     &           ,ieta0fablo0,ieta0fablo1,ieta0fablo2
     &           ,ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
     &           ,eta1fab
     &           ,ieta1fablo0,ieta1fablo1,ieta1fablo2
     &           ,ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
     &           ,eta2fab
     &           ,ieta2fablo0,ieta2fablo1,ieta2fablo2
     &           ,ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
     &           ,lam0fab
     &           ,ilam0fablo0,ilam0fablo1,ilam0fablo2
     &           ,ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
     &           ,lam1fab
     &           ,ilam1fablo0,ilam1fablo1,ilam1fablo2
     &           ,ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
     &           ,lam2fab
     &           ,ilam2fablo0,ilam2fablo1,ilam2fablo2
     &           ,ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      REAL_T lphi(0:2)
      integer iv(0:2)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iacofablo0,iacofablo1,iacofablo2
      integer iacofabhi0,iacofabhi1,iacofabhi2
      REAL_T acofab(
     &           iacofablo0:iacofabhi0,
     &           iacofablo1:iacofabhi1,
     &           iacofablo2:iacofabhi2)
      integer ieta0fablo0,ieta0fablo1,ieta0fablo2
      integer ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
      REAL_T eta0fab(
     &           ieta0fablo0:ieta0fabhi0,
     &           ieta0fablo1:ieta0fabhi1,
     &           ieta0fablo2:ieta0fabhi2)
      integer ieta1fablo0,ieta1fablo1,ieta1fablo2
      integer ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
      REAL_T eta1fab(
     &           ieta1fablo0:ieta1fabhi0,
     &           ieta1fablo1:ieta1fabhi1,
     &           ieta1fablo2:ieta1fabhi2)
      integer ieta2fablo0,ieta2fablo1,ieta2fablo2
      integer ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
      REAL_T eta2fab(
     &           ieta2fablo0:ieta2fabhi0,
     &           ieta2fablo1:ieta2fabhi1,
     &           ieta2fablo2:ieta2fabhi2)
      integer ilam0fablo0,ilam0fablo1,ilam0fablo2
      integer ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
      REAL_T lam0fab(
     &           ilam0fablo0:ilam0fabhi0,
     &           ilam0fablo1:ilam0fabhi1,
     &           ilam0fablo2:ilam0fabhi2)
      integer ilam1fablo0,ilam1fablo1,ilam1fablo2
      integer ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
      REAL_T lam1fab(
     &           ilam1fablo0:ilam1fabhi0,
     &           ilam1fablo1:ilam1fabhi1,
     &           ilam1fablo2:ilam1fabhi2)
      integer ilam2fablo0,ilam2fablo1,ilam2fablo2
      integer ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
      REAL_T lam2fab(
     &           ilam2fablo0:ilam2fabhi0,
     &           ilam2fablo1:ilam2fabhi1,
     &           ilam2fablo2:ilam2fabhi2)
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      real_t divf(0:CH_SPACEDIM-1)
      real_t aphi(0:CH_SPACEDIM-1)
      real_t  etaL(0:CH_SPACEDIM-1)
      real_t  etaH(0:CH_SPACEDIM-1)
      real_t  lamL(0:CH_SPACEDIM-1)
      real_t  lamH(0:CH_SPACEDIM-1)
      real_t  fluxL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  fluxH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divuL(0:CH_SPACEDIM-1)
      real_t  divuH(0:CH_SPACEDIM-1)
      integer i,j,k, facedir, derivdir, veldir
      integer iif,jjf,kkf
      integer iid,jjd,kkd
      
      i = iv(0)
      j = iv(1)
      k = iv(2)
      
      etaL(0) =eta0fab(i   ,j   ,k   )
      etaL(1) =eta1fab(i   ,j   ,k   )
      etaL(2) =eta2fab(i   ,j   ,k   )
                                  
      etaH(0) =eta0fab(i+1 ,j   ,k   )
      etaH(1) =eta1fab(i   ,j+1 ,k   )
      etaH(2) =eta2fab(i   ,j   ,k+1 )
      
      lamL(0) =lam0fab(i   ,j   ,k   )
      lamL(1) =lam1fab(i   ,j   ,k   )
      lamL(2) =lam2fab(i   ,j   ,k   )
                                  
      lamH(0) =lam0fab(i+1 ,j   ,k   )
      lamH(1) =lam1fab(i   ,j+1 ,k   )
      lamH(2) =lam2fab(i   ,j   ,k+1 )
      do facedir = 0, CH_SPACEDIM-1
         
         iif = chf_id(facedir, 0)
         jjf = chf_id(facedir, 1)
         kkf = chf_id(facedir, 2)
         do derivdir = 0, CH_SPACEDIM-1
            
            iid = chf_id(derivdir, 0)
            jjd = chf_id(derivdir, 1)
            kkd = chf_id(derivdir, 2)
            do veldir = 0, CH_SPACEDIM-1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phi(i+iid,j+jjd,k+kkd,veldir) - phi(i    ,j    ,k    ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phi(i    ,j    ,k    ,veldir) - phi(i-iid,j-jjd,k-kkd,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid+iif,j+jjd+jjf,k+kkd+kkf,veldir) - phi(i-iid+iif,j-jjd+jjf,k-kkd+kkf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,k+kkd    ,veldir) - phi(i-iid    ,j-jjd    ,k-kkd    ,veldir)  )
                  gphiL(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid-iif,j+jjd-jjf,k+kkd-kkf,veldir) - phi(i-iid-iif,j-jjd-jjf,k-kkd-kkf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,k+kkd    ,veldir) - phi(i-iid    ,j-jjd    ,k-kkd    ,veldir)  )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         divuL(facedir) = zero
         divuH(facedir) = zero
         do veldir = 0, CH_SPACEDIM-1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, facedir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, facedir)
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         do veldir = 0, CH_SPACEDIM-1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldir, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldir, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(facedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(facedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, CH_SPACEDIM-1
         divf(veldir) = zero
         do facedir = 0, CH_SPACEDIM-1
            divf(veldir) = divf(veldir)
     $           + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      
      aphi(0) = acofab(i,j,k)*phi(i,j,k, 0)
      aphi(1) = acofab(i,j,k)*phi(i,j,k, 1)
      aphi(2) = acofab(i,j,k)*phi(i,j,k, 2)
      
      lphi(0) =  alpha*aphi(0) + beta*divf(0)
      lphi(1) =  alpha*aphi(1) + beta*divf(1)
      lphi(2) =  alpha*aphi(2) + beta*divf(2)
      ch_flops=ch_flops + 4 + 6*CH_SPACEDIM + 12*(CH_SPACEDIM*CH_SPACEDIM)
      return
      end
      subroutine GSRBNWOEBVTOP(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1,ilambdalo2
     &           ,ilambdahi0,ilambdahi1,ilambdahi2
     &           ,nlambdacomp
     &           ,acofab
     &           ,iacofablo0,iacofablo1,iacofablo2
     &           ,iacofabhi0,iacofabhi1,iacofabhi2
     &           ,eta0fab
     &           ,ieta0fablo0,ieta0fablo1,ieta0fablo2
     &           ,ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
     &           ,eta1fab
     &           ,ieta1fablo0,ieta1fablo1,ieta1fablo2
     &           ,ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
     &           ,eta2fab
     &           ,ieta2fablo0,ieta2fablo1,ieta2fablo2
     &           ,ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
     &           ,lam0fab
     &           ,ilam0fablo0,ilam0fablo1,ilam0fablo2
     &           ,ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
     &           ,lam1fab
     &           ,ilam1fablo0,ilam1fablo1,ilam1fablo2
     &           ,ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
     &           ,lam2fab
     &           ,ilam2fablo0,ilam2fablo1,ilam2fablo2
     &           ,ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
     &           ,dx
     &           ,alpha
     &           ,beta
     &           ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     &           ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           ilambdalo2:ilambdahi2,
     &           0:nlambdacomp-1)
      integer iacofablo0,iacofablo1,iacofablo2
      integer iacofabhi0,iacofabhi1,iacofabhi2
      REAL_T acofab(
     &           iacofablo0:iacofabhi0,
     &           iacofablo1:iacofabhi1,
     &           iacofablo2:iacofabhi2)
      integer ieta0fablo0,ieta0fablo1,ieta0fablo2
      integer ieta0fabhi0,ieta0fabhi1,ieta0fabhi2
      REAL_T eta0fab(
     &           ieta0fablo0:ieta0fabhi0,
     &           ieta0fablo1:ieta0fabhi1,
     &           ieta0fablo2:ieta0fabhi2)
      integer ieta1fablo0,ieta1fablo1,ieta1fablo2
      integer ieta1fabhi0,ieta1fabhi1,ieta1fabhi2
      REAL_T eta1fab(
     &           ieta1fablo0:ieta1fabhi0,
     &           ieta1fablo1:ieta1fabhi1,
     &           ieta1fablo2:ieta1fabhi2)
      integer ieta2fablo0,ieta2fablo1,ieta2fablo2
      integer ieta2fabhi0,ieta2fabhi1,ieta2fabhi2
      REAL_T eta2fab(
     &           ieta2fablo0:ieta2fabhi0,
     &           ieta2fablo1:ieta2fabhi1,
     &           ieta2fablo2:ieta2fabhi2)
      integer ilam0fablo0,ilam0fablo1,ilam0fablo2
      integer ilam0fabhi0,ilam0fabhi1,ilam0fabhi2
      REAL_T lam0fab(
     &           ilam0fablo0:ilam0fabhi0,
     &           ilam0fablo1:ilam0fabhi1,
     &           ilam0fablo2:ilam0fabhi2)
      integer ilam1fablo0,ilam1fablo1,ilam1fablo2
      integer ilam1fabhi0,ilam1fabhi1,ilam1fabhi2
      REAL_T lam1fab(
     &           ilam1fablo0:ilam1fabhi0,
     &           ilam1fablo1:ilam1fabhi1,
     &           ilam1fablo2:ilam1fabhi2)
      integer ilam2fablo0,ilam2fablo1,ilam2fablo2
      integer ilam2fabhi0,ilam2fabhi1,ilam2fabhi2
      REAL_T lam2fab(
     &           ilam2fablo0:ilam2fabhi0,
     &           ilam2fablo1:ilam2fabhi1,
     &           ilam2fablo2:ilam2fabhi2)
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      REAL_T lphi(0:CH_SPACEDIM-1)
      real_t divf(0:CH_SPACEDIM-1)
      real_t aphi(0:CH_SPACEDIM-1)
      real_t  etaL(0:CH_SPACEDIM-1)
      real_t  etaH(0:CH_SPACEDIM-1)
      real_t  lamL(0:CH_SPACEDIM-1)
      real_t  lamH(0:CH_SPACEDIM-1)
      real_t  fluxL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  fluxH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divuL(0:CH_SPACEDIM-1)
      real_t  divuH(0:CH_SPACEDIM-1)
      integer i,j,k, facedir, derivdir, veldir
      integer iif,jjf,kkf
      integer iid,jjd,kkd, ncolors
      
      do k = icoloredboxlo2,icoloredboxhi2,2
      do j = icoloredboxlo1,icoloredboxhi1,2
      do i = icoloredboxlo0,icoloredboxhi0,2

      
      etaL(0) =eta0fab(i   ,j   ,k   )
      etaL(1) =eta1fab(i   ,j   ,k   )
      etaL(2) =eta2fab(i   ,j   ,k   )
                                  
      etaH(0) =eta0fab(i+1 ,j   ,k   )
      etaH(1) =eta1fab(i   ,j+1 ,k   )
      etaH(2) =eta2fab(i   ,j   ,k+1 )
      
      lamL(0) =lam0fab(i   ,j   ,k   )
      lamL(1) =lam1fab(i   ,j   ,k   )
      lamL(2) =lam2fab(i   ,j   ,k   )
                                  
      lamH(0) =lam0fab(i+1 ,j   ,k   )
      lamH(1) =lam1fab(i   ,j+1 ,k   )
      lamH(2) =lam2fab(i   ,j   ,k+1 )
      do facedir = 0, CH_SPACEDIM-1
         
         iif = chf_id(facedir, 0)
         jjf = chf_id(facedir, 1)
         kkf = chf_id(facedir, 2)
         do derivdir = 0, CH_SPACEDIM-1
            
            iid = chf_id(derivdir, 0)
            jjd = chf_id(derivdir, 1)
            kkd = chf_id(derivdir, 2)
            do veldir = 0, CH_SPACEDIM-1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phi(i+iid,j+jjd,k+kkd,veldir) - phi(i    ,j    ,k    ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phi(i    ,j    ,k    ,veldir) - phi(i-iid,j-jjd,k-kkd,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid+iif,j+jjd+jjf,k+kkd+kkf,veldir) - phi(i-iid+iif,j-jjd+jjf,k-kkd+kkf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,k+kkd    ,veldir) - phi(i-iid    ,j-jjd    ,k-kkd    ,veldir)  )
                  gphiL(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phi(i+iid-iif,j+jjd-jjf,k+kkd-kkf,veldir) - phi(i-iid-iif,j-jjd-jjf,k-kkd-kkf,veldir) + 
     $                 phi(i+iid    ,j+jjd    ,k+kkd    ,veldir) - phi(i-iid    ,j-jjd    ,k-kkd    ,veldir)  )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         divuL(facedir) = zero
         divuH(facedir) = zero
         do veldir = 0, CH_SPACEDIM-1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, facedir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, facedir)
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         do veldir = 0, CH_SPACEDIM-1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldir, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldir, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(facedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(facedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, CH_SPACEDIM-1
         divf(veldir) = zero
         do facedir = 0, CH_SPACEDIM-1
            divf(veldir) = divf(veldir)
     $           + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      
      aphi(0) = acofab(i,j,k)*phi(i,j,k, 0)
      aphi(1) = acofab(i,j,k)*phi(i,j,k, 1)
      aphi(2) = acofab(i,j,k)*phi(i,j,k, 2)
      
      lphi(0) =  alpha*aphi(0) + beta*divf(0)
      lphi(1) =  alpha*aphi(1) + beta*divf(1)
      lphi(2) =  alpha*aphi(2) + beta*divf(2)
      
      phi(i,j,k, 0) =  phi(i,j,k,0) + lambda(i,j,k,0)*(rhs(i,j,k,0) - lphi(0))
      phi(i,j,k, 1) =  phi(i,j,k,1) + lambda(i,j,k,1)*(rhs(i,j,k,1) - lphi(1))
      phi(i,j,k, 2) =  phi(i,j,k,2) + lambda(i,j,k,2)*(rhs(i,j,k,2) - lphi(2))
      
      enddo
      enddo
      enddo
      ncolors = 2 *2 *2
      ch_flops=ch_flops+(icoloredboxhi0- icoloredboxlo0+1)*(icoloredboxhi1- icoloredboxlo1+1)*(icoloredboxhi2- icoloredboxlo2+1)*(4 + 9*CH_SPACEDIM + 12*(CH_SPACEDIM*CH_SPACEDIM))/ncolors
      return
      end
      subroutine CELLGRADEBVTOP(
     &           grad
     &           ,igradlo0,igradlo1,igradlo2
     &           ,igradhi0,igradhi1,igradhi2
     &           ,vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,dx
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           ,haslo
     &           ,hashi
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           igradlo1:igradhi1,
     &           igradlo2:igradhi2)
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2)
      REAL_T dx
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer haslo
      integer hashi
      integer divdir
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      kk = chf_id(divdir, 2)
      
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

      grad(i,j,k) =
     $     (    vel(i+ii,j+jj,k+kk)
     $     -    vel(i-ii,j-jj,k-kk) )/(two*dx)
      
      enddo
      enddo
      enddo
      if(haslo.eq.1) then
         
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         grad(i,j,k) =
     $        (    vel(i+ii,j+jj,k+kk)
     $        -    vel(i   ,j   ,k   ) )/(dx)
         
      enddo
      enddo
      enddo
      endif
      if(hashi.eq.1) then
         
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         grad(i,j,k) =
     $        (    vel(i   ,j   ,k   )
     $        -    vel(i-ii,j-jj,k-kk) )/(dx)
         
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine INCRAPPLYEBVTOP(
     &           lhs
     &           ,ilhslo0,ilhslo1,ilhslo2
     &           ,ilhshi0,ilhshi1,ilhshi2
     &           ,interiorflux
     &           ,iinteriorfluxlo0,iinteriorfluxlo1,iinteriorfluxlo2
     &           ,iinteriorfluxhi0,iinteriorfluxhi1,iinteriorfluxhi2
     &           ,domainfluxlo
     &           ,idomainfluxlolo0,idomainfluxlolo1,idomainfluxlolo2
     &           ,idomainfluxlohi0,idomainfluxlohi1,idomainfluxlohi2
     &           ,domainfluxhi
     &           ,idomainfluxhilo0,idomainfluxhilo1,idomainfluxhilo2
     &           ,idomainfluxhihi0,idomainfluxhihi1,idomainfluxhihi2
     &           ,beta
     &           ,dx
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           ,haslo
     &           ,hashi
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL_T lhs(
     &           ilhslo0:ilhshi0,
     &           ilhslo1:ilhshi1,
     &           ilhslo2:ilhshi2)
      integer iinteriorfluxlo0,iinteriorfluxlo1,iinteriorfluxlo2
      integer iinteriorfluxhi0,iinteriorfluxhi1,iinteriorfluxhi2
      REAL_T interiorflux(
     &           iinteriorfluxlo0:iinteriorfluxhi0,
     &           iinteriorfluxlo1:iinteriorfluxhi1,
     &           iinteriorfluxlo2:iinteriorfluxhi2)
      integer idomainfluxlolo0,idomainfluxlolo1,idomainfluxlolo2
      integer idomainfluxlohi0,idomainfluxlohi1,idomainfluxlohi2
      REAL_T domainfluxlo(
     &           idomainfluxlolo0:idomainfluxlohi0,
     &           idomainfluxlolo1:idomainfluxlohi1,
     &           idomainfluxlolo2:idomainfluxlohi2)
      integer idomainfluxhilo0,idomainfluxhilo1,idomainfluxhilo2
      integer idomainfluxhihi0,idomainfluxhihi1,idomainfluxhihi2
      REAL_T domainfluxhi(
     &           idomainfluxhilo0:idomainfluxhihi0,
     &           idomainfluxhilo1:idomainfluxhihi1,
     &           idomainfluxhilo2:idomainfluxhihi2)
      REAL_T beta
      REAL_T dx
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
     $     +beta*
     $     (interiorflux(i+ii,j+jj,k+kk)
     $     -interiorflux(i   ,j   ,k   ))/dx
      
      enddo
      enddo
      enddo
      if(haslo .eq. 1) then
         
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

         lhs(i,j,k) = lhs(i,j,k)
     $        + beta*
     $        (interiorflux(i+ii,j+jj,k+kk)
     $        -domainfluxlo(i   ,j   ,k   ))/dx
         
      enddo
      enddo
      enddo
      endif
      if(hashi .eq. 1) then
         
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

         lhs(i,j,k) = lhs(i,j,k)
     $        + beta*
     $        (domainfluxhi(i+ii,j+jj,k+kk)
     $        -interiorflux(i   ,j   ,k   ))/dx
         
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine NORMALGRADVISCDIRCH(
     &           gradvelface
     &           ,igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
     &           ,igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
     &           ,velcomp
     &           ,ivelcomplo0,ivelcomplo1,ivelcomplo2
     &           ,ivelcomphi0,ivelcomphi1,ivelcomphi2
     &           ,velvalu
     &           ,ivelvalulo0,ivelvalulo1,ivelvalulo2
     &           ,ivelvaluhi0,ivelvaluhi1,ivelvaluhi2
     &           ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     &           ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     &           ,dx
     &           ,iside
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
      integer igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0,
     &           igradvelfacelo1:igradvelfacehi1,
     &           igradvelfacelo2:igradvelfacehi2)
      integer ivelcomplo0,ivelcomplo1,ivelcomplo2
      integer ivelcomphi0,ivelcomphi1,ivelcomphi2
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0,
     &           ivelcomplo1:ivelcomphi1,
     &           ivelcomplo2:ivelcomphi2)
      integer ivelvalulo0,ivelvalulo1,ivelvalulo2
      integer ivelvaluhi0,ivelvaluhi1,ivelvaluhi2
      REAL_T velvalu(
     &           ivelvalulo0:ivelvaluhi0,
     &           ivelvalulo1:ivelvaluhi1,
     &           ivelvalulo2:ivelvaluhi2)
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      REAL_T dx
      integer iside
      integer divdir
      integer i,id,j,jd,k,kd
      real_t val0, val1, val2, denom
      
      id = chf_id(divdir, 0)
      jd = chf_id(divdir, 1)
      kd = chf_id(divdir, 2)
      if(iside.eq.-1) then
         
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

         val0 = velvalu(i   ,j   ,k   )
         val1 = velcomp(i   ,j   ,k   )
         val2 = velcomp(i+id,j+jd,k+kd)
         denom = half*dx
         gradvelface(i,j,k) = (val1-val0)/denom
         
      enddo
      enddo
      enddo
      else
         
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

         val0 = velvalu(i     ,j     ,k     )
         val1 = velcomp(i-  id,j-  jd,k-  kd)
         val2 = velcomp(i-2*id,j-2*jd,k-2*kd)
         denom = -half*dx
         gradvelface(i,j,k) = (val1-val0)/denom
         
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine SLIPWALLGRAD(
     &           gradvelface
     &           ,igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
     &           ,igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
     &           ,velcomp
     &           ,ivelcomplo0,ivelcomplo1,ivelcomplo2
     &           ,ivelcomphi0,ivelcomphi1,ivelcomphi2
     &           ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     &           ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     &           ,dx
     &           ,veldir
     &           ,divdir
     &           ,iside
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           ,haslo
     &           ,hashi
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
      integer igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0,
     &           igradvelfacelo1:igradvelfacehi1,
     &           igradvelfacelo2:igradvelfacehi2)
      integer ivelcomplo0,ivelcomplo1,ivelcomplo2
      integer ivelcomphi0,ivelcomphi1,ivelcomphi2
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0,
     &           ivelcomplo1:ivelcomphi1,
     &           ivelcomplo2:ivelcomphi2)
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      REAL_T dx
      integer veldir
      integer divdir
      integer iside
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer haslo
      integer hashi
      integer facedir
      integer i,id,iv,j,jd,jv,k,kd,kv
      real_t  val1, val2
      
      id = chf_id(facedir, 0)
      jd = chf_id(facedir, 1)
      kd = chf_id(facedir, 2)
      
      iv = chf_id(divdir, 0)
      jv = chf_id(divdir, 1)
      kv = chf_id(divdir, 2)
      if(facedir.eq.veldir) then
         if(iside.eq.-1) then
            
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

            val1 = velcomp(i   ,j   ,k   )
            gradvelface(i,j,k) = two*val1/dx
            
      enddo
      enddo
      enddo
         else
            
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

            val1 = velcomp(i-  id,j-  jd,k-  kd)
            gradvelface(i,j,k) = -two*val1/dx
            
      enddo
      enddo
      enddo
         endif
      else
         if(iside.eq.-1) then
            
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

            val1 = velcomp(i+iv,j+jv,k+kv)
            val2 = velcomp(i-iv,j-jv,k-kv)
            gradvelface(i,j,k) = (val1-val2)/(two*dx)
            
      enddo
      enddo
      enddo
            if(haslo.eq.1) then
               
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

               val1 = velcomp(i+iv,j+jv,k+kv)
               val2 = velcomp(i   ,j   ,k   )
               gradvelface(i,j,k) = (val1-val2)/(dx)
               
      enddo
      enddo
      enddo
            endif
            if(hashi.eq.1) then
               
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

               val1 = velcomp(i   ,j   ,k   )
               val2 = velcomp(i-iv,j-jv,k-kv)
               gradvelface(i,j,k) = (val1-val2)/(dx)
               
      enddo
      enddo
      enddo
            endif
         else
         
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

            val1 = velcomp(i-id+iv,j-jd+jv,k-kd+kv)
            val2 = velcomp(i-id-iv,j-jd-jv,k-kd-kv)
            gradvelface(i,j,k) = (val1-val2)/(two*dx)
         
      enddo
      enddo
      enddo
         if(haslo.eq.1) then
            
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

            val1 = velcomp(i-id+iv,j-jd+jv,k-kd+kv)
            val2 = velcomp(i-id   ,j-jd   ,k-kd   )
            gradvelface(i,j,k) = (val1-val2)/(dx)
            
      enddo
      enddo
      enddo
         endif
         if(hashi.eq.1) then
            
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

            val1 = velcomp(i-id   ,j-jd   ,k-kd   )
            val2 = velcomp(i-id-iv,j-jd-jv,k-kd-kv)
            gradvelface(i,j,k) = (val1-val2)/(dx)
            
      enddo
      enddo
      enddo
         endif
         endif
      endif
      return
      end
      subroutine VELDOTSIGMA(
     &           veldotsig
     &           ,iveldotsiglo0,iveldotsiglo1,iveldotsiglo2
     &           ,iveldotsighi0,iveldotsighi1,iveldotsighi2
     &           ,nveldotsigcomp
     &           ,vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,nvelcomp
     &           ,sig
     &           ,isiglo0,isiglo1,isiglo2
     &           ,isighi0,isighi1,isighi2
     &           ,nsigcomp
     &           ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     &           ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nveldotsigcomp
      integer iveldotsiglo0,iveldotsiglo1,iveldotsiglo2
      integer iveldotsighi0,iveldotsighi1,iveldotsighi2
      REAL_T veldotsig(
     &           iveldotsiglo0:iveldotsighi0,
     &           iveldotsiglo1:iveldotsighi1,
     &           iveldotsiglo2:iveldotsighi2,
     &           0:nveldotsigcomp-1)
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2,
     &           0:nvelcomp-1)
      integer nsigcomp
      integer isiglo0,isiglo1,isiglo2
      integer isighi0,isighi1,isighi2
      REAL_T sig(
     &           isiglo0:isighi0,
     &           isiglo1:isighi1,
     &           isiglo2:isighi2,
     &           0:nsigcomp-1)
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k
      integer icomp
      real_t sum
      
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

      sum = zero
      do icomp = 0, CH_SPACEDIM-1
         sum = sum + vel(i,j,k, icomp)*sig(i,j,k, icomp)
      enddo
      veldotsig(i,j,k, 0) = sum
      
      enddo
      enddo
      enddo
      return
      end
      subroutine INCRPOINTDOTPROD(
     &           sigmadotgradu
     &           ,isigmadotgradulo0,isigmadotgradulo1,isigmadotgradulo2
     &           ,isigmadotgraduhi0,isigmadotgraduhi1,isigmadotgraduhi2
     &           ,gradu
     &           ,igradulo0,igradulo1,igradulo2
     &           ,igraduhi0,igraduhi1,igraduhi2
     &           ,sigma
     &           ,isigmalo0,isigmalo1,isigmalo2
     &           ,isigmahi0,isigmahi1,isigmahi2
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmadotgradulo0,isigmadotgradulo1,isigmadotgradulo2
      integer isigmadotgraduhi0,isigmadotgraduhi1,isigmadotgraduhi2
      REAL_T sigmadotgradu(
     &           isigmadotgradulo0:isigmadotgraduhi0,
     &           isigmadotgradulo1:isigmadotgraduhi1,
     &           isigmadotgradulo2:isigmadotgraduhi2)
      integer igradulo0,igradulo1,igradulo2
      integer igraduhi0,igraduhi1,igraduhi2
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1,
     &           igradulo2:igraduhi2)
      integer isigmalo0,isigmalo1,isigmalo2
      integer isigmahi0,isigmahi1,isigmahi2
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1,
     &           isigmalo2:isigmahi2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      real_t incr
      
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

      incr =  gradu(i,j,k)*sigma(i,j,k)
      sigmadotgradu(i,j,k) = sigmadotgradu(i,j,k) + incr
      
      enddo
      enddo
      enddo
      return
      end
      subroutine INCRCCSIGMALAMBDA(
     &           sigma
     &           ,isigmalo0,isigmalo1,isigmalo2
     &           ,isigmahi0,isigmahi1,isigmahi2
     &           ,gradu
     &           ,igradulo0,igradulo1,igradulo2
     &           ,igraduhi0,igraduhi1,igraduhi2
     &           ,ngraducomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1,ilambdalo2
     &           ,ilambdahi0,ilambdahi1,ilambdahi2
     &           ,igridlo0,igridlo1,igridlo2
     &           ,igridhi0,igridhi1,igridhi2
     &           ,diagcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmalo0,isigmalo1,isigmalo2
      integer isigmahi0,isigmahi1,isigmahi2
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1,
     &           isigmalo2:isigmahi2)
      integer ngraducomp
      integer igradulo0,igradulo1,igradulo2
      integer igraduhi0,igraduhi1,igraduhi2
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1,
     &           igradulo2:igraduhi2,
     &           0:ngraducomp-1)
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           ilambdalo2:ilambdahi2)
      integer igridlo0,igridlo1,igridlo2
      integer igridhi0,igridhi1,igridhi2
      integer diagcomp
      integer i,j,k
      
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

        sigma(i,j,k) = sigma(i,j,k)
     $     + lambda(i,j,k)*gradu(i,j,k, diagcomp)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine INCRCCSIGMAETA(
     &           sigma
     &           ,isigmalo0,isigmalo1,isigmalo2
     &           ,isigmahi0,isigmahi1,isigmahi2
     &           ,gradu
     &           ,igradulo0,igradulo1,igradulo2
     &           ,igraduhi0,igraduhi1,igraduhi2
     &           ,ngraducomp
     &           ,eta
     &           ,ietalo0,ietalo1,ietalo2
     &           ,ietahi0,ietahi1,ietahi2
     &           ,igridlo0,igridlo1,igridlo2
     &           ,igridhi0,igridhi1,igridhi2
     &           ,gradcomp
     &           ,trancomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isigmalo0,isigmalo1,isigmalo2
      integer isigmahi0,isigmahi1,isigmahi2
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1,
     &           isigmalo2:isigmahi2)
      integer ngraducomp
      integer igradulo0,igradulo1,igradulo2
      integer igraduhi0,igraduhi1,igraduhi2
      REAL_T gradu(
     &           igradulo0:igraduhi0,
     &           igradulo1:igraduhi1,
     &           igradulo2:igraduhi2,
     &           0:ngraducomp-1)
      integer ietalo0,ietalo1,ietalo2
      integer ietahi0,ietahi1,ietahi2
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1,
     &           ietalo2:ietahi2)
      integer igridlo0,igridlo1,igridlo2
      integer igridhi0,igridhi1,igridhi2
      integer gradcomp
      integer trancomp
      integer i,j,k
      
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

      sigma(i,j,k) = sigma(i,j,k)
     $     +  eta(i,j,k)*
     $     (gradu(i,j,k, gradcomp)
     $     +gradu(i,j,k, trancomp))
      
      enddo
      enddo
      enddo
      return
      end
      subroutine INCRFACETOCELLAVERAGE(
     &           cellcoef
     &           ,icellcoeflo0,icellcoeflo1,icellcoeflo2
     &           ,icellcoefhi0,icellcoefhi1,icellcoefhi2
     &           ,facecoef
     &           ,ifacecoeflo0,ifacecoeflo1,ifacecoeflo2
     &           ,ifacecoefhi0,ifacecoefhi1,ifacecoefhi2
     &           ,facedir
     &           ,igridlo0,igridlo1,igridlo2
     &           ,igridhi0,igridhi1,igridhi2
     &           ,nfacesper
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icellcoeflo0,icellcoeflo1,icellcoeflo2
      integer icellcoefhi0,icellcoefhi1,icellcoefhi2
      REAL_T cellcoef(
     &           icellcoeflo0:icellcoefhi0,
     &           icellcoeflo1:icellcoefhi1,
     &           icellcoeflo2:icellcoefhi2)
      integer ifacecoeflo0,ifacecoeflo1,ifacecoeflo2
      integer ifacecoefhi0,ifacecoefhi1,ifacecoefhi2
      REAL_T facecoef(
     &           ifacecoeflo0:ifacecoefhi0,
     &           ifacecoeflo1:ifacecoefhi1,
     &           ifacecoeflo2:ifacecoefhi2)
      integer facedir
      integer igridlo0,igridlo1,igridlo2
      integer igridhi0,igridhi1,igridhi2
      REAL_T nfacesper
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

        cellcoef(i,j,k) = cellcoef(i,j,k)   +
     $     (facecoef(i+ii,j+jj,k+kk)
     $     +facecoef(i   ,j   ,k   ))/nfacesper
      
      enddo
      enddo
      enddo
      return
      end
