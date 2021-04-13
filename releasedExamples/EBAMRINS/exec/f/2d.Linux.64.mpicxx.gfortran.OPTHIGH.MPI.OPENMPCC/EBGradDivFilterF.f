      subroutine EBGDFFACEDIVINCR(
     & divvel
     & ,idivvello0,idivvello1
     & ,idivvelhi0,idivvelhi1
     & ,vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,nvelcomp
     & ,gradvel
     & ,igradvello0,igradvello1
     & ,igradvelhi0,igradvelhi1
     & ,ngradvelcomp
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
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
      integer idivvello0,idivvello1
      integer idivvelhi0,idivvelhi1
      REAL*8 divvel(
     & idivvello0:idivvelhi0,
     & idivvello1:idivvelhi1)
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & 0:nvelcomp-1)
      integer ngradvelcomp
      integer igradvello0,igradvello1
      integer igradvelhi0,igradvelhi1
      REAL*8 gradvel(
     & igradvello0:igradvelhi0,
     & igradvello1:igradvelhi1,
     & 0:ngradvelcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      integer facedir
      integer divdir
      integer ii,i,jj,j, gradcomp, sdim
      sdim = 2
          ii = chf_id(facedir, 0)
          jj = chf_id(facedir, 1)
          if (facedir .eq. divdir) then
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
             divvel(i,j) = divvel(i,j) +
     $ ( vel(i ,j ,facedir)
     $ - vel(i-ii,j-jj,facedir) )/dx
      enddo
      enddo
          else
             gradcomp = divdir + sdim*divdir
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
             divvel(i,j) = divvel(i,j) +
     $ ( gradvel(i ,j , gradcomp)
     $ + gradvel(i-ii,j-jj, gradcomp) )/(2.0d0)
      enddo
      enddo
          endif
      return
      end
      subroutine EBGDFGRADVEL(
     & grad
     & ,igradlo0,igradlo1
     & ,igradhi0,igradhi1
     & ,vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,iloboxlo0,iloboxlo1
     & ,iloboxhi0,iloboxhi1
     & ,haslo
     & ,ihiboxlo0,ihiboxlo1
     & ,ihiboxhi0,ihiboxhi1
     & ,hashi
     & ,icenterboxlo0,icenterboxlo1
     & ,icenterboxhi0,icenterboxhi1
     & ,dx
     & ,divdir
     & ,loworderoneside
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradlo0,igradlo1
      integer igradhi0,igradhi1
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1)
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1)
      integer iloboxlo0,iloboxlo1
      integer iloboxhi0,iloboxhi1
      integer haslo
      integer ihiboxlo0,ihiboxlo1
      integer ihiboxhi0,ihiboxhi1
      integer hashi
      integer icenterboxlo0,icenterboxlo1
      integer icenterboxhi0,icenterboxhi1
      REAL*8 dx
      integer divdir
      integer loworderoneside
      integer ii,i,jj,j
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0
      grad(i,j) =
     $ ( vel(i+ii,j+jj)
     $ - vel(i-ii,j-jj) )/((2.0d0)*dx)
      enddo
      enddo
      if(haslo.eq.1) then
         if(loworderoneside.eq.1) then
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            grad(i,j) =
     $ ( vel(i+ii,j+jj)
     $ - vel(i ,j ))/(dx)
      enddo
      enddo
         else
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0
            grad(i,j) =
     $ ((4.0d0)*(vel(i+ ii,j+ jj) - vel(i,j))
     $ - (vel(i+2*ii,j+2*jj) - vel(i,j)))/((2.0d0)*dx)
      enddo
      enddo
         endif
      endif
      if(hashi.eq.1) then
         if(loworderoneside.eq.1) then
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            grad(i,j) =
     $ ( vel(i-ii,j-jj)
     $ - vel(i ,j ))/(-dx)
      enddo
      enddo
         else
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0
            grad(i,j) =
     $ ((4.0d0)*(vel(i- ii,j- jj) - vel(i,j))
     $ - (vel(i-2*ii,j-2*jj) - vel(i,j)))/(-(2.0d0)*dx)
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine EBGDFCELLGRAD(
     & graddiv
     & ,igraddivlo0,igraddivlo1
     & ,igraddivhi0,igraddivhi1
     & ,div
     & ,idivlo0,idivlo1
     & ,idivhi0,idivhi1
     & ,lambda
     & ,ilambdalo0,ilambdalo1
     & ,ilambdahi0,ilambdahi1
     & ,iregionlo0,iregionlo1
     & ,iregionhi0,iregionhi1
     & ,dx
     & ,facedir
     & ,imultbylambda
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igraddivlo0,igraddivlo1
      integer igraddivhi0,igraddivhi1
      REAL*8 graddiv(
     & igraddivlo0:igraddivhi0,
     & igraddivlo1:igraddivhi1)
      integer idivlo0,idivlo1
      integer idivhi0,idivhi1
      REAL*8 div(
     & idivlo0:idivhi0,
     & idivlo1:idivhi1)
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      integer facedir
      integer imultbylambda
      integer ii,i,jj,j
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
      graddiv(i,j) =
     $ ( div(i+ii,j+jj)
     $ - div(i ,j ) )/dx
      if(imultbylambda .eq. 1) then
         graddiv(i,j) = lambda(i,j)*graddiv(i,j)
      endif
      enddo
      enddo
      return
      end
