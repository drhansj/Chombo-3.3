#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine VOLWGTSUM(
     &           src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,volfrac
     &           ,ivolfraclo0,ivolfraclo1,ivolfraclo2
     &           ,ivolfrachi0,ivolfrachi1,ivolfrachi2
     &           ,nvolfraccomp
     &           ,norm
     &           ,volume
     &           ,comp
     &           ,pval
     &           ,idoreg
     &           ,idoirr
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      integer nvolfraccomp
      integer ivolfraclo0,ivolfraclo1,ivolfraclo2
      integer ivolfrachi0,ivolfrachi1,ivolfrachi2
      REAL_T volfrac(
     &           ivolfraclo0:ivolfrachi0,
     &           ivolfraclo1:ivolfrachi1,
     &           ivolfraclo2:ivolfrachi2,
     &           0:nvolfraccomp-1)
      REAL_T norm
      REAL_T volume
      integer comp
      integer pval
      integer idoreg
      integer idoirr
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer i,j,k
      integer ncompsrc,counter
      REAL_T eps, rabs, vfrac
      logical isreg, iscov, isirr, usethis, useirr, usereg
      ncompsrc = nsrccomp
      if(pval .lt. 0) then
         call MAYDAY_ERROR
      endif
      if(ncompsrc .le. comp) then
         call MAYDAY_ERROR()
      endif
      eps = 1.0e-9
      usereg = (idoreg .eq. 1)
      useirr = (idoirr .eq. 1)
      counter=0
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

      vfrac = volfrac(i,j,k,0)
      iscov = (vfrac .lt. eps)
      isreg = (abs(one-vfrac) .lt. eps)
      isirr = (.not. isreg).and.(.not. iscov)
      usethis = ((isreg .and. usereg).or.(isirr .and. useirr))
      if(usethis) then
         counter=counter+1
         if(pval .eq. 0) then
            rabs  = abs(src(i,j,k,comp))
            norm  = max(norm,  rabs)
         elseif(pval .eq. 1) then
            rabs  = abs(src(i,j,k,comp))
            norm  = norm + vfrac*rabs
         elseif(pval .eq. 2) then
            rabs  = src(i,j,k,comp)*src(i,j,k,comp)
            norm  = norm + vfrac*rabs
         else
            call MAYDAY_ERROR()
         endif
         volume = volume + vfrac
      endif
      
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1) + counter*(pval+1)
      return
      end
      subroutine ADDTWOFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer D_DECL(i,j,k)
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &   ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp+numcomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        dst(i,j,k,ndest) +
     &        src(i,j,k,nsrc)
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*numcomp
      return
      end
      subroutine SCALEADDTWOFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,scale
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      REAL_T scale
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer D_DECL(i,j,k)
      integer nsrc, ndest, ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &   ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        dst(i,j,k,ndest) +
     &        src(i,j,k,nsrc) * scale
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*2*numcomp
      return
      end
      subroutine AXBYFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,x
     &           ,ixlo0,ixlo1,ixlo2
     &           ,ixhi0,ixhi1,ixhi2
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1,iylo2
     &           ,iyhi0,iyhi1,iyhi2
     &           ,nycomp
     &           ,a
     &           ,b
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nxcomp
      integer ixlo0,ixlo1,ixlo2
      integer ixhi0,ixhi1,ixhi2
      REAL_T x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           ixlo2:ixhi2,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1,iylo2
      integer iyhi0,iyhi1,iyhi2
      REAL_T y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           iylo2:iyhi2,
     &           0:nycomp-1)
      REAL_T a
      REAL_T b
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer D_DECL(i,j,k)
      integer nsrc, ndest, ncompsrc, ncompdst
      ncompsrc = nxcomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &   ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        x(i,j,k,nsrc) * a
     &      + y(i,j,k,nsrc) * b
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*3*numcomp
      return
      end
      subroutine AXBYFABCOMP(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,x
     &           ,ixlo0,ixlo1,ixlo2
     &           ,ixhi0,ixhi1,ixhi2
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1,iylo2
     &           ,iyhi0,iyhi1,iyhi2
     &           ,nycomp
     &           ,a
     &           ,b
     &           ,destcomp
     &           ,xcomp
     &           ,ycomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nxcomp
      integer ixlo0,ixlo1,ixlo2
      integer ixhi0,ixhi1,ixhi2
      REAL_T x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           ixlo2:ixhi2,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1,iylo2
      integer iyhi0,iyhi1,iyhi2
      REAL_T y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           iylo2:iyhi2,
     &           0:nycomp-1)
      REAL_T a
      REAL_T b
      integer destcomp
      integer xcomp
      integer ycomp
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer D_DECL(i,j,k)
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

      dst(i,j,k,destcomp) =
     &     x(i,j,k,xcomp) * a
     &     + y(i,j,k,ycomp) * b
      
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*3
      return
      end
      subroutine SUBTRACTTWOFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j,k
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &     ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        dst(i,j,k,ndest) -
     &        src(i,j,k,nsrc)
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*numcomp
      return
      end
      subroutine MULTIPLYTWOFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j,k
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &   ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        dst(i,j,k,ndest)*
     &        src(i,j,k,nsrc)
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*numcomp
      return
      end
      subroutine DIVIDETWOFAB(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,nsrccomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,srccomp
     &           ,destcomp
     &           ,numcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      integer nsrccomp
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2,
     &           0:nsrccomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer srccomp
      integer destcomp
      integer numcomp
      integer i,j,k
      integer nsrc, ndest,ncompsrc, ncompdst
      ncompsrc = nsrccomp
      ncompdst = ndstcomp
      if(((srccomp+numcomp).gt.ncompsrc).or.
     &   ((destcomp+numcomp).gt. ncompdst)) then
         call MAYDAY_ERROR()
      endif
      ndest = destcomp
      do nsrc = srccomp, srccomp + numcomp - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,ndest) =
     &        dst(i,j,k,ndest)/
     &        src(i,j,k,nsrc)
         
      enddo
      enddo
      enddo
         ndest = ndest + 1
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*numcomp
      return
      end
      subroutine SUBTRACTFABR(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      REAL_T src
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer i,j,k
      integer n, ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,n) =
     &        dst(i,j,k,n) - src
         
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*ncompdst
      return
      end
      subroutine ADDFABR(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      REAL_T src
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer i,j,k
      integer n, ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,n) =
     &        dst(i,j,k,n) + src
         
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*ncompdst
      return
      end
      subroutine MULTIPLYFABR(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      REAL_T src
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer i,j,k
      integer n,ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,n) =
     &        dst(i,j,k,n) * src
         
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*ncompdst
      return
      end
      subroutine DIVIDEFABR(
     &           dst
     &           ,idstlo0,idstlo1,idstlo2
     &           ,idsthi0,idsthi1,idsthi2
     &           ,ndstcomp
     &           ,src
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ndstcomp
      integer idstlo0,idstlo1,idstlo2
      integer idsthi0,idsthi1,idsthi2
      REAL_T dst(
     &           idstlo0:idsthi0,
     &           idstlo1:idsthi1,
     &           idstlo2:idsthi2,
     &           0:ndstcomp-1)
      REAL_T src
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer i,j,k
      integer n,ncompdst
      ncompdst = ndstcomp
      do n = 0, ncompdst - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         dst(i,j,k,n) =
     &        dst(i,j,k,n)/src
         
      enddo
      enddo
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(iregionhi1- iregionlo1+1)*(iregionhi2- iregionlo2+1)*ncompdst
      return
      end
      subroutine EBDOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0,ibfablo1,ibfablo2
     &           ,ibfabhi0,ibfabhi1,ibfabhi2
     &           ,nbfabcomp
     &           ,volfrac
     &           ,ivolfraclo0,ivolfraclo1,ivolfraclo2
     &           ,ivolfrachi0,ivolfrachi1,ivolfrachi2
     &           ,nvolfraccomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,icomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0,ibfablo1,ibfablo2
      integer ibfabhi0,ibfabhi1,ibfabhi2
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           ibfablo1:ibfabhi1,
     &           ibfablo2:ibfabhi2,
     &           0:nbfabcomp-1)
      integer nvolfraccomp
      integer ivolfraclo0,ivolfraclo1,ivolfraclo2
      integer ivolfrachi0,ivolfrachi1,ivolfrachi2
      REAL_T volfrac(
     &           ivolfraclo0:ivolfrachi0,
     &           ivolfraclo1:ivolfrachi1,
     &           ivolfraclo2:ivolfrachi2,
     &           0:nvolfraccomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer icomp
      integer i,j,k, ncomp, counter
      REAL_T eps
      if(icomp .lt. 0) then
         call MAYDAY_ERROR()
      endif
      ncomp = nafabcomp
      if(ncomp .le. icomp) then
         call MAYDAY_ERROR()
      endif
      ncomp = nbfabcomp
      if(ncomp .le. icomp) then
         call MAYDAY_ERROR()
      endif
      dotprodout = zero
      eps = 1.0e-10
      counter=0
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

      if(volfrac(i,j,k,0) .gt. eps) then
         counter=counter+1
         dotprodout = dotprodout +
     &        afab(i,j,k,icomp)*bfab(i,j,k,icomp)
      endif
      
      enddo
      enddo
      enddo
      ch_flops=ch_flops+counter*2
      return
      end
      subroutine MAXFAB(
     &           aval
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,acomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T aval
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer acomp
      integer i,j,k
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

        aval = max(aval, afab(i,j,k,acomp))
      
      enddo
      enddo
      enddo
      return
      end
      subroutine MINFAB(
     &           aval
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,acomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T aval
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer acomp
      integer i,j,k
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

        aval = min(aval, afab(i,j,k,acomp))
      
      enddo
      enddo
      enddo
      return
      end
