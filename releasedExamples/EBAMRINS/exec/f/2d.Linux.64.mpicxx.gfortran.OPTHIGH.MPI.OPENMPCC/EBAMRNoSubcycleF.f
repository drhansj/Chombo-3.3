      subroutine MULTVELINTOFLUX(
     & flux
     & ,ifluxlo0,ifluxlo1
     & ,ifluxhi0,ifluxhi1
     & ,nfluxcomp
     & ,vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,nvelcomp
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & ,ncomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & 0:nfluxcomp-1)
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & 0:nvelcomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer ncomp
      integer i,j
      integer ivar
      do ivar=0, ncomp-1
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
         flux(i,j, ivar) =
     & flux(i,j, ivar)*vel(i,j, 0)
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTEVORT(
     & vort
     & ,ivortlo0,ivortlo1
     & ,ivorthi0,ivorthi1
     & ,vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,nvelcomp
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & ,dx
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivortlo0,ivortlo1
      integer ivorthi0,ivorthi1
      REAL*8 vort(
     & ivortlo0:ivorthi0,
     & ivortlo1:ivorthi1)
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & 0:nvelcomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      REAL*8 dx
      integer dir
      integer i,j
      integer dir1, dir2
      REAL*8 oneOnDx
      REAL*8 vortTemp
      integer ii1, jj1,kk1
      integer ii2,jj2,kk2
      oneOnDx = (0.500d0)/dx
      dir1 = -1
      dir2 = -1
      if ((2 .eq.2).or.(dir.eq.2)) then
         dir1 = 0
         dir2 = 1
      else if (dir.eq.0) then
         dir1= 1
         dir2= 2
      else if (dir.eq.1) then
         dir1=2
         dir2=0
      else
         write(*,*) 'COMPUTEVORT: bad direction = ', dir
         call MAYDAY_ERROR()
      endif
      ii1 = CHF_ID(dir1,0)
      jj1 = CHF_ID(dir1,1)
      kk1 = CHF_ID(dir1,2)
      ii2 = CHF_ID(dir2,0)
      jj2 = CHF_ID(dir2,1)
      kk2 = CHF_ID(dir2,2)
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
      vortTemp = oneOnDx*((vel(i+ii1,j+jj1,dir2)
     & -vel(i-ii1,j-jj1,dir2))
     & -(vel(i+ii2,j+jj2,dir1)
     & -vel(i-ii2,j-jj2,dir1)))
      vort(i,j) = oneOnDx*((vel(i+ii1,j+jj1,dir2)
     & -vel(i-ii1,j-jj1,dir2))
     & -(vel(i+ii2,j+jj2,dir1)
     & -vel(i-ii2,j-jj2,dir1)))
      enddo
      enddo
      return
      end
      subroutine KINETICENERGY(
     & energy
     & ,ienergylo0,ienergylo1
     & ,ienergyhi0,ienergyhi1
     & ,vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,nvelcomp
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ienergylo0,ienergylo1
      integer ienergyhi0,ienergyhi1
      REAL*8 energy(
     & ienergylo0:ienergyhi0,
     & ienergylo1:ienergyhi1)
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & 0:nvelcomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer i, j
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
      energy(i,j) = (0.500d0)*(vel(i,j,0)**2
     & +vel(i,j,1)**2)
      enddo
      enddo
      return
      end
      subroutine MAGVECT(
     & magVector
     & ,imagVectorlo0,imagVectorlo1
     & ,imagVectorhi0,imagVectorhi1
     & ,vect
     & ,ivectlo0,ivectlo1
     & ,ivecthi0,ivecthi1
     & ,nvectcomp
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imagVectorlo0,imagVectorlo1
      integer imagVectorhi0,imagVectorhi1
      REAL*8 magVector(
     & imagVectorlo0:imagVectorhi0,
     & imagVectorlo1:imagVectorhi1)
      integer nvectcomp
      integer ivectlo0,ivectlo1
      integer ivecthi0,ivecthi1
      REAL*8 vect(
     & ivectlo0:ivecthi0,
     & ivectlo1:ivecthi1,
     & 0:nvectcomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer i,j
      integer dir
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
         magVector(i,j) = 0
      enddo
      enddo
      do dir=0, nvectcomp-1
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
            magVector(i,j) = magVector(i,j)
     & +vect(i,j,dir)*vect(i,j,dir)
      enddo
      enddo
      enddo
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0
        magVector(i,j) = sqrt(magVector(i,j))
      enddo
      enddo
      return
      end
      subroutine CENTERED_DERIV(
     & deriv
     & ,iderivlo0,iderivlo1
     & ,iderivhi0,iderivhi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,igridlo0,igridlo1
     & ,igridhi0,igridhi1
     & ,dx
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iderivlo0,iderivlo1
      integer iderivhi0,iderivhi1
      REAL*8 deriv(
     & iderivlo0:iderivhi0,
     & iderivlo1:iderivhi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer igridlo0,igridlo1
      integer igridhi0,igridhi1
      REAL*8 dx
      integer dir
      integer i,j
      REAL oneOnDx
       oneOnDx = 1.0d0/dx
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0
         deriv(i,j) = oneOnDx*(phi(i+CHF_ID(dir,0),
     & j+CHF_ID(dir,1))
     & -phi(i-CHF_ID(dir,0),
     & j-CHF_ID(dir,1)))
      enddo
      enddo
       return
       end
