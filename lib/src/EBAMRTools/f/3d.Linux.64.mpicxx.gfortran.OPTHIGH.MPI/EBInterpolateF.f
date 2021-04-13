      subroutine EBINTERPCONSTANT(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,refratio
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2)
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer refratio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      do kcc = iblo2,ibhi2
      do jcc = iblo1,ibhi1
      do icc = iblo0,ibhi0
      do kk = ibreflo2,ibrefhi2
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0
      iff = icc*refratio + ii
      jff = jcc*refratio + jj
      kff = kcc*refratio + kk
      fine(iff,jff,kff) = coarse(icc,jcc,kcc)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EBCENTRALSLOPE(
     & slope
     & ,islopelo0,islopelo1,islopelo2
     & ,islopehi0,islopehi1,islopehi2
     & ,state
     & ,istatelo0,istatelo1,istatelo2
     & ,istatehi0,istatehi1,istatehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL*8 slope(
     & islopelo0:islopehi0,
     & islopelo1:islopehi1,
     & islopelo2:islopehi2)
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & istatelo2:istatehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer dir
      integer i,ii, j,jj, k,kk
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      do k = iblo2,ibhi2
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0
      slope(i,j,k) = (0.500d0) *(
     & state(i+ii, j+jj, k+kk) -
     & state(i-ii, j-jj, k-kk))
      enddo
      enddo
      enddo
      return
      end
      subroutine EBHISIDESLOPE(
     & slope
     & ,islopelo0,islopelo1,islopelo2
     & ,islopehi0,islopehi1,islopehi2
     & ,state
     & ,istatelo0,istatelo1,istatelo2
     & ,istatehi0,istatehi1,istatehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL*8 slope(
     & islopelo0:islopehi0,
     & islopelo1:islopehi1,
     & islopelo2:islopehi2)
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & istatelo2:istatehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer dir
      integer i,ii, j,jj, k,kk
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      do k = iblo2,ibhi2
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0
      slope(i,j,k) =
     & state(i+ii,j+jj,k+kk) -
     & state(i ,j ,k )
      enddo
      enddo
      enddo
      return
      end
      subroutine EBLOSIDESLOPE(
     & slope
     & ,islopelo0,islopelo1,islopelo2
     & ,islopehi0,islopehi1,islopehi2
     & ,state
     & ,istatelo0,istatelo1,istatelo2
     & ,istatehi0,istatehi1,istatehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL*8 slope(
     & islopelo0:islopehi0,
     & islopelo1:islopehi1,
     & islopelo2:islopehi2)
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & istatelo2:istatehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer dir
      integer i,ii, j,jj, k,kk
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      do k = iblo2,ibhi2
      do j = iblo1,ibhi1
      do i = iblo0,ibhi0
      slope(i,j,k) =
     & state(i ,j ,k ) -
     & state(i-ii,j-jj,k-kk)
      enddo
      enddo
      enddo
      return
      end
      subroutine EBMAXMINMOD(
     & mmslope
     & ,immslopelo0,immslopelo1,immslopelo2
     & ,immslopehi0,immslopehi1,immslopehi2
     & ,loslope
     & ,iloslopelo0,iloslopelo1,iloslopelo2
     & ,iloslopehi0,iloslopehi1,iloslopehi2
     & ,hislope
     & ,ihislopelo0,ihislopelo1,ihislopelo2
     & ,ihislopehi0,ihislopehi1,ihislopehi2
     & ,islopeboxlo0,islopeboxlo1,islopeboxlo2
     & ,islopeboxhi0,islopeboxhi1,islopeboxhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer immslopelo0,immslopelo1,immslopelo2
      integer immslopehi0,immslopehi1,immslopehi2
      REAL*8 mmslope(
     & immslopelo0:immslopehi0,
     & immslopelo1:immslopehi1,
     & immslopelo2:immslopehi2)
      integer iloslopelo0,iloslopelo1,iloslopelo2
      integer iloslopehi0,iloslopehi1,iloslopehi2
      REAL*8 loslope(
     & iloslopelo0:iloslopehi0,
     & iloslopelo1:iloslopehi1,
     & iloslopelo2:iloslopehi2)
      integer ihislopelo0,ihislopelo1,ihislopelo2
      integer ihislopehi0,ihislopehi1,ihislopehi2
      REAL*8 hislope(
     & ihislopelo0:ihislopehi0,
     & ihislopelo1:ihislopehi1,
     & ihislopelo2:ihislopehi2)
      integer islopeboxlo0,islopeboxlo1,islopeboxlo2
      integer islopeboxhi0,islopeboxhi1,islopeboxhi2
      integer i,j,k
      REAL*8 deltal, deltar, mono, rsign, finslope
      do k = islopeboxlo2,islopeboxhi2
      do j = islopeboxlo1,islopeboxhi1
      do i = islopeboxlo0,islopeboxhi0
      deltal = loslope(i,j,k)
      deltar = hislope(i,j,k)
      mono = deltal*deltar
      if(mono .gt. (0.0d0)) then
         rsign = sign((1.0d0), deltal + deltar)
         finslope = rsign*(min(abs(deltal), abs(deltar)))
      else
         finslope = (0.0d0)
      endif
      mmslope(i,j,k) = finslope
      enddo
      enddo
      enddo
      return
      end
      subroutine EBINTERPLINEAR(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,slope
     & ,islopelo0,islopelo1,islopelo2
     & ,islopehi0,islopehi1,islopehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,dir
     & ,refratio
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2)
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL*8 slope(
     & islopelo0:islopehi0,
     & islopelo1:islopehi1,
     & islopelo2:islopehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer dir
      integer refratio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer icc, jcc, kcc
      integer iff, jff, kff
      integer ii, jj, kk
      integer id
      REAL*8 dxf
      do kcc = iblo2,ibhi2
      do jcc = iblo1,ibhi1
      do icc = iblo0,ibhi0
      do kk = ibreflo2,ibrefhi2
      do jj = ibreflo1,ibrefhi1
      do ii = ibreflo0,ibrefhi0
      iff = icc*refratio + ii
      jff = jcc*refratio + jj
      kff = kcc*refratio + kk
      if(dir .eq. 0) then
         id = ii
      else if(dir .eq. 1) then
         id = jj
      else if(dir .eq. 2) then
         id = kk
      endif
      dxf = -(0.500d0) +((id+(0.500d0)) / refratio)
      fine(iff,jff,kff) =
     & fine(iff,jff,kff) +
     & dxf*slope(icc, jcc, kcc)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine EBINTERPSMOOTHERLINEAR(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,refratio
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2)
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer refratio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      return
      end
      subroutine EBINTERPQUADRATIC(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,refratio
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2)
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer refratio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      return
      end
      subroutine EBINTERPQUADSHIFT(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,coarse
     & ,icoarselo0,icoarselo1,icoarselo2
     & ,icoarsehi0,icoarsehi1,icoarsehi2
     & ,iblo0,iblo1,iblo2
     & ,ibhi0,ibhi1,ibhi2
     & ,refratio
     & ,ishift
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2)
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & icoarselo2:icoarsehi2)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer refratio
      integer ishift(0:2)
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      return
      end
