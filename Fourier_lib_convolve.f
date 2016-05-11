      subroutine toFourier
      use fourier_lib
      implicit none
c
c   Do Fourier transform integrations needed to convert data on a
c   theta, zeta grid [stored in array f(i=1,nznt)] to a set
c   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
c   Typically, the number of grid points in each direction needs
c   to be > 3*number of modes used in each direction to avoid
c   aliasing errors(implies nznt > 9*mnmx).
c   
c      do mn=1,mnmx     !loop over Fourier modes
c       fnm(mn) = 0.
c      do i=1,nznt      !loop over theta,zeta grid
c       if(sin_type .eq. 1 .and. cos_type .eq. 0) then
c         fnm(mn) = fnm(mn) + f(i)*sin_toF(i,mn)
c       else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
c         fnm(mn) = fnm(mn) + f(i)*cos_toF(i,mn)
c       endif
c      end do
c      end do
c
       if(sin_type .eq. 1 .and. cos_type .eq. 0) then
         fnm = matmul(f,sin_toF)
       else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
         fnm = matmul(f,cos_toF)
       endif
c
      return
      end
c
      subroutine old_toFourier
      use fourier_lib
      implicit none
      real(kind=rprec) dum,dnorm
c
c   Do Fourier transform integrations needed to convert data on a
c   theta, zeta grid [stored in array f(i=1,nznt)] to a set
c   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
c   Typically, the number of grid points in each direction needs
c   to be > 3*number of modes used in each direction to avoid
c   aliasing errors(implies nznt > 9*mnmx).
c   
      do mn=1,mnmx     !loop over Fourier modes
       fnm(mn) = 0.
c       dnorm = 2.*real(nfp)/real(nznt)
       dnorm = 2./real(nznt)
       dum = abs(rn(mn)) + abs(rm(mn))
       if(nint(dum) .eq. 0) dnorm = .5*dnorm
      do i=1,nznt      !loop over theta,zeta grid
c       arg = -rn(mn)*ztgrd(i) + rm(mn)*thtgrd(i)
       if(sin_type .eq. 1 .and. cos_type .eq. 0) then
         fnm(mn) = fnm(mn) + f(i)*sin_ar(i,mn)*dnorm
       else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
         fnm(mn) = fnm(mn) + f(i)*cos_ar(i,mn)*dnorm
       endif
      end do
      end do
c
      return
      end
c
      subroutine toReal
c
c    Convert Fourier mode representation [stored in array anm(mn=1,mnmx)]
c    to values of function on a regularly spaced 2D grid
c    [stored in array f(i=1,nznt)].
c
      use fourier_lib
      implicit none
      do i=1,nznt
      f(i) = 0.
        do mn=1,mnmx
         if(sin_type .eq. 1 .and. cos_type .eq. 0) then
          f(i) = f(i) + anm(mn)*sin_ar(i,mn)
         else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
          f(i) = f(i) + anm(mn)*cos_ar(i,mn)
         endif
        end do
       end do
      end
c
c
c
      subroutine dbydth
c
c    Take the theta derivative of the input Fourier amplitude array, fnm
c    and place the result in the output Fourier amplitude array, anm.
c    Changes to the sin/cos parity are reflected through the sin_type and
c    cos_type variables.
c
      use fourier_lib
      implicit none
      do i=1,mnmx
         if(sin_type .eq. 1 .and. cos_type .eq. 0) then
          anm(i) = rm(i)*fnm(i)
         else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
          anm(i) = -rm(i)*fnm(i)
         endif
        end do
        if(sin_type .eq. 1 .and. cos_type .eq. 0) then
          sin_type = 0; cos_type = 1
         else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
          sin_type = 1; cos_type = 0
        endif
       end
c
c
c
      subroutine dbydzt
c
c    Take the zeta derivative of the input Fourier amplitude array, fnm
c    and place the result in the output Fourier amplitude array, anm.
c    Changes to the sin/cos parity are reflected through the sin_type and
c    cos_type variables.
c
      use fourier_lib
      implicit none
      do i=1,mnmx
         if(sin_type .eq. 1 .and. cos_type .eq. 0) then
          anm(i) = -rn(i)*fnm(i)
          else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
          anm(i) = rn(i)*fnm(i)
         endif
        end do
        if(sin_type .eq. 1 .and. cos_type .eq. 0) then
          sin_type = 0; cos_type = 1
         else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
          sin_type = 1; cos_type = 0
        endif
       end
