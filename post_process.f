      real*8, allocatable :: omega(:,:)
      real*8, allocatable :: alpha_r(:,:), alpha_i(:,:), beta(:,:)
      real*8, allocatable :: omega_r(:,:), omega_i(:,:)
      complex*16, allocatable :: omega2(:,:)
      real*8, allocatable :: temp(:)
      real*8 :: r_pt, lam_min,
     1   lam_max, cond_no, lambda_real, pi
      integer :: nang2, irads, iopt, m_emax, n_emax,
     1    itot
      pi = 4.*atan(1.d0)
      open(unit=4,file="alfven_spec",status="old")
      open(unit=7,file="alfven_post",status="unknown")
      open(unit=9,file="cond_no",status="unknown")
      open(unit=8,file="data_post",status="old")
      read(8,*) iopt, nang2, irads, isym_pos
c     if positive-definite matrix option is used in stellgap
c      (subroutine DSYGV) then isym_pos = 1
c     if general matrix option is used in stellgap
c      (subroutine DGEGV) then isym_pos = 0
c      isym_pos = 0
      write(*,*) iopt,nang2,irads, isym_pos
      allocate(omega(irads,nang2),stat=istat)
      allocate(omega2(irads,nang2),stat=istat)
      allocate(omega_r(irads,nang2),stat=istat)
      allocate(omega_i(irads,nang2),stat=istat)
      allocate(beta(irads,nang2),stat=istat)
      allocate(alpha_r(irads,nang2),stat=istat)
      allocate(alpha_i(irads,nang2),stat=istat)
      allocate(temp(nang2),stat=istat)
      itot = irads*nang2
c      fnorm = 1./(2000.*pi)    !converts omega from radians/sec to kHz
      fnorm = 1.      !conversion to kHz is now done in stellgap
c      fnorm = 1./224.
      write(7,'(i10)') itot
      lam_min = 1.d+10
      lam_max = -1.d+10
      do ir=1,irads
       do i=1,nang2
        if(isym_pos .eq. 1) then
        read(4,'(2(e15.7,2x),i4,2x,i4)') r_pt, omega(ir,i),
     1   m_emax, n_emax
       temp(i) = omega(ir,i)**2
       write(7,'(2(e15.7,2x),i4,2x,i4)') r_pt, omega(ir,i)*fnorm,
     1   m_emax, n_emax
     
       else if(isym_pos .eq. 0) then
        read(4,'(4(e15.7,2x),i4,2x,i4)') r_pt, alpha_r(ir,i),
     1     alpha_i(ir,i),beta(ir,i),m_emax, n_emax
c      alpha_r(ir,i) = abs(alpha_r(ir,i))
      if(beta(ir,i) .gt. 0.d0) then
         lambda_real = alpha_r(ir,i)/beta(ir,i)
         if(lambda_real .gt. 0.) then
           lam_min = min(lam_min,abs(lambda_real))
           lam_max = max(lam_max,abs(lambda_real))
         endif
         lambda_imag = alpha_i(ir,i)/beta(ir,i)
         omega2(ir,i) = dcmplx(lambda_real,lambda_imag)
       else if(beta(ir,i) .le. 0.d0) then
         omega2(ir,i) = 1.e+10
       endif
       omega_r(ir,i) = real(cdsqrt(omega2(ir,i)))
       omega_i(ir,i) = imag(cdsqrt(omega2(ir,i)))
       if(omega_i(ir,i) .eq. 0.d0) then
        if(r_pt .lt. 0.99) then
         write(7,'(2(e15.7,2x),i4,2x,i4)') r_pt, omega_r(ir,i)*fnorm,
     1     m_emax, n_emax
        end if
       end if
     
       endif
       
       end do
       if(isym_pos .eq. 1) then
         lam_min = minval(temp)
         lam_max = maxval(temp)
       endif
       if(lam_min .ne. 0.) cond_no = lam_max/lam_min
       write(9,'(e15.7,2x,e15.7)') r_pt, cond_no
c       write(*,*) r_pt, cond_no
       end do
       end
c
