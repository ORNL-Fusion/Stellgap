*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program stellgap.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
c   Before running, check that the irads variable is set consistently with the
c   number of surfaces in the tae_boozer_data file. irads should be the number
c   of surfaces in the original VMEC run minus 1 or 2 to avoid edge problems.
c
c  10/21/21
c   The adiabatic index, gamma, is set to 1.35, as is more appropriate to AE simulations
c   (For example, see VAN ZEELAND, M. A., HEIDBRINK, W. W., SHARAPOV, S. E., et al,
c   “Electron cyclotron heating can drastically alter reversed shear Alfvén eigenmode
c   activity in DIII-D through finite pressure effects,” Nucl. Fusion 56 112007 (2016).)
c   gamma is currently set in the source code, but can be varied and the code
c   recompiled, if needed.
c
c  10/28/16 modifications: added ion_profile = 4 option. This uses a polynomial
c   fit similar to ion_profile = 1, but the independent variable is sqrt(toroidal flux)
c   rather than toroidal flux. Also, the same option of using sqrt(toroidal flux) is
c   added for the electron temperature profile. A control variable te_profile is added.
c   For te_profile 1 the usual poly fit with toroidal flux is used (default). For
c   te_profile = 2, a poly fit to sqrt(toroidal flux) is used.
c
c  Modified 5/22/2013 making irads and ir_fine_scl to be command-line arguments
c   (if not given, default values of 128 and 41 are used)
c            e.g., xstgap_snd 159 2048
c
c  12/12/2011: Slow-sound option included. By setting slow_sound = .true. the
c   dense sound continua at low frequencies can be removed so that one only
c   sees the low-frequency gap of the Alfven continua. This is based on a
c   low beta approximation and was first suggested in M.S. Chu, et al.,
c   Phys. Fluids B: 4(11), 3713 (1992). The basic idea is to remove the sound
c   continua since these continua are rather dense and the sound waves
c   will experience strong continuum damping. Also, see W. Deng, et al.,
c   NF (2012), Appendix A.
c  Modified 12/2008 so that ion density profile, ion mass and several other
c   parameters are read in through the plasma.dat file rather than being
c   hardwired into the code. (see explanation below for variables contained
c   in plasma.dat file).
c
c  Modified on 9/18/2008 to offer both parallel or serial options, depending on
c   use of the precompiler.
c
c This version implements the coupled Alfven/soundwave continuum equations
c  that were derived from a Langrangian starting point.
c
c compile/load using:
c
c   cp stellgap_soundwave.f temp.c
c   cpp -E -P -C -DSERIAL temp.c > temp.f
c   ifort $OPT -c temp.f
c   mv temp.o stellgap_soundwave.o
c   rm temp.c temp.f
c   ifort -c $OPT Fourier_lib_convolve.f
c   ifort -c $OPT fitpack.f
c   ifort $OPT -o xstgap_snd stellgap_soundwave.o Fourier_lib_convolve.o fitpack.o -I. \
c    -L$MKLPATH -I$MKLINCLUDE -lmkl_intel -lmkl_lapack -lmkl_core -lguide -lmkl_intel_thread
c
c  This version of stellgap.f has been modified to be compatible with ae_solve.f.
c  The frequencies are output in kHz. The mode selection is done using the fourier.dat
c  file. The equilibrium data comes from the tae_data_boozer file. Following the
c  stellgap run, a post-processing code, post_process.f must be run. The gap plot
c  is then made using the results in the alfven_post file.
c  The central ion density, ion density profile, and ion mass are set from the
c  plasma.dat input file. 
c  Two other variables that the user should set are irads and ir_fine_scl. These
c  are set througth paramter statements in the beginning of program tae_continua:
c      integer, parameter :: ir_fine_scl = 1000
c      integer, parameter :: irads = 39, irad3 = 3*irads
c  irads = # of flux surfaces in the initial equilibrium (wout file) - 2.
c  The outermost and innermost surfaces from the initial equilibrium
c  are left off to avoid inaccuracies that are sometimes present near
c  these regions.
c  ir_fine_scl = # of surfaces desired in the continuum output. Generally, to
c  resolve fine scale features in the Alfven continua, more surfaces are
c  desireable than in the original equilibrium data. The code performs
c  interpolations to allow this.
c
c
c  Input parameters:   (fourier.dat file).
c
c   ith, izt = number of theta and zeta grid points in the original mesh data
c     (in file tae_data_vmec) calculated by xcobra_vmec_tae.
c    This file contains data
c    for iota, phipc, |B|, gssup (the psi-psi contra-variant metric element,
c    the Jacobian, and the contravariant poloidal and toroidal components of
c    the magnetic field)
c   nfp = number of field periods (for a tokamak run, set nfp = 1)
c   mpol, ntor = number of poloidal and toroidal modes used in the Fourier
c    expansion of the above gridded quantities suplied by the tae_data_vmec
c    data file.  To avoid anti-aliasing errors, mpol < 0.5*ith and
c    ntor < 0.5*izt.
c   mp_col = number of m,n mode pairs used in the representation of the
c    eigebfunction
c   nt_col = number of toroidal modes used in the representation of the
c    eigenfunction (should be an integral multiple of nfp)
c   mode_family = toroidal mode number about which the toroidal eigenfunction
c     spectrum is built. i.e., for m = 0,
c     n = mode_family, mode_family + nfp, mode_family + 2*nfp,
c         ..., mode_family + nt_col
c      while for m .ne. 0,
c      n = -nt_col + mode_family, -nt_col + nfp + mode_family, -nt_col
c            + 2*nfp + mode_family, ..., nt_col + mode_family
c
c   Input parameters:   (plasma.dat file).
c
c   ion_to_proton_mass = m_ion/m_proton
c   ion_density_0 = n_ion(0) = ion density (in m**-3) at magnetic axis
c   ion_profile = integer parameter that determines form of n_ion(r)/n_ion(0) fit
c      for ion_profile = 0 ion_density = [iota(rho)/iota(0)]**2
c      for ion_profile = 1 ion_density = polynomial fit
c                          = nion(1) + nion(2)*rho + nion(3)*(rho**2)
c                           + nion(4)*(rho**3) + ... + nion(9)*(rho**8)
c               (here rho = normalized toroidal flux)
c      for ion_profile = 2 ion_density = ion_density = constant = 1.
c      for ion_profile = 3 ion_density = [1 - aion*(rho**bion)]**cion
c
c   nion(9) = array of polynomial fit coefficients for ion_profile = 1
c   aion, bion, cion = parameters for ion_profile = 3 fit
c   jdqz_data = logical variable that is used in ae3d, but not in stellgap.f
c                (included here only so that same plasma.dat file can be used)
c   egnout_form = "binr" or "asci" - like jdqz_data this variable is used in
c                  ae3d, but not stellgap.f
c
c compile: sh bld_stellgap_soundwave_ver1
c
c
c
      module kind_spec
c
c     Double precision
c
      integer, parameter :: rprec = selected_real_kind(12,100)
      integer, parameter :: iprec = selected_int_kind(8)
      end module kind_spec
c
c Variables that need to be set for each time the configuration is changed:
c   ith, izt, nfp, mpol, ntor, mp_col, nt_col, mode_family, tokamak
c
      module fourier_lib
      use kind_spec
      implicit none
      integer::ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
      real(kind=rprec), parameter :: parity_gss = 1.,
     1     parity_bsupth = 1., parity_bsupzt = 1.
      integer :: i, j, lg, nu, nl, m, n, mn, istat, sin_type, cos_type
c   Equilibrium coefficient arrays
      real(kind=rprec), allocatable :: thtgrd(:),ztgrd(:),rn(:),rm(:)
      real(kind=rprec), allocatable :: fnm(:),f(:),anm(:)
      real(kind=rprec), allocatable :: cos_ar(:,:),sin_ar(:,:)
      real(kind=rprec), allocatable :: cos_toF(:,:),sin_toF(:,:)
      real(kind=rprec) twopi, arg
      logical, parameter :: ipos_def_sym=.false.
c    The subset_eq flag allows one to just request a bracketed subset of
c     eigenvalues, rather than all eigenvalues. This option has not been 
c     developed beyond the call. Initial tests have not indicated it
c     speeds things up.
      logical, parameter :: subset_eq = .false.
c   Eigenfunction arrays,variables
      integer :: count, mn_col, ith_col, izt_col
      real(kind=rprec), allocatable :: rn_col(:),rm_col(:)
      real(kind=rprec), allocatable :: rn2_col(:),rm2_col(:),rnm_col(:)
      integer, allocatable :: in_col(:), im_col(:)
      integer, allocatable :: in_col_aug(:), im_col_aug(:)
      integer, allocatable :: nw(:), mwl(:), mwu(:)
      contains
c      
      subroutine readin
       open(unit=20,file="fourier.dat",status="old")
       read(20,*) nfp, ith, izt, mode_family
       mpol = ith*2/5; ntor = izt*2/5
       nznt=ith*izt; mnmx=(2*ntor+1)*mpol-ntor
       read(20,*) ntors
       allocate(nw(ntors), stat=istat)
       allocate(mwl(ntors), stat=istat)
       allocate(mwu(ntors), stat=istat)
       do i=1,ntors
        read(20,*) nw(i), mwl(i), mwu(i)
       end do       
       close(unit=20)
      end subroutine readin
c            
c      
      subroutine trig_array
      use kind_spec
      implicit none
      real(kind=rprec) :: dum,dnorm
      allocate(thtgrd(nznt), stat=istat)
      allocate(ztgrd(nznt), stat=istat)
      allocate(rm(mnmx), stat=istat)
      allocate(rn(mnmx), stat=istat)
      allocate(fnm(mnmx), stat=istat)
      allocate(f(nznt), stat=istat)
      allocate(anm(mnmx), stat=istat)
      allocate(cos_ar(nznt,mnmx), stat=istat)
      allocate(sin_ar(nznt,mnmx), stat=istat)
      allocate(cos_toF(nznt,mnmx), stat=istat)
      allocate(sin_toF(nznt,mnmx), stat=istat)
      twopi = 8.*atan(1.)
c    Generate theta, zeta grid
      lg = 0
      do i=1,izt
      do j=1,ith
       lg = lg + 1
       ztgrd(lg) = twopi*real(i-1)/(real(nfp*izt))
       thtgrd(lg) = twopi*real(j-1)/real(ith)
      end do
      end do
c    Generate Fourier mode distribution
      mn = 0
      nu = ntor
      do m=0,mpol-1
       nl = -ntor
       if(m .eq. 0) nl = 0
        do n = nl,nu
         mn = mn + 1
         rm(mn) = real(m)
         rn(mn) = real(n*nfp)
        enddo
      enddo
      do i=1,nznt
       do mn=1,mnmx
        arg = -rn(mn)*ztgrd(i) + rm(mn)*thtgrd(i)
        cos_ar(i,mn) = cos(arg)
        sin_ar(i,mn) = sin(arg)
        dnorm = 2./real(nznt)
        dum = abs(rn(mn)) + abs(rm(mn))
        if(nint(dum) .eq. 0) dnorm = .5*dnorm
        cos_toF(i,mn) = cos_ar(i,mn)*dnorm
        sin_toF(i,mn) = sin_ar(i,mn)*dnorm
       end do
      end do
      end subroutine trig_array
      
      subroutine convolution_array
      use kind_spec
      implicit none
c     First, count the number of modes to be used
      count = 0
      do i=1,ntors
       count = count + mwu(i) - mwl(i) + 1
      end do
      mn_col = count
c     Allocate memory for arrays needed in eigenfunction calculations
      allocate(rm_col(mn_col), stat=istat)
      allocate(rn_col(mn_col), stat=istat)
      allocate(im_col(mn_col), stat=istat)
      allocate(in_col(mn_col), stat=istat)
      allocate(rm2_col(mn_col), stat=istat)
      allocate(rn2_col(mn_col), stat=istat)
      allocate(rnm_col(mn_col), stat=istat)
c     Create mode number arrays needed for eigenfunction
       count = 0
       do n = 1,ntors
        do m = mwl(n), mwu(n)
         count = count + 1
         rm_col(count) = real(m)
         rn_col(count) = real(nw(n))
         im_col(count) = m
         in_col(count) = nw(n)
         rm2_col(count) = rm_col(count)*rm_col(count)
         rn2_col(count) = rn_col(count)*rn_col(count)
         rnm_col(count) = rm_col(count)*rn_col(count)
        end do
      end do         
      mn_col = count
c      write(*,*) mn_col, count
      end subroutine convolution_array
c
c
      subroutine scs_convolve(ans,m1,n1,m2,n2,meq,neq)
      use kind_spec
      implicit none
      real(kind=rprec) :: tht_int1, tht_int2, tht_int3, tht_int4,
     1  zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
      integer :: m1, m2, n1, n2, meq, neq
      integer :: sm1, sm2, sn1, sn2, smeq, sneq
c
c     This subroutine calculates the 2D integral for tht and zeta running
c     from 0 to 2*PI of:
c     sin(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*sin(m2*tht - n2*zeta)
c     The value returned is actually (2/PI) times this integral.
c
      sm1=1;sm2=1;sn1=1;sn2=1;smeq=1;sneq=1
      if(m1 .ne. 0) sm1 = m1/abs(m1)
      if(m2 .ne. 0) sm2 = m2/abs(m2)
      if(n1 .ne. 0) sn1 = n1/abs(n1)
      if(n2 .ne. 0) sn2 = n2/abs(n2)
      if(meq .ne. 0) smeq = meq/abs(meq)
      if(neq .ne. 0) sneq = neq/abs(neq)
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      call css(tht_int1, meq, m1, m2)
      call ccc(zeta_int1, n1, n2, neq)
      call ccc(tht_int2, m1, m2, meq)
      call css(zeta_int2, neq, n1, n2)
      call css(tht_int3, m2, m1, meq)
      call css(zeta_int3, n1, n2, neq)
      call css(tht_int4, m1, m2, meq)
      call css(zeta_int4, n2, n1, neq)
      ans = tht_int1*zeta_int1*sm1*sm2
     1    + tht_int2*zeta_int2*sn1*sn2
     2    - tht_int3*zeta_int3*sm1*smeq*sn2*sneq
     3    - tht_int4*zeta_int4*sm2*smeq*sn1*sneq
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      end subroutine scs_convolve
c
c
c
      subroutine ccc_convolve(ans,m1,n1,m2,n2,meq,neq)
      use kind_spec
      implicit none
      real(kind=rprec) :: tht_int1, tht_int2, tht_int3, tht_int4,
     1  zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
      integer :: m1, m2, n1, n2, meq, neq
      integer :: sm1, sm2, sn1, sn2, smeq, sneq
c
c     This subroutine calculates the 2D integral for tht and zeta running
c     from 0 to 2*PI of:
c     cos(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*cos(m2*tht - n2*zeta)
c     The value returned is actually (2/PI) times this integral.
c
      sm1=1;sm2=1;sn1=1;sn2=1;smeq=1;sneq=1
      if(m1 .ne. 0) sm1 = m1/abs(m1)
      if(m2 .ne. 0) sm2 = m2/abs(m2)
      if(n1 .ne. 0) sn1 = n1/abs(n1)
      if(n2 .ne. 0) sn2 = n2/abs(n2)
      if(meq .ne. 0) smeq = meq/abs(meq)
      if(neq .ne. 0) sneq = neq/abs(neq)
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      call ccc(tht_int1, m1, m2, meq)
      call ccc(zeta_int1, n1, n2, neq)
      call css(tht_int2, meq, m1, m2)
      call css(zeta_int2, neq, n1, n2)
      call css(tht_int3, m1, m2, meq)
      call css(zeta_int3, n1, n2, neq)
      call css(tht_int4, m2, m1, meq)
      call css(zeta_int4, n2, n1, neq)
      ans = tht_int1*zeta_int1
     1    + tht_int2*zeta_int2*sm1*sm2*sn1*sn2
     2    + tht_int3*zeta_int3*sm2*smeq*sn2*sneq
     3    + tht_int4*zeta_int4*sm1*smeq*sn1*sneq
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      end subroutine ccc_convolve
c
c
c
      subroutine ccc(result, i, j, k)
      use kind_spec
      real(kind=rprec), parameter :: zero = 0, one = 1,
     1    two = 2, four = 4
      real(kind=rprec) :: result
      integer :: i, j, k, izeros
c
c     This subroutine calculates the 1D integral of
c     cos(i*x)*cos(j*x)*cos(k*x)
c     for x running from 0 to 2*PI.
c     The value returned is actually (2/PI) times this integral.
c
      result = zero
      izeros = 0
      if(i .eq. 0) izeros = izeros + 1
      if(j .eq. 0) izeros = izeros + 1
      if(k .eq. 0) izeros = izeros + 1
      if(izeros .ne. 0) then
       if(izeros .eq. 3) then
        result = four
        return
       else if(izeros .eq. 2) then
        result = zero
        return
       else if(izeros .eq. 1) then
        if(i .eq. 0) then
         result = zero
         if(j*j .eq. k*k) then
          result = two
         endif
         return
        else if(j .eq. 0) then
          result = zero
          if(i*i .eq. k*k) then
          result = two
         endif
         return
        else if(k .eq. 0) then
          result = zero
          if(i*i .eq. j*j) then
          result = two
         endif
         return
        endif
       endif
      else if (izeros .eq. 0) then
       if(k .eq. (i-j)) result = one
       if(k .eq. (i+j)) result = one
       if(k .eq. (j-i)) result = one
       if(k .eq. -(i+j)) result = one
      endif
      end subroutine ccc
c
c
c
      subroutine css(result, k, i, j)
      use kind_spec
      real(kind=rprec), parameter :: one = 1, neg_one = -1, zero = 0,
     1  two = 2, neg_two = -2
      integer :: i, j, k
      real(kind=rprec) :: result
c
c     This subroutine calculates the 1D integral of
c     cos(k*x)*sin(i*x)*sin(j*x)
c     for x running from 0 to 2*PI.
c     The value returned is actually (2/PI) times this integral.
c
      result = zero
      if(i .eq. 0 .or. j .eq. 0) return
      if(k .eq. 0) then
        if(i*i .ne. j*j) then
         return
        else if(i*i .eq. j*j) then
         result = 2
         return
        endif
      endif
      if(k .eq. (i-j)) result = one
      if(k .eq. (i+j)) result = neg_one
      if(k .eq. (j-i)) result = one
      if(k .eq. -(i+j)) result = neg_one
      end subroutine css
c
c
      subroutine trg_deallocate
      deallocate(thtgrd,ztgrd,rm,rn,cos_ar,sin_ar,f,fnm,anm,
     1   rm_col,rn_col,im_col,in_col,
     2   rm2_col,rn2_col,rnm_col)
      end subroutine trg_deallocate
      
      end module fourier_lib
c-----------------------------------------------------------------
      program tae_continua
      use fourier_lib
      use kind_spec
      implicit none
#if defined (PARALLEL)      
      include 'mpif.h'
#endif
c
c
c
      integer :: ir_fine_scl, irads, irad3
c      integer, parameter :: ir_fine_scl = 256
      integer, parameter :: ir_out = 12800
      logical, parameter ::  mat_out = .false.
c      integer, parameter :: irads = 89, irad3 = 3*irads
      integer, parameter :: iopt = 1, iradt = 300
      real(kind=rprec), parameter :: R0 = 1.0
      real(kind=rprec), parameter :: mass_proton = 1.67d-27
      real(kind=rprec) :: ion_density_0, mass_ion, ion_to_proton_mass,
     1   aion, bion, cion
      real(kind=rprec) :: nion(10), telec(10)
      integer :: ir, irr, il, iu, num_eq, jj, ion_profile,
     >  te_profile
      character*4 :: egnout_form
      logical :: jdqz_data
      real(kind=rprec), parameter :: delta_iota = 0.0d0
      real(kind=rprec) :: temp_elec_0
      real(kind=rprec), parameter :: ev_to_jou = 1.6d-19
      real(kind=rprec) :: gamma
      real(kind=rprec), parameter :: zion = 1.d0
      real(kind=rprec), dimension(:,:,:), allocatable :: bfield,
     1  gsssup,rjacob, bsupth, bsupzt, bfield_lrg,
     2  gsssup_lrg, bsupth_lrg, bsupzt_lrg, rjacob_lrg
      real(kind=rprec), dimension(:), allocatable :: theta_tae,
     1 zeta_tae, f1_nm, f3a_nm, f3b_nm, f3c_nm, bsupth_tmp, bsupzt_tmp,
     2 g1a_nm, g1b_nm, g1c_nm, g2_nm, g3_nm, c2b_nm, c1_nm,
     3 g4_nm, g5_nm, eta_sq, zeta_sq
      real(kind=rprec), dimension(:,:), allocatable :: f1, f3a,
     1 f3b, f3c, g, g1a, g1b, g1c, g2, g3, c2b, c1, dB_dtht, dB_dzet,
     2 g4, g5, dlg_dzet, dlg_dtht, dsg_dzet, dsg_dtht
     
c      real(kind=rprec), dimension(ir_fine_scl) :: bavg, cs2, te
c      real(kind=rprec), dimension(ir_fine_scl) :: ion_density, iota_r,
c     1  jtor_lrg, jpol_lrg, phipc_lrg, rho_ion
c      real(kind=rprec), dimension(irads) :: nsurf,iotac,phipc,rho
c      real(kind=rprec), dimension(irads) :: sp_fit, jpol, jtor, bav
c      real(kind=rprec), dimension(irad3) :: yp, temp, ypi, tempi,
c     1  ypjt, ypjp, ypphp

      real(kind=rprec), allocatable, dimension(:) :: bavg, cs2, te
      real(kind=rprec), allocatable, dimension(:) :: ion_density,
     1   iota_r, jtor_lrg, jpol_lrg, phipc_lrg, rho_ion
      real(kind=rprec), allocatable, dimension(:) :: nsurf,iotac,
     1   phipc,rho, sp_fit, jpol, jtor, bav, yp, temp, ypi, tempi,
     1   ypjt, ypjp, ypphp
     
      real(kind=rprec), allocatable :: beta(:), aux(:), eig_vect(:),
     1  gss0(:), gss1(:), gss2(:), gss3(:), lg0(:), lg1(:),
     2  lg2(:), lg3(:), b0(:), b1(:), b2(:), b3(:)
      real(kind=rprec), allocatable :: amat(:,:),temp_amat(:,:),
     1           temp_bmat(:,:),bmat(:,:), z_r(:,:),
     2           vr(:,:), vl(:,:), vr_keep(:,:), b11(:,:),
     3           b22(:,:)
      real(kind=rprec), allocatable :: bgrad2(:,:,:)
      real(kind=rprec), allocatable :: omega(:), work(:), alfr(:),
     1  alfi(:), om_keep(:), phi_temp0(:), phi_temp1(:),
     2  psi_temp0(:), psi_temp1(:), phi_norm(:), psi_norm(:)
      complex*16, allocatable :: alpha(:)
      real(kind=rprec) :: sp1, sp2, sigma_spl, r_pt,
     1  iotac_r, curv2, ra, ra2, ra3, ra4, ra5, ra6
      real(kind=rprec) :: amat1, amat2, amat3, amat4, amat5
      real(kind=rprec) :: f1_coef, f2a_coef, f3a_coef,
     1  f3b_coef, f3c_coef, ftest
      real(kind=rprec) :: dm1,dm2,dm3,dm4,dm5
      real(kind=rprec) :: f1max, eig_max, mu0,
     1  f3amax, f3bmax, f3cmax, ccci, scsi, test,
     2  egl, egu, abstol, ssci, cssi, alf_max, snd_max,
     3  alf_rms, snd_rms, scale_khz, bfavg, va2, f1_avg
      integer :: ispl_opt, ierr_spl, naux, irmv, jrmv
      integer :: ier, npes, mype, numrads
      integer :: ni, nj, mi, mj, ieq, meq, neq, ii
      integer :: j_max_index, irdc
      integer :: m_emax, n_emax, isym_opt, in, im, jcnt
      integer :: m_emax_phi, n_emax_phi, m_emax_psi, n_emax_psi,
     1           mx, kx, incx, incy
      real(kind=rprec) :: eig_max_phi, eig_max_psi, alfx, betx, ddot
      integer :: ldvl, ldvr, lwork, info, nn, m_red,
     1   mn_col2, mn_col2_rdcd
      integer, allocatable :: ifail(:), iwork(:), itht(:), izet(:)
      character*20 outfile, gamfile
      character*3 procnum
      character*1 jobz
      character*10 date, time, zone
      character*10 arg1, arg2
      integer values(8)
      logical :: slow_sound, cyl
      namelist /plasma_input/ ion_to_proton_mass, ion_density_0,
     >    ion_profile, jdqz_data, egnout_form, gamma, nion,
     >    aion, bion, cion, temp_elec_0, slow_sound, telec,
     >    te_profile
c  Defaults
      temp_elec_0 = 636.57d0; slow_sound = .true.; cyl = .false.
      te_profile = 1; ion_profile = 3; gamma = 1.35d0
      call getarg(1,arg1)
      call getarg(2,arg2)
      irads = 41; ir_fine_scl = 128  !default values
      read(arg1,'(i3)') irads
      read(arg2,'(i4)') ir_fine_scl
      irad3 = 3*irads
      write(*,'("irads = ",i4," irads3 = ",i4," ir_fine_scl = ",i5)')
     >       irads, irad3, ir_fine_scl
          
      allocate(bavg(ir_fine_scl),cs2(ir_fine_scl),
     >  te(ir_fine_scl), rho_ion(ir_fine_scl),
     >  ion_density(ir_fine_scl),iota_r(ir_fine_scl),
     >  jtor_lrg(ir_fine_scl),jpol_lrg(ir_fine_scl),
     >  phipc_lrg(ir_fine_scl), stat=istat)
      allocate(nsurf(irads),iotac(irads),
     >  phipc(irads),rho(irads),sp_fit(irads),
     >  jpol(irads), jtor(irads), bav(irads),
     >  yp(irads), temp(irads), ypi(irads), tempi(irads),
     >  ypjt(irads), ypjp(irads), ypphp(irads), stat=istat)
      
      telec(1) = 1.
      telec(2:10) = 0.
      open(unit=4,file="plasma.dat",status="old")
      read(4,plasma_input)
      close(unit=4)
      mass_ion = mass_proton*ion_to_proton_mass
#if defined (PARALLEL)      
C Get NPES and MYPE.  Requires initialization of MPI.
      call MPI_INIT(ier)
      if (ier .ne. 0) then
        write(*,'("mpi_init returned ier =")') ier
        stop
        endif
      call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ier)
      if (ier .ne. 0) then
        write(*,'("mpi_comm_size returned ier =")') ier
        stop
        endif
      call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ier)
      if (ier .ne. 0) then
        write(*,'("mpi_comm_rank returned ier =")') ier
        stop
        endif
#endif
#if defined (SERIAL)      
       mype = 0
       npes = 1
#endif
c
       jobz = 'V'

#if defined (PARALLEL)      
       write(procnum,'(i3)') mype
       write(*,*) procnum
       outfile = "alfven_spec" // trim(adjustl(procnum))
       gamfile = "gam_spec" // trim(adjustl(procnum))
#endif
#if defined (SERIAL)
       outfile = "alfven_spec"
       gamfile = "gam_spec"      
#endif

c
       if(mype .eq. 0) then
        write(*,*) mype, npes
        write(*,*) procnum
        write(*,*) outfile, gamfile
        call date_and_time(date, time, zone, values)
        write(*,*) time
       endif
       if(ipos_def_sym) isym_opt = 1
       if(.NOT.ipos_def_sym) isym_opt = 0
c
c    Generate Fourier arrays
c
      call readin
      call trig_array
      call convolution_array
      irdc = 0
      mn_col2_rdcd = 2*mn_col
      if(irdc .eq. 1) mn_col2_rdcd = 2*mn_col - 1
      do i=1,mn_col
       if(mype .eq. 0) then
         write(*,'(i2,2(2x,i2),2(2x,f6.3))') i, im_col(i),
     >    in_col(i), rm_col(i), rn_col(i)
       end if
       if(im_col(i) .eq. 0 .and. in_col(i) .eq. 0) irdc = i
      end do
      irdc = 0       !initially used for Lagrangian version
      if(mype .eq. 0) write(*,*) irdc
      if(irdc .gt. 0) then
       allocate(temp_amat(mn_col2_rdcd, mn_col2_rdcd), stat=istat)
       allocate(temp_bmat(mn_col2_rdcd, mn_col2_rdcd), stat=istat)
       allocate(im_col_aug(mn_col2_rdcd), stat=istat)
       allocate(in_col_aug(mn_col2_rdcd), stat=istat)
       allocate(eta_sq(mn_col), stat=istat)
       allocate(zeta_sq(mn_col2_rdcd-mn_col), stat=istat)
       do i=1,mn_col
        im_col_aug(i) = im_col(i)
        in_col_aug(i) = in_col(i)
       end do
       do i=mn_col+1,mn_col+irdc-1
        im_col_aug(i) = im_col(i-mn_col)
        in_col_aug(i) = in_col(i-mn_col)
       end do
       do i=mn_col+irdc+1,2*mn_col
        im_col_aug(i-1) = im_col(i-mn_col)
        in_col_aug(i-1) = in_col(i-mn_col)
       end do
      else if(irdc .eq. 0) then
       allocate(temp_amat(2*mn_col, 2*mn_col), stat=istat)
       allocate(temp_bmat(2*mn_col, 2*mn_col), stat=istat)
       allocate(im_col_aug(2*mn_col), stat=istat)
       allocate(in_col_aug(2*mn_col), stat=istat)
       allocate(eta_sq(mn_col), stat=istat)
       allocate(zeta_sq(mn_col), stat=istat)
       do i=1,mn_col
        im_col_aug(i) = im_col(i)
	im_col_aug(i+mn_col) = im_col(i)
        in_col_aug(i) = in_col(i)
	in_col_aug(i+mn_col) = in_col(i)
       end do
      end if

c   check augmented mode number arrays:      
c      do i=1,mn_col
c       write(*,*) im_col(i), in_col(i)
c      end do
c      write(*,'(///)')
c      do i=1,mn_col2_rdcd
c       write(*,*) im_col_aug(i), in_col_aug(i)
c      end do
      
      num_eq = mnmx
      naux = 20*mn_col
      mu0 = 2.d-7*twopi
      scale_khz = (1.d+3*twopi)**2
      jcnt = 0
c      write(*,*) mype, npes, mu0
c      write(*,*) procnum
c      write(*,*) outfile
      allocate(aux(naux), stat=istat)
      allocate(alpha(2*mn_col), stat=istat)
      allocate(beta(2*mn_col), stat=istat)
      allocate(eig_vect(2*mn_col), stat=istat)
      allocate(z_r(2*mn_col,2*mn_col), stat=istat)
      allocate(omega(2*mn_col), stat=istat)
      allocate(amat(2*mn_col,2*mn_col), stat=istat)
      allocate(bmat(2*mn_col,2*mn_col), stat=istat)
      allocate(ifail(2*mn_col), stat=istat)
      allocate(om_keep(2*mn_col), stat=istat)
      allocate(vr_keep(2*mn_col,2*mn_col), stat=istat)
      ldvl = 2*mn_col; ldvr = 2*mn_col; lwork = 40*mn_col
      allocate(work(lwork), stat=istat)
      allocate(iwork(lwork), stat=istat)
      allocate(alfr(2*mn_col), stat=istat)
      allocate(alfi(2*mn_col), stat=istat)
      allocate(vl(2*mn_col,2*mn_col), stat=istat)
      allocate(vr(2*mn_col,2*mn_col), stat=istat)
      allocate(b11(mn_col,mn_col), stat=istat)
      allocate(b22(mn_col,mn_col), stat=istat)
      allocate(phi_temp0(mn_col), stat=istat)
      allocate(phi_temp1(mn_col), stat=istat)
      allocate(psi_temp0(mn_col), stat=istat)
      allocate(psi_temp1(mn_col), stat=istat)
      allocate(phi_norm(2*mn_col), stat=istat)
      allocate(psi_norm(2*mn_col), stat=istat)
      if(mype .eq. 0) then
       write(*,*) izt, ith, irads, mnmx, nznt, mn_col
      end if
      allocate(bfield(izt,ith,irads), stat=istat)
      allocate(gsssup(izt,ith,irads), stat=istat)
      allocate(rjacob(izt,ith,irads), stat=istat)
      allocate(bsupth(izt,ith,irads), stat=istat)
      allocate(bsupzt(izt,ith,irads), stat=istat)
      allocate(bfield_lrg(izt,ith,ir_fine_scl), stat=istat)
      allocate(gsssup_lrg(izt,ith,ir_fine_scl), stat=istat)
      allocate(bsupth_lrg(izt,ith,ir_fine_scl), stat=istat)
      allocate(bsupzt_lrg(izt,ith,ir_fine_scl), stat=istat)
      allocate(rjacob_lrg(izt,ith,ir_fine_scl), stat=istat)
      allocate(theta_tae(ith), stat=istat)
      allocate(zeta_tae(izt), stat=istat)
      allocate(f1_nm(mnmx), stat=istat)
      allocate(f3a_nm(mnmx), stat=istat)
      allocate(f3b_nm(mnmx), stat=istat)
      allocate(f3c_nm(mnmx), stat=istat)
      allocate(g1a_nm(mnmx), stat=istat)
      allocate(g1b_nm(mnmx), stat=istat)
      allocate(g1c_nm(mnmx), stat=istat)
      allocate(g2_nm(mnmx), stat=istat)
      allocate(g3_nm(mnmx), stat=istat)
      allocate(g4_nm(mnmx), stat=istat)
      allocate(g5_nm(mnmx), stat=istat)
      allocate(c1_nm(mnmx), stat=istat)
      allocate(c2b_nm(mnmx), stat=istat)
      allocate(bsupth_tmp(nznt), stat=istat)
      allocate(bsupzt_tmp(nznt), stat=istat)
      allocate(itht(nznt), stat=istat)
      allocate(izet(nznt), stat=istat)
      allocate(f1(izt,ith), stat=istat)
      allocate(f3a(izt,ith), stat=istat)
      allocate(f3b(izt,ith), stat=istat)
      allocate(f3c(izt,ith), stat=istat)
      allocate(g(izt,ith), stat=istat)
      allocate(g1a(izt,ith), stat=istat)
      allocate(g1b(izt,ith), stat=istat)
      allocate(g1c(izt,ith), stat=istat)
      allocate(g2(izt,ith), stat=istat)
      allocate(g3(izt,ith), stat=istat)
      allocate(g4(izt,ith), stat=istat)
      allocate(g5(izt,ith), stat=istat)
      allocate(c1(izt,ith), stat=istat)
      allocate(c2b(izt,ith), stat=istat)
      allocate(dB_dtht(izt,ith), stat=istat)
      allocate(dB_dzet(izt,ith), stat=istat)
      allocate(dlg_dtht(izt,ith), stat=istat)
      allocate(dlg_dzet(izt,ith), stat=istat)
      allocate(dsg_dtht(izt,ith), stat=istat)
      allocate(dsg_dzet(izt,ith), stat=istat)
      allocate(gss0(izt), stat=istat)
      allocate(gss1(izt), stat=istat)
      allocate(gss2(izt), stat=istat)
      allocate(gss3(izt), stat=istat)
      allocate(lg0(izt), stat=istat)
      allocate(lg1(izt), stat=istat)
      allocate(lg2(izt), stat=istat)
      allocate(lg3(izt), stat=istat)
      allocate(b0(izt), stat=istat)
      allocate(b1(izt), stat=istat)
      allocate(b2(izt), stat=istat)
      allocate(b3(izt), stat=istat)
c
c     files for output
c
c      write(*,*) trim(adjustl(outfile))
      
      open(unit=21,file=outfile,status="unknown")
      open(unit=23,file=gamfile,status="unknown")
      if(mype .eq. 0) then
       open(unit=7,file="data_post",status="unknown")
       open(unit=9,file="ion_profile",status="unknown")
       open(unit=31,file="alfven_profiles",status="unknown")
       open(unit=11,file="gss_G.dat",status="unknown")
      end if
c
c    Boozer coordinates input - new ae-mode-structure input
c
      open(unit=20,file="tae_data_boozer",status="old")
      do ir = 1,irads
       read(20,'(1x,i3,4(2x,e15.7))') nn,
     1   iotac(ir),phipc(ir),jtor(ir),jpol(ir)
       nsurf(ir) = dble(nn)
       iotac(ir) = iotac(ir) + delta_iota
       write(0,*) ir,(1./iotac(ir))
c       write(31,'(1x,i3,4(2x,e15.7))') nn,
c     1   iotac(ir),phipc(ir),jtor(ir),jpol(ir)
       do i = 1,izt
        do j = 1,ith
       read(20,'(1x,4(e24.12,2x),e24.12)') theta_tae(j),
     1   zeta_tae(i),
     1   bfield(i,j,ir), gsssup(i,j,ir),
     2   dm1
       bav(ir)  = bav(ir) + bfield(i,j,ir)/(dble(izt*ith))

       read(20,'(1x,4(e24.12,2x),e24.12)') dm2,
     1   dm3, dm4, dm5,
     2   rjacob(i,j,ir)
c      Remove phipc factor since rjacob always occurs here in conjunction with grad-parallel
c       or as an overall multiplier prior to integration
       rjacob(i,j,ir) = rjacob(i,j,ir)/phipc(ir)
       end do
       end do
       end do
       bfavg = sum(bav)/(dble(irads))
c       do i=1,ith
c        write(*,*) i, ith, theta_tae(i)
c       end do
c       stop
c
c    Record equilibrium and eigenfunction Fourier mode list
c
      if(mype .eq. 0) then
       open(unit=22,file="modes",status="unknown")
       write(22,'("Equilibrium modes:",/)')
       write(22,'("meq    neq")')
       do i=1,mnmx
        meq = rm(i)
        neq = rn(i)
        write(22,'(i3,3x,i3)') meq, neq
       end do
       write(22,'(///,"Eigenvector modes:",/)')
       write(22,'("m    n")')
       do i=1,mn_col
        write(22,'(i3,3x,i3)') im_col(i), in_col(i)
       end do
       close(unit=22)
      endif
c
c   Make spline fits and fill in fine_radius_scale arrays
c
      do ir=1,irads
        rho(ir) = nsurf(ir)/nsurf(irads)
c	 write(*,*) rho(ir)
      end do
      
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,iotac,sp1,sp2,ispl_opt,ypi,
     >   tempi,sigma_spl,ierr_spl)     
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,jtor,sp1,sp2,ispl_opt,ypjt,
     >   tempi,sigma_spl,ierr_spl)
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,jpol,sp1,sp2,ispl_opt,ypjp,
     >   tempi,sigma_spl,ierr_spl)
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,phipc,sp1,sp2,ispl_opt,ypphp,
     >   tempi,sigma_spl,ierr_spl)
     
      do irr = 1,ir_fine_scl
       r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
        ra = sqrt(r_pt); ra2 = ra**2; ra3 = ra**3; ra4 = ra**4
	ra5 = ra**5; ra6 = ra**6
        iota_r(irr) = curv2(r_pt,irads,rho,iotac,ypi,sigma_spl)
        jtor_lrg(irr) = curv2(r_pt,irads,rho,jtor,ypjt,sigma_spl)
        jpol_lrg(irr) = curv2(r_pt,irads,rho,jpol,ypjp,sigma_spl)
        phipc_lrg(irr) = curv2(r_pt,irads,rho,phipc,ypphp,sigma_spl)
c        te(irr) = temp_elec_0*ev_to_jou*(1.0821 - 2.0243*ra
c     >            + 0.95204*ra2 + 0.54071*ra3 - 0.50144*ra4)    !HSX
c        te(irr) = temp_elec_0*ev_to_jou*(1.0 - 0.98*ra2)         !LHD RSAE
c	te(irr) = temp_elec_0*ev_to_jou                          !MST
	if(te_profile .eq. 1) then
	te(irr) = (telec(1) + r_pt*telec(2) + telec(3)*(r_pt**2)
     1   + telec(4)*(r_pt**3) + telec(5)*(r_pt**4) + telec(6)*(r_pt**5)
     2   + telec(7)*(r_pt**6) + telec(8)*(r_pt**7)
     3   + telec(9)*(r_pt**8))*ev_to_jou*temp_elec_0
	else if(te_profile .eq. 2) then
	  te(irr) = (telec(1) + ra*telec(2) + telec(3)*ra2
     1     + telec(4)*ra3 + telec(5)*ra4 + telec(6)*ra5
     2     + telec(7)*ra6 + telec(8)*ra6*ra + telec(9)*ra6*ra2
     3     + telec(10)*ra6*ra3)*ev_to_jou*temp_elec_0
        end if
     
     	if(ion_profile .eq. 0) then
	  ion_density(irr) = (iota_r(irr)/iota_r(1))**2   !profile that lines up gaps
	else if(ion_profile .eq. 1) then
	  ion_density(irr) = nion(1) + r_pt*nion(2) + nion(3)*(r_pt**2)
     1     + nion(4)*(r_pt**3) + nion(5)*(r_pt**4) + nion(6)*(r_pt**5)
     2     + nion(7)*(r_pt**6) + nion(8)*(r_pt**7) + nion(9)*(r_pt**8)
     3     + nion(10)*(r_pt**9)
        else if(ion_profile .eq. 2) then
	  ion_density(irr) = 1.d0
	else if(ion_profile .eq. 3) then
	  ion_density(irr) = ((1. - (r_pt**bion))**cion + aion)
     1                       /(1. + aion)
	else if(ion_profile .eq. 4) then
	  ion_density(irr) = nion(1) + ra*nion(2) + nion(3)*ra2
     1     + nion(4)*ra3 + nion(5)*ra4 + nion(6)*ra5
     2     + nion(7)*ra6 + nion(8)*ra6*ra + nion(9)*ra6*ra2
     3     + nion(10)*ra6*ra3
	end if

c       For standard HSX_QHS case, Te0 = 636.57eV; ne0 = 1.5551e+18
 
       if(te(irr) .le. 0.) te(irr) = temp_elec_0*ev_to_jou*0.05
       if(mype .eq. 0) then
        write(31,'(e14.5,3(2x,e14.5))') ra,
     >      ion_density(irr)*ion_density_0, te(irr)/ev_to_jou,
     >      iota_r(irr)
       end if
c
        rho_ion(irr) = ion_density(irr)*mass_ion*ion_density_0
	cs2(irr) = gamma*zion*te(irr)/mass_ion
      end do     !irr = 1,ir_fine_scl
      if(mype .eq. 0) close(unit=31)
      if(ierr_spl .ne. 0) write(*,'("spline error 1",i3)') ierr_spl
      do i=1,izt
       do j=1,ith
c
c       Bfield array
c
        do ir=1,irads
         sp_fit(ir) = bfield(i,j,ir)
        end do      !ir=1,irads
        ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
        call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     >        temp,sigma_spl,ierr_spl)
      if(ierr_spl .ne. 0) write(*,'("spline error 2",i3)') ierr_spl
      do irr = 1,ir_fine_scl
       r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
       bfield_lrg(i,j,irr) = curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
       if(irr .eq. ir_fine_scl/2) then
        if(j .eq. 1) b0(i) = bfield_lrg(i,j,irr)
        if(j .eq. 1+(ith-1)/4) b1(i) = bfield_lrg(i,j,irr)
        if(j .eq. 2+(ith-1)/2) b2(i) = bfield_lrg(i,j,irr)
        if(j .eq. 2+3*(ith-1)/4) b3(i) = bfield_lrg(i,j,irr)
       end if
      end do     !irr = 1,ir_fine_scl
c
c
c       Jacobian array
c
        do ir=1,irads
         sp_fit(ir) = rjacob(i,j,ir)
        end do      !ir=1,irads
        ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
        call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     >        temp,sigma_spl,ierr_spl)
      if(ierr_spl .ne. 0) write(*,'("spline error 2",i3)') ierr_spl
      do irr = 1,ir_fine_scl
       r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
       rjacob_lrg(i,j,irr) = curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
      end do     !irr = 1,ir_fine_scl
c
c       Gsssup array
c
        do ir=1,irads
         sp_fit(ir) = gsssup(i,j,ir)
        end do
        ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
        call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     >        temp,sigma_spl,ierr_spl)
      if(ierr_spl .ne. 0) write(*,'("spline error 3",i3)') ierr_spl
      do irr = 1,ir_fine_scl
       r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
       gsssup_lrg(i,j,irr) = curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
c    Save diagnostic info.
       if(irr .eq. ir_fine_scl/2) then
        if(j .eq. 1) gss0(i) = gsssup_lrg(i,j,irr)
        if(j .eq. 1+(ith-1)/4) gss1(i) = gsssup_lrg(i,j,irr)
        if(j .eq. 2+(ith-1)/2) gss2(i) = gsssup_lrg(i,j,irr)
        if(j .eq. 2+3*(ith-1)/4) gss3(i) = gsssup_lrg(i,j,irr)
       end if
      end do     !irr = 1,ir_fine_scl
      end do     !j=1,ith
      end do     !i=1,izt
c
c
c     Fine_scale arrays and spline fits finished
c
       numrads = ir_fine_scl/npes
c       write(*,*) numrads
       do ir = (mype*numrads+1),(mype*numrads+numrads)
c      do ir=1,ir_fine_scl
        r_pt = rho(1)+real(ir-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
c
c    Form d|B|/dtheta and d|B|/dzeta
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = bfield_lrg(i,j,ir)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      call dbydth
      sin_type = 1; cos_type = 0
      call toReal
      do lg=1,nznt
       dB_dtht(izet(lg),itht(lg)) = f(lg)
      end do
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = bfield_lrg(i,j,ir)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      call dbydzt
      sin_type = 1; cos_type = 0
      call toReal
      do lg=1,nznt
       dB_dzet(izet(lg),itht(lg)) = f(lg)
      end do      
c
c    Form coefficients involving G in real space
c
      do i=1,izt
       do j=1,ith
       
        g(i,j) = (jtor_lrg(ir)*dB_dzet(i,j) + jpol_lrg(ir)*dB_dtht(i,j))   !latest test change
     1      /(phipc_lrg(ir)*rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
c     1      *phipc_lrg(ir)/(rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
        g1c(i,j) = cs2(ir)*g(i,j)/(bfield_lrg(i,j,ir)**2)
	g1b(i,j) = iota_r(ir)*g1c(i,j)
	g1a(i,j) = (iota_r(ir)**2)*g1c(i,j)
	g2(i,j) = mu0*rho_ion(ir)*cs2(ir)*rjacob_lrg(i,j,ir)*(g(i,j)**2)
	g3(i,j) = mu0*rho_ion(ir)*cs2(ir)*g(i,j)
	c1(i,j) = rjacob_lrg(i,j,ir)
	c2b(i,j) = cs2(ir)/(rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
c    Save diagnostic info.
       if(ir .eq. ir_fine_scl/2) then
        if(j .eq. 1) lg0(i) = g(i,j)
        if(j .eq. 1+(ith-1)/4) lg1(i) = g(i,j)
        if(j .eq. 2+(ith-1)/2) lg2(i) = g(i,j)
        if(j .eq. 2+3*(ith-1)/4) lg3(i) = g(i,j)
       end if
	
       end do
      end do
c
c      go to 109      !Eventually remove following uneeded section
c
c     Form grad_prl of G and gss:
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = g(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      call dbydth
      sin_type = 0; cos_type = 1
      call toReal
      do lg=1,nznt
       dlg_dtht(izet(lg),itht(lg)) = f(lg)
      end do
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = g(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      call dbydzt
      sin_type = 0; cos_type = 1
      call toReal
      do lg=1,nznt
       dlg_dzet(izet(lg),itht(lg)) = f(lg)
      end do      
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = gsssup_lrg(i,j,ir)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      call dbydth
      sin_type = 1; cos_type = 0
      call toReal
      do lg=1,nznt
       dsg_dtht(izet(lg),itht(lg)) = f(lg)
      end do
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
	itht(lg) = j; izet(lg) = i
        f(lg) = gsssup_lrg(i,j,ir)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      call dbydzt
      sin_type = 1; cos_type = 0
      call toReal
      do lg=1,nznt
       dsg_dzet(izet(lg),itht(lg)) = f(lg)
      end do      

c
c    Form coefficients involving grad_prl of g, G in real space
c
      do i=1,izt
       do j=1,ith
       
        g4(i,j) = (dsg_dzet(i,j) + iota_r(ir)*dsg_dtht(i,j))
     1            /(rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
        g5(i,j) = cs2(ir)*(dlg_dzet(i,j) + iota_r(ir)*dlg_dtht(i,j))
     1            /((bfield_lrg(i,j,ir)**2))
	
       end do
      end do
c
c
 109  continue

c
c     Generate Fourier spectra of above gridded coefficients
c
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g1a(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      g1a_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g1b(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      g1b_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g1c(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      g1c_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g2(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      g2_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g3(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      g3_nm(:) = fnm(:)
c
c      go to 110      !Eventually remove following uneeded section
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g4(i,j)
       end do
      end do
      sin_type = 1; cos_type = 0
      call toFourier
      g4_nm(:) = fnm(:)
c
 110  continue
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = g5(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      g5_nm(:) = fnm(:)
      
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = c1(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      c1_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = c2b(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      c2b_nm(:) = fnm(:)
c
c    Form coefficients involving gss
c
      f1_avg = 0.
      do i=1,izt
       do j=1,ith
         f1(i,j) = mu0*rho_ion(ir)*gsssup_lrg(i,j,ir)
     1            *rjacob_lrg(i,j,ir)/(bfield_lrg(i,j,ir)**2)
         f1_avg = f1_avg + f1(i,j)/real(izt*ith)
         f3c(i,j) = gsssup_lrg(i,j,ir)/
     1              (rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
         f3b(i,j) = iota_r(ir)*f3c(i,j)
         f3a(i,j) = iota_r(ir)*f3b(i,j)
       end do
      end do
      if(cyl) then
       f1(:,:) = f1_avg
      endif
c
c
c     Generate Fourier spectra of above gridded coefficients
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = f1(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f1_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = f3a(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3a_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = f3b(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3b_nm(:) = fnm(:)
c
      lg = 0
      do i=1,izt
       do j=1,ith
        lg = lg + 1
        f(lg) = f3c(i,j)
       end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3c_nm(:) = fnm(:)
c
c     Build A and B matrices
c
      amat(:,:) = 0.d0; bmat(:,:) = 0.d0
c    (1,1) Block
c
      do i=1,mn_col
        do j=1,mn_col
         ni = in_col(i)
         nj = in_col(j)
         mi = im_col(i)
         mj = im_col(j)
          do ieq = 1,mnmx
            meq = rm(ieq)
            neq = rn(ieq)
            call ccc_convolve(ccci,mi,ni,mj,nj,meq,neq)
            call scs_convolve(scsi,mi,ni,mj,nj,meq,neq)
	    
            bmat(i,j) = bmat(i,j) - ccci*f1_nm(ieq)*scale_khz
              amat(i,j) = amat(i,j)
     1         - (scsi*(f3a_nm(ieq)*rm_col(i)*rm_col(j)
     2               - f3b_nm(ieq)*rm_col(j)*rn_col(i)
     3               - f3b_nm(ieq)*rn_col(j)*rm_col(i)
     4               + f3c_nm(ieq)*rn_col(j)*rn_col(i)))
     5               - g2_nm(ieq)*ccci
     
           end do              !ieq = 1,mnmx
        end do                 !j=1,mn_col
       end do                  !i=1,mn_col
c
c    (2,1) Block
c
       do i=mn_col+1,2*mn_col
        do j=1,mn_col
	 ii = i-mn_col
         ni = in_col(ii)
         nj = in_col(j)
         mi = im_col(ii)
         mj = im_col(j)
          do ieq = 1,mnmx
            meq = rm(ieq)
            neq = rn(ieq)
            call ccc_convolve(ccci,mi,ni,mj,nj,meq,neq)
            call scs_convolve(cssi,mj,nj,meq,neq,mi,ni)
	    
	      bmat(i,j) = 0.d0
              amat(i,j) = amat(i,j)
     1                    - ccci*g5_nm(ieq)
     3                    - cssi*(g1c_nm(ieq)*rn_col(j)
     4                            - g1b_nm(ieq)*rm_col(j))

           end do              !ieq = 1,mnmx
        end do                 !j=1,mn_col
       end do                  !i=1,mn_col
c
c    (1,2) Block
c
       do i=1,mn_col
        do j=mn_col+1,2*mn_col
	 jj = j-mn_col
         ni = in_col(i)
         nj = in_col(jj)
         mi = im_col(i)
         mj = im_col(jj)
          do ieq = 1,mnmx
            meq = rm(ieq)
            neq = rn(ieq)
            call scs_convolve(ssci,mi,ni,meq,neq,mj,nj)
            call scs_convolve(cssi,mj,nj,meq,neq,mi,ni)
	    
	      bmat(i,j) = 0.d0
              amat(i,j) = amat(i,j)
     1                    - ssci*g3_nm(ieq)*(rn_col(i)
     2                            - iota_r(ir)*rm_col(i))
     3                    + cssi*g3_nm(ieq)*(rn_col(jj)
     4                            - iota_r(ir)*rm_col(jj))
     
           end do              !ieq = 1,mnmx
        end do                 !j=1,mn_col
       end do                  !i=1,mn_col
c
c    (2,2) Block
c
       do i=mn_col+1,2*mn_col
        do j=mn_col+1,2*mn_col
	 ii = i-mn_col
	 jj = j-mn_col
         ni = in_col(ii)
         nj = in_col(jj)
         mi = im_col(ii)
         mj = im_col(jj)
          do ieq = 1,mnmx
            meq = rm(ieq)
            neq = rn(ieq)
            call ccc_convolve(ccci,mi,ni,mj,nj,meq,neq)
	    
	     bmat(i,j) = bmat(i,j) - c1_nm(ieq)*ccci*scale_khz
	     if(slow_sound) then
              amat(i,j) = 0.d0
	     else
              amat(i,j) = amat(i,j)
     1        - (ccci*(c2b_nm(ieq)*rm_col(jj)*rm_col(jj)*(iota_r(ir)**2)
     2            - 2.d0*c2b_nm(ieq)*rm_col(jj)*rn_col(jj)*iota_r(ir)
     3            + c2b_nm(ieq)*rn_col(jj)*rn_col(jj)))
             end if
     
           end do              !ieq = 1,mnmx
        end do                 !j=1,mn_col
       end do                  !i=1,mn_col
	 
c
       mn_col2 = 2*mn_col
c
c
c
c     If irdc > 0 (i.e., m=0,n=0 mode is present) remove the psi(0,0) components
c      since the psi equation does not have a constant solution. These components are at
c      row = mn_col + irdc. Also, remove column at col = mn_col + irdc.
c
      if(mype .eq. 0) write(*,*) irdc, mn_col, mn_col2
c
      irdc = 0
      if(irdc .eq. 0) then
       do i=1,mn_col2
        do j=1,mn_col2
	 temp_amat(i,j) = amat(i,j)
	 temp_bmat(i,j) = bmat(i,j)
	end do
       end do       
      else if(irdc .gt. 0) then
       irmv = mn_col+irdc; jrmv=mn_col+irdc
       do i=1,mn_col2
        do j=1,mn_col2
	 if(i .eq. irmv) cycle
	 if(j .eq. jrmv) cycle
	 if(i .lt. irmv .and. j .lt. jrmv) then
	  temp_amat(i,j) = amat(i,j)
	  temp_bmat(i,j) = bmat(i,j)
	 else if(i .gt. irmv .and. j .lt. jrmv) then
	  temp_amat(i-1,j) = amat(i,j)
	  temp_bmat(i-1,j) = bmat(i,j)
	 else if(i .lt. irmv .and. j .gt. jrmv) then
	  temp_amat(i,j-1) = amat(i,j)
	  temp_bmat(i,j-1) = bmat(i,j)
	 else if(i .gt. irmv .and. j .gt. jrmv) then
	  temp_amat(i-1,j-1) = amat(i,j)
	  temp_bmat(i-1,j-1) = bmat(i,j)
	 endif
	 
	end do
       end do       
       mn_col2 = mn_col2_rdcd
       end if
c
c
c   Test output of matrices. Only useful for smaller size problems.
c
       if(mat_out) then
       if(ir .eq. ir_fine_scl/2) then
        open(unit=32,file="amat.dat",status="unknown")
        open(unit=33,file="bmat.dat",status="unknown")
	do i=1,mn_col2
	 write(32,44) (temp_amat(i,j),j=1,mn_col2)
	 write(33,44) (temp_bmat(i,j),j=1,mn_col2)
c	 do j=1,2*mn_col
c	  write(32,*) i,j,amat(i,j)
c	  write(33,*) i,j,bmat(i,j)
c	 end do
	end do
 44     format(e12.4,8(2x,e12.4),///)   !Need to change format depending on matrix size
	close(unit=32)
	close(unit=33)
       end if  !ir .eq. ir_fine_scl/2
       end if  !mat_out

c
c     Call matrix eigenvalue solver. This solves the system A*x = lambda*B*x
c
c      where:         [ A11   A12 ]           [ B11   B12 ]        [ Phi ]
c                A =  [           ]      B =  [ 	  ]    x = [     ]
c                     [ A21   A22 ]	      [ B21   B22 ]        [ Psi ]
c

       egl = 1.d-2; egu = 0.6d0; abstol = 1.d-8
       il = 0; iu = 0

       if(.NOT.ipos_def_sym) then
        call dggev('N','V',mn_col2,temp_amat,mn_col2,temp_bmat,mn_col2,
     1     alfr,alfi,beta,vl,ldvl,vr,ldvr,work,lwork,info)
       else if(ipos_def_sym .and. .not.subset_eq) then
        call dsygv(iopt,jobz,'L',mn_col2,temp_amat,mn_col2,temp_bmat,
     1     mn_col2,omega, work, lwork, info)
       
       else if(ipos_def_sym .and. subset_eq) then
        call dsygvx(iopt,'V','V','L',mn_col2,temp_amat,mn_col2,
     1     temp_bmat,mn_col2,egl, egu, il, iu, abstol, m_red,
     >     omega, vr, mn_col2,work, lwork, iwork, ifail, info)
       
       endif
c       if(info .ne. 0) write(*,'("info = ",i8)') info
c       write(*,*) ir, info, phipc_lrg(ir)
c
      do i=1,mn_col2
       if(iopt .eq. 1) then
       do j=1,mn_col2
        if(ipos_def_sym .and. jobz .eq. 'V' .and. .not.subset_eq)
     1          eig_vect(j) = abs(temp_amat(j,i))
        if(.NOT.ipos_def_sym) eig_vect(j) = abs(vr(j,i))
c	if(ir .eq. 5 .and. i .eq. mn_col/2)
c     >   write(*,'(e15.7,2x,i4,2x,i4)') eig_vect(j),
c     >   im_col(j), in_col(j)

       end do   !j = 1,mn_col2
c
c     Save data from one radial position and for f < 300 kHz
c       
c       if(ir .eq. ir_out .and. beta(i) .gt. 0.) then
c	ftest = sqrt(alfr(i)/beta(i))
c	write(*,'(e12.5,3(2x,e12.5))') ftest, alfr(i), alfi(i), beta(i)
c	  if(ftest .lt. 3.d+2) then
c	   jcnt = jcnt + 1
c	   vr_keep(:,jcnt) = vr(:,i)
c	   om_keep(jcnt) = ftest
c	  endif
c        endif  !if(ir .eq. ir_out)
       
c       eig_max = -1.d+30
c       do j=1,mn_col2
c        if(eig_vect(j) .gt. eig_max) then
c	 eig_max = eig_vect(j)
c	 j_max_index = j
c	end if
c       end do
c        m_emax = im_col_aug(j_max_index)
c        n_emax = in_col_aug(j_max_index)
c
c     Compute norms (phi_transpose*B11*phi and psi_transpose*B22*psi)
c      of the Alfven and acoustic components as a way to separate the
c      continua.
c
       do kx = 1,mn_col
        phi_temp0(kx) = vr(kx,i)
	psi_temp0(kx) = vr(kx+mn_col,i)
	do mx = 1,mn_col
	 b11(kx,mx) = bmat(kx,mx)
	 b22(kx,mx) = bmat(kx+mn_col,mx+mn_col)
	end do
       end do
       incx = 1; incy = 1; alfx = 1.; betx = 0.
       call dgemv('N',mn_col,mn_col,alfx,b11,mn_col,phi_temp0,
     >   incx,betx,phi_temp1,incy)
       call dgemv('N',mn_col,mn_col,alfx,b22,mn_col,psi_temp0,
     >   incx,betx,psi_temp1,incy)
       incx = 1; incy = 1
       va2 = bfavg**2/(mu0*rho_ion(ir)/scale_khz)
       phi_norm(i) = va2*ddot(mn_col,phi_temp0,incx,phi_temp1,incy)
       incx = 1; incy = 1
       psi_norm(i) = ddot(mn_col,psi_temp0,incx,psi_temp1,incy)

       eig_max_phi = -1.d+30
       eig_max_psi = -1.d+30
       do j=1,mn_col
        if(eig_vect(j) .gt. eig_max_phi) then
	 eig_max_phi = eig_vect(j)
	 j_max_index = j
	end if
       end do
        m_emax_phi = im_col_aug(j_max_index)
        n_emax_phi = in_col_aug(j_max_index)

       do j=mn_col+1,mn_col2
        if(eig_vect(j) .gt. eig_max_psi) then
	 eig_max_psi = eig_vect(j)
	 j_max_index = j
	end if
       end do
        m_emax_psi = im_col_aug(j_max_index)
        n_emax_psi = in_col_aug(j_max_index)
	
       alf_max = -1.d+30
        do j=1,mn_col
         eta_sq(j) = eig_vect(j)**2
	 if(eig_vect(j) .gt. alf_max) alf_max = eig_vect(j)
	end do
       alf_rms = sqrt(sum(eta_sq))
       snd_max = -1.d+30
        do j=mn_col+1, mn_col2
         zeta_sq(j-mn_col) = eig_vect(j)**2
	 if(eig_vect(j) .gt. snd_max) snd_max = eig_vect(j)
	end do
       snd_rms = sqrt(sum(zeta_sq))
        
       if(ipos_def_sym) then
         write(21,'(2(e15.7,2x),i4,2x,i4)') r_pt,
     1     sqrt(abs(omega(i)))/(1000.*twopi),
     1     m_emax, n_emax
       else if(.NOT.ipos_def_sym) then
        write(21,'(4(e15.7,2x),i4,3(2x,i4),2(2x,e15.7))')
     1   r_pt, alfr(i), alfi(i),beta(i),m_emax_phi, n_emax_phi,
     2   m_emax_psi, n_emax_psi, phi_norm(i), psi_norm(i)
c        if(alf_max .gt. 1.e-6*snd_max)then
         write(23,'(e15.7,5(2x,e15.7))') r_pt,
     1        (sqrt(alfr(i)/beta(i)))/(1000.*twopi),
     2        alf_max, snd_max, alf_rms, snd_rms
c	endif
       endif
      else if(iopt .eq. 0) then
       if(ipos_def_sym) then
        write(21,'(e15.7,2x,e15.7)') r_pt,
     1    sqrt(abs(omega(i)))/(1000.*twopi)
       else if(.NOT.ipos_def_sym) then
        if(r_pt .lt. 0.94) then
         write(21,'(e15.7,3(2x,e15.7))') r_pt,alfr(i),
     1     alfi(i),beta(i)
        endif
       endif
      endif
      end do          !do i=1,mn_col2
      if(mype .eq. 0) then
       write(*,82) ir,alfr(mn_col2/10),alfi(mn_col2/10),
     1 beta(mn_col2/10),alfr(mn_col2/2),alfi(mn_col2/2),beta(mn_col2/2)
      end if
 
      end do       !ir=1,ir_fine_scl
c
 82   format(i4,3(3x,e12.4),/,6x,e12.5,2(3x,e12.4))
c      if(ir .eq. ir_out) then
       open(unit=77,file="mode_structures.dat",status="unknown")
       write(77,*) mn_col, mn_col2
       write(77,*) jcnt
       write(*,*) jcnt
       do j = 1,jcnt
        write(*,*) om_keep(j)
        write(77,*) om_keep(j)
	do i=1,mn_col2
          im = im_col_aug(i)
          in = in_col_aug(i)
	 write(77,*) im,in,vr_keep(i,j)
	end do
	write(77,'(///)')
       end do
       close(unit=77)
c      end if    !if(ir .eq. ir_out)
      if(mype .eq. 0) then
       write(11,'(e14.4,3(2x,e14.4))') theta_tae(1),
     >    theta_tae(1+(ith-1)/4), theta_tae(2+(ith-1)/2),
     >    theta_tae(2+3*(ith-1)/4)
       do i=1,izt
        write(11,'(e14.5,12(2x,e14.5))') zeta_tae(i), gss0(i),
     >    gss1(i), gss2(i), gss3(i), lg0(i), lg1(i), lg2(i),
     >    lg3(i), b0(i), b1(i), b2(i), b3(i)
       end do
      end if
c
      close(unit=20)
      close(unit=21)
      if(mype .eq. 0) write(*,
     >  '("modes = ",i5,2x,"no. of radial points = ",i5)')
     >  mn_col,ir_fine_scl
      if(mype .eq. 0) write(7,'(i5,3(2x,i5))') iopt,mn_col2,
     1    ir_fine_scl,isym_opt
      if(mype .eq. 0) close(unit=7)
      call trg_deallocate
      deallocate(alpha, beta, aux)
       if(mype .eq. 0) then
         call date_and_time(date, time, zone, values)
         write(*,*) time
       endif
#if defined (PARALLEL)      
        call MPI_FINALIZE(ier)
#endif
      stop
      end program tae_continua
