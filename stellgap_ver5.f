*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program stellgap.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
c
c  Note: this code uses normalized toroidal flux (variables: rho, r_pt) as
c        an independent variable. If a radius-like variable is desired
c        (i.e., sqrt(tor_flux) or radius based on sqrt(area/pi)), this must be
c        constructed in subsequent codes that use the AE3D eigenfunction data
c
c   Before running, check that the irads variable is set consistently with the
c   number of surfaces in the tae_boozer_data file. irads should be the number
c   of surfaces in the original VMEC run minus 1 or 2 to avoid edge problems.
c
c  Modified 2/23/2013 to properly treat RFPs where there is a reversal surface
c    (q -> 0, iota -> infinity). This is handled through options controlled by
c    the logical variable lrfp. lrfp is read from the first line of tae_data_boozer,
c    which is written by xmetric. For analysis of systems with lrfp = .true. both
c    updated versions of metric_element_create.f and stellagap.f are required.
c    If using an older tae_data_boozer file, a read error will result since the
c    first line will not have the lrfp logical. One can edit tae_data_boozer
c    and add "F" in the first line to resolve this. However, to treat cases with
c    lrfp = .true. it is necessary to have also used the upgraded xmetric. By
c    default lrfp is set to .false.
c
c  Modified 7/6/2012 making ir_fine_scl and irads to be command-line arguments
c   (if not given, default values of 128 and 41 are used)
c
c  Modified 12/2008 so that ion density profile, ion mass and several other
c   parameters are read in through the plasma.dat file rather than being
c   hardwired into the code. (see explanation below for variables contained
c   in plasma.dat file).
c
c  Modified on 9/18/2008 to offer both parallel or serial options, depending on
c   use of the precompiler.
c
c  This is the stellgap code for calculating Alfven contina in stellarators. The
c  associated publication that describes this calculation is:
c
c   D. A. Spong, R. Sanchez, A. Weller, "Shear Alfven Continua in Stellarators,"
c     Phys. of Plasmas, Vol. 10, pg. 3217, August, 2003.
c
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
c     (in file tae_data_vmec) calculated by the xmetric code.
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
c                           + nion(4)*(rho**3) + ... + nion(9)*(rho**8) + nion(10)*(rho**9)
c               (here rho = normalized toroidal flux)
c      for ion_profile = 2 ion_density = ion_density = constant = 1.
c      for ion_profile = 3 ion_density = [1 - aion*(rho**bion)]**cion
c
c   nion(10) = array of polynomial fit coefficients for ion_profile = 1
c   aion, bion, cion = parameters for ion_profile = 3 fit
c   jdqz_data = logical variable that is used in ae3d, but not in stellgap.f
c                (included here only so that same plasma.dat file can be used)
c   egnout_form = "binr" or "asci" - like jdqz_data this variable is used in
c                  ae3d, but not stellgap.f
c
c compile: sh bld_stellgap
c
c  This code is can be run in serial or parallel mode, depending on whether
c   is it run through the precompiler with SERIAL or PARALLEL set.
c
c

C-----------------------------------------------------------------------
C     KIND SPEC MODULE
C-----------------------------------------------------------------------
      module kind_spec
      integer, parameter :: rprec = selected_real_kind(12,100)
      integer, parameter :: iprec = selected_int_kind(8)
      end module kind_spec

C-----------------------------------------------------------------------
C     FOURIER_LIB
C-----------------------------------------------------------------------
      module fourier_lib
      use kind_spec
      implicit none
      integer::ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
      real(kind=rprec), parameter :: parity_gss = 1.,
     1     parity_bsupth = 1., parity_bsupzt = 1.
      integer :: i, j, lg, nu, nl, m, n, mn, istat, sin_type, cos_type
c     Equilibrium coefficient arrays
      real(kind=rprec), allocatable :: thtgrd(:),ztgrd(:),rn(:),rm(:)
      real(kind=rprec), allocatable :: fnm(:),f(:),anm(:)
      real(kind=rprec), allocatable :: cos_ar(:,:),sin_ar(:,:)
      real(kind=rprec), allocatable :: cos_toF(:,:),sin_toF(:,:)
      real(kind=rprec) twopi, arg
      logical, parameter :: ipos_def_sym=.false.
c     The subset_eq flag allows one to just request a bracketed subset of
c     eigenvalues, rather than all eigenvalues. This option has not been 
c     developed beyond the call. Initial tests have not indicated it
c     speeds things up.
      logical, parameter :: subset_eq = .false.
c     Eigenfunction arrays,variables
      integer :: count, mn_col, ith_col, izt_col
      real(kind=rprec), allocatable :: rn_col(:),rm_col(:)
      real(kind=rprec), allocatable :: rn2_col(:),rm2_col(:),rnm_col(:)
      integer, allocatable :: in_col(:), im_col(:)
      integer, allocatable :: nw(:), mwl(:), mwu(:)             !change when modes are added
      contains

c     Process the Fourier.dat file  
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

c     Trigonometric transforms 
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
c     Generate theta, zeta grid
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
      do n = 1,ntors                           !change when modes are added
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
      end subroutine convolution_array
c
c
      subroutine scs_convolve(ans,m1,n1,m2,n2,meq,neq)
      use kind_spec
      implicit none
      real(kind=rprec) :: tht_int1, tht_int2, tht_int3, tht_int4,
     1          zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
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
            if(j*j .eq. k*k) result = two
            return
          else if(j .eq. 0) then
            result = zero
            if(i*i .eq. k*k) result = two
            return
          else if(k .eq. 0) then
            result = zero
            if(i*i .eq. j*j) result = two
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
!-----------------------------------------------------------------------
!   M o d u l e s
!-----------------------------------------------------------------------
      use fourier_lib
      use kind_spec

!-----------------------------------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------------------------------
      implicit none
#if defined (PARALLEL)      
      include 'mpif.h'
#endif
      integer, parameter :: iopt = 1
      real(kind=rprec), parameter :: R0 = 1.0
      real(kind=rprec), parameter :: mass_proton = 1.67d-27

C-----------------------------------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------------------------------
      integer :: ir_fine_scl, irads, irad3
      real(kind=rprec) :: ion_density_0, mass_ion, ion_to_proton_mass,
     1   aion, bion, cion
      real(kind=rprec) :: nion(10)
      integer :: ir, irr, il, iu, ion_profile
      character*4 :: egnout_form
      logical :: jdqz_data
      real(kind=rprec), dimension(:,:,:), allocatable :: bfield,
     1  gsssup,rjacob, bsupth, bsupzt, bfield_lrg,
     2  gsssup_lrg, bsupth_lrg, bsupzt_lrg, rjacob_lrg
      real(kind=rprec), dimension(:), allocatable :: theta_tae,
     1 zeta_tae, f1_nm, f3a_nm, f3b_nm, f3c_nm, bsupth_tmp, bsupzt_tmp
      real(kind=rprec), dimension(:,:), allocatable :: f1, f3a,
     1 f3b, f3c
      real(kind=rprec), dimension(:), allocatable :: bavg, mu0_rho_ion,
     1                  ion_density, iota_r, iota_r_inv 
      real(kind=rprec), allocatable :: beta(:), aux(:), eig_vect(:)
      real(kind=rprec), allocatable :: amat(:,:),
     1           bmat(:,:), z_r(:,:), vr(:,:), vl(:,:)
      real(kind=rprec), allocatable :: bgrad2(:,:,:)
      real(kind=rprec), allocatable :: omega(:), work(:), alfr(:),
     1  alfi(:)
      complex*16, allocatable :: alpha(:)      
      real(kind=rprec), dimension(:), allocatable :: nsurf, iotac,
     1                   phipc, rho, sp_fit, b2, iotac_inv      
      real(kind=rprec) :: sp1, sp2, sigma_spl, r_pt,
     1  iotac_r, curv2, dum1, dum2, dum3, dum4     
      real(kind=rprec), dimension(:), allocatable :: yp, temp,
     1                  ypi, tempi      
      real(kind=rprec) :: amat1, amat2, amat3, amat4, amat5
      real(kind=rprec) :: f1_coef, f2a_coef, f3a_coef,
     1  f3b_coef, f3c_coef
      real(kind=rprec) :: dm1,dm2,dm3,dm4,dm5, bfavg, va
      real(kind=rprec) :: f1max, eig_max, mu0, scale_khz,
     1  f3amax, f3bmax, f3cmax, ccci, scsi, test,
     2  egl, egu, abstol, f1_avg, f3_avg
      real(kind=rprec) :: ra, ra2, ra3, ra4, ra5, ra6
      integer :: ispl_opt, ierr_spl, num_eq, naux
      integer :: ier, npes, mype, numrads
      integer :: ni, nj, mi, mj, ieq, meq, neq, ii, jj
      integer :: j_max_index, ios
      integer :: m_emax, n_emax, isym_opt, in, im
      integer :: ldvl, ldvr, lwork, info, nn, m_red
      integer, allocatable :: ifail(:), iwork(:)
      logical :: cyl, lrfp
      character*20 outfile
      character*3 procnum
      character*1 jobz
      character*10 date, time, zone
      character*10 arg1, arg2
      integer values(8)

C-----------------------------------------------------------------------
C     Fortran Namelists
C-----------------------------------------------------------------------
      namelist /plasma_input/ ion_to_proton_mass, ion_density_0,
     >    ion_profile, jdqz_data, egnout_form, nion,
     >    aion, bion, cion

C-----------------------------------------------------------------------
C     Local runtime settings
C-----------------------------------------------------------------------
      cyl = .false.; lrfp = .false.
      irads = 41; ir_fine_scl = 128  !default values
      jobz = 'V'

C-----------------------------------------------------------------------
C     Get command line input
C-----------------------------------------------------------------------
      call getarg(1,arg1)
      call getarg(2,arg2)
      read(arg1,'(i3)') irads
      read(arg2,'(i4)') ir_fine_scl
      irad3 = 3*irads
      write(*,'("irads = ",i4," irads3 = ",i4," ir_fine_scl = ",i5)')
     >       irads, irad3, ir_fine_scl

C-----------------------------------------------------------------------
C     Allocate quantities
C----------------------------------------------------------------------- 
      allocate(bavg(ir_fine_scl),mu0_rho_ion(ir_fine_scl),
     >  ion_density(ir_fine_scl),iota_r(ir_fine_scl),
     >  iota_r_inv(ir_fine_scl),stat=istat)
      allocate(nsurf(irads),iotac(irads),iotac_inv(irads),
     >  phipc(irads),rho(irads),sp_fit(irads),
     >  b2(irads), stat=istat)
      allocate(yp(irad3),temp(irad3),ypi(irad3),
     >  tempi(irad3), stat=istat) 

C-----------------------------------------------------------------------
C     Read the plasma.dat file
C----------------------------------------------------------------------- 
      open(unit=4,file="plasma.dat",status="old")
      read(4,plasma_input)
      close(unit=4)
      mass_ion = mass_proton*ion_to_proton_mass

C-----------------------------------------------------------------------
C     Create MPI Communicators if needed
C----------------------------------------------------------------------- 
      mype = 0; npes = 1
      outfile = "alfven_spec"  
#if defined (PARALLEL)      
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
       write(procnum,'(i3)') mype
       write(*,*) procnum
       outfile = "alfven_spec" // trim(adjustl(procnum))
#endif

C-----------------------------------------------------------------------
C     Output runtime information
C----------------------------------------------------------------------- 
      if(mype .eq. 0) then
        write(*,*) mype, npes
        write(*,*) procnum
        write(*,*) outfile
        call date_and_time(date, time, zone, values)
        write(*,*) time
      endif
      if(ipos_def_sym) isym_opt = 1
      if(.NOT.ipos_def_sym) isym_opt = 0

C-----------------------------------------------------------------------
C     Generate Fourier arrays
C----------------------------------------------------------------------- 
      call readin
      call trig_array
      call convolution_array
      num_eq = mnmx
      naux = 10*mn_col
      mu0 = 2.d-7*twopi
      scale_khz = (1.d+3*twopi)**2

C-----------------------------------------------------------------------
C     Allocate arrays
C-----------------------------------------------------------------------
      allocate(aux(naux), stat=istat)
      allocate(alpha(mn_col), stat=istat)
      allocate(beta(mn_col), stat=istat)
      allocate(eig_vect(mn_col), stat=istat)
      allocate(z_r(mn_col,mn_col), stat=istat)
      allocate(omega(mn_col), stat=istat)
      allocate(amat(mn_col,mn_col), stat=istat)
      allocate(bmat(mn_col,mn_col), stat=istat)
      allocate(ifail(mn_col), stat=istat)
      ldvl = mn_col; ldvr = mn_col; lwork = 20*mn_col
      allocate(work(lwork), stat=istat)
      allocate(iwork(lwork), stat=istat)
      allocate(alfr(mn_col), stat=istat)
      allocate(alfi(mn_col), stat=istat)
      allocate(vl(mn_col,mn_col), stat=istat)
      allocate(vr(mn_col,mn_col), stat=istat)
      if(mype .eq. 0)
     >    write(*,*) izt, ith, irads, mnmx, nznt, mn_col
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
      allocate(bsupth_tmp(nznt), stat=istat)
      allocate(bsupzt_tmp(nznt), stat=istat)
      allocate(f1(izt,ith), stat=istat)
      allocate(f3a(izt,ith), stat=istat)
      allocate(f3b(izt,ith), stat=istat)
      allocate(f3c(izt,ith), stat=istat)

C-----------------------------------------------------------------------
C     Open files for output
C----------------------------------------------------------------------- 
      open(unit=21,file=trim(adjustl(outfile)),status="unknown")
      open(unit=8,file="coef_arrays",status="unknown")
      if(mype .eq. 0) then
        open(unit=7,file="data_post",status="unknown")
        open(unit=9,file="ion_profile",status="unknown")
      endif

C-----------------------------------------------------------------------
C     Boozer coordinates input - new ae-mode-structure input
C----------------------------------------------------------------------- 
      open(unit=20,file="tae_data_boozer",status="old")
c     LRFP logic is all screwy
c      read(20,'(L)') lrfp    !uncomment and remove next line if lrfp added to tae_data_boozer
      lrfp = .false.
      if(lrfp) write(*,'("Using RFP settings")')
      if(.not.lrfp) write(*,'("Using tokamak/stellarator settings")')
      do ir = 1,irads
        read(20,'(1x,i3,4(2x,e15.7))') nn,iotac(ir),phipc(ir),dum1,dum2
        iotac_inv(ir) = 1.d0/iotac(ir)
        nsurf(ir) = dble(nn)
        b2(ir) = 0.d0
        do i = 1,izt
          do j = 1,ith
            read(20,'(1x,4(e24.12,2x),e24.12)') theta_tae(j), 
     1          zeta_tae(i), bfield(i,j,ir), gsssup(i,j,ir), dm1
            b2(ir)  = b2(ir) + bfield(i,j,ir)/(dble(izt*ith))
            read(20,'(1x,4(e24.12,2x),e24.12)') dm2, dm3, dm4, dm5,
     1          rjacob(i,j,ir)
            rjacob(i,j,ir) = rjacob(i,j,ir)
            if(.not.lrfp) rjacob(i,j,ir) = rjacob(i,j,ir)/phipc(ir)
          end do
        end do
      end do
      bfavg = sum(b2)/(dble(irads))

C-----------------------------------------------------------------------
C     Record equilibrium and eigenfunction Fourier mode list
C-----------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C     Make RHO array
C-----------------------------------------------------------------------
      do ir=1,irads
        rho(ir) = nsurf(ir)/nsurf(irads)
      end do

C-----------------------------------------------------------------------
C     Make spline fits of iota
C-----------------------------------------------------------------------
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,iotac,sp1,sp2,ispl_opt,ypi,
     >   tempi,sigma_spl,ierr_spl)
      if(ierr_spl .ne. 0) write(*,'("spline error iota",i3)') ierr_spl

C-----------------------------------------------------------------------
C     Compute ion density and Alfven velocity
C-----------------------------------------------------------------------
      do irr = 1,ir_fine_scl
        r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
        ra = sqrt(r_pt); ra2 = ra**2; ra3 = ra**3; ra4 = ra**4
        ra5 = ra**5; ra6 = ra**6
        iota_r(irr) = curv2(r_pt,irads,rho,iotac,ypi,sigma_spl)
        if(ion_profile .eq. 0) then
          ion_density(irr) = (iota_r(irr)/iotac(1))**2   !profile that lines up gaps
        else if(ion_profile .eq. 1) then
          ion_density(irr) = nion(1) + r_pt*nion(2) + nion(3)*(r_pt**2)
     1     + nion(4)*(r_pt**3) + nion(5)*(r_pt**4) + nion(6)*(r_pt**5)
     2     + nion(7)*(r_pt**6) + nion(8)*(r_pt**7) + nion(9)*(r_pt**8)
     3     + nion(10)*(r_pt**9)
        else if(ion_profile .eq. 2) then
          ion_density(irr) = 1.d0
        else if(ion_profile .eq. 3) then
          ion_density(irr) = (1. - aion*(r_pt**bion))**cion
        end if
        mu0_rho_ion(irr)=mu0*mass_ion*ion_density_0
     1    *ion_density(irr)*scale_khz
        va = sqrt(bfavg**2/(mu0_rho_ion(irr)/scale_khz))
        if(mype .eq. 0)
     1    write(9,67) r_pt,ion_density_0*ion_density(irr),
     2                iota_r(irr),va
      end do     !irr = 1,ir_fine_scl
 67   format(e12.5,3(3x,e12.5))

C-----------------------------------------------------------------------
C     Make spline fits of 1/iota
C-----------------------------------------------------------------------
      ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
      call curv1(irads,rho,iotac_inv,sp1,sp2,ispl_opt,ypi,
     1   tempi,sigma_spl,ierr_spl)
      if(ierr_spl .ne. 0) 
     1 write(*,'("spline error 1/iota",i3)') ierr_spl
      do irr = 1,ir_fine_scl
        r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
        iota_r_inv(irr) = curv2(r_pt,irads,rho,iotac_inv,ypi,sigma_spl)
      end do     !irr = 1,ir_fine_scl


C-----------------------------------------------------------------------
C     Make spline fits of 3D Arrays
C-----------------------------------------------------------------------
      do i=1,izt
        do j=1,ith

C-----------------------------------------------------------------------
C         Make spline fits of B-Field
C-----------------------------------------------------------------------
          do ir=1,irads
            sp_fit(ir) = bfield(i,j,ir)
          end do      !ir=1,irads
          ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
          call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     1               temp,sigma_spl,ierr_spl)
          if(ierr_spl .ne. 0) write(*,'("spline error B",i3)') ierr_spl
          do irr = 1,ir_fine_scl
            r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1             /real(ir_fine_scl-1)
            bfield_lrg(i,j,irr) = 
     1          curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
          end do     !irr = 1,ir_fine_scl

C-----------------------------------------------------------------------
C         Make spline fits of Jacobian
C-----------------------------------------------------------------------
          do ir=1,irads
            sp_fit(ir) = rjacob(i,j,ir)
          end do      !ir=1,irads
          ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
          call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     1               temp,sigma_spl,ierr_spl)
          if(ierr_spl .ne. 0) 
     1        write(*,'("spline error Jacobian",i3)') ierr_spl
          do irr = 1,ir_fine_scl
            r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
            rjacob_lrg(i,j,irr) = 
     1           curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
          end do     !irr = 1,ir_fine_scl

C-----------------------------------------------------------------------
C         Make spline fits of Gsssup
C-----------------------------------------------------------------------
          do ir=1,irads
            sp_fit(ir) = gsssup(i,j,ir)
          end do
          ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
          call curv1(irads,rho,sp_fit,sp1,sp2,ispl_opt,yp,
     1        temp,sigma_spl,ierr_spl)
          if(ierr_spl .ne. 0) 
     1        write(*,'("spline error GSSSUP",i3)') ierr_spl
          do irr = 1,ir_fine_scl
            r_pt = rho(1)+real(irr-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
            gsssup_lrg(i,j,irr) = 
     1           curv2(r_pt,irads,rho,sp_fit,yp,sigma_spl)
          end do     !irr = 1,ir_fine_scl
        end do     !j=1,ith
      end do     !i=1,izt
      close(unit=9)

C-----------------------------------------------------------------------
C     Compute F Coefficients
C-----------------------------------------------------------------------
      numrads = ir_fine_scl/npes
      do ir = (mype*numrads+1),(mype*numrads+numrads)
        r_pt = rho(1)+real(ir-1)*(rho(irads) - rho(1))
     1           /real(ir_fine_scl-1)
        f1_avg = 0.;f3_avg = 0.
        do i=1,izt
          do j=1,ith
            f1(i,j) = gsssup_lrg(i,j,ir)*rjacob_lrg(i,j,ir)/
     1                 (bfield_lrg(i,j,ir)**2)
            f1_avg = f1_avg + f1(i,j)/real(izt*ith)
            if(.not.lrfp) then
              f3c(i,j) = gsssup_lrg(i,j,ir)/
     1              (rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
              f3_avg = f3_avg + f3c(i,j)/real(izt*ith)
              f3b(i,j) = iota_r(ir)*f3c(i,j)
              f3a(i,j) = iota_r(ir)*f3b(i,j)
            else if(lrfp) then
              f3a(i,j) = gsssup_lrg(i,j,ir)/
     1              (rjacob_lrg(i,j,ir)*(bfield_lrg(i,j,ir)**2))
              f3b(i,j) = f3a(i,j)*iota_r_inv(ir)
              f3c(i,j) = f3b(i,j)*iota_r_inv(ir)
              f3_avg = f3_avg + f3a(i,j)/real(izt*ith)
            end if
          end do
        end do
        if(cyl) then
          if(.not.lrfp) then
            do i=1,izt
              do j=1,ith
                f1(i,j) = f1_avg
              end do
            end do
          else if(lrfp) then
            do i=1,izt
              do j=1,ith
                f1(i,j) = f1_avg
              end do
            end do
          endif
        endif

C-----------------------------------------------------------------------
C       Generate Fourier spectra of above gridded coefficients
C-----------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C       Write out coefficient spectra if half way out in flux
C-----------------------------------------------------------------------
        if(ir .eq. ir_fine_scl/2) then
          do mn = 1,mnmx
            write(8,'(f6.1,2x,f6.1,4(2x,e15.7))')rm(mn),rn(mn),
     1          f1_nm(mn),f3a_nm(mn),
     2          f3b_nm(mn),f3c_nm(mn)
          end do
        end if

C-----------------------------------------------------------------------
C       Build A and B matrices
C-----------------------------------------------------------------------
        do i=1,mn_col
          do j=1,mn_col
            bmat(i,j) = 0
            amat(i,j) = 0
            if(i .lt. j .and. ipos_def_sym) cycle     !i.e. symmetric storage mode: only need bottom half
            ni = in_col(i)
            nj = in_col(j)
            mi = im_col(i)
            mj = im_col(j)
            do ieq = 1,mnmx
              meq = rm(ieq)
              neq = rn(ieq)
              call ccc_convolve(ccci,mi,ni,mj,nj,meq,neq)
              call scs_convolve(scsi,mi,ni,mj,nj,meq,neq)
              bmat(i,j) = bmat(i,j) - ccci*f1_nm(ieq)*mu0_rho_ion(ir)
              amat(i,j) = amat(i,j)
     1                  - (scsi*(f3a_nm(ieq)*rm_col(i)*rm_col(j)
     2                  - f3b_nm(ieq)*rm_col(j)*rn_col(i)
     3                  - f3b_nm(ieq)*rn_col(j)*rm_col(i)
     4                  + f3c_nm(ieq)*rn_col(j)*rn_col(i)))
            end do              !ieq = 1,mnmx
          end do                 !j=1,mn_col
        end do                  !i=1,mn_col

C-----------------------------------------------------------------------
C       Call matrix eigenvalue solver
C-----------------------------------------------------------------------
        egl = 1.d-2; egu = 0.6d0; abstol = 1.d-8; info = 0
        il = 0; iu = 0
        if(.NOT.ipos_def_sym) then
          call dggev('N','V',mn_col,amat,mn_col,bmat,mn_col,alfr,alfi,
     1               beta,vl,ldvl,vr,ldvr,work,lwork,info)
        else if(ipos_def_sym .and. .not.subset_eq) then
          call dsygv(iopt,jobz,'L',mn_col,amat,mn_col,bmat,mn_col,
     1               omega, work, lwork, info)
        else if(ipos_def_sym .and. subset_eq) then
          call dsygvx(iopt,'V','V','L',mn_col,amat,mn_col,bmat,mn_col,
     1              egl, egu, il, iu, abstol, m_red, omega, vr, mn_col,
     2              work, lwork, iwork, ifail, info)
        endif
        if(info .ne. 0) write(*,'("info = ",i8)') info
        do i=1,mn_col
          if(iopt .eq. 1) then
            do j=1,mn_col
              if(ipos_def_sym .and. jobz .eq. 'V' .and. .not.subset_eq)
     1           eig_vect(j) = abs(amat(j,i))
              if(.NOT.ipos_def_sym) eig_vect(j) = abs(vr(j,i))
            end do
       
            eig_max = -1.d+30
            do j=1,mn_col
              if(eig_vect(j) .gt. eig_max) then
                eig_max = eig_vect(j)
                j_max_index = j
              end if
            end do
            m_emax = im_col(j_max_index)
            n_emax = in_col(j_max_index)

            if(ipos_def_sym) then
              write(21,'(2(e15.7,2x),i4,2x,i4)')
     1               r_pt, sqrt(abs(omega(i))), m_emax, n_emax
            else if(.NOT.ipos_def_sym) then
              write(21,'(4(e15.7,2x),i4,2x,i4)')
     1               r_pt, alfr(i), alfi(i), beta(i), m_emax, n_emax
            endif
          else if(iopt .eq. 0) then
            if(ipos_def_sym) then
              write(21,'(e15.7,2x,e15.7)') r_pt,sqrt(abs(omega(i)))
            else if(.NOT.ipos_def_sym) then
              write(21,'(e15.7,3(2x,e15.7))') r_pt,alfr(i),
     1           alfi(i),beta(i)
            endif
          endif
        end do          !do i=1,mn_col
      end do       !ir=1,ir_fine_scl

C-----------------------------------------------------------------------
C     Close Open files
C-----------------------------------------------------------------------
      close(unit=20)
      close(unit=21)
      close(unit=8)
      if(mype .eq. 0) write(*,
     > '("modes = ",i5,2x,"no. of radial points = ",i5)')
     >  mn_col,ir_fine_scl
      if(mype .eq. 0) write(7,'(i5,3(2x,i5))') iopt,mn_col,ir_fine_scl,
     >    isym_opt
      if(mype .eq. 0) close(unit=7)

C-----------------------------------------------------------------------
C     Deallocate
C-----------------------------------------------------------------------
      call trg_deallocate
      deallocate(alpha, beta, aux)

C-----------------------------------------------------------------------
C     Write time
C-----------------------------------------------------------------------
      if(mype .eq. 0) then
        call date_and_time(date, time, zone, values)
        write(*,*) time
      endif

C-----------------------------------------------------------------------
C     Finalize MPI
C-----------------------------------------------------------------------
#if defined (PARALLEL)      
        call MPI_FINALIZE(ier)
#endif
      end program tae_continua
