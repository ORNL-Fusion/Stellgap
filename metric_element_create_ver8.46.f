*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program metric_element_create.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
c   7/29/2010: A test was included to see if all the surfaces to the edge were included in
c              the xbooz_xform calculation. i.e., if(rmncbh(1,nsd) .eq. 0.)... This can also
c              be seen from running ncdump on the boozmn file and comparing ns_b to comput_surfs.
c              If the last surface is not included, then there can be problems with the outer
c              surface in boozmn since averages between surfaces are made. This particularly
c              shows up when xmetric is used for the outer surface metric calculation. In this
c              case, the outer surface area can be too low by a factor of ~4 due to the outer R
c              and z being down by 1/2 (from averaging over a 0 plus the correct value).
c
      PROGRAM XMETRIC

!-----------------------------------------------------------------------
!   M o d u l e s
!-----------------------------------------------------------------------
      use read_boozer_mod, iota_bw=>iota_b, bmnc_bw=>bmnc_b,
     &    rmnc_bw=>rmnc_b,zmns_bw=>zmns_b,gmn_bw=>gmnc_b,
     &    pmns_bw=>pmns_b,pres_bw=>pres_b,phip_bw=>phip_b,
     &    bvco_bw=>bvco_b,buco_bw=>buco_b,ixm_bw=>ixm_b,
     &    ixn_bw=>ixn_b
      use stel_kinds

!-----------------------------------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------------------------------
      implicit none
      real(kind=rprec), parameter :: p5=0.5_dp, one=1.0_dp
      REAL(kind=rprec), PARAMETER :: twopi  = 6.28318530717958623_dp
      REAL(kind=rprec), PARAMETER :: mu_0   = 2.0e-7_dp*twopi
      real(kind=rprec), parameter :: zero = 0.0_dp, two = 2.0_dp

C-----------------------------------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------------------------------
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: hiota, hpres,
     1   hphip, hjtor, hjpol
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: iotaf, jpolf,
     1   jtorf, phipf, presf
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: iotapf, jpolpf,
     1   jtorpf, phippf, prespf
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: xmb,xnb,xm,xn
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: radii, radii_flux
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmncbf, zmnsbf,
     1   pmnsbf, bmncbf
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmncpbf,zmnspbf,
     1   pmnspbf, bmncpbf
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmncbh,zmnsbh,
     1   pmnsbh, bmncbh
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: bfield,
     1  bfields, bfieldze, bfieldth, rjacob, rjacobze, rjacobth,
     2  rjacobs,gsssup, gttsup, gzzsup, gstsup,
     3  gszsup, gtzsup, brho
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: jprl_coef0,
     1  jprl_coef1, jprl_coef2
      REAL(kind=rprec), DIMENSION(3,3) :: mat_lowr
      REAL(kind=rprec), DIMENSION(3,3) :: mat_upr
      REAL(kind=rprec), DIMENSION(3,3) :: mat_diag
       integer :: nsurf
       real(kind=rprec):: phipc, iotac, jtorc, jpolc,
     1  presc, prespc, jtorpc, jpolpc, phippc, iotapc 
       real(kind=rprec) :: phis, phize, phith, rboo, rs,
     2  rze, rth, zs, zze, zth, arg, ccosi, ssine, num1,
     3  num2, ibf2, ibf3, rboo2, rjac2, rjac2i, gtssub, gstsub,
     5  gzssub, gszsub, gtzsub, gztsub, gttsub, gsssub,
     6  gzzsub, gaasup, gsasup, gtssup, gzssup, gztsup,
     8  den1, cka, beta, t1, t2, t3, t4, t5, cks,
     9  thetang, zetang, det, error, rjac2_avg, det_avg,
     1  zboo, phiboo, xr, yr, zr
      INTEGER :: nfp, nsd, mnboz, i, j, k, mn, kz, kt, ks
      REAL(kind=rprec) :: r0,b0,amin,beta0,r0max,r0min,
     1  aspect, betaxis, ohs2, mat_test_diag, mat_test_offdiag
      character arg1*40,warg1*45,bozout*45,arg2*40
      integer nargs, numargs, numchars, itheta, izeta, nznt
      integer iargc, unit_no, istat, ierr, ig, lf, is
      real viz_flux   !for plotting interior flux surfaces in AVS
      real surf_area_element, surf_area_total
      logical lasym, viz, test_jacob, test_upr_lowr,
     1  make_stellgap_data,make_full_torus,surf_compute

C-----------------------------------------------------------------------
C     Local runtime settings
C-----------------------------------------------------------------------
      itheta = 80                 ! # of theta gridpoints
      izeta = 80                  ! # of zeta gridpoints
      make_stellgap_data = .true. ! Output tae_data_boozer file
      viz = .false.                ! Output visualization (AVS) data
      viz_flux = 0.5              ! Surface to output AVS data
      make_full_torus = .false.   ! Output full_torus_coords file
      surf_compute = .false.      ! Output surf_area_elements file
      test_jacob = .false.        ! Testing of Jacobian for consistency
      test_upr_lowr = .false.     ! Testing of Boozer metrics
C-----------------------------------------------------------------------
C     Get command line input
C-----------------------------------------------------------------------
      numargs = iargc()
      CALL getarg(1,arg1)
      if( numargs.lt.1 )then
        print *,' Try calling with -help for options.'
        stop
      endif
      i = 2
      IF (arg1 == '-help') i=1
      DO WHILE (i<=numargs)
        CALL GETARG(i,arg2)
        SELECT CASE (arg2)
          CASE ("-surf_compute")
            surf_compute = .TRUE.
          CASE ("-test_jacob")
            test_jacob = .TRUE.
          CASE ("-test_upr_lowr")
            test_upr_lowr = .TRUE.
          CASE ("-make_full_torus")
            make_full_torus = .TRUE.
          CASE ("-viz")
            viz = .TRUE.
            i = i + 1
            CALL GETARG(i,arg2)
            READ(arg2,*) viz_flux
          CASE ("-itheta")
            i = i + 1
            CALL GETARG(i,arg2)
            READ(arg2,*) itheta
          CASE ("-izeta")
            i = i + 1
            CALL GETARG(i,arg2)
            READ(arg2,*) izeta
          CASE ("-help","-h")
            WRITE(6, '(A)') ' Initialization code for STELLGAP/AE3D'
            WRITE(6, '(A)') ' Usage: xmetric EXT <options>'
            WRITE(6, '(A)') '    EXT: Boozmn file extension'
            WRITE(6, '(A)') '    <options>'
            WRITE(6, '(A)') 
     1         '     -itheta 90:       Number of theta gridpoints'
            WRITE(6, '(A)') 
     1         '     -izeta 90:        Number of zeta gridpoints'
            WRITE(6, '(A)') 
     1         '     -surf_compute:    Produce surf_area_elements file'
            WRITE(6, '(A)') '     -test_jacob:      Jacobian test'
            WRITE(6, '(A)') '     -test_upr_lowr:   Metric Test'
            WRITE(6, '(A)') 
     1         '     -make_full_torus: Produce full_torus_coords file'
            WRITE(6, '(A,A)') 
     1         '     -viz 0.5:          Produce cart_coords and',
     2         ' metrics file at s'
            WRITE(6, '(A)') 
     1         '     -help:             Produce this help message'
            STOP
        END SELECT
        i = i + 1
      END DO

C-----------------------------------------------------------------------
C     Open output files
C-----------------------------------------------------------------------
      if(surf_compute) then
        open(unit=47,file="surf_area_elements",status="unknown")
      endif
      if(make_stellgap_data) then
        open(unit=20,file="tae_data_boozer",status="unknown")
      endif
      if(viz) then
        WRITE(6,'(A,F5.3)') "Using viz_flux = ",viz_flux 
        open(unit=12,file="cart_coords",status="unknown")
        open(unit=14,file="metrics",status="unknown")
      endif
      open(unit=15,file="ae_metric.dat",status="unknown")

C-----------------------------------------------------------------------
C     Read Boozmn file
C-----------------------------------------------------------------------
      warg1 = arg1
      WRITE(6,*) 'READING: ',warg1
      call read_boozer_file(warg1,ierr)
      if (istat.ne.0) stop 22

C-----------------------------------------------------------------------
C     Process Boozmn data
C-----------------------------------------------------------------------
      nfp = nfp_b
      nsd = ns_b
      aspect = aspect_b
      r0max = rmax_b 
      r0min = rmin_b
      betaxis = betaxis_b
      mnboz = mnboz_b
      if(surf_compute) write(47,*) nfp, izeta, itheta

C-----------------------------------------------------------------------
C     Process IOTA, PRES, PHIP, JTOR, JPOL
C        Pressure is given in Pascals for versions > 6.00
C        PHIP = -PHIP_BOOZER
C        JTOR = I = -BUCO_BOOZER
C        JPOL = J = BVCO_BOOZER
C        All on half-mesh
C-----------------------------------------------------------------------
      allocate (hiota(nsd), hpres(nsd), hjpol(nsd), hjtor(nsd), 
     1          hphip(nsd), jprl_coef0(nsd), jprl_coef1(nsd),
     2          jprl_coef2(nsd), stat=istat)
      if (istat .ne. 0) stop 23
      do k=1,nsd
        hiota(k) = iota_bw(k)
        hpres(k) = mu_0*pres_bw(k)
        hphip(k) = -phip_bw(k)
        hjpol(k) = bvco_bw(k)
        hjtor(k) = -buco_bw(k) 
      end do
      r0 = (r0max+r0min)/two
      amin = r0/aspect
      beta0 = betaxis
      write(*,'(///)')
      WRITE(6,'(4X,A,I3)') 'ITHETA: ',itheta
      WRITE(6,'(4X,A,I3)') 'IZETA: ',izeta
      WRITE(6,'(4X,A,I3)') 'FIELD PERIODS: ',nfp
      WRITE(6,'(4X,A,I3)') 'NS: ',nsd
      WRITE(6,'(4X,A,F6.3)') 'ASPECT: ',ASPECT
      !write(*,*) nfp,nsd,aspect
      write(*,'(/)')
      WRITE(6,'(4X,A,F6.3)') 'R0MAX: ',r0max
      WRITE(6,'(4X,A,F6.3)') 'R0MIN: ',r0min
      WRITE(6,'(4X,A,F6.4)') 'BETAAXIS: ',betaxis
      WRITE(6,'(4X,A,I6)') 'MNBOZ: ',mnboz
      !write(*,*) r0max,r0min,betaxis,mnboz
      write(*,'(/)')
      WRITE(6,'(4X,A,F6.3)') 'R0: ',r0
      WRITE(6,'(4X,A,F6.3)') 'AMIN: ',AMIN
      WRITE(6,'(4X,A,F6.3)') 'BETA0: ',beta0
      !write(*,*) r0,amin,beta0

C-----------------------------------------------------------------------
C     Process Boozmn Fourier data (all on half mesh)
C-----------------------------------------------------------------------
      allocate (xm(mnboz), xn(mnboz), stat=istat)
      if (istat .ne. 0) stop 24
      do mn=1,mnboz
      xm(mn) = ixm_bw(mn)
      xn(mn) = ixn_bw(mn)            
      end do

      allocate (rmncbh(mnboz,nsd), zmnsbh(mnboz,nsd), 
     1           pmnsbh(mnboz,nsd), bmncbh(mnboz,nsd), stat = istat)
      if (istat .ne. 0) stop 25

      zmnsbh=zero; rmncbh=zero; pmnsbh=zero; bmncbh=zero

      do k=1, nsd
        do mn = 1,mnboz
          bmncbh(mn,k) = bmnc_bw(mn,k)
          rmncbh(mn,k) = rmnc_bw(mn,k)
          zmnsbh(mn,k) = zmns_bw(mn,k)
          pmnsbh(mn,k) = pmns_bw(mn,k)                  
        enddo
      enddo

C-----------------------------------------------------------------------
C     Deallocate Boozer data
C-----------------------------------------------------------------------
      call read_boozer_deallocate

C-----------------------------------------------------------------------
C     Compute radial helper variables
C        r/a= (radii_flux)**1/2
C        s = radii_flux
C-----------------------------------------------------------------------
      allocate (radii(nsd), radii_flux(nsd), stat = istat)
      do i=1,nsd
        radii(i)=sqrt(real(i-1,kind=rprec)/(nsd-1))
        radii_flux(i) = real(i-1,kind=rprec)/(nsd-1)
      enddo

C-----------------------------------------------------------------------
C     Warn if missing last surface from Boozmn file
C-----------------------------------------------------------------------
      if(rmncbh(1,nsd) .eq. 0.) then
        WRITE(6,'(A)') "WARNING: Last surface not found in boozmn file"
      end if
   
C-----------------------------------------------------------------------
C     Interpolate all quantities from half to full radial mesh
C-----------------------------------------------------------------------
      allocate(iotaf(nsd), presf(nsd), jtorf(nsd), jpolf(nsd),
     1         phipf(nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error 01'

      iotaf=zero ; presf =zero; jtorf =zero; jpolf=zero; phipf=zero

      iotaf(2:nsd-1) = p5*(hiota(2:nsd-1) + hiota(3:nsd))
      presf(2:nsd-1) = p5*(hpres(2:nsd-1) + hpres(3:nsd))
      jtorf(2:nsd-1) = p5*(hjtor(2:nsd-1) + hjtor(3:nsd))
      jpolf(2:nsd-1) = p5*(hjpol(2:nsd-1) + hjpol(3:nsd))
      phipf(2:nsd-1) = p5*(hphip(2:nsd-1) + hphip(3:nsd))
 
C-----------------------------------------------------------------------
C     Compute derivative quantities on radial full mesh
C-----------------------------------------------------------------------       
      allocate(iotapf(nsd), jpolpf(nsd), jtorpf(nsd),
     1         phippf(nsd), prespf(nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error 02'

      iotapf=zero; prespf =zero; jtorpf =zero; jpolpf=zero; phippf=zero

      ohs2 = dble(nsd-1)
      iotapf(2:nsd-1) = ohs2*(hiota(3:nsd) - hiota(2:nsd-1))
      prespf(2:nsd-1) = ohs2*(hpres(3:nsd) - hpres(2:nsd-1))
      jtorpf(2:nsd-1) = ohs2*(hjtor(3:nsd) - hjtor(2:nsd-1))
      jpolpf(2:nsd-1) = ohs2*(hjpol(3:nsd) - hjpol(2:nsd-1))
      phippf(2:nsd-1) = ohs2*(hphip(3:nsd) - hphip(2:nsd-1))
 
C-----------------------------------------------------------------------
C     Deallocate half mesh quantities
C-----------------------------------------------------------------------        
      deallocate(hiota,hpres,hjtor,hjpol,hphip)

C-----------------------------------------------------------------------
C     store (m,n)-descriptors
C-----------------------------------------------------------------------  
      allocate (xnb(mnboz), xmb(mnboz), stat=k)
      if (k .ne. 0) stop 'Allocation error 03'
      xnb=xn
      xmb=xm

C-----------------------------------------------------------------------
C     Boozer fourier coefs and derivatives on full radial mesh
C-----------------------------------------------------------------------  
      allocate (rmncbf(mnboz,nsd), zmnsbf(mnboz,nsd), pmnsbf(mnboz,nsd),
     1          bmncbf(mnboz,nsd), rmncpbf(mnboz,nsd), 
     2          zmnspbf(mnboz,nsd),pmnspbf(mnboz,nsd),
     3          bmncpbf(mnboz,nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error 04'
 
      rmncbf=zero; zmnsbf=zero; pmnsbf=zero; bmncbf=zero
      rmncpbf=zero; zmnspbf=zero; pmnspbf=zero; bmncpbf=zero
      do mn = 1,mnboz
        do k = 1,nsd-1
          rmncbf(mn,k) = p5*(rmncbh(mn,k+1)+rmncbh(mn,k))
          zmnsbf(mn,k) = p5*(zmnsbh(mn,k+1)+zmnsbh(mn,k))
          pmnsbf(mn,k) = p5*(pmnsbh(mn,k+1)+pmnsbh(mn,k))
          bmncbf(mn,k) = p5*(bmncbh(mn,k+1)+bmncbh(mn,k))       
          rmncpbf(mn,k) = ohs2*(rmncbh(mn,k+1)-rmncbh(mn,k))
          zmnspbf(mn,k) = ohs2*(zmnsbh(mn,k+1)-zmnsbh(mn,k))
          pmnspbf(mn,k) = ohs2*(pmnsbh(mn,k+1)-pmnsbh(mn,k))
          bmncpbf(mn,k) = ohs2*(bmncbh(mn,k+1)-bmncbh(mn,k))
        enddo
      enddo

C-----------------------------------------------------------------------
C     Deallocate half mesh quantities
C-----------------------------------------------------------------------
      deallocate(xn, xm, rmncbh, zmnsbh, pmnsbh, bmncbh, stat=k)

C-----------------------------------------------------------------------
C     Allocate realspace variables (ig = itheta*izeta)
C-----------------------------------------------------------------------
      ig = itheta*izeta
      allocate(bfield(ig,nsd), bfields(ig,nsd), bfieldze(ig,nsd),
     1         bfieldth(ig,nsd), rjacob(ig,nsd),
     1         rjacobze(ig,nsd), rjacobth(ig,nsd),
     2         rjacobs(ig,nsd), gsssup(ig,nsd), gttsup(ig,nsd),
     3         gzzsup(ig,nsd), gstsup(ig,nsd), gszsup(ig,nsd),
     4         gtzsup(ig,nsd), brho(ig,nsd), stat=istat)
      if (istat .ne. 0) stop 'Allocation error 05'

C-----------------------------------------------------------------------
C     Begin main computation loop over surfaces
C-----------------------------------------------------------------------     
      flux_surface: do ks = 2,nsd-1 
C-----------------------------------------------------------------------
C       Radial quantities
C----------------------------------------------------------------------- 
        phipc=phipf(ks)      ! toroidal magnetic flux
        iotac=iotaf(ks)      ! iota
        iotapc=iotapf(ks)    ! radial iota derivative
        jtorc=jtorf(ks)      ! torodial current density
        jpolc=jpolf(ks)      ! poloidal current density
        presc=presf(ks)      ! pressure
        prespc=prespf(ks)    ! radial pressure gradient
        jtorpc=jtorpf(ks)    ! radial derivative toroidal current density
        jpolpc=jpolpf(ks)    ! radial derivative poloidal current density
        phippc=phippf(ks)    ! radial derivative toroidal flux
        jprl_coef0(ks) = (jpolc*jtorpc - jpolpc*jtorc)
     1                   /((jpolc + jtorc*iotac)*phipc)
        jprl_coef1(ks) = jtorc/((jpolc + jtorc*iotac)*phipc)
        jprl_coef2(ks) = -jpolc/((jpolc + jtorc*iotac)*phipc)

        if(make_stellgap_data) then
          write(20,21) ks,iotac,phipc,jtorc,jpolc
  21      format(1x,i3,4(2x,e15.7))
        endif
C-----------------------------------------------------------------------
C       Fourier Inversion
C----------------------------------------------------------------------- 
        lf = 0
        error = 0.; rjac2_avg = 0.; det_avg = 0.
        surf_area_total = 0.
        toroidal: do kz = 1,izeta
          poloidal: do kt = 1,itheta
            lf = lf + 1
            thetang = twopi*dble(kt-1)/dble(itheta-1)
            zetang = twopi*dble(kz-1)/(dble(nfp)*dble(izeta-1))
!...        initialize for Fourier inversion
            bfield(lf,ks) = zero
            bfieldze(lf,ks) = zero
            bfields(lf,ks) = zero
            bfieldth(lf,ks) = zero
            phis = zero
            phize = one  ! since Phi_cyl  =  Phi_VMEC  = Phi_BOOZ + Pmn
            phith = zero
            rboo = zero; zboo = zero; phiboo = zero
            rs = zero
            rze = zero
            rth = zero
            zs = zero
            zze = zero
            zth = zero  
!...        Fourier invert B, R, Z, Phi and derivatives
            fourier: do j = 1,mnboz
              arg = xmb(j)*thetang-zetang*xnb(j)
              ccosi = cos(arg)
              ssine = sin(arg)
              bfield(lf,ks) = bfield(lf,ks)+bmncbf(j,ks)*ccosi    !magnetic field magnitude
              bfields(lf,ks) = bfields(lf,ks)+bmncpbf(j,ks)*ccosi ! ..... radial derivative
              bfieldze(lf,ks) = bfieldze(lf,ks)
     1                          +xnb(j)*bmncbf(j,ks)*ssine        ! ..... zeta derivative
              bfieldth(lf,ks) = bfieldth(lf,ks)
     1                          -xmb(j)*bmncbf(j,ks)*ssine        ! ..... theta derivative
              rboo = rboo+rmncbf(j,ks)*ccosi                      ! cylindrical R 
              zboo = zboo+zmnsbf(j,ks)*ssine                      ! cylindrical Z
              phiboo = phiboo+pmnsbf(j,ks)*ssine                  ! cylindrical Phi
              rth = rth-rmncbf(j,ks)*xmb(j)*ssine                 ! ..... theta derivative
              rze = rze+rmncbf(j,ks)*xnb(j)*ssine                 ! ..... zeta derivative
              rs = rs+rmncpbf(j,ks)*ccosi                         ! ..... radial derivative
              zth = zth+zmnsbf(j,ks)*xmb(j)*ccosi                 ! cylindrical Z: theta derivative
              zze = zze-zmnsbf(j,ks)*xnb(j)*ccosi                 ! ..... zeta derivative
              zs = zs+zmnspbf(j,ks)*ssine                         ! ..... radial derivative
              phith = phith+pmnsbf(j,ks)*xmb(j)*ccosi             ! cylindrical Phi: theta derivative
              phize = phize-pmnsbf(j,ks)*xnb(j)*ccosi             ! ..... zeta derivative
              phis = phis+pmnspbf(j,ks)*ssine                     ! ..... radial derivative
            enddo fourier
            phiboo = phiboo + zetang

!...        Jacobian from cylindrical to Boozer  !note: this is sqrt(g), not 1/sqrt(g)
            num1 = (iotac*jtorc-jpolc)*phipc    
            num2 = (iotac*jtorc-jpolc)*phippc+(iotapc*jtorc+
     1              iotac*jtorpc-jpolpc)*phipc 
            ibf2 = one/bfield(lf,ks)**2    
            ibf3 = -2.0_dp/bfield(lf,ks)**3
            rjacob(lf,ks) = num1*ibf2                           ! jacobian from cyl. to Boozer
            rjacobze(lf,ks) = num1*bfieldze(lf,ks)*ibf3         ! ..... zeta derivative
            rjacobth(lf,ks) = num1*bfieldth(lf,ks)*ibf3         ! ..... theta derivative
            rjacobs(lf,ks) = num1*bfields(lf,ks)*ibf3+num2*ibf2 ! ..... radial derivative 
            rboo2 = rboo**2
            rjac2 = rjacob(lf,ks)**2
            rjac2i = one/rjac2
 
!...        Boozer lower metric elements
            gtssub = rth*rs+zth*zs+rboo2*phith*phis
            gstsub = gtssub
            gzssub = rze*rs+zze*zs+rboo2*phize*phis
            gszsub = gzssub                                   
            gtzsub = rth*rze+zth*zze+rboo2*phith*phize
            gztsub = gtzsub
            gttsub = (rth**2)+(zth**2)+rboo2*(phith**2)
            gsssub = (rs**2)+(zs**2)+rboo2*(phis**2)
            gzzsub = (rze**2)+(zze**2)+rboo2*(phize**2)
            if (test_upr_lowr) then
              mat_lowr(1,1)=gsssub
              mat_lowr(1,2)=gstsub
              mat_lowr(1,3)=gszsub
              mat_lowr(2,1)=gtssub
              mat_lowr(2,2)=gttsub
              mat_lowr(2,3)=gtzsub
              mat_lowr(3,1)=gzssub
              mat_lowr(3,2)=gztsub
              mat_lowr(3,3)=gzzsub
            endif

!...        Compute surface area
            if (ks .eq. nsd-1 .and. surf_compute) then
              surf_area_element = sqrt(gttsub*gzzsub - gtzsub*gtzsub)
              surf_area_total = surf_area_total + surf_area_element
     1                    *(twopi**2)/(dble(itheta-1)*dble(izeta-1))
              write(47,'(3e16.8)') zetang,thetang,surf_area_element
            endif
       
!...        Consistency test:
            if(test_jacob) then
              det = gsssub*(gttsub*gzzsub - gztsub*gtzsub)
     1            + gstsub*(gtzsub*gzssub - gtssub*gzzsub)
     2            + gszsub*(gtssub*gztsub - gttsub*gzssub)
              det_avg = det_avg + det/dble(ig)
              rjac2_avg = rjac2_avg + rjac2/dble(ig)
              error = error + (det - rjac2)/(rjac2*dble(ig))
            end if

!...        Boozer upper metric elements
            gsssup(lf,ks) = (gttsub*gzzsub-gtzsub*gztsub)*rjac2i
            gttsup(lf,ks) = (gsssub*gzzsub-gszsub*gzssub)*rjac2i
            gzzsup(lf,ks) = (gsssub*gttsub-gstsub*gtssub)*rjac2i  
            gstsup(lf,ks) = (gztsub*gzssub-gstsub*gzzsub)*rjac2i       
            gtssup = gstsup(lf,ks)
            gszsup(lf,ks) = (gtssub*gtzsub-gttsub*gzssub)*rjac2i
            gzssup = gszsup(lf,ks)
            gtzsup(lf,ks) = (gstsub*gzssub-gsssub*gztsub)*rjac2i
            gztsup = gtzsup(lf,ks)
            brho(lf,ks) = -(gstsup(lf,ks)*jtorc + gszsup(lf,ks)*jpolc)
     1                     /gsssup(lf,ks)

!...        Consistency test:
            if (test_upr_lowr) then
              mat_upr(1,1)=gsssup(lf,ks);mat_upr(1,2)=gstsup(lf,ks)
              mat_upr(1,3)=gszsup(lf,ks)
              mat_upr(2,1)=gtssup;mat_upr(2,2)=gttsup(lf,ks)
              mat_upr(2,3)=gtzsup(lf,ks)
              mat_upr(3,1)=gzssup;mat_upr(3,2)=gztsup
              mat_upr(3,3)=gzzsup(lf,ks)
              mat_diag = matmul(mat_lowr,mat_upr)
              mat_test_diag = (mat_diag(1,1)+mat_diag(2,2)
     1                        +mat_diag(3,3))/3.
              mat_test_offdiag = (mat_diag(1,2)+mat_diag(1,3)
     1                           +mat_diag(2,1)+mat_diag(2,3)
     2                           +mat_diag(3,1)+mat_diag(3,2))/6.
              if(abs(mat_test_diag-1.) .gt. .04) then
                write(*,'(e15.7,2x,e15.7)') mat_test_diag, 
     1                                      mat_test_offdiag
              end if
            endif

!...        Upper metric elements in (s,alpha,phi)-coord.
            gaasup = gttsup(lf,ks)+(iotac**2)*gzzsup(lf,ks)
     1              -2.*iotac*gtzsup(lf,ks)
            gsasup = gstsup(lf,ks)-iotac*gszsup(lf,ks)              

!...        Boozer geodesic curvature
            den1 = 2.0_dp*(iotac*jtorc-jpolc)*rjacob(lf,ks)
            cka = (jpolc*rjacobth(lf,ks)+jtorc*rjacobze(lf,ks))/den1        

!...        Boozer normal curvature
            beta = (gstsup(lf,ks)*jtorc-gszsup(lf,ks)*jpolc)
     1             /gsssup(lf,ks)    
            t1 = iotapc*jtorc+iotac*jtorpc-jpolpc+(phippc/phipc)
     1          *(iotac*jtorc-jpolc)
            t2 = 2.0_dp*prespc/phipc
            t3 = -(iotac*jtorc-jpolc)*rjacobs(lf,ks)
            t4 = -beta*iotac*rjacobth(lf,ks)
            t5 = -beta*rjacobze(lf,ks)
            cks = (rjacob(lf,ks)*t1+rjac2*t2+t3+t4+t5)/den1

!...        Section to write out 3D data for AVS:
            if(ks .eq. (nsd-1) .and. viz) then
               xr = rboo*cos(phiboo)
               yr = rboo*sin(phiboo)
               zr = zboo
               write(12,'(e15.7,2(2x,e15.7))') xr,yr,zr
               write(14,'(e15.7,7(2x,e15.7))') gsssup(lf,ks),
     1            gttsup(lf,ks), gzzsup(lf,ks), gstsup(lf,ks),
     2            gszsup(lf,ks), gtzsup(lf,ks), bfield(lf,ks)
            endif    !ks .eq. nsd-1 

!...        Output the STELLGAP information
            if(make_stellgap_data) then
              write(20,'(1x,4(e24.12,2x),e24.12)') thetang,
     1            zetang, bfield(lf,ks), gsssup(lf,ks),
     2            gtzsup(lf,ks)

              write(20,'(1x,4(e24.12,2x),e24.12)') gttsup(lf,ks),
     1            gzzsup(lf,ks), gstsup(lf,ks), gszsup(lf,ks),
     2            rjacob(lf,ks)
            endif

          end do poloidal
        end do toroidal
        if(ks .eq. 3) nznt = lf
        if(test_jacob) then
          write(*,'(i3,3(2x,e15.7))') ks, error, det_avg, rjac2_avg
        endif
      end do flux_surface

C-----------------------------------------------------------------------
C     Output Metric data (for AE3D)
C----------------------------------------------------------------------- 
      if(surf_compute)
     1   write(*,'("Total surface area = ",e15.7,"(m**2)")')
     2           surf_area_total
      write(15,*) (nsd-2), izeta, itheta, nznt
      do ks = 2,nsd-1
        write(15,19) iotaf(ks), iotapf(ks), jpolf(ks), jpolpf(ks),
     1               jtorf(ks), jtorpf(ks), phipf(ks), phippf(ks)
      end do
      do ks = 2,nsd-1
        do lf=1,nznt
          write(15,18) rjacob(lf,ks), bfield(lf,ks),
     1                 gsssup(lf,ks),
     2                 gttsup(lf,ks), gzzsup(lf,ks),
     3                 gstsup(lf,ks),
     4                 gszsup(lf,ks), gtzsup(lf,ks),
     5                 bfields(lf,ks),
     6                 bfieldth(lf,ks), bfieldze(lf,ks)
        end do
      end do

C-----------------------------------------------------------------------
C     Output J_prl/B coefficients (for AE3D)
C----------------------------------------------------------------------- 
      do ks = 2,nsd-1
        write(15,48) jprl_coef0(ks),jprl_coef1(ks),
     1               jprl_coef2(ks), prespf(ks)
      end do
      do ks = 2,nsd-1
        do lf=1,nznt
          write(15,49) brho(lf,ks)
        end do
      end do
   18  format(e15.7,10(2x,e15.7))
   19  format(e15.7,7(2x,e15.7))
   48  format(e15.7,3(2x,e15.7))
   49  format(e15.7)

C-----------------------------------------------------------------------
C     Output full torus data (full_torus_coords)
C       Note we redefine itheta and izeta here
C----------------------------------------------------------------------- 
      if(make_full_torus) then
        itheta = 300; izeta = 300
        ks = (nsd-1)*viz_flux
        open(unit=17,file="full_torus_coords",status="unknown")
        do kz = 1,izeta
          do kt = 1,itheta
            thetang = twopi*dble(kt-1)/dble(itheta-1)
            zetang = 0.9*twopi*dble(kz-1)/dble(izeta-1)
            rboo = 0.d0; zboo = 0.d0; phiboo = 0.d0
            do j = 1,mnboz             ! Fourier invert B, R, Z, Phi and derivatives
              arg = xmb(j)*thetang-zetang*xnb(j)
              ccosi = cos(arg)
              ssine = sin(arg)
              rboo = rboo+rmncbf(j,ks)*ccosi     ! cylindrical R
              zboo = zboo+zmnsbf(j,ks)*ssine     ! cylindrical Z
              phiboo = phiboo+pmnsbf(j,ks)*ssine ! cylindrical Phi
            enddo
            phiboo = phiboo + zetang
            xr = rboo*cos(phiboo)
            yr = rboo*sin(phiboo)
            zr = zboo
            write(17,'(e15.7,2(2x,e15.7))') xr,yr,zr
          end do
        end do

        open(unit=18,file="torus_slice_coords",status="unknown")
        write(*,*) nsd-2, itheta
        zetang = 0.9*twopi
        do ks = 2,int((nsd-1)*viz_flux)
          do kt = 1,itheta
            thetang = twopi*dble(kt-1)/dble(itheta-1)
            rboo = 0.d0; zboo = 0.d0; phiboo = 0.d0
            do j = 1,mnboz             ! Fourier invert B, R, Z, Phi and derivatives
              arg = xmb(j)*thetang-zetang*xnb(j)
              ccosi = cos(arg)
              ssine = sin(arg)
              rboo = rboo+rmncbf(j,ks)*ccosi                    ! cylindrical R
              zboo = zboo+zmnsbf(j,ks)*ssine                    ! cylindrical Z
              phiboo = phiboo+pmnsbf(j,ks)*ssine                ! cylindrical Phi
            enddo
            phiboo = phiboo + zetang
            xr = rboo*cos(phiboo)
            yr = rboo*sin(phiboo)
            zr = zboo
            write(18,'(e16.7,2(2x,e15.7))') xr,yr,zr
          end do
        end do
      endif
      END PROGRAM XMETRIC 
