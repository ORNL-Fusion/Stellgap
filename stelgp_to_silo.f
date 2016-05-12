c
c  this code reads data from alfven_post and converts it
c   to a .silo file for use in Visit. Note below that independent
c   variable axis can be either tor_flux or sqrt(tor_flux). Data
c   coming in from stellgap, post_process is in tor_flux.
c
      implicit none
      include
     >"/Users/dsp/fortran_code_development/Silo-4.6.2/include/silo.inc"
       real, allocatable :: r(:), omega(:), rr(:), w(:)
       integer, allocatable :: m(:), n(:), nn(:), mm(:)
       real :: freq_norm, freq_sort_max, dum1, dum2, dum3
       integer :: npts, i, ierr, dbfile, ndims, istat,
     >     err, ic, idum, jdum
c       freq_norm = 1./584.8
       freq_norm = 1.
       freq_sort_max = 300.
       open(unit=4,file="alfven_post",status="old")
       open(unit=5,file="freq_sort",status="unknown")
       open(unit=6,file="m11.dat",status="unknown")
       open(unit=7,file="m31.dat",status="unknown")
       open(unit=9,file="m51.dat",status="unknown")
       read(4,*) npts
       ic = 0
       do i=1,9000000
        read(4,*,END=71) dum1,dum2,dum3
	ic = ic + 1
       end do
  71   continue
       npts = ic
       close(unit=4)
       write(5,*) freq_sort_max
       allocate(r(npts), stat=istat)
       allocate(omega(npts), stat=istat)
       allocate(rr(npts), stat=istat)
       allocate(w(npts), stat=istat)
       allocate(m(npts), stat=istat)
       allocate(n(npts), stat=istat)
       allocate(nn(npts), stat=istat)
       allocate(mm(npts), stat=istat)
       open(unit=4,file="alfven_post",status="old")
       read(4,*) idum
       ic = 0
       do i=1,npts
       read(4,*) r(i),omega(i),m(i),n(i)
c       read(4,*) r(i),omega(i),idum,jdum,m(i),n(i)
        omega(i) = omega(i)*freq_norm
c	if(omega(i) .lt. 300.) then
	 ic = ic + 1
	 w(ic) = omega(i)
c	 rr(ic) = r(i)        !plot vs. psi
	 rr(ic) = sqrt(r(i))  !plot vs. sqrt(psi)
c	 rr(ic) = -0.00087078 + 1.1526*sqrt(r(i)) - 0.11225*r(i)
c     >            + 0.21002*r(i)*sqrt(r(i)) -0.24732*r(i)*r(i)  !plot vs. sqrt(psi-pol)
	 nn(ic) = n(i)
	 mm(ic) = m(i)
	 if(omega(i) .lt. freq_sort_max) then
	  write(5,*) m(i), n(i)
	end if
	if(mm(ic) .eq. 1 .and. abs(nn(ic)) .eq. 1) then
	  write(6,*) r(i), omega(i)
	else if(mm(ic) .eq. 3 .and. abs(nn(ic)) .eq. 1) then
	  write(7,*) r(i), omega(i)
	else if(mm(ic) .eq. 5 .and. abs(nn(ic)) .eq. 1) then
	  write(9,*) r(i), omega(i)
	endif
       end do
       npts = ic
       write(*,*) rr(1), rr(npts)


c
c     open and write silo formatted file for Visit
c        
      ierr = dbcreate("stellgap.silo", 13, DB_CLOBBER, DB_LOCAL,
     >  "Comment about the data", 22, DB_PDB, dbfile)
      if(dbfile.eq.-1) then
      write (6,*) 'Could not create Silo file!\n'
      goto 10000
      endif
c Add other Silo calls here.
       ndims = 2
       err = dbputpm(dbfile, "pointmesh",9, ndims, rr, w,
     >  DB_F77NULL, npts, DB_FLOAT, DB_F77NULL, ierr)
       write(*,*) ierr
       err = dbputpv1(dbfile, "m", 1, "pointmesh", 9,
     >  mm, npts, DB_INT, DB_F77NULL, ierr)
       write(*,*) ierr
       err = dbputpv1(dbfile, "n", 1, "pointmesh", 9,
     >  nn, npts, DB_INT, DB_F77NULL, ierr)
       write(*,*) ierr
c Close the Silo file. 
      ierr = dbclose(dbfile)
10000  stop
       stop
       end
