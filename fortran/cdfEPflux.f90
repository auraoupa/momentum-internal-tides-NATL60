PROGRAM cdfEPflux
  !!======================================================================
  !!                     ***  PROGRAM  cdfEPflux  ***
  !!=====================================================================
  !!  ** Purpose : Compute the zonal and merdidional component of the EP flux
  !!
  !!  ** Method  : Formula from 2019, Shakespear & Hogg article 
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE eos

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: nlev               ! number of output levels
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: ijarg              !
  INTEGER(KIND=4)                           :: ivar               !
  INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
  INTEGER(KIND=4), DIMENSION(2)             :: ipk                ! output variable properties
  INTEGER(KIND=4), DIMENSION(2)             :: id_varout          ! output variable properties
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nilev             ! level to be processed

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f ! horizontql metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn             ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn             ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: wn, bn             ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: n2                 ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ff                 ! mask
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim, zdep, gdep    ! time counter
  REAL(KIND=4)                              :: zmask              ! mask at T point for -T option

  REAL(KIND=4)                                 :: rau0  = 1000.   ! density of water
  REAL(KIND=4)                                 :: grav  = 9.81    ! Gravity

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: depx, depy              ! curl 

  CHARACTER(LEN=256)                        :: cf_ufil,cf_vfil    ! file names
  CHARACTER(LEN=256)                        :: cf_wfil,cf_bfil    ! file names
  CHARACTER(LEN=256)                        :: cf_n2fil            ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'epflux.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v         ! variable names
  CHARACTER(LEN=256)                        :: cv_w, cv_b         ! variable names
  CHARACTER(LEN=256)                        :: cv_n2              ! variable names
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(2)             :: stypvar            ! structure for attibutes

  LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                   :: lnc4=.FALSE.       ! flag for netcdf4 output with chunking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfEPflux -u U-file U-var -v V-file V-var -w W-file W-var -b buoy-file buoy-var '
     PRINT *,'         -n n2-file n2-var ...  [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the EPflux, at a specified level. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file U-var : file and variable name for zonal velocity'
     PRINT *,'       -v V-file V-var : file and variable name for meridional velocity'
     PRINT *,'       -w W-file W-var : file and variable name for vertical velocity'
     PRINT *,'       -b buoy-file buoy-var : file and variable name for buoyancy'
     PRINT *,'       -n n2-file n2-var : file and variable name for brunt vaisala frequency'
     PRINT * 
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out) 
     PRINT *,'       [-nc4] : use netcdf4 output with chunking and deflation 1'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : epx and epy , units : Nm^-2'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-u'    ) ; CALL getarg(ijarg, cf_ufil) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_u   ) ; ijarg=ijarg+1
     CASE ('-v'    ) ; CALL getarg(ijarg, cf_vfil) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_v   ) ; ijarg=ijarg+1
     CASE ('-w'    ) ; CALL getarg(ijarg, cf_wfil) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_w   ) ; ijarg=ijarg+1
     CASE ('-b'    ) ; CALL getarg(ijarg, cf_bfil) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_b   ) ; ijarg=ijarg+1
     CASE ('-n'    ) ; CALL getarg(ijarg, cf_n2fil) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_n2   ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4    = .TRUE.
     CASE ('-o'    ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  ENDDO


  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  lchk = chkfile(cf_wfil ) .OR. lchk
  lchk = chkfile(cf_bfil ) .OR. lchk
  lchk = chkfile(cf_n2fil ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  npiglo = getdim(cf_bfil,cn_x)
  npjglo = getdim(cf_bfil,cn_y)
  npk    = getdim(cf_bfil,cn_z)
  npt    = getdim(cf_bfil,cn_t) 

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt
  PRINT *, 'nlev   = ',nlev

  IF ( nlev==0 ) THEN
     PRINT *, 'nlev=0, assume 1'
     nlev=1
  END IF

  IF ( npt==0 ) THEN
     PRINT *, 'npt=0, assume 1'
     npt=1
  END IF
  ! 

  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( zun(npiglo,npjglo)  , zvn(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( wn(npiglo,npjglo)  , bn(npiglo,npjglo)  )
  ALLOCATE ( n2(npiglo,npjglo)  )
  ALLOCATE ( ff(npiglo,npjglo) )
  ALLOCATE ( depx(npiglo,npjglo) , depy(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( gdep(nlev) , zdep(npk))
  ALLOCATE ( nilev(nlev) )

  ! fills in gdep
  zdep(:) = getvar1d(cf_bfil, cn_vdeptht, npk )
  DO jk=1,nlev
     gdep(jk) = zdep(jk)
  ENDDO

  zun(:,:) = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
  zvn(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)

  CALL CreateOutput

  ff(:,:) =  getvar(cn_fhgr, 'ff',1, npiglo,npjglo)
  PRINT *, 'ff ok '

  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt
     DO jk = 1, nlev
        un(:,:) =  getvar(cf_ufil, cv_u, jk ,npiglo,npjglo, ktime=jt)
        vn(:,:) =  getvar(cf_vfil, cv_v, jk ,npiglo,npjglo, ktime=jt)
        wn(:,:) =  getvar(cf_wfil, cv_w, jk ,npiglo,npjglo, ktime=jt)
        bn(:,:) =  getvar(cf_bfil, cv_b, jk ,npiglo,npjglo, ktime=jt)
        n2(:,:) =  getvar(cf_n2fil, cv_n2, jk ,npiglo,npjglo, ktime=jt)

        depx(:,:) =  un * wn - ff * vn * bn / n2 
        depy(:,:) =  vn * wn + ff * un * bn / n2 
        ! write drotn on file at level k and at time jt
        ierr = putvar(ncout, id_varout(1), depx, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(2), depy, jk, npiglo, npjglo, ktime=jt)
     ENDDO
  END DO
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------

    ivar=0
    ! Zonal EP flux
    ivar                            = ivar+1
    ipk(ivar) = nlev  !  nlevel so far
    stypvar(ivar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname = 'EPx'
    stypvar(ivar)%cunits = 'Nm-2'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = -1000.
    stypvar(ivar)%valid_max         =  1000.
    stypvar(ivar)%clong_name        = 'Zonal Ep flux'
    stypvar(ivar)%cshort_name       = 'EPx'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TYX'
    ! Merdional EP flux
    ivar                            = ivar+1
    ipk(ivar) = nlev  !  nlevel so far
    stypvar(ivar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname = 'EPy'
    stypvar(ivar)%cunits = 'Nm-2'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = -1000.
    stypvar(ivar)%valid_max         =  1000.
    stypvar(ivar)%clong_name        = 'Meridional Ep flux'
    stypvar(ivar)%cshort_name       = 'EPy'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_bfil, npiglo, npjglo, nlev           , ld_nc4=lnc4)
    ierr  = createvar   (ncout , stypvar, 2,      ipk,    id_varout      , ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_bfil, npiglo, npjglo, nlev, pnavlon=zun, pnavlat=zvn, pdep=gdep)

    tim  = getvar1d(cf_bfil, cn_vtimec, npt      )
    ierr = putvar1d(ncout,   tim,       npt,  'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfEPflux


