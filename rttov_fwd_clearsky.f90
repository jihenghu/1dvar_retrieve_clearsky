
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing input cloud data
!!     the RTTOV optical depth coefficient file
!!     the RTTOV hydrotable file
!!
!!   The subroutine returns the retrieved emissivity  

subroutine rttov_fwd_clearsky( &
		imonth,nlevel,nchannel,incident, &
		channels,& 					!! [nchannel] channel list
		longitude,latitude,&        !!  Geolications
		plevel,phalf,&				!!  pressure levels in hpa
		vapor,&						!!  vapor mixing ratio in kg/kg
		temps,&						!!  air temperature profiles in K
		skt,&						!!  skin temperature in K
		psrf,&						!!  surface pressure in hpa
		t2m,&						!!  2m air temperature in K
		snowc,&						!!  snow cover [0~1]
		smc,&						!!  soil top layer moisture in m3/m3
		TBin,&					    !!  GMI top of atmos. brigtness temperature in K
		emiss_ret,&  				!!  to retrieve emissivities
		emiss_telsem)  				!!  1stguess emissivities


  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
		 surftype_sea,        &
         sensor_id_mw,        &
         sensor_id_po,		  &
         hydro_index_rain,    &
         hydro_index_snow,    &
         hydro_index_graupel, &
         hydro_index_clw,     &
         hydro_index_ciw

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_radiance,      &
		 rttov_transmission,  &
         rttov_chanprof,      &
         rttov_emissivity,	  &
         rttov_scatt_emis_retrieval_type

  ! The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY : &
        rttov_emis_atlas_data, &
        atlas_type_ir, atlas_type_mw

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit


  IMPLICIT NONE

  include "rttov_scatt.interface"
  include "rttov_parallel_scatt.interface"
  include "rttov_read_scattcoeffs.interface"
  include "rttov_dealloc_scattcoeffs.interface"
  include "rttov_scatt_setupindex.interface"

  include "rttov_read_coefs.interface"
  include "rttov_dealloc_coefs.interface"
  include "rttov_alloc_direct.interface"
  include "rttov_init_emis_refl.interface"
  include "rttov_print_opts_scatt.interface"
  include "rttov_print_profile.interface"
  include "rttov_print_cld_profile.interface"
  include "rttov_skipcommentline.interface"

!!!!!!!!!
  include "rttov_parallel_direct.interface"
  include "rttov_user_options_checkinput.interface"
  include "rttov_print_opts.interface"

! Use emissivity atlas
  include "rttov_setup_emis_atlas.interface"
  include "rttov_get_emis.interface"
  include "rttov_deallocate_emis_atlas.interface"

! Retrieve emissivity interface
  include "rttov_alloc_emis_ret_terms.interface"
  include "rttov_scatt_emis_terms.interface"
  include "rttov_scatt_emis_retrieval.interface"


!! ===========================================================================
  INTEGER ,intent(in) 	:: nchannel,nlevel

  REAL 	,intent(in)								:: incident
  INTEGER, dimension(nchannel) ,intent(in) 		:: channels
  real*4,intent(in) 			:: longitude,latitude
  real*4,intent(in)  		:: skt,psrf,t2m,snowc,smc
  real*4, dimension(nlevel) ,intent(in)			:: plevel,phalf
  real*4, dimension(nlevel),intent(in)  	:: vapor,temps
  real*4, dimension(nchannel),intent(in)	:: TBin
  real*4, dimension(nchannel),intent(out)	:: emiss_ret,emiss_telsem
!! ===========================================================================

  INTEGER(KIND=jpim), PARAMETER :: nhydro_frac = 1 ! number of hydrometeor fractions is one, i.e. a single profile of cloud cover

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                :: opts                     ! Options structure - leave this set to defaults
  TYPE(rttov_options_scatt)          :: opts_scatt               ! RTTOV-SCATT options structure
  TYPE(rttov_coefs)                  :: coefs                    ! Coefficients structure
  TYPE(rttov_scatt_coef)             :: coef_scatt               ! RTTOV-SCATT coefficients structure
  TYPE(rttov_chanprof),      POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  INTEGER(KIND=jpim),        POINTER :: frequencies(:) => NULL() ! Channel indexes for hydrotable lookup
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)  => NULL() ! Flags to specify channels to simulate
  LOGICAL(KIND=jplm),        POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),    POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),       POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)=> NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(rttov_radiance)               :: radiance                 ! Output radiances

  TYPE(rttov_emis_atlas_data)     	 :: emis_atlas               ! Data structure for emissivity atlas
  !!!   
  INTEGER(KIND=jpim) 			  	 :: err_alloc
  Type(rttov_scatt_emis_retrieval_type) :: emis_retrieval_terms  ! radiance terms for emissivity retrieve
  !!! 
  INTEGER(KIND=jpim)                 :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: atlas_type
  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=22)  :: NameOfRoutine = 'rttov_retrieve_gmi_emiss'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: hydrotable_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim),intent(in) :: imonth
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: ilev
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios


  REAL(kind=jprb), ALLOCATABLE  :: obs_tb (:)     ! Observed TB 
  REAL(kind=jprb), ALLOCATABLE  :: land_emis (:)  ! Retrieved emissivity

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV-SCATT are as follows:
  !   1. Specify required RTTOV-SCATT options
  !   2. Read coefficients and hydrotable file
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof and frequencies arrays by calling rttov_scatt_setupindex
  !   5. Read input profile(s)
  !   6. Set up surface emissivity
  !   7. Call rttov_scatt and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV-SCATT interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  hydrotable_filename='/home/jihenghu/rttov13/rtcoef_rttov13/hydrotable/hydrotable_gpm_gmi.dat'
  coef_filename='/home/jihenghu/rttov13/rtcoef_rttov13/rttov13pred54L/rtcoef_gpm_1_gmi.dat'
  

  
  nprof=1
  nlevels=nlevel
  nchannels=nchannel

  ALLOCATE(channel_list(nchannel))
  channel_list=channels!(/1,2,3,4,5,6,7,8,9,10/)
  nthreads=1
  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV-SCATT options structure
  ! --------------------------------------------------------------------------

  ! The rttov_options structure (opts) should be left with its default values.
  ! RTTOV-SCATT only allows access to a limited number of RTTOV options: these
  ! are set in the rttov_options_scatt structure (opts_scatt).

  ! For example:
  opts_scatt % interp_mode = 1                    ! Set interpolation method
  opts_scatt % config % verbose =.FALSE.          ! Enable printing of warnings
  opts_scatt % config % do_checkinput            =  .FALSE. !.TRUE. !
  ! See user guide for full list of RTTOV-SCATT options

  !========== Interactive inputs == end ==============
  !===================================================

  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)

  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Read the RTTOV-SCATT hydrotable file
  CALL rttov_read_scattcoeffs(errorstatus, opts_scatt, coefs, coef_scatt, file_coef=hydrotable_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading RTTOV-SCATT coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof
! print*,nchanprof
  !! Allocate tb_obs
  ALLOCATE(obs_tb(nchanprof))
  ALLOCATE(land_emis(nchanprof))

  ! Allocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,                 &
        1_jpim,                      &  ! 1 => allocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        coefs,                       &
        radiance = radiance,         &
        calcemis = calcemis,         &
        emissivity = emissivity,     &
        frequencies = frequencies,   &
        coef_scatt = coef_scatt,     &
        nhydro_frac = nhydro_frac,   &
        cld_profiles = cld_profiles, &
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  !! Allocate emis_retrieval_terms
	!! @param[out]    err                     status on exit
	!! @param[in]     nchanprof               size of the chanprof array (total number of channels being simulated)
	!! @param[in,out] emis_retrieval_terms    emissivity retrieval terms structure to allocate/deallocate
	!! @param[in]     asw                     1_jpim => allocate; 0_jpim => deallocate
  Call rttov_alloc_emis_ret_terms(err_alloc, nchanprof, emis_retrieval_terms, 1_jpim)
  IF (err_alloc /= errorstatus_success) THEN
	WRITE(*,*) 'allocation error for emis_retrieval_terms structures'
	CALL rttov_exit(errorstatus)
  ENDIF


  ! Initialise the RTTOV emissivity atlas
  ! (this loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases)
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. &
      coefs%coef%id_sensor == sensor_id_po) THEN
    atlas_type = atlas_type_mw ! MW atlas
  ELSE
    atlas_type = atlas_type_ir ! IR atlas
  ENDIF
  ! use TELSEM2 ATLAS
  CALL rttov_setup_emis_atlas(          &
              errorstatus,              &
              opts,                     &
              imonth,                   &
              atlas_type,               & ! Selects MW (1) or IR (2)
              emis_atlas,               &
              path = '/home/jihenghu/rttov13/emis_data', & ! The default path to atlas data
              coefs = coefs) ! This is mandatory for the CNRM MW atlas, ignored by TELSEM2;
                             ! if supplied for IR atlases they are initialised for this sensor
                             ! and this makes the atlas much faster to access
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error initialising emissivity atlas'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 4. Populate chanprof and frequencies arrays
  ! --------------------------------------------------------------------------

  ! RTTOV-SCATT requires the frequencies array to be populated by a call to
  ! rttov_scatt_setupindex. This also populates the chanprof array. To specify
  ! only a subset of channels (i.e. those in channel_list) an array of flags is
  ! passed in (use_chan).

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,coefs%coef%fmv_chn))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO j = 1, nprof
    use_chan(j,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
        errorstatus,        &
        nprof,              &
        coefs%coef%fmv_chn, &
        coefs,              &
        coef_scatt,         & 
        nchanprof,          &
        chanprof,           &
        frequencies,        &
        use_chan)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error finding channels, frequencies and polarisations'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

 
  ! Read gas units for profiles
  profiles(:) % gas_units =1

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

	joff = (iprof-1_jpim) * nchannels
	obs_tb(1+joff:nchannels+joff)=TBin
	
	profiles    (iprof) % p= plevel                         ! full level pressure (hPa)
	profiles    (iprof) % t= temps                       ! temperature (K)
	profiles    (iprof) % q= vapor                         ! specific humidity (ppmv or kg/kg - as read above)
   ! cfr,ciwc,clwc,crwc,cswc
    cld_profiles(iprof) % ph= phalf                        ! half level pressure (hPa)
	cld_profiles(iprof) % hydro_frac(:,1)=0.0_jprb             ! cloud cover (0-1)
	cld_profiles(iprof) % hydro(:,hydro_index_clw)= 0.0_jprb          ! cloud liquid water (kg/kg)
	cld_profiles(iprof) % hydro(:,hydro_index_ciw)= 0.0_jprb          ! cloud ice water (kg/kg)
	cld_profiles(iprof) % hydro(:,hydro_index_rain)= 0.0_jprb         ! rain (kg/kg)
	cld_profiles(iprof) % hydro(:,hydro_index_snow)= 0.0_jprb           ! snow (kg/kg)

    ! graupel (kg/kg) (not in data file)
    cld_profiles(iprof) % hydro(:,hydro_index_graupel) = 0.0_jprb

    ! 2 meter air variables
    profiles(iprof) % s2m % t=t2m
	profiles(iprof) % s2m % q= vapor(nlevels)
	profiles(iprof) % s2m % p= plevel(nlevels)
	profiles(iprof) % s2m % u= 0.0_jprb 
	profiles(iprof) % s2m % v= 0.0_jprb 
	
	where(profiles(iprof)%p(:)>psrf)profiles(iprof)%q=1.E-9
	where(profiles(iprof)%p(:)>psrf)cld_profiles(iprof)%hydro_frac(:,1)=0.0
	where(profiles(iprof)%p(:)>psrf)cld_profiles(iprof)%hydro(:,hydro_index_clw)=0.0
	where(profiles(iprof)%p(:)>psrf)cld_profiles(iprof)%hydro(:,hydro_index_ciw)=0.0
	where(profiles(iprof)%p(:)>psrf)cld_profiles(iprof)%hydro(:,hydro_index_rain)=0.0
	where(profiles(iprof)%p(:)>psrf)cld_profiles(iprof)%hydro(:,hydro_index_snow)=0.0
	
	
	! print*,profiles(iprof) % s2m % q
    ! The bottom-most half pressure level is taken as the 2m pressure
    cld_profiles(iprof) % ph(nlevels+1) = profiles(iprof) % s2m % p
    ! cld_profiles(iprof) % ph(nlevels+1) = plevel(nlevels)
	! where(cld_profiles(iprof)%ph>profiles(iprof)%s2m%p) cld_profiles(iprof)%ph=profiles(iprof)%s2m%p-0.1

    ! Skin variables
    profiles(iprof) % skin % t= skt
	profiles(iprof) % skin % salinity= 35.0  !! 3.5%, typical 'cause only valid over ocean
	profiles(iprof) % skin % fastem= (/3.0, 5.0, 15.0, 0.1, 0.3/)  !! actually no need, 'cause a ATLAS will be used.
    profiles(iprof) % skin % snow_fraction=snowc   !! actually no need, 'cause a ATLAS will be used.
    profiles(iprof) % skin % soil_moisture=smc     !! actually no need, 'cause a ATLAS will be used.

    ! Surface type and water type
    profiles(iprof) % skin % surftype= 0! land only
    profiles(iprof) % skin % watertype=1
   
    ! Elevation, latitude and longitude
    profiles(iprof) % elevation= 0
	profiles(iprof) % latitude= latitude
	profiles(iprof) % longitude= longitude

    ! Satellite angles
    profiles(iprof) % zenangle= incident!52.8
    profiles(iprof) % azangle = 0
    profiles(iprof) % sunzenangle= 0
    profiles(iprof) % sunazangle = 0

  ENDDO
  ! CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------
  ! Use emissivity atlas
  CALL rttov_get_emis(             &
            errorstatus,           &
            opts,                  &
            chanprof,              &
            profiles,              &
            coefs,                 &
            emis_atlas,            &
            emissivity(:) % emis_in)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error reading emissivity atlas'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Calculate emissivity within RTTOV where the atlas emissivity value is
  ! zero or less
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! Use the RTTOV emissivity and BRDF calculations over sea surfaces
  DO j = 1, SIZE(chanprof)
    IF (profiles(chanprof(j)%prof) % skin % surftype == surftype_sea) THEN
      calcemis(j) = .TRUE.
    ENDIF
  ENDDO

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV-SCATT forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_scatt ( &
        errorstatus,         				&! out   error flag
        opts_scatt,          				&! in    RTTOV-SCATT options structure
        nlevels,             				&! in    number of profile levels
        chanprof,            				&! in    channel and profile index structure
        frequencies,         				&! in    channel indexes for hydrotable lookup
        profiles,            				&! in    profile array
        cld_profiles,        				&! in    cloud/hydrometeor profile array
        coefs,               				&! in    coefficients structure
        coef_scatt,          				&! in    hydrotable structure
        calcemis,             &! in    flag for internal emissivity calcs
        emissivity,           &! inout input/output emissivities per channel
        radiance,            &! inout computed radiances
		emis_retrieval_terms=emis_retrieval_terms)             	
  ELSE
    CALL rttov_parallel_scatt ( &
        errorstatus,                    &! out   error flag
        opts_scatt,                     &! in    RTTOV-SCATT options structure
        nlevels,                        &! in    number of profile levels
        chanprof,                       &! in    channel and profile index structure
        frequencies,                    &! in    channel indexes for hydrotable lookup
        profiles,                       &! in    profile array
        cld_profiles,                   &! in    cloud/hydrometeor profile array
        coefs,                          &! in    coefficients structure
        coef_scatt,                     &! in    hydrotable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,         &! inout input/output emissivities per channel
        radiance,            &! inout computed radiances
		emis_retrieval_terms=emis_retrieval_terms, &
        nthreads = nthreads)  ! in    number of threads to use
  ENDIF
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_scatt error'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 8. Call Emssivity retrieval interface
  ! --------------------------------------------------------------------------
  !! @param[in]     chanprof      channels and profiles simulated by RTTOV-SCATT
  !! @param[in]     coefs         RTTOV coefficients structure
  !! @param[in]     emis_terms    output radiances and corresponding BTs
  !! @param[in]     obs_tb        observed BTs corresponding to simulated BTs
  !! @param[out]    land_emis     output retrieved emissivities
  call rttov_scatt_emis_retrieval(&
								  chanprof,             &
								  coefs,                &
								  emis_retrieval_terms, &
								  obs_tb,               &
								  land_emis)

  !=====================================================
  !============== Output results == start ==============

  DO iprof = 1, nprof
    joff = (iprof-1_jpim) * nchannels
    ! TBsimu(:)=radiance%bt((1+joff):(nchannels+joff))
    emiss_telsem=emissivity((1+joff):(nchannels+joff)) % emis_out
	emiss_ret= land_emis((1+joff):(nchannels+joff))
  ENDDO

  !============== Output results == end ==============
  !=====================================================

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE(use_chan, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,                 &
        0_jpim,                      &  ! 0 => deallocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        coefs,                       &
        radiance = radiance,         &
        calcemis = calcemis,         &
        emissivity = emissivity,     &
        frequencies = frequencies,   &
        coef_scatt = coef_scatt,     &
        nhydro_frac = nhydro_frac,   &
        cld_profiles = cld_profiles)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_scattcoeffs(coef_scatt)

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


  Call rttov_alloc_emis_ret_terms(err_alloc, nchanprof, emis_retrieval_terms, 0_jpim)
  IF (err_alloc /= errorstatus_success) THEN
	WRITE(*,*) 'deallocation error for emis_retrieval_terms structures'
	CALL rttov_exit(errorstatus)
  ENDIF
  
  Deallocate(channel_list)
  Deallocate(obs_tb)
  Deallocate(land_emis)

end subroutine rttov_fwd_clearsky
