!!  Adapted from RTTOV13.2 example_k.f90
!!  As a subroutine to compute the jacobian matrices of multichannel 
!!           TBs versus profiles, skt and emissivity.
!!
!!	By, Jiheng Hu, 2024/1/29, University of Michigan.

!! @Inputs
!! imonth,longitude,latitude,&
		! TBobs,Tatm,QWatm,LST,T2m,Snowc,Smc,SfcPress,EmissAnalyt,Emiss1st,
		
subroutine rttov_fwd_jacobian(nlevel,nchannel,incident,plevel,&
				vapor,temps,skt,t2m,Emissin,TBout,Emss_K,Ta_K,Qw_K,LST_K)

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         gas_unit_specconc

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

include "rttov_k.interface"
include "rttov_parallel_k.interface"
include "rttov_read_coefs.interface"
include "rttov_dealloc_coefs.interface"
include "rttov_alloc_k.interface"
include "rttov_init_emis_refl.interface"
include "rttov_init_prof.interface"
include "rttov_init_rad.interface"
include "rttov_init_transmission.interface"
include "rttov_user_options_checkinput.interface"
include "rttov_print_opts.interface"
include "rttov_print_profile.interface"
include "rttov_skipcommentline.interface"

!! ===========================================================================
  INTEGER ,intent(in) 								:: nchannel,nlevel
  REAL 	,intent(in)									:: incident
  real*4,intent(in)  								:: skt,t2m
  real*4, dimension(nlevel) ,intent(in)				:: plevel
  real*4, dimension(nlevel),intent(in)  			:: vapor,temps
  
  real*4, dimension(nchannel),intent(in)			:: Emissin
  real*4, dimension(nchannel)			:: Emissin_opt
  real*4, dimension(nchannel),intent(out)			:: TBout
  real*4, dimension(nchannel),intent(out)			:: Emss_K
  real*4, dimension(nchannel),intent(out)			:: LST_K
  real*4, dimension(nchannel,nlevel),intent(out) 	:: Ta_K,Qw_K
  INTEGER :: i
!! ===========================================================================



  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER	:: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER	:: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                       ! Options structure
  TYPE(rttov_coefs)                :: coefs                      ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)      => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)      => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)    => NULL() ! Input/output surface emissivity
  TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)  => NULL() ! Emissivity Jacobians
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)      => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:)   => NULL() ! Input/output surface BRDF
  TYPE(rttov_reflectance), POINTER :: reflectance_k(:) => NULL() ! Reflectance Jacobians
  TYPE(rttov_profile),     POINTER :: profiles(:)      => NULL() ! Input profiles
  TYPE(rttov_profile),     POINTER :: profiles_k(:)    => NULL() ! Output Jacobians
  TYPE(rttov_transmission)         :: transmission               ! Output transmittances
  TYPE(rttov_transmission)         :: transmission_k             ! Transmittance Jacobians
  TYPE(rttov_radiance)             :: radiance                   ! Output radiances
  TYPE(rttov_radiance)             :: radiance_k                 ! Radiance Jacobians

  INTEGER(KIND=jpim)               :: errorstatus                ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=9)  :: NameOfRoutine = 'rttov_fwd_jacobian'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  CHARACTER(LEN=11)  :: gas_unit
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch, l
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============
  coef_filename='/home/jihenghu/rttov13/rtcoef_rttov13/rttov13pred54L/rtcoef_gpm_1_gmi.dat'
  nprof=1
  nlevels=nlevel
  dosolar=0
  nchannels=nchannel
  ALLOCATE(channel_list(nchannels))
  channel_list=[(i,i=1,nchannel)]
  nthreads=1

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_all % ozone_data         = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_all % co2_data           = .FALSE. !   when supplying a profile of the
  opts % rt_all % n2o_data           = .FALSE. !   given trace gas (ensure the
  opts % rt_all % ch4_data           = .FALSE. !   coef file supports the gas)
  opts % rt_all % co_data            = .FALSE. !
  opts % rt_all % so2_data           = .FALSE. !
  opts % rt_mw  % clw_data           = .FALSE. !
  opts % rt_all % switchrad          = .TRUE.  ! Input K perturbation in BT

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings
  opts % config % do_checkinput      =  .FALSE. !.TRUE. !

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

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_k
  CALL rttov_alloc_k( &
        errorstatus,                 &
        1_jpim,                      &  ! 1 => allocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        profiles_k,                  &
        coefs,                       &
        transmission,                &
        transmission_k,              &
        radiance,                    &
        radiance_k,                  &
        calcemis=calcemis,           &
        emissivity=emissivity,       &
        emissivity_k=emissivity_k,   &
        calcrefl=calcrefl,           &
        reflectance=reflectance,     &
        reflectance_k=reflectance_k, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_k structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  ! Read gas units for profiles
  profiles(:) % gas_units =1

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg - as read above)
    profiles(iprof) % p =	plevel

    profiles(iprof) % t =	temps

    profiles(iprof) % q =	vapor

    ! 2 meter air variables
    profiles(iprof) % s2m % t 		  =	t2m
    profiles(iprof) % s2m % q 		  =	vapor(nlevels)
    profiles(iprof) % s2m % p 		  =	plevel(nlevels)
    profiles(iprof) % s2m % u 		  =	0.0_jprb 
    profiles(iprof) % s2m % v 		  =	0.0_jprb 

    ! Skin variables
    profiles(iprof) % skin % t        = skt
    profiles(iprof) % skin % salinity = 35.0  ! Salinity only applies to FASTEM over sea
    profiles(iprof) % skin % fastem   =  (/3.0, 5.0, 15.0, 0.1, 0.3/)      ! FASTEM only applies to MW instruments
    profiles(iprof) % skin % snow_fraction=0.     !! actually no need, 'cause a ATLAS will be used.
    profiles(iprof) % skin % soil_moisture=0.     !! actually no need, 'cause a ATLAS will be used.


    ! Surface type and water type
    profiles(iprof) % skin % surftype = 0 ! land only
    profiles(iprof) % skin % watertype= 1

    ! Elevation, latitude and longitude
    profiles(iprof) % elevation = 0
    profiles(iprof) % latitude  = 0.
    profiles(iprof) % longitude = 0.

    ! Satellite and solar angles
    profiles(iprof) % zenangle =incident
    profiles(iprof) % azangle =0.
    profiles(iprof) % sunzenangle =0.
    profiles(iprof) % sunazangle =0.

    ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
    profiles(iprof) % ctp = 0.
    profiles(iprof) % cfraction =0.0

  ENDDO

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities or reflectances
  ! so we initialise all inputs to zero
  CALL rttov_init_emis_refl(emissivity, reflectance)

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  ! calcemis(:) =(emissivity(:) % emis_in <= 0._jprb)
  calcemis(:) =.False.  
	Emissin_opt=Emissin
  ! where(Emissin.gt.1.0) Emissin_opt=1.0
  ! where(Emissin.lt.0.0) Emissin_opt=0.0
  emissivity(:) % emis_in = Emissin_opt
  ! Calculate reflectances within RTTOV where the input BRDF value is zero or
  ! less (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV K model
  ! --------------------------------------------------------------------------

  ! The input/output K variables must be initialised to zero before every call to rttov_k:

  ! Initialise RTTOV Jacobian structures to zero
  CALL rttov_init_prof(profiles_k(:))
  CALL rttov_init_rad(radiance_k)
  CALL rttov_init_transmission(transmission_k)
  CALL rttov_init_emis_refl(emissivity_k, reflectance_k)

  ! Set input perturbation in radiance_k:
  !   If switchrad is TRUE the perturbation in bt(:) is used for "thermal" channels
  !   If switchrad is FALSE the perturbation in total(:) is used
  !   For solar-only channels the perturbation in total(:) is *always* used
  !   It is harmless to specify inputs in both bt/total in any case: RTTOV will use the
  !     appropriate input for each channel
  radiance_k % total(:) = 1._jprb
  radiance_k % bt(:) = 1._jprb

  IF (nthreads <= 1) THEN
    CALL rttov_k(                          &
            errorstatus,                   &! out   error flag
            chanprof,                      &! in    channel and profile index structure
            opts,                          &! in    options structure
            profiles,                      &! in    profile array
            profiles_k,                    &! inout Jacobian array
            coefs,                         &! in    coefficients structure
            transmission,                  &! inout computed transmittances
            transmission_k,                &! inout transmittance Jacobians
            radiance,                      &! inout computed radiances
            radiance_k,                    &! inout input radiance/BT perturbation
            calcemis      = calcemis,      &! in    flag for internal emissivity calcs
            emissivity    = emissivity,    &! inout input/output emissivities per channel
            emissivity_k  = emissivity_k,  &! inout emissivity Jacobians
            calcrefl      = calcrefl,      &! in    flag for internal BRDF calcs
            reflectance   = reflectance,   &! inout input/output BRDFs per channel
            reflectance_k = reflectance_k)  ! inout BRDF Jacobians
  ELSE
    CALL rttov_parallel_k(                 &
            errorstatus,                   &! out   error flag
            chanprof,                      &! in    channel and profile index structure
            opts,                          &! in    options structure
            profiles,                      &! in    profile array
            profiles_k,                    &! inout Jacobian array
            coefs,                         &! in    coefficients structure
            transmission,                  &! inout computed transmittances
            transmission_k,                &! inout transmittance Jacobians
            radiance,                      &! inout computed radiances
            radiance_k,                    &! inout input radiance/BT perturbation
            calcemis      = calcemis,      &! in    flag for internal emissivity calcs
            emissivity    = emissivity,    &! inout input/output emissivities per channel
            emissivity_k  = emissivity_k,  &! inout emissivity Jacobians
            calcrefl      = calcrefl,      &! in    flag for internal BRDF calcs
            reflectance   = reflectance,   &! inout input/output BRDFs per channel
            reflectance_k = reflectance_k, &! inout BRDF Jacobians
            nthreads      = nthreads)       ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_k error'
    CALL rttov_exit(errorstatus)
  ENDIF

  !=====================================================
  !============== Output results == start ==============

  ! --- Output the DIRECT results -------------------------------------------

	! CALCULATED BRIGHTNESS TEMPERATURES (K):
	TBout = radiance % bt

	! Jacobians:'
	LST_K = profiles_k%skin%t   !!  K/K 
	Emss_K = emissivity_k%emis_out  !! K/1
	
	DO iprof=1,nchannel
		Ta_K(iprof,:)=profiles_k(iprof)%t !! K/K  
		Qw_K(iprof,:)=profiles_k(iprof)%q !! K/(kg/kg)
	END DO

  !============== Output results == end ==============
  !=====================================================

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_k
  CALL rttov_alloc_k( &
        errorstatus,                 &
        0_jpim,                      &  ! 0 => deallocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        profiles_k,                  &
        coefs,                       &
        transmission,                &
        transmission_k,              &
        radiance,                    &
        radiance_k,                  &
        calcemis=calcemis,           &
        emissivity=emissivity,       &
        emissivity_k=emissivity_k,   &
        calcrefl=calcrefl,           &
        reflectance=reflectance,     &
        reflectance_k=reflectance_k)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_k structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF

END SUBROUTINE rttov_fwd_jacobian
