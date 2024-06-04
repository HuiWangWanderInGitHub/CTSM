module VOCEmissionMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Volatile organic compound emission
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use clm_varctl         , only : iulog
  use clm_varpar         , only : maxveg, nlevcan
  use pftconMod          , only : ndllf_evr_tmp_tree,  ndllf_evr_brl_tree
  use pftconMod          , only : ndllf_dcd_brl_tree,  nbrdlf_evr_trp_tree
  use pftconMod          , only : nbrdlf_evr_tmp_tree, nbrdlf_dcd_brl_shrub
  use pftconMod          , only : nbrdlf_dcd_trp_tree, nbrdlf_dcd_tmp_tree
  use pftconMod          , only : nbrdlf_dcd_brl_tree, nbrdlf_evr_shrub
  use pftconMod          , only : nc3_arctic_grass   , nc3crop
  use pftconMod          , only : nc4_grass,           noveg
  use shr_megan_mod      , only : shr_megan_megcomps_n, shr_megan_megcomp_t, shr_megan_linkedlist
  use shr_megan_mod      , only : shr_megan_mechcomps_n, shr_megan_mechcomps, shr_megan_mapped_emisfctrs
  use MEGANFactorsMod    , only : Agro, Amat, Anew, Aold, betaT, ct1, ct2, LDF, Ceo
  use decompMod          , only : bounds_type
  use abortutils         , only : endrun
  use fileutils          , only : getfil
  use clm_varcon         , only : grlnd
  use atm2lndType        , only : atm2lnd_type
  use CanopyStateType    , only : canopystate_type
  use PhotosynthesisMod  , only : photosyns_type
  use WaterStateBulkType , only : waterstatebulk_type
  use SoilStateType      , only : soilstate_type
  use SolarAbsorbedType  , only : solarabs_type
  use TemperatureType    , only : temperature_type
  use PatchType          , only : patch                
  !*********************Stressed VOC*******************************************
  use shr_megan_mod      , only : shr_megan_ht_option,shr_megan_lt_option,shr_megan_hw_option
  use shr_megan_mod      , only : shr_megan_aq_option,shr_megan_sm_option
  use MEGANFactorsMod    , only : CAQ,CHT,CLT,CHW 
  use MEGANFactorsMod    , only : TAQ,THT,TLT,THW 
  use MEGANFactorsMod    , only : DTAQ,DTHT,DTLT,DTHW 
  use ColumnType         , only : col                
!2020-11-05,by Hui,Ref: Jiang et al. 2018.
  use EnergyFluxType  , only : energyflux_type

  !
  implicit none
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
  !
  ! !PUBLIC TYPES:
  type, public :: vocemis_type
     real(r8) , pointer, private :: Eopt_out_patch    (:)   ! Eopt coefficient
     real(r8) , pointer, private :: topt_out_patch    (:)   ! topt coefficient
     real(r8) , pointer, private :: alpha_out_patch   (:)   ! alpha coefficient
     real(r8) , pointer, private :: cp_out_patch      (:)   ! cp coefficient
     real(r8) , pointer, private :: paru_out_patch    (:)   ! 
     real(r8) , pointer, private :: par24u_out_patch  (:)   ! 
     real(r8) , pointer, private :: par240u_out_patch (:)   !  
     real(r8) , pointer, private :: para_out_patch    (:)   ! 
     real(r8) , pointer, private :: par24a_out_patch  (:)   ! 
     real(r8) , pointer, private :: par240a_out_patch (:)   ! 
     real(r8) , pointer, private :: gamma_out_patch   (:)   ! 
     real(r8) , pointer, private :: gammaL_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaT_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaP_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaA_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaS_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaC_out_patch  (:)   ! 
     real(r8) , pointer, private :: vocflx_tot_patch  (:)   ! total VOC flux into atmosphere [moles/m2/sec] 
     real(r8) , pointer, PUBLIC  :: vocflx_patch      (:,:) ! (num_mech_comps) MEGAN flux [moles/m2/sec] 
     real(r8) , pointer, private :: efisop_grc        (:,:) ! gridcell isoprene emission factors


     !************stressed VOC********************
     real(r8) , pointer, private :: gamma_ht_out_patch  (:)   ! 
     real(r8) , pointer, private :: gamma_lt_out_patch  (:)   ! 
     real(r8) , pointer, private :: gamma_hw_out_patch  (:)   ! 
     real(r8) , pointer, private :: gamma_aq_out_patch  (:)   ! 
     !high temperature lasting time
     real(r8) , pointer, private :: gamma_ht_all_patch  (:,:)   ! 
     real(r8) , pointer, private :: gamma_lt_all_patch  (:,:)   ! 
     real(r8) , pointer, private :: gamma_hw_all_patch  (:,:)   ! 
     real(r8) , pointer, private :: gamma_aq_all_patch  (:,:)   ! 
     !high temperature lasting time
     integer  , pointer, private :: ht_lst_time_patch       (:,:) 
     !low temperature lasting time
     integer  , pointer, private :: lt_lst_time_patch       (:,:) 
     !high wind speed lasting time
     integer  , pointer, private :: hw_lst_time_patch       (:,:)
 
     logical  , pointer, private :: ht_trigger_patch        (:,:) 
     logical  , pointer, private :: lt_trigger_patch        (:,:) 
     logical  , pointer, private :: hw_trigger_patch        (:,:) 


   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
  end type vocemis_type
  !
  ! !PRIVATE TYPES:
  type :: megan_out_type
     ! VOC fluxes structure for CLM history output
     real(r8), pointer, private  :: flux_out(:)   ! patch MEGAN flux [ug C m-2 h-1]
  end type megan_out_type
  type(megan_out_type), private, pointer :: meg_out(:) ! (n_megan_comps) points to output fluxes
  !
  logical, parameter :: debug = .false.
  !logical, parameter :: debug = .true.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    use clm_varctl     , only : use_fates, use_fates_sp
    class(vocemis_type) :: this
    type(bounds_type), intent(in)    :: bounds  

    if ( shr_megan_mechcomps_n > 0) then
       if ( use_fates .and. (.not. use_fates_sp) ) then
           call endrun( msg='ERROR: MEGAN currently does NOT work with FATES outside of FATES-SP mode (see github issue #115)'//&
                     errMsg(sourcefile, __LINE__))
       end if
       call this%InitAllocate(bounds) 
       call this%InitHistory(bounds)
       call this%InitCold(bounds)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate memory for module datatypes
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use shr_megan_mod   , only : shr_megan_factors_file
    use MEGANFactorsMod , only : megan_factors_init, megan_factors_get
    use clm_varpar      , only : mxpft
    !
    ! !ARGUMENTS:
    class(vocemis_type) :: this
    type(bounds_type)  , intent(in)  :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer            :: i, imeg
    integer            :: class_num
    real(r8)           :: factors(mxpft+1)
    real(r8)           :: molec_wght
    integer            :: begg, endg
    integer            :: begp, endp
    type(shr_megan_megcomp_t), pointer :: meg_cmp
    !-----------------------------------------------------------------------


    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endp = bounds%endp

    call megan_factors_init( shr_megan_factors_file )
    
    meg_cmp => shr_megan_linkedlist
    do while(associated(meg_cmp))
       allocate(meg_cmp%emis_factors(maxveg))
       call megan_factors_get( trim(meg_cmp%name), factors, class_num, molec_wght )
       meg_cmp%emis_factors(1:maxveg) = factors(1:maxveg)
       meg_cmp%class_number = class_num
       meg_cmp%molec_weight = molec_wght
       meg_cmp => meg_cmp%next_megcomp
    enddo

    allocate(this%Eopt_out_patch    (begp:endp)) ; this%EOPT_out_patch    (:)   = nan
    allocate(this%topt_out_patch    (begp:endp)) ; this%topt_out_patch    (:)   = nan
    allocate(this%alpha_out_patch   (begp:endp)) ; this%alpha_out_patch   (:)   = nan
    allocate(this%cp_out_patch      (begp:endp)) ; this%cp_out_patch      (:)   = nan   
    allocate(this%para_out_patch    (begp:endp)) ; this%para_out_patch    (:)   = nan   
    allocate(this%par24a_out_patch  (begp:endp)) ; this%par24a_out_patch  (:)   = nan 
    allocate(this%par240a_out_patch (begp:endp)) ; this%par240a_out_patch (:)   = nan
    allocate(this%paru_out_patch    (begp:endp)) ; this%paru_out_patch    (:)   = nan   
    allocate(this%par24u_out_patch  (begp:endp)) ; this%par24u_out_patch  (:)   = nan 
    allocate(this%par240u_out_patch (begp:endp)) ; this%par240u_out_patch (:)   = nan
    allocate(this%gamma_out_patch   (begp:endp)) ; this%gamma_out_patch   (:)   = nan
    allocate(this%gammaL_out_patch  (begp:endp)) ; this%gammaL_out_patch  (:)   = nan
    allocate(this%gammaT_out_patch  (begp:endp)) ; this%gammaT_out_patch  (:)   = nan
    allocate(this%gammaP_out_patch  (begp:endp)) ; this%gammaP_out_patch  (:)   = nan
    allocate(this%gammaA_out_patch  (begp:endp)) ; this%gammaA_out_patch  (:)   = nan
    allocate(this%gammaS_out_patch  (begp:endp)) ; this%gammaS_out_patch  (:)   = nan
    allocate(this%gammaC_out_patch  (begp:endp)) ; this%gammaC_out_patch  (:)   = nan

    allocate(this%vocflx_tot_patch  (begp:endp));  this%vocflx_tot_patch  (:)   = nan
    allocate(this%efisop_grc      (6,begg:endg));  this%efisop_grc        (:,:) = nan
    !stressed VOC
    allocate(this%gamma_ht_all_patch  (begp:endp,201)) ; this%gamma_ht_all_patch  (:,:)   = nan
    allocate(this%gamma_lt_all_patch  (begp:endp,201)) ; this%gamma_lt_all_patch  (:,:)   = nan
    allocate(this%gamma_hw_all_patch  (begp:endp,201)) ; this%gamma_hw_all_patch  (:,:)   = nan
    allocate(this%gamma_aq_all_patch  (begp:endp,201)) ; this%gamma_aq_all_patch  (:,:)   = nan
    
    allocate(this%gamma_ht_out_patch  (begp:endp)) ; this%gamma_ht_out_patch  (:)   = nan
    allocate(this%gamma_lt_out_patch  (begp:endp)) ; this%gamma_lt_out_patch  (:)   = nan
    allocate(this%gamma_hw_out_patch  (begp:endp)) ; this%gamma_hw_out_patch  (:)   = nan
    allocate(this%gamma_aq_out_patch  (begp:endp)) ; this%gamma_aq_out_patch  (:)   = nan
    
    allocate(this%ht_lst_time_patch  (begp:endp,201));  this%ht_lst_time_patch  (:,:)   = 0
    allocate(this%lt_lst_time_patch  (begp:endp,201));  this%lt_lst_time_patch  (:,:)   = 0
    allocate(this%hw_lst_time_patch  (begp:endp,201));  this%hw_lst_time_patch  (:,:)   = 0

    allocate(this%ht_trigger_patch   (begp:endp,201));  this%ht_trigger_patch  (:,:)   = .False.
    allocate(this%lt_trigger_patch   (begp:endp,201));  this%lt_trigger_patch  (:,:)   = .False.
    allocate(this%hw_trigger_patch   (begp:endp,201));  this%hw_trigger_patch  (:,:)   = .False.

    allocate(meg_out(shr_megan_megcomps_n)) 
    do i=1,shr_megan_megcomps_n
       allocate(meg_out(i)%flux_out(begp:endp))
       meg_out(i)%flux_out(:) = 0._r8
    end do

    allocate(this%vocflx_patch(begp:endp,1:shr_megan_mechcomps_n)) 
    this%vocflx_patch(:,1:shr_megan_mechcomps_n)= nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize history output fields for MEGAN emissions diagnositics
    !
    ! !USES 
    use clm_varcon  , only : spval
    use histFileMod , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(vocemis_type) :: this
    type(bounds_type), intent(in) :: bounds  
    real(r8), pointer :: ptr_1d(:)
    !
    ! !LOCAL VARIABLES
    integer :: imeg, ii
    integer :: begp, endp
    type(shr_megan_megcomp_t), pointer :: meg_cmp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    if (shr_megan_megcomps_n>0) then

       ! loop over megan compounds
       meg_cmp => shr_megan_linkedlist
       do while(associated(meg_cmp))
          imeg = meg_cmp%index

          call hist_addfld1d ( fname='MEG_'//trim(meg_cmp%name), units='kg/m2/sec',  &
               avgflag='A', long_name='MEGAN flux', &
               ptr_patch=meg_out(imeg)%flux_out, set_lake=0._r8, set_urb=0._r8 )

          meg_cmp => meg_cmp%next_megcomp
       enddo
       
       this%vocflx_tot_patch(begp:endp)= spval
       call hist_addfld1d (fname='VOCFLXT', units='moles/m2/sec',  &
            avgflag='A', long_name='total VOC flux into atmosphere', &
            ptr_patch=this%vocflx_tot_patch, set_lake=0._r8, set_urb=0._r8, default='inactive')

       this%gamma_out_patch(begp:endp)   = spval
       call hist_addfld1d (fname='GAMMA', units='non',  &
            avgflag='A', long_name='total gamma for VOC calc', &
            ptr_patch=this%gamma_out_patch, set_lake=0._r8, default='inactive')

       this%gammaL_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAL', units='non',  &
            avgflag='A', long_name='gamma L for VOC calc', &
            ptr_patch=this%gammaL_out_patch, set_lake=0._r8, default='inactive')

       this%gammaT_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAT', units='non',  &
            avgflag='A', long_name='gamma T for VOC calc', &
            ptr_patch=this%gammaT_out_patch, set_lake=0._r8, default='inactive')

       this%gammaP_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAP', units='non',  &
            avgflag='A', long_name='gamma P for VOC calc', &
            ptr_patch=this%gammaP_out_patch, set_lake=0._r8, default='inactive')

       this%gammaA_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAA', units='non',  &
            avgflag='A', long_name='gamma A for VOC calc', &
            ptr_patch=this%gammaA_out_patch, set_lake=0._r8, default='inactive')

       this%gammaS_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAS', units='non',  &
            avgflag='A', long_name='gamma S for VOC calc', &
            ptr_patch=this%gammaS_out_patch, set_lake=0._r8, default='inactive')

       this%gammaC_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMAC', units='non',  &
            avgflag='A', long_name='gamma C for VOC calc', &
            ptr_patch=this%gammaC_out_patch, set_lake=0._r8, default='inactive')

       this%EOPT_out_patch(begp:endp) = spval
       call hist_addfld1d (fname='EOPT', units='non',  &
            avgflag='A', long_name='Eopt coefficient for VOC calc', &
            ptr_patch=this%Eopt_out_patch, set_lake=0._r8, default='inactive')

       this%topt_out_patch(begp:endp)    = spval
       call hist_addfld1d (fname='TOPT', units='non',  &
            avgflag='A', long_name='topt coefficient for VOC calc', &
            ptr_patch=this%topt_out_patch, set_lake=0._r8, default='inactive')

       this%alpha_out_patch(begp:endp)    = spval
       call hist_addfld1d (fname='ALPHA', units='non',  &
            avgflag='A', long_name='alpha coefficient for VOC calc', &
            ptr_patch=this%alpha_out_patch, set_lake=0._r8, default='inactive')

       this%cp_out_patch(begp:endp)      = spval   
       call hist_addfld1d (fname='currentPatch', units='non',  &
            avgflag='A', long_name='currentPatch coefficient for VOC calc', &
            ptr_patch=this%cp_out_patch, set_lake=0._r8, default='inactive')

       this%paru_out_patch(begp:endp)    = spval   
       call hist_addfld1d (fname='PAR_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR', &
            ptr_patch=this%paru_out_patch, set_lake=0._r8, default='inactive')

       this%par24u_out_patch(begp:endp)  = spval 
       call hist_addfld1d (fname='PAR24_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (24 hrs)', &
            ptr_patch=this%par24u_out_patch, set_lake=0._r8, default='inactive')

       this%par240u_out_patch(begp:endp) = spval
       call hist_addfld1d (fname='PAR240_sun', units='umol/m2/s', &
            avgflag='A', long_name='sunlit PAR (240 hrs)', &
            ptr_patch=this%par240u_out_patch, set_lake=0._r8, default='inactive')

       this%para_out_patch(begp:endp)    = spval   
       call hist_addfld1d (fname='PAR_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR', &
            ptr_patch=this%para_out_patch, set_lake=0._r8, default='inactive')

       this%par24a_out_patch(begp:endp)  = spval 
       call hist_addfld1d (fname='PAR24_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (24 hrs)', &
            ptr_patch=this%par24a_out_patch, set_lake=0._r8, default='inactive')

       this%par240a_out_patch(begp:endp) = spval
       call hist_addfld1d (fname='PAR240_shade', units='umol/m2/s', &
            avgflag='A', long_name='shade PAR (240 hrs)', &
            ptr_patch=this%par240a_out_patch, set_lake=0._r8, default='inactive')
       !stressed VOC
       this%gamma_ht_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMA_ht', units='non',  &
            avgflag='A', long_name='gamma ht for VOC calc', &
            ptr_patch=this%gamma_ht_out_patch, set_lake=0._r8, default='inactive')
       this%gamma_lt_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMA_lt', units='non',  &
            avgflag='A', long_name='gamma lt for VOC calc', &
            ptr_patch=this%gamma_lt_out_patch, set_lake=0._r8, default='inactive')
       this%gamma_hw_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMA_hw', units='non',  &
            avgflag='A', long_name='gamma hw for VOC calc', &
            ptr_patch=this%gamma_hw_out_patch, set_lake=0._r8, default='inactive')
       this%gamma_aq_out_patch(begp:endp)  = spval
       call hist_addfld1d (fname='GAMMA_aq', units='non',  &
            avgflag='A', long_name='gamma aq for VOC calc', &
            ptr_patch=this%gamma_aq_out_patch, set_lake=0._r8, default='inactive')

    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES
    use ncdio_pio
    use clm_varctl, only : fsurdat
    !
    ! !ARGUMENTS:
    class(vocemis_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar 
    integer            :: begg, endg
    type(file_desc_t)  :: ncid       ! netcdf id
    character(len=256) :: locfn      ! local filename
    real(r8) ,pointer  :: temp_ef(:) ! read in - temporary EFs 
    !-----------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    ! Time constant

    allocate(temp_ef(begg:endg))

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_io(ncid=ncid, varname='EF1_BTR', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_BTR'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(1,begg:endg)=temp_ef(begg:endg)

    call ncd_io(ncid=ncid, varname='EF1_FET', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_FET'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(2,begg:endg)=temp_ef(begg:endg)

    call ncd_io(ncid=ncid, varname='EF1_FDT', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_FDT'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(3,begg:endg)=temp_ef(begg:endg)

    call ncd_io(ncid=ncid, varname='EF1_SHR', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_SHR'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(4,begg:endg)=temp_ef(begg:endg)

    call ncd_io(ncid=ncid, varname='EF1_GRS', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_GRS'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(5,begg:endg)=temp_ef(begg:endg)

    call ncd_io(ncid=ncid, varname='EF1_CRP', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg='iniTimeConst: errror reading EF1_CRP'//errMsg(sourcefile, __LINE__))
    end if
    this%efisop_grc(6,begg:endg)=temp_ef(begg:endg)

    deallocate(temp_ef)

    call ncd_pio_closefile(ncid)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine VOCEmission (bounds, num_soilp, filter_soilp, &
       atm2lnd_inst, canopystate_inst, photosyns_inst, temperature_inst, &
       vocemis_inst, energyflux_inst, soilstate_inst, waterstatebulk_inst)
    !
    ! ! NEW DESCRIPTION
    ! Volatile organic compound emission
    ! This code simulates volatile organic compound emissions following
    ! MEGAN (Model of Emissions of Gases and Aerosols from Nature) v2.1 
    ! for 20 compound classes. The original description of this
    ! algorithm (for isoprene only) can be found in Guenther et al., 2006
    ! (we follow equations 2-9, 16-17, 20 for explicit canopy).
    ! The model scheme came be described as:
    !    E= epsilon * gamma * rho
    ! VOC flux (E) [ug m-2 h-1] is calculated from baseline emission
    ! factors (epsilon) [ug m-2 h-1] which are specified for each of the 16
    ! CLM Patches (in input file) OR in the case of isoprene, from
    ! mapped EFs for each PATCH which reflect species divergence of emissions,
    ! particularly in North America. 
    ! The emission activity factor (gamma) [unitless] for includes 
    ! dependence on PPFT, temperature, LAI, leaf age and soil moisture.
    ! For isoprene only we also include the effect of CO2 inhibition as
    ! described by Heald et al., 2009. 
    ! The canopy environment constant was calculated offline for CLM+CAM at 
    ! standard conditions.
    ! We assume that the escape efficiency (rho) here is unity following
    ! Guenther et al., 2006.
    ! A manuscript describing MEGAN 2.1 and the implementation in CLM is
    ! in preparation: Guenther, Heald et al., 2012
    ! Subroutine written to operate at the patch level.
    !
    ! Input: <filename> to be read in with EFs and some parameters.  
    !        Currently these are set in procedure init_EF_params
    ! Output: vocflx(shr_megan_mechcomps_n) !VOC flux [moles/m2/sec]
    !
    ! !USES:
    use subgridAveMod        , only : p2g
    use shr_const_mod        , only : SHR_CONST_CDAY
    use clm_time_manager   , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                  
    integer                , intent(in)    :: num_soilp               ! number of columns in soil patch filter
    integer                , intent(in)    :: filter_soilp(num_soilp) ! patch filter for soil
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(photosyns_type)   , intent(in)    :: photosyns_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(vocemis_type)     , intent(inout) :: vocemis_inst
    !by Hui,
    type(energyflux_type)  , intent(in)    :: energyflux_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    !
    ! !REVISION HISTORY:
    ! 4/29/11: Colette L. Heald: expand MEGAN to 20 compound classes
    ! 7 Feb 2012: Francis Vitt: Implemented capability to specify MEGAN emissions in namelist
    !                           and read in MEGAN factors from file.
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,p,g,c                ! indices
    real(r8) :: epsilon                 ! emission factor [ug m-2 h-1]
    real(r8) :: gamma                   ! activity factor (accounting for light, T, age, LAI conditions)
    real(r8) :: gamma_p                 ! activity factor for PPFD
    real(r8) :: gamma_l                 ! activity factor for LAI (just LAI in MEGAN v3)
    real(r8) :: gamma_t                 ! activity factor for temperature
    real(r8) :: gamma_tp                ! activity factor for temperature and light (MEGAN v3 )
    real(r8) :: gamma_tld               ! activity factor for temperature (light dependent)
    real(r8) :: gamma_tli               ! activity factor for temperature (light independent)
    real(r8) :: gamma_a                 ! activity factor for leaf age
    real(r8) :: gamma_sm                ! activity factor for soil moisture
    real(r8) :: gamma_c                 ! activity factor for CO2 (only isoprene)
    real(r8) :: par_sun                 ! temporary
    real(r8) :: par24_sun               ! temporary
    real(r8) :: par240_sun              ! temporary
    real(r8) :: par_sha                 ! temporary
    real(r8) :: par24_sha               ! temporary
    real(r8) :: par240_sha              ! temporary
    !**********************stressed VOC**************************
    real(r8) :: gamma_ht
    real(r8) :: gamma_lt
    real(r8) :: gamma_hw
    real(r8) :: gamma_aq
    integer  :: stress_period
    ! value, obtained from ATM
    real(r8), parameter :: forc_ozone = 100._r8 * 1.e-9_r8  ! ozone partial pressure [mol/mol]
 
    integer                            :: class_num, n_meg_comps, imech, imeg, ii
    character(len=16)                  :: mech_name
    type(shr_megan_megcomp_t), pointer :: meg_cmp
    real(r8)                           :: cp, alpha,  Eopt, topt  ! for history output
    real(r8)                           :: co2_ppmv

    real(r8)                           :: vocflx_meg(shr_megan_megcomps_n)

    ! factor used convert MEGAN units [micro-gra/m2/hr] to CAM srf emis units [g/m2/sec]
    real(r8), parameter :: megemis_units_factor_old = 1._r8/3600._r8/1.e6_r8
    
    ! factor used convert MEGAN units [nmol/m2/sec] to [mol/m2/sec]
    ! for CAM srf emis units [g/m2/sec]
    real(r8), parameter :: megemis_units_factor_new = 1.e-9_r8

    character(len=32), parameter :: subname = "VOCEmission"
    !-----------------------------------------------------------------------
     real(r8) :: root_depth(0:maxveg)    ! Root depth [m]

        ! root depth (m) (defined based on Zeng et al., 2001, cf Guenther 2006)
        root_depth(noveg)                                     = 0._r8   ! bare-soil
        root_depth(ndllf_evr_tmp_tree:ndllf_evr_brl_tree)     = 1.8_r8  ! evergreen tree
        root_depth(ndllf_dcd_brl_tree)                        = 2.0_r8  ! needleleaf deciduous boreal tree
        root_depth(nbrdlf_evr_trp_tree:nbrdlf_evr_tmp_tree)   = 3.0_r8  ! broadleaf evergreen tree
        root_depth(nbrdlf_dcd_trp_tree:nbrdlf_dcd_brl_tree)   = 2.0_r8  ! broadleaf deciduous tree
        root_depth(nbrdlf_evr_shrub:nbrdlf_dcd_brl_shrub)     = 2.5_r8  ! shrub
        root_depth(nc3_arctic_grass:maxveg)                   = 1.5_r8  ! grass/crop
    !
    if ( shr_megan_mechcomps_n < 1) return

    if ( nlevcan /= 1 )then
       call endrun( subname//' error: can NOT work without nlevcan == 1' )
    end if

    associate(                                                    & 
         btran_min     => energyflux_inst%btran_min_patch           , & ! Input: [real(r8) (:)    ]  transpiration wetness factor (0 to 1)
         btran         => energyflux_inst%btran_patch           , & ! Input: [real(r8) (:)    ]  transpiration wetness factor (0 to 1)
         vcmax_z       => photosyns_inst%vcmax_z_patch          , & ! Input: [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)

         !************************************************************************
         !*************************for old drought algorithm**********************
         !************************************************************************
         dz           => col%dz                                , & ! Input:  [real(r8) (:,:) ]  depth of layer (m)                              
         bsw          => soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b" (nlevgrnd)             
         clayfrac     => soilstate_inst%clayfrac_patch           , & ! Input:  [real(r8) (:)   ]  fraction of soil that is clay                     
         sandfrac     => soilstate_inst%sandfrac_patch           , & ! Input:  [real(r8) (:)   ]  fraction of soil that is sand                     
         watsat       => soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity) (nlevgrnd)
         sucsat       => soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm) (nlevgrnd)            
         h2osoi_vol   => waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (m3/m3)                   
         h2osoi_ice   => waterstatebulk_inst%h2osoi_ice_col        , & ! Input:  [real(r8) (:,:) ]  ice soil content (kg/m3)                        
         !**********************************************************************
         !**********************************************************************
         
         forc_solad    => atm2lnd_inst%forc_solad_downscaled_col, & ! Input:  [real(r8) (:,:) ]  direct beam radiation (visible only)            
         forc_solai    => atm2lnd_inst%forc_solai_grc           , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (visible only)            
         forc_pbot     => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  downscaled atmospheric pressure (Pa)                          
         forc_pco2     => atm2lnd_inst%forc_pco2_grc            , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             
         forc_solad24  => atm2lnd_inst%fsd24_patch              , & ! Input:  [real(r8) (:)   ]  direct beam radiation last 24hrs  (visible only)  
         forc_solad240 => atm2lnd_inst%fsd240_patch             , & ! Input:  [real(r8) (:)   ]  direct beam radiation last 240hrs (visible only)  
         forc_solai24  => atm2lnd_inst%fsi24_patch              , & ! Input:  [real(r8) (:)   ]  diffuse radiation  last 24hrs     (visible only)  
         forc_solai240 => atm2lnd_inst%fsi240_patch             , & ! Input:  [real(r8) (:)   ]  diffuse radiation  last 240hrs    (visible only)  

         fsun          => canopystate_inst%fsun_patch           , & ! Input:  [real(r8) (:)   ]  sunlit fraction of canopy                         
         fsun24        => canopystate_inst%fsun24_patch         , & ! Input:  [real(r8) (:)   ]  sunlit fraction of canopy last 24 hrs             
         fsun240       => canopystate_inst%fsun240_patch        , & ! Input:  [real(r8) (:)   ]  sunlit fraction of canopy last 240 hrs            
         elai          => canopystate_inst%elai_patch           , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         elai240       => canopystate_inst%elai240_patch        , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow last 240 hrs

         cisun_z       => photosyns_inst%cisun_z_patch          , & ! Input:  [real(r8) (:,:) ]  sunlit intracellular CO2 (Pa)
         cisha_z       => photosyns_inst%cisha_z_patch          , & ! Input:  [real(r8) (:,:) ]  shaded intracellular CO2 (Pa)
         
         t_veg         => temperature_inst%t_veg_patch          , & ! Input:  [real(r8) (:)   ]  patch vegetation temperature (Kelvin)
         t_veg24       => temperature_inst%t_veg24_patch        , & ! Input:  [real(r8) (:)   ]  avg patch vegetation temperature for last 24 hrs
         t_veg240      => temperature_inst%t_veg240_patch       , & ! Input:  [real(r8) (:)   ]  avg patch vegetation temperature for last 240 hrs
         
         Eopt_out      => vocemis_inst%Eopt_out_patch           , & ! Output: [real(r8) (:)   ]                                                    
         topt_out      => vocemis_inst%topt_out_patch           , & ! Output: [real(r8) (:)   ]                                                    
         alpha_out     => vocemis_inst%alpha_out_patch          , & ! Output: [real(r8) (:)   ]                                                    
         cp_out        => vocemis_inst%cp_out_patch             , & ! Output: [real(r8) (:)   ]                                                    
         paru_out      => vocemis_inst%paru_out_patch           , & ! Output: [real(r8) (:)   ]                                                    
         par24u_out    => vocemis_inst%par24u_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         par240u_out   => vocemis_inst%par240u_out_patch        , & ! Output: [real(r8) (:)   ]                                                    
         para_out      => vocemis_inst%para_out_patch           , & ! Output: [real(r8) (:)   ]                                                    
         par24a_out    => vocemis_inst%par24a_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         par240a_out   => vocemis_inst%par240a_out_patch        , & ! Output: [real(r8) (:)   ]                                                    
         gammaL_out    => vocemis_inst%gammaL_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gammaT_out    => vocemis_inst%gammaT_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gammaP_out    => vocemis_inst%gammaP_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gammaA_out    => vocemis_inst%gammaA_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gammaS_out    => vocemis_inst%gammaS_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gammaC_out    => vocemis_inst%gammaC_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_out     => vocemis_inst%gamma_out_patch          , & ! Output: [real(r8) (:)   ]                                                    
         vocflx        => vocemis_inst%vocflx_patch             , & ! Output: [real(r8) (:,:) ]  VOC flux [moles/m2/sec]
         vocflx_tot    => vocemis_inst%vocflx_tot_patch         , & ! Output: [real(r8) (:)   ]  VOC flux [moles/m2/sec]
         ! stressed VOC
         forc_wind     => atm2lnd_inst%forc_wind_grc           , & ! Input:  [real(r8) (:,:) ]  atmospheric wind speed (m/s)            
         !forc_wind_grc
         ! stress lasting time
         ht_lst_time   => vocemis_inst%ht_lst_time_patch       , &
         lt_lst_time   => vocemis_inst%lt_lst_time_patch       , &
         hw_lst_time   => vocemis_inst%hw_lst_time_patch       , &
         ht_trigger    => vocemis_inst%ht_trigger_patch        , &
         lt_trigger    => vocemis_inst%lt_trigger_patch        , &
         hw_trigger    => vocemis_inst%hw_trigger_patch        , &
         ! stressed VOC out
         gamma_ht_out    => vocemis_inst%gamma_ht_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_lt_out    => vocemis_inst%gamma_lt_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_hw_out    => vocemis_inst%gamma_hw_out_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_aq_out    => vocemis_inst%gamma_aq_out_patch         , & ! Output: [real(r8) (:)   ]                                                    

         gamma_ht_all    => vocemis_inst%gamma_ht_all_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_lt_all    => vocemis_inst%gamma_lt_all_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_hw_all    => vocemis_inst%gamma_hw_all_patch         , & ! Output: [real(r8) (:)   ]                                                    
         gamma_aq_all    => vocemis_inst%gamma_aq_all_patch          & ! Output: [real(r8) (:)   ]                                                    
         )

    ! initialize variables which get passed to the atmosphere
    vocflx(bounds%begp:bounds%endp,:)   = 0._r8
    vocflx_tot(bounds%begp:bounds%endp) = 0._r8

    ! stressed period
    stress_period = nint(SHR_CONST_CDAY) / get_step_size()

    do imeg=1,shr_megan_megcomps_n
      meg_out(imeg)%flux_out(bounds%begp:bounds%endp) = 0._r8
    enddo
          
    ! Begin loop over points
    !_______________________________________________________________________________
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       g = patch%gridcell(p)
       c = patch%column(p)

       ! initialize EF
       epsilon=0._r8
       
       ! initalize to zero since this might not alway get set
       ! this needs to be within the fp loop ... 
       vocflx_meg(:) = 0._r8

       ! calculate VOC emissions for non-bare ground Patches
       if (patch%itype(p) > 0) then 
          gamma=0._r8

          ! Calculate PAR: multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)
          !------------------------
          ! SUN:
          par_sun    = (forc_solad(c,1)  + fsun(p)    * forc_solai(g,1))  * 4.6_r8
          par24_sun  = (forc_solad24(p)  + fsun24(p)  * forc_solai24(p))  * 4.6_r8
          par240_sun = (forc_solad240(p) + fsun240(p) * forc_solai240(p)) * 4.6_r8

          ! SHADE:
          par_sha    = ((1._r8 - fsun(p))    * forc_solai(g,1))  * 4.6_r8
          par24_sha  = ((1._r8 - fsun24(p))  * forc_solai24(p))  * 4.6_r8
          par240_sha = ((1._r8 - fsun240(p)) * forc_solai240(p)) * 4.6_r8

          ! Activity factor for LAI (MEGAN v3) 
          gamma_l =  get_gamma_L(fsun240(p),elai(p))

          

          ! Activity factor for soil moisture: all species (commented out for now)
          !**********Soil Moisture Stress******************8
          if (shr_megan_sm_option .eq. 0) then
                gamma_sm = 1.0_r8 
          else if (shr_megan_sm_option .eq. 1) then
                    gamma_sm = get_gamma_SM_test(btran(p))
          else if (shr_megan_sm_option .eq. 3) then
                    gamma_sm = get_gamma_SM_test2(btran(p))
          else if (shr_megan_sm_option .eq. 4) then
                    gamma_sm = get_gamma_SM_test3(btran(p))
          else if (shr_megan_sm_option .eq. 2) then

             if (btran(p)>=0.6) then
                gamma_sm = 1.0_r8
             else
                gamma_sm = vcmax_z(p,1)/37
             endif
             if(gamma_sm>1.0)then
                gamma_sm = 1.0_r8
             endif 
          
          else
                gamma_sm = 1.0_r8 
          end if


          ! Loop through VOCs for light, temperature and leaf age activity factor & apply
          ! all final activity factors to baseline emission factors
          !_______________________________________________________________________________

          ! loop over megan compounds
          meg_cmp => shr_megan_linkedlist
          meg_cmp_loop: do while(associated(meg_cmp))
             imeg = meg_cmp%index

             ! set emis factor
             ! if specified, set EF for isoprene with mapped values
             if ( trim(meg_cmp%name) == 'isoprene' .and. shr_megan_mapped_emisfctrs) then
                epsilon = get_map_EF(patch%itype(p),g, vocemis_inst)
             else
                epsilon = meg_cmp%emis_factors(patch%itype(p))
             end if

             class_num = meg_cmp%class_number

             ! Activity factor for PPFD
             !gamma_p = get_gamma_P(par_sun, par_sha, fsun(p), cp, alpha)
             gamma_p = get_gamma_P(par_sun, par24_sun, par240_sun, par_sha, par24_sha, par240_sha, &
                  fsun(p), fsun240(p), forc_solad240(p),forc_solai240(p), LDF(class_num), cp, alpha)

             ! Activity factor for T
             gamma_tld = get_gamma_tld(t_veg240(p), t_veg24(p),t_veg(p), ct1(class_num), ct2(class_num),&
                                  Ceo(class_num), Eopt, topt)
             gamma_tli = get_gamma_tli(t_veg(p), betaT(class_num))
             
             !gamma_t = LDF(class_num)*gamma_tld + (1. - LDF(class_num))*gamma_tli
             !gamma_tp = LDF(class_num)*gamma_tld*gamma_p + (1. - LDF(class_num))*gamma_tli
             gamma_t = get_gamma_T(t_veg240(p), t_veg24(p),t_veg(p), ct1(class_num), ct2(class_num),&
                                   betaT(class_num),LDF(class_num), Ceo(class_num), Eopt, topt)

             ! Activity factor for Leaf Age
             gamma_a = get_gamma_A(patch%itype(p), elai240(p),elai(p),class_num)

             ! Activity factor for CO2 (only for isoprene)
             if (trim(meg_cmp%name) == 'isoprene') then 
                co2_ppmv = 1.e6_r8*forc_pco2(g)/forc_pbot(c)
                gamma_c = get_gamma_C(cisun_z(p,1),cisha_z(p,1),forc_pbot(c),fsun(p), co2_ppmv)
             else
                gamma_c = 1._r8
             end if
             !stresssed VOC
             !**********High Temperature Stress*****************
             if (shr_megan_ht_option) then
                  gamma_ht = get_gamma_ht(t_veg(p),CHT(class_num),&
                             THT(class_num),DTHT(class_num)) 
                  if( ht_trigger(p,imeg) )then
                     if( gamma_ht .gt. gamma_ht_all(p,imeg))then
                        !update ht_lst_time
                        ht_lst_time(p,imeg) = 0
                     else
                        gamma_ht = gamma_ht_all(p,imeg) 
                     endif

                     if( ht_lst_time(p,imeg) .le. stress_period ) then! nstep <= 24h
                        !Don't change ht_lst_time
                        ht_lst_time(p,imeg) = ht_lst_time(p,imeg) + 1
                     else
                        ht_lst_time(p,imeg) = 0
                        ht_trigger(p,imeg) = .False.
                     end if
                  else
                     if (gamma_ht .gt. 1._r8)then
                        ht_trigger(p,imeg) = .True.
                        ht_lst_time(p,imeg) = 1
                     end if
                  end if
             else
                  gamma_ht = 1.0_r8
             end if 
             !**********Low Temperature Stress******************8
             if (shr_megan_lt_option) then
                  gamma_lt = get_gamma_lt(t_veg(p),CLT(class_num),&
                             TLT(class_num),DTLT(class_num))
                  if( lt_trigger(p,imeg) )then
                     if( gamma_lt .gt. gamma_lt_all(p,imeg))then
                        !update lt_lst_time
                        lt_lst_time(p,imeg) = 0 
                     else
                        gamma_lt = gamma_lt_all(p,imeg)
                     endif

                     if( lt_lst_time(p,imeg) .le. stress_period ) then! nstep <= 24h
                        !Don't change lt_lst_time
                        lt_lst_time(p,imeg) = lt_lst_time(p,imeg) + 1
                     else
                        lt_lst_time(p,imeg) = 0
                        lt_trigger(p,imeg) = .False.
                     end if
                  else
                     if(gamma_lt .gt. 1._r8)then
                        lt_trigger(p,imeg) = .True.
                        lt_lst_time(p,imeg) = 1
                     end if
                  end if
             else
                  gamma_lt = 1.0_r8
             end if 
             !**********High Windspeed Stress******************8
             if (shr_megan_hw_option) then
                  gamma_hw = get_gamma_hw(forc_wind(p),CHW(class_num),&
                             THW(class_num),DTHW(class_num))
                  if( hw_trigger(p,imeg) )then
                     if( gamma_hw .gt. gamma_hw_all(p,imeg))then
                        !update hw_lst_time
                        hw_lst_time(p,imeg) = 0 
                     else
                        gamma_hw = gamma_hw_all(p,imeg)
                     endif

                     if( hw_lst_time(p,imeg) .le. stress_period ) then! nstep <= 24h
                       !Don't change ht_lst_time
                       hw_lst_time(p,imeg) = hw_lst_time(p,imeg) + 1
                     else
                       hw_lst_time(p,imeg) = 0
                       hw_trigger(p,imeg) = .False.
                     end if
                  else
                     if (gamma_hw .gt. 1._r8)then
                     hw_trigger(p,imeg) = .True.
                     hw_lst_time(p,imeg) = 1
                     end if
                  end if
             else
                  gamma_hw = 1.0_r8
             end if 
             !**********Air Quality (ozone) Stress******************8
             if (shr_megan_aq_option) then
                  gamma_aq = get_gamma_aq(forc_ozone,CAQ(class_num),&
                             TAQ(class_num),DTAQ(class_num))
             else
                  gamma_aq = 1.0_r8
             end if 

             ! Calculate total scaling factor
             if ( trim(meg_cmp%name) == 'isoprene') then
             gamma = gamma_l * gamma_sm * gamma_a * gamma_p * gamma_t * gamma_c * &
                     gamma_ht*gamma_lt*gamma_hw*gamma_aq
             ! gamma_tp replace gamma_T*gamma_p in MEGAN v3
             !gamma = gamma_l * gamma_sm * gamma_a * gamma_tp * gamma_c * &
             !        gamma_ht*gamma_lt*gamma_hw*gamma_aq
             else
             gamma = gamma_l * gamma_a * gamma_p * gamma_T * gamma_c * &
                     gamma_ht*gamma_lt*gamma_hw*gamma_aq
             !gamma = gamma_l * gamma_a * gamma_tp * gamma_c * &
             !        gamma_ht*gamma_lt*gamma_hw*gamma_aq
             endif

             if ( (gamma >=0.0_r8) .and. (gamma< 100._r8) ) then

                vocflx_meg(imeg) =  meg_cmp%coeff * epsilon * gamma * megemis_units_factor_old / meg_cmp%molec_weight ! moles/m2/sec
        
                meg_out(imeg)%flux_out(p) = meg_out(imeg)%flux_out(p) &
                                          + epsilon * gamma * megemis_units_factor_old*1.e-3_r8 ! Kg/m2/sec

                   !stressed VOC
                   gamma_ht_all(p,imeg) = gamma_ht
                   gamma_lt_all(p,imeg) = gamma_lt
                   gamma_hw_all(p,imeg) = gamma_hw
                   gamma_aq_all(p,imeg) = gamma_aq
                if (class_num == 4 )then
                   gamma_ht_out(p) = gamma_ht
                   gamma_lt_out(p) = gamma_lt
                   gamma_hw_out(p) = gamma_hw
                   gamma_aq_out(p) = gamma_aq
                end if

                if (imeg==1) then 
                   !stress species
                   gamma_out(p)=gamma
                   gammaP_out(p)=gamma_p
                   gammaT_out(p)=gamma_t
                   gammaA_out(p)=gamma_a
                   gammaS_out(p)=gamma_sm
                   gammaL_out(p)=gamma_l
                   gammaC_out(p)=gamma_c

                   paru_out(p)=par_sun
                   par24u_out(p)=par24_sun
                   par240u_out(p)=par240_sun

                   para_out(p)=par_sha
                   par24a_out(p)=par24_sha
                   par240a_out(p)=par240_sha

                   alpha_out(p)=alpha
                   cp_out(p)=cp

                   topt_out(p)=topt
                   Eopt_out(p)=Eopt

                end if
             endif

             !if (debug .and. gamma > 0.0_r8) then
             if (debug ) then
                write(iulog,*) 'MEGAN: n, megan name, epsilon, gamma, vocflx: ', &
                     imeg, meg_cmp%name, epsilon, gamma, vocflx_meg(imeg), gamma_p,gamma_t,gamma_a,gamma_sm,gamma_l
             endif

             meg_cmp => meg_cmp%next_megcomp
          enddo meg_cmp_loop

          ! sum up the megan compound fluxes for the fluxes of chem mechanism compounds 
          do imech = 1,shr_megan_mechcomps_n
             n_meg_comps = shr_megan_mechcomps(imech)%n_megan_comps
             do imeg = 1,n_meg_comps ! loop over number of megan compounds that make up the nth mechanism compoud
                ii = shr_megan_mechcomps(imech)%megan_comps(imeg)%ptr%index
                vocflx(p,imech) = vocflx(p,imech) + vocflx_meg(ii)
             enddo
             vocflx_tot(p) = vocflx_tot(p) + vocflx(p,imech) ! moles/m2/sec
          enddo

       end if ! patch%itype(1:15 only)

    enddo ! fp 


  end associate
  end subroutine VOCEmission

  !-----------------------------------------------------------------------
  function get_map_EF(ivt_in, g_in, vocemis_inst)
    !
    ! Get mapped EF for isoprene
    ! Use gridded values for 6 Patches specified by MEGAN following
    ! Guenther et al. (2006).  Map the maxveg CLM Patches to these 6.
    ! Units: [ug m-2 h-1] 
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ivt_in
    integer, intent(in) :: g_in
    type(vocemis_type), intent(in) :: vocemis_inst
    !
    ! !LOCAL VARIABLES:
    real(r8)            :: get_map_EF
    !-----------------------------------------------------------------------

    ! vocemis_inst%efisop_patch ! Output: [real(r8) (:,:)]  emission factors for isoprene for each patch [ug m-2 h-1]

    get_map_EF = 0._r8
    
    if (     ivt_in == ndllf_evr_tmp_tree  &
         .or.     ivt_in == ndllf_evr_brl_tree) then   !fineleaf evergreen
       get_map_EF = vocemis_inst%efisop_grc(2,g_in)
    else if (ivt_in == ndllf_dcd_brl_tree) then        !fineleaf deciduous
       get_map_EF = vocemis_inst%efisop_grc(3,g_in)
    else if (ivt_in >= nbrdlf_evr_trp_tree &
         .and.    ivt_in <= nbrdlf_dcd_brl_tree) then  !broadleaf trees
       get_map_EF = vocemis_inst%efisop_grc(1,g_in)
    else if (ivt_in >= nbrdlf_evr_shrub &
         .and.    ivt_in <= nbrdlf_dcd_brl_shrub) then !shrubs
       get_map_EF = vocemis_inst%efisop_grc(4,g_in)
    else if (ivt_in >= nc3_arctic_grass &
         .and.    ivt_in <= nc4_grass) then            !grass
       get_map_EF = vocemis_inst%efisop_grc(5,g_in)
    else if (ivt_in >= nc3crop) then                   !crops
       get_map_EF = vocemis_inst%efisop_grc(6,g_in)
    end if

  end function get_map_EF

!  !-----------------------------------------------------------------------
!  function get_gamma_P(par_sun_in, par_sha_in, fsun_in, cp, alpha) 
!    !
!    ! Activity factor for PPFD (MEGAN3): all light dependent species
!    !-------------------------
!    ! Currently gamma_p doesn't rely on the long-term light conditions in MEGAN3
!    ! With distinction between sunlit and shaded leafs, weight scalings by
!    ! fsun and fshade 
!    ! !ARGUMENTS:
!    implicit none
!    real(r8),intent(in) :: par_sun_in
!    real(r8),intent(in) :: par_sha_in
!    real(r8),intent(in) :: fsun_in
!    real(r8),intent(out):: cp                      ! temporary
!    real(r8),intent(out):: alpha                   ! temporary
!    !
!    ! !LOCAL VARIABLES:
!    real(r8) :: gamma_p_LDF             ! activity factor for PPFD
!    real(r8) :: get_gamma_P             ! return value
!    real(r8), parameter :: alpha_fix = 0.004_r8  ! empirical coefficient
!    real(r8), parameter :: cp_fix = 1.03_r8      ! empirical coefficient
!    !-----------------------------------------------------------------------
!    alpha = alpha_fix
!    cp = cp_fix
!    ! SUN:
!    gamma_p_LDF = fsun_in * ( cp * alpha * par_sun_in * (1._r8 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5_r8) )
!    ! SHADE:
!    gamma_p_LDF = gamma_p_LDF + (1._r8-fsun_in) * (cp*alpha*par_sha_in*(1._r8 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5_r8))
!
!    !the LDF will be used for gamma_tp
!    get_gamma_P =  gamma_p_LDF
!
!  end function get_gamma_P
  !-----------------------------------------------------------------------
  function get_gamma_P(par_sun_in, par24_sun_in, par240_sun_in, par_sha_in, par24_sha_in, par240_sha_in, &
       fsun_in, fsun240_in, forc_solad240_in,forc_solai240_in, LDF_in, cp, alpha) 
    !
    ! Activity factor for PPFD (Guenther et al., 2006): all light dependent species
    !-------------------------
    ! With distinction between sunlit and shaded leafs, weight scalings by
    ! fsun and fshade 
    ! Scale total incident par by fraction of sunlit leaves (added on 1/2002)

    ! fvitt -- forc_solad240, forc_solai240 can be zero when CLM finidat is specified
    !          which will cause par240 to be zero and produce NaNs via log(par240)
    ! dml   -- fsun240 can be equal to or greater than one before 10 day averages are
    !           set on startup or if a new patch comes online during land cover change.
    !           Avoid this problem by only doing calculations with fsun240 when fsun240 is
    !           between 0 and 1
    !
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: par_sun_in
    real(r8),intent(in) :: par24_sun_in
    real(r8),intent(in) :: par240_sun_in
    real(r8),intent(in) :: par_sha_in
    real(r8),intent(in) :: par24_sha_in
    real(r8),intent(in) :: par240_sha_in
    real(r8),intent(in) :: fsun_in
    real(r8),intent(in) :: fsun240_in
    real(r8),intent(in) :: forc_solad240_in
    real(r8),intent(in) :: forc_solai240_in
    real(r8),intent(in) :: LDF_in
    real(r8),intent(out):: cp                      ! temporary
    real(r8),intent(out):: alpha                   ! temporary
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gamma_p_LDF             ! activity factor for PPFD
    real(r8) :: get_gamma_P             ! return value
    real(r8), parameter :: ca1 = 0.004_r8        ! empirical coefficent for alpha
    real(r8), parameter :: ca2 = 0.0005_r8       ! empirical coefficent for alpha
    real(r8), parameter :: ca3 = 0.0468_r8       ! empirical coefficent for cp
    real(r8), parameter :: par0_sun = 200._r8    ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: par0_shade = 50._r8   ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: alpha_fix = 0.001_r8  ! empirical coefficient
    real(r8), parameter :: cp_fix = 1.21_r8      ! empirical coefficient
    !-----------------------------------------------------------------------

    if ( (fsun240_in > 0._r8) .and. (fsun240_in < 1._r8) .and.  (forc_solad240_in > 0._r8) &
         .and. (forc_solai240_in > 0._r8)) then
       ! With alpha and cp calculated based on eq 6 and 7:
       ! Note indexing for accumulated variables is all at patch level
       ! SUN:
       alpha = ca1 - ca2 * log(par240_sun_in)
       cp = ca3 * exp(ca2 * (par24_sun_in-par0_sun))*par240_sun_in**(0.6_r8)
       gamma_p_LDF = fsun_in * ( cp * alpha * par_sun_in * (1._r8 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5_r8) )
       ! SHADE:
       alpha = ca1 - ca2 * log(par240_sha_in)
       cp = ca3 * exp(ca2 * (par_sha_in-par0_shade))*par240_sha_in**(0.6_r8)
       gamma_p_LDF = gamma_p_LDF + (1._r8-fsun_in) * (cp*alpha*par_sha_in*(1._r8 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5_r8))
    else
       ! With fixed alpha and cp (from MEGAN User's Guide):
       ! SUN: direct + diffuse  
       alpha = alpha_fix
       cp = cp_fix
       gamma_p_LDF = fsun_in * ( cp * alpha*par_sun_in * (1._r8 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5_r8) )
       ! SHADE: diffuse 
       gamma_p_LDF = gamma_p_LDF + (1._r8-fsun_in) * (cp*alpha*par_sha_in*(1._r8 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5_r8))
    end if

    ! Calculate total activity factor for PPFD accounting for light-dependent fraction
    get_gamma_P = (1._r8 - LDF_in) + LDF_in * gamma_p_LDF

  end function get_gamma_P

  !-----------------------------------------------------------------------
  function get_gamma_L(fsun240_in,elai_in)
    !
    ! Activity factor for LAI (Guenther et al., 2006): all species
    ! Guenther et al., 2006 eq 3
    !
    ! !USES:
    use clm_varcon   , only : denice
    use clm_varpar   , only : nlevsoi
    !
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: fsun240_in
    real(r8),intent(in) :: elai_in
    real(r8)            :: get_gamma_L             ! return value
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: cce = 0.30_r8                   ! factor to set emissions to unity @ std
    real(r8), parameter :: cce1 = 0.24_r8                  ! same as Cce but for non-accumulated vars
    real(r8), parameter :: CDEA = 0.925_r8                 ! CDEA value in MEGANv3 
    !-----------------------------------------------------------------------
    !   get_gamma_L = CDEA* elai_in
    if ( (fsun240_in > 0.0_r8) .and. (fsun240_in < 1.e30_r8) ) then 
       get_gamma_L = cce * elai_in
    else
       get_gamma_L = cce1 * elai_in
    end if

  end function get_gamma_L

  !-----------------------------------------------------------------------

  function get_gamma_SM_test(btran_in)
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: btran_in
    !!!------- the drought algorithm by Kc, 20210820--------
    !real(r8), parameter :: frac = 1.4          
    !real(r8), parameter :: k1 = -5.86        
    !real(r8), parameter :: b1 = 3.50        
    !real(r8), parameter :: k2 = -18.25           
    !real(r8), parameter :: b2 = 5e6
    !-------------------------------------------------------
    
    !!!------- the drought algorithm by Kc, 20211020--------
    real(r8), parameter :: frac = 1.4          
    real(r8), parameter :: k1 = -7.4463      
    real(r8), parameter :: b1 = 3.2552       
    real(r8), parameter :: k2 = -28.7629         
    real(r8), parameter :: b2 = 2.35e6
    real(r8)            :: get_gamma_SM_test
             if (btran_in >= 1.) then
                get_gamma_SM_test = 1!1/(1+b*exp(k*(theta1-h2osoi_vol_in(j)))) 
             !else if (btran_in > 0.78.and. btran_in < 1) then
             !   get_gamma_SM_test = 1.18
             else
                !get_gamma_SM_test = frac*(1/(1+b1*exp(k1*(btran_in))))*((1-1/frac)/(1+b2*exp(k2*((1.6-btran_in))))+1/frac) 
                get_gamma_SM_test = frac*(1/(1+b1*exp(k1*(btran_in-0.2))))*((1-1/frac)/(1+b2*exp(k2*((1.5-btran_in))))+1/frac) 
             endif 
 
  end function get_gamma_SM_test
  !-----------------------------------------------------------------------

  function get_gamma_SM_test2(btran_in)
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: btran_in
    !!!------- the drought algorithm by Kc, 20210820--------
    
    !!!------- the drought algorithm by Kc, 20211020--------
    real(r8), parameter :: frac = 1.4          
    real(r8), parameter :: k1 = -7.4463      
    real(r8), parameter :: b1 = 3.2552       
    real(r8), parameter :: k2 = -28.7629         
    real(r8), parameter :: b2 = 2.35e6
    real(r8)            :: get_gamma_SM_test2
             if (btran_in >= 1.) then
                get_gamma_SM_test2 = 1!1/(1+b*exp(k*(theta1-h2osoi_vol_in(j)))) 
             else
                get_gamma_SM_test2 = 1/(1+b1*exp(k1*(btran_in-0.2)))
             endif 
 
  end function get_gamma_SM_test2
  
!-----------------------------------------------------------------------
  function get_gamma_SM_test3(btran_in)
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: btran_in
    !!!------- the drought algorithm by Kc, 20210820--------
    
    !!!------- the drought algorithm by Kc, 20211020--------
    real(r8), parameter :: k1 = -15.34      
    real(r8), parameter :: b1 = 3.14       
    real(r8)            :: get_gamma_SM_test3
             if (btran_in >= 1.) then
                get_gamma_SM_test3 = 1!1/(1+b*exp(k*(theta1-h2osoi_vol_in(j)))) 
             else
                get_gamma_SM_test3 = 1/(1+b1*exp(k1*(btran_in-0.2)))
             endif 
 
  end function get_gamma_SM_test3

  !-----------------------------------------------------------------------
  function get_gamma_tld(t_veg240_in, t_veg24_in,t_veg_in, ct1_in, ct2_in, Ceo_in, Eopt, topt)

    ! Activity factor for temperature (light dependent) 
    !--------------------------------
    ! Calculate both a light-dependent fraction as in Guenther et al., 2006 for isoprene
    ! of a max saturation type form. Also caculate a light-independent fraction of the
    ! form of an exponential. Final activity factor depends on light dependent fraction
    ! of compound type.
    !
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: t_veg240_in
    real(r8),intent(in) :: t_veg24_in
    real(r8),intent(in) :: t_veg_in
    real(r8),intent(in) :: ct1_in
    real(r8),intent(in) :: ct2_in
    real(r8),intent(in) :: Ceo_in
    real(r8),intent(out) :: Eopt                    ! temporary 
    real(r8),intent(out) :: topt                    ! temporary 
    !
    ! !LOCAL VARIABLES:
    real(r8) :: get_gamma_tld
    real(r8) :: x                       ! temporary 
    real(r8), parameter :: co1 = 313._r8                   ! empirical coefficient
    real(r8), parameter :: co2 = 0.6_r8                    ! empirical coefficient
    real(r8), parameter :: co4 = 0.05_r8                   ! empirical coefficient
    real(r8), parameter :: tstd0 = 297_r8                  ! std temperature [K]
    real(r8), parameter :: topt_fix = 317._r8              ! std temperature [K]
    real(r8), parameter :: Eopt_fix = 2.26_r8              ! empirical coefficient
    real(r8), parameter :: ct3 = 0.00831_r8                ! empirical coefficient (0.0083 in User's Guide)
    !-----------------------------------------------------------------------

    ! Light dependent fraction (Guenther et al., 2006)
    if ( (t_veg240_in > 0.0_r8) .and. (t_veg240_in < 1.e30_r8) ) then 
       ! topt and Eopt from eq 8 and 9:
       topt = co1 + (co2 * (t_veg240_in-tstd0))
       Eopt = Ceo_in * exp (co4 * (t_veg24_in-tstd0)) * exp(co4 * (t_veg240_in -tstd0))
    else
       topt = topt_fix
       Eopt = Eopt_fix
    endif
    x = ( (1._r8/topt) - (1._r8/(t_veg_in)) ) / ct3
    get_gamma_tld = Eopt * ( ct2_in * exp(ct1_in * x)/(ct2_in - ct1_in * (1._r8 - exp(ct2_in * x))) )
    

  end function get_gamma_tld
  
  function get_gamma_tli(t_veg_in,betaT_in)

    ! Activity factor for temperature (light independent) 
    !--------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: t_veg_in
    real(r8),intent(in) :: betaT_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: get_gamma_tli
    real(r8), parameter :: tstd = 303.15_r8                ! std temperature [K]
    !-----------------------------------------------------------------------

    ! Light independent fraction (of exp(beta T) form)
    get_gamma_tli = exp(betaT_in * (t_veg_in - tstd))
  end function get_gamma_tli

  !-----------------------------------------------------------------------
  function get_gamma_T(t_veg240_in, t_veg24_in,t_veg_in, ct1_in, ct2_in, betaT_in, LDF_in, Ceo_in, Eopt, topt)

    ! Activity factor for temperature 
    !--------------------------------
    ! Calculate both a light-dependent fraction as in Guenther et al., 2006 for isoprene
    ! of a max saturation type form. Also caculate a light-independent fraction of the
    ! form of an exponential. Final activity factor depends on light dependent fraction
    ! of compound type.
    !
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: t_veg240_in
    real(r8),intent(in) :: t_veg24_in
    real(r8),intent(in) :: t_veg_in
    real(r8),intent(in) :: ct1_in
    real(r8),intent(in) :: ct2_in
    real(r8),intent(in) :: betaT_in
    real(r8),intent(in) :: LDF_in
    real(r8),intent(in) :: Ceo_in
    real(r8),intent(out) :: Eopt                    ! temporary 
    real(r8),intent(out) :: topt                    ! temporary 
    !
    ! !LOCAL VARIABLES:
    real(r8) :: get_gamma_T
    real(r8) :: gamma_t_LDF             ! activity factor for temperature
    real(r8) :: gamma_t_LIF             ! activity factor for temperature
    real(r8) :: x                       ! temporary 
    real(r8), parameter :: co1 = 313._r8                   ! empirical coefficient
    real(r8), parameter :: co2 = 0.6_r8                    ! empirical coefficient
    real(r8), parameter :: co4 = 0.05_r8                   ! empirical coefficient
    real(r8), parameter :: tstd0 = 297_r8                  ! std temperature [K]
    real(r8), parameter :: topt_fix = 317._r8              ! std temperature [K]
    real(r8), parameter :: Eopt_fix = 2.26_r8              ! empirical coefficient
    real(r8), parameter :: ct3 = 0.00831_r8                ! empirical coefficient (0.0083 in User's Guide)
    real(r8), parameter :: tstd = 303.15_r8                ! std temperature [K]
    real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
    !-----------------------------------------------------------------------

    ! Light dependent fraction (Guenther et al., 2006)
    if ( (t_veg240_in > 0.0_r8) .and. (t_veg240_in < 1.e30_r8) ) then 
       ! topt and Eopt from eq 8 and 9:
       topt = co1 + (co2 * (t_veg240_in-tstd0))
       Eopt = Ceo_in * exp (co4 * (t_veg24_in-tstd0)) * exp(co4 * (t_veg240_in -tstd0))
    else
       topt = topt_fix
       Eopt = Eopt_fix
    endif
    x = ( (1._r8/topt) - (1._r8/(t_veg_in)) ) / ct3
    gamma_t_LDF = Eopt * ( ct2_in * exp(ct1_in * x)/(ct2_in - ct1_in * (1._r8 - exp(ct2_in * x))) )
    
    
    ! Light independent fraction (of exp(beta T) form)
    gamma_t_LIF = exp(betaT_in * (t_veg_in - tstd))
    
    ! Calculate total activity factor for light as a function of light-dependent fraction
    !--------------------------------
    get_gamma_T = (1-LDF_in)*gamma_T_LIF + LDF_in*gamma_T_LDF 

  end function get_gamma_T



  !-----------------------------------------------------------------------
  function get_gamma_A(ivt_in, elai240_in, elai_in, nclass_in)

    ! Activity factor for leaf age (Guenther et al., 2006)
    !-----------------------------
    ! If not CNDV elai is constant therefore gamma_a=1.0
    ! gamma_a set to unity for evergreens (Patches 1, 2, 4, 5)
    ! Note that we assume here that the time step is shorter than the number of 
    !days after budbreak required to induce isoprene emissions (ti=12 days) and 
    ! the number of days after budbreak to reach peak emission (tm=28 days)
    !
    ! !ARGUMENTS:
    implicit none
    integer,intent(in)  :: ivt_in
    integer,intent(in)  :: nclass_in
    real(r8),intent(in) :: elai240_in
    real(r8),intent(in) :: elai_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: get_gamma_A
    real(r8) :: elai_prev               ! lai for previous timestep
    real(r8) :: fnew, fgro, fmat, fold  ! fractions of leaves at different phenological stages
    !-----------------------------------------------------------------------
    if ( (ivt_in == ndllf_dcd_brl_tree) .or. (ivt_in >= nbrdlf_dcd_trp_tree) ) then  ! non-evergreen

       if ( (elai240_in > 0.0_r8) .and. (elai240_in < 1.e30_r8) )then 
          elai_prev = 2._r8*elai240_in-elai_in  ! have accumulated average lai over last 10 days
          if (elai_prev == elai_in) then
             fnew = 0.0_r8
             fgro = 0.1_r8
             fmat = 0.8_r8
             fold = 0.1_r8
          else if (elai_prev > elai_in) then
             fnew = 0.0_r8
             fgro = 0.0_r8
             fmat = 1.0_r8 - (elai_prev - elai_in)/elai_prev
             fold = (elai_prev - elai_in)/elai_prev
          else if (elai_prev < elai_in) then
             fnew = 1 - (elai_prev / elai_in)
             fgro = 0.0_r8
             fmat = (elai_prev / elai_in)
             fold = 0.0_r8
          end if
          
          get_gamma_A = fnew*Anew(nclass_in) + fgro*Agro(nclass_in) + fmat*Amat(nclass_in) + fold*Aold(nclass_in)

       else
          get_gamma_A = 1.0_r8
       end if
       
    else
       get_gamma_A = 1.0_r8
    end if
    

  end function get_gamma_A

  !-----------------------------------------------------------------------
  function get_gamma_C(cisun_in,cisha_in,forc_pbot_in,fsun_in, co2_ppmv)
    
    ! Activity factor for instantaneous CO2 changes (Heald et al., 2009)
    !-------------------------
    ! With distinction between sunlit and shaded leaves, weight scalings by
    ! fsun and fshade 
    !
    ! !CALLED FROM: VOCEmission
    !
    ! !REVISION HISTORY:
    ! Author: Colette L. Heald (11/30/11)
    !         Louisa K. Emmons (16/03/2015) - implement Colette's intended code
    !                                         and use atmosphere CO2 (not nml setting)
    !
    ! !USES:
    !    use clm_varctl,    only : co2_ppmv      ! corresponds to CCSM_CO2_PPMV set in env_conf.xml
    !
    ! !ARGUMENTS:
    implicit none
    ! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: cisun_in
    real(r8),intent(in) :: cisha_in
    real(r8),intent(in) :: forc_pbot_in
    real(r8),intent(in) :: fsun_in
    real(r8),intent(in) :: co2_ppmv

    real(r8)            :: get_gamma_C

    ! local variables
    real(r8)            :: Ismax            ! empirical coeff for CO2 
    real(r8)            :: h                ! empirical coeff for CO2 
    real(r8)            :: Cstar            ! empirical coeff for CO2 
    real(r8)            :: fint             ! interpolation fraction for CO2
    real(r8)            :: ci               ! temporary sunlight/shade weighted cisun & cisha (umolCO2/mol)
    real(r8)            :: gamma_ci         ! short-term exposure gamma
    real(r8)            :: gamma_ca         ! long-term exposure gamma
    real(r8), parameter :: Ismax_ca   = 1.344_r8  ! Estimated asymptote at which further decreases in intercellular CO2 have a negligible effect on isoprene emission
    real(r8), parameter :: h_ca       = 1.4614_r8 ! Exponential scalar
    real(r8), parameter :: Cstar_ca   = 585._r8   ! Scaling coefficient
    real(r8), parameter :: CiCa_ratio = 0.7_r8    ! Ratio of intercellular CO2 to atmospheric CO2
    !-----------------------------------------------------------------------


    ! LONG-TERM EXPOSURE (based on ambient CO2, Ca)
    !-----------------------------------------------------------------------------
    gamma_ca = Ismax_ca - ((Ismax_ca * (CiCa_ratio*co2_ppmv)**h_ca) / (Cstar_ca**h_ca + (CiCa_ratio*co2_ppmv)**h_ca) )


    ! SHORT-TERM EXPOSURE (based on intercellular CO2, Ci)
    !-----------------------------------------------------------------------------
    ! Determine long-term CO2 growth environment (ie. ambient CO2) and interpolate
    ! parameters
    if ( co2_ppmv < 400._r8 ) then
       Ismax   = 1.072_r8
       h       = 1.70_r8
       Cstar   = 1218._r8
    else if ( (co2_ppmv > 400._r8) .and. (co2_ppmv < 600._r8) ) then
       fint    = (co2_ppmv - 400._r8)/200._r8
       Ismax   = fint*1.036_r8 + (1.- fint)*1.072_r8
       h       = fint*2.0125_r8 + (1.- fint)*1.70_r8
       Cstar   = fint*1150._r8 + (1.- fint)*1218._r8
    else if ( (co2_ppmv > 600._r8) .and. (co2_ppmv < 800._r8) ) then
       fint    = (co2_ppmv - 600._r8)/200._r8
       Ismax   = fint*1.046_r8 + (1.- fint)*1.036_r8
       h       = fint*1.5380_r8 + (1.- fint)*2.0125_r8
       Cstar   = fint*2025._r8 + (1.- fint)*1150._r8
    else if ( co2_ppmv > 800._r8 ) then
       Ismax   = 1.014_r8
       h       = 2.861_r8
       Cstar   = 1525._r8
    end if

    ! Intercellular CO2 concentrations (ci) given in Pa, divide by atmos
    ! pressure to get mixing ratio (umolCO2/mol)
    if ( (cisun_in .eq. cisun_in) .and. (cisha_in .eq. cisha_in) .and. (forc_pbot_in > 0._r8) .and. (fsun_in > 0._r8) ) then
       ci = ( fsun_in*cisun_in + (1._r8-fsun_in)*cisha_in )/forc_pbot_in * 1.e6_r8
       gamma_ci = Ismax - ( (Ismax*ci**h)/(Cstar**h+ci**h) ) 
    else if ( (cisun_in > 0.0_r8) .and. (cisun_in < 1.e30_r8) .and. (forc_pbot_in > 0._r8) .and. (fsun_in .eq. 1._r8) ) then
       ci = cisun_in/forc_pbot_in * 1.e6_r8
       gamma_ci = Ismax - ( (Ismax*ci**h)/(Cstar**h+ci**h) ) 
    else if ( (cisha_in > 0.0_r8) .and. (cisha_in < 1.e30_r8)  .and. (forc_pbot_in > 0._r8) .and. (fsun_in .eq. 0._r8) ) then
       ci = cisha_in/forc_pbot_in * 1.e6_r8
       gamma_ci = Ismax - ( (Ismax*ci**h)/(Cstar**h+ci**h) ) 
    else
       gamma_ci = 1._r8
    end if

    get_gamma_C = gamma_ci * gamma_ca

  end function get_gamma_C
!===================stressed VOC=============================
!========================High Temperature====================
  !-----------------------------------------------------------------------
  function get_gamma_ht(maxT_in,CHT_in,THT_in,DTHT_in)

    !--------------------------------
    !--------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: CHT_in
    real(r8),intent(in) :: THT_in
    real(r8),intent(in) :: DTHT_in
    real(r8),intent(in) :: maxT_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: THTK,t1
    real(r8) :: get_gamma_ht

    !real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
    if (THT_in .lt. 173.15)then !Temperature < -100 Degree C
    THTK = 273.15+THT_in
    else
    THTK = THT_in
    endif 
    t1 = THTK + DTHT_in
    !init_value
    get_gamma_ht = 1._r8
 
    if ( maxT_in .le. THTK ) then 
    get_gamma_ht = 1._r8 
    else if ( (maxT_in .gt. THTK) .and. (maxT_in .lt. t1 ) ) then
    get_gamma_ht = 1._r8 + (CHT_in - 1._r8)* (maxT_in - THTK)/DTHT_in
    else
    get_gamma_ht = CHT_in
    endif
  end function get_gamma_ht

!========================Low Temperature====================
  function get_gamma_lt(minT_in,CLT_in,TLT_in,DTLT_in)

    !--------------------------------
    !--------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: CLT_in
    real(r8),intent(in) :: TLT_in
    real(r8),intent(in) :: DTLT_in
    real(r8),intent(in) :: minT_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: TLTK,t1
    real(r8) :: get_gamma_lt

    !real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]

    if (TLT_in .lt. 173.15)then !Temperature < -100 Degree C
    TLTK = 273.15+TLT_in
    else
    TLTK = TLT_in
    endif 
    t1 = TLTK - DTLT_in
    !init_value
    get_gamma_lt = 1._r8 
    if ( minT_in .ge. TLTK ) then 
    get_gamma_lt = 1._r8 
    else if ( (minT_in .lt. TLTK) .and. (minT_in .gt. t1 ) ) then
    get_gamma_lt = 1._r8 + (CLT_in - 1.0)* (TLTK - minT_in)/DTLT_in
    else
    get_gamma_lt = CLT_in
    endif
  end function get_gamma_lt
!========================High Wind Speed====================
  function get_gamma_hw(maxWS_in,CHW_in,THW_in,DTHW_in)

    !--------------------------------
    !--------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: CHW_in
    real(r8),intent(in) :: THW_in
    real(r8),intent(in) :: DTHW_in
    real(r8),intent(in) :: maxWS_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: t1
    real(r8) :: get_gamma_hw

    !real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
    t1 = THW_in + DTHW_in
    !init_value
    get_gamma_hw = 1._r8 
    if ( maxWS_in .le. THW_in ) then 
    get_gamma_hw = 1._r8 
    else if ( (maxWS_in .gt. THW_in) .and. (maxWS_in .lt. t1 ) ) then
    get_gamma_hw = 1._r8 + (CHW_in - 1.0)* (maxWS_in - THW_in)/DTHW_in
    else
    get_gamma_hw = CHW_in
    end if
  end function get_gamma_hw
!========================Air Quality====================
  function get_gamma_aq(AQ_in,CAQ_in,TAQ_in,DTAQ_in)

    !--------------------------------
    !--------------------------------
    ! !ARGUMENTS:
    implicit none
    real(r8),intent(in) :: CAQ_in
    real(r8),intent(in) :: TAQ_in
    real(r8),intent(in) :: DTAQ_in
    real(r8),intent(in) :: AQ_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: t1
    real(r8) :: get_gamma_aq

    !real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
    t1 = TAQ_in + DTAQ_in
    !init_value
    get_gamma_aq = 1._r8 
    if ( AQ_in .le. TAQ_in ) then 
    get_gamma_aq = 1._r8 
    else if ( (AQ_in .gt. TAQ_in) .and. (AQ_in .lt. t1 ) ) then
    get_gamma_aq = 1._r8 + (CAQ_in - 1.0)* (AQ_in - TAQ_in)/DTAQ_in
    else
    get_gamma_aq = CAQ_in
    end if
  end function get_gamma_aq
!============================================================


end module VOCEmissionMod


