module condtend

!
! script adapted from OsloAeroSec condtend.F90 and aeronucl.F90 by smb
!
! TODO: make copy of state for substepping
! TODO: currently stuck in an infinite while loop
! TODO: get rid of: use aero_sectional,     only: secConstIndex index for aerosol in q array
! TODO: get rid of chemistryindex -> find physics index, then - imozart
! TODO: re-implement organics.. somehow :P
! TODO: add new mass for SO4 to aero_state%aero_range_state(1)%mmr(:,:,SO4)
! TODO: add condensate to aero_state%aero_range_state(:)%mmr(:,:,SO4)
! TODO: make sure unit of massDistrib is correct, in kg/kg, not number concentration
! TODO: make diagnostic output, addfld and outfld for nucleation and condensation tendencies
! TODO: verify all the calculations, this is just a quick cleanup of the osloaerosec version, and I have not checked the physics yet.
! TODO: check for unused variables, and remove them
! TODO: make an attribute in aerosol properties for the precursor gas, and de CAM-ify the code, so that it is not hardcoded to H2SO4
! GOOD LUCK!
  ! CAM/NorESM specific stuff
   use phys_control, only: phys_getopts
   use chem_mods,    only: gas_pcnst
   use mo_tracname,  only: solsym
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld, addfld, add_default, fieldname_len, horiz_only
   use physconst,    only: rair, gravit, pi
   use chem_mods,    only: adv_mass !molecular weights from mozart
   use ppgrid,       only: pcols, pver, pverp
   use constituents, only: cnst_get_ind
   use spmd_utils,   only: masterproc

   ! Sectional aerosol specific stuff
   use aerosol_properties_mod,           only: aerosol_properties
   use sectional_aerosol_properties_mod, only: sectional_aerosol_properties
   use aerosol_state_mod,                only: aerosol_state, ptr2d_t
   use sectional_aerosol_state_mod,      only: sectional_aerosol_state

   implicit none

   ! public subroutines, used in aero_model
   public :: condensation_init
   public :: condtend_sub_super

   private

   ! module variables
   real(r8), allocatable :: cond_sink_norm(:)       ![m3/#/s] condensation sink per particle in bin i
   real(r8), allocatable :: bin_centers(:)          ![m] bin centers
   integer  :: l_h2so4_chem
   integer  :: imozart
   integer  :: l_h2so4
   integer  :: sulfate_specprop_ndx ! index for sulfate species in aero_props

contains

!==========================================================================================
! Initialize Condensation
!==========================================================================================

   subroutine condensation_init(aero_props)

      !condensation coefficients:
      !Theory: Poling et al, "The properties of gases and liquids"
      !5th edition, eqn 11-4-4

      ! dummy arguments
      type(sectional_aerosol_properties), intent(in) :: aero_props

      ! local parameters
      real(r8), parameter :: aunit   = 1.6606e-27_r8  ! [kg] Atomic mass unit
      real(r8), parameter :: boltz   = 1.3806e-23_r8  ! [J/K/molec]
      real(r8), parameter :: t0      = 273.15_r8      ! [K] standard temperature
      real(r8), parameter :: p0      = 101325.0_r8    ! [Pa] Standard pressure
      real(r8), parameter :: radair  = 1.73e-10_r8    ! [m] Typical air molecule collision radius
      real(r8), parameter :: Mair    = 28.97_r8       ! [amu/molec] Molecular weight for dry air

      ! Diffusion volumes for simple molecules [Poling et al], table 11-1
! TODO: vad seems like some kind of property for h2so4?
      real(r8), parameter :: vad     = 51.96_r8       ![cm3/mol]
      real(r8), parameter :: vadAir  = 19.7_r8        ![cm3/mol]
      real(r8), parameter :: aThird  = 1.0_r8/3.0_r8
      real(r8), parameter :: cm2Tom2 = 1.e-4_r8       ! convert from cm2 ==> m2

      ! local variables
      real(r8), allocatable :: diff_coeff(:)          ! [m2/s] Diffusion coefficient sectional

      real(r8) :: mfv             ![m] mean free path
      real(r8) :: diff            ![m2/s] diffusion coefficient for cond. vap
      real(r8) :: Mdual           ![molec/amu] 1/M_1 + 1/M_2
      real(r8) :: radmol          ![m] radius molecule
      real(r8) :: th              !thermal velocity
      integer  :: ibin, ispecprop ! indexes
      character(len=16) :: type   ! type of species in aero_props

      ! local variables for output
      logical  :: history_aerosol
!      character(len=fieldname_len+3) :: fieldname_donor
!      character(len=fieldname_len+3) :: fieldname_receiver
!      character(128)                 :: long_name
!      character(8)                   :: unit

      !-----------------------------------------------------------------------------------

      allocate(cond_sink_norm(aero_props%nbins()))
      allocate(diff_coeff(aero_props%nbins()))
      allocate(bin_centers(aero_props%nbins()))

      ! Couple the condenseable vapours to chemical species for properties and indexes
      ! add dimension for several species
! TODO: get this from properties?
      ! gas phase h2so4
      call cnst_get_ind('H2SO4',l_h2so4, abort=.true.)

      ! find first Mozart tracer imozart
      call cnst_get_ind(trim(solsym(1)), imozart, abort=.true.)

      ! Find the chemistry index for gas H2SO4
      l_h2so4_chem = l_h2so4 - imozart + 1

      do ispecprop = 1, aero_props%nspecies_tot() !TODO: make nicer to have several species
        type = aero_props%spectype(ispecprop)
        if (trim(type) == 'sulfate') then
            sulfate_specprop_ndx = ispecprop
            exit
        end if
      end do

      bin_centers = aero_props%bin_centers(aero_props%nbins())

      ! pick up densities and weights from aerosol properties
      radmol = ( 3.0_r8*aero_props%molecular_weight(sulfate_specprop_ndx)*aunit &
                 / ( 4.0_r8*pi*aero_props%density(sulfate_specprop_ndx) ) )**aThird          ! Radius of molecul (straight forward assuming spherical)
      Mdual  = 2.0_r8/(1.0_r8/Mair+1.0_r8/aero_props%molecular_weight(sulfate_specprop_ndx)) ! factor of [1/m_1 + 1_m2]

      ! thermal velocity for H2SO4 in air (m/s)
      ! https://en.wikipedia.org/wiki/Thermal_velocity
      th = sqrt(8.0_r8*boltz*t0/ ( pi*aero_props%molecular_weight(sulfate_specprop_ndx) * aunit ) )

      ! calculating microphysical parameters from equations in Ch. 8 of Seinfeld & Pandis (1998):
      ! mean free path for molec in air (m)
      mfv = 1.0_r8 / ( pi*sqrt( 1.0_r8+aero_props%molecular_weight(sulfate_specprop_ndx)/Mair ) &
                    * (radair+radmol)**2 * p0/(boltz*t0) )

      ! Solve eqn 11-4.4 in Poling et al
      ! (A bit hard to follow units here, but result in the book is in cm2/s)..
      ! so scale by "cm2Tom2" to get m2/sec
! TODO: what is Vad?
      diff = cm2Tom2   &
         *0.00143_r8*t0**1.75_r8     &
         /((p0/1.0e5_r8)*sqrt(Mdual)   &
         *(((Vad)**aThird+(Vadair)**aThird)**2))

      do ibin = 1, aero_props%nbins()         !all bins receive condensation
          ! Correct for non-continuum effects, formula is from
          ! Chuang and Penner, Tellus, 1995, sticking coeffient from
          ! Vignati et al, JGR, 2004
          diff_coeff(ibin) = diff  &    !original diffusion coefficient
               / ( bin_centers(ibin) / (bin_centers(ibin) + mfv )  &  ! non-continuum correction factor
               +4.0_r8*diff/ (1._r8*th*bin_centers(ibin) ) )
      enddo

      !Find sink per particle in mode "imode"
      !Eqn 13 in Kulmala et al, Tellus 53B, 2001, pp 479
      !http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.530411.x/abstract

      cond_sink_norm = 0.0_r8
      do ibin = 1, aero_props%nbins()
         ! Since we do not sum over bins, it is slightly different than above (normnk=1, no summing)
! TODO: what?
         cond_sink_norm(ibin) =   + 4.0_r8*pi * diff_coeff(ibin) &    ![m2/s] diffusion coefficient
                                                * bin_centers(ibin)     ![m] radius of bin
      end do

      !Initialize output
      call phys_getopts(history_aerosol_out = history_aerosol)

! TODO: addflds
! condensation tendency for SO4_R1, maybe num_1, ...

   end subroutine condensation_init

!==========================================================================================
! Condensation + Nucleation called from aero_model
! manages time sub-stepping by calling condtend_sub
!==========================================================================================

   subroutine condtend_sub_super(lchnk, q, cond_vap_gasprod, temperature, &
               pmid, pdel, dt, ncol, pblh, zm, qh20, aero_props, aero_state)
      ! Calculates nucleation rate and condensation rate of aerosols
      !
      ! This method calls the condtend method. If the timestep needs to be split in two,
      ! this will be done here and the condtend method will be called several times.
      ! This method also writes output once condend is done.
      !

      use cam_history,  only: outfld,fieldname_len
      use constituents, only: pcnst  ! h2so4 and soa nucleation (cka)

      ! dummy arguments
      type(sectional_aerosol_properties), intent(in) :: aero_props
      type(sectional_aerosol_state), intent(inout) :: aero_state

      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: ncol                  ! number of columns
      real(r8), intent(in)    :: temperature(:,:)      ! Temperature (K)
      real(r8), intent(in)    :: pmid(:,:)             ! [Pa] pressure at mid point
      real(r8), intent(in)    :: pdel(:,:)             ! [Pa] difference in grid cell
      real(r8), intent(in)    :: cond_vap_gasprod(:,:) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in)    :: dt                    ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(:)               ! pbl height (m)
      real(r8), intent(in)    :: zm(:,:)               ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(:,:)             ! specific humidity (kg/kg)
      real(r8), intent(inout) :: q(:,:,:) ! TMR [kg/kg] including moisture

      real(r8) :: q_t0(pcols,pver,gas_pcnst) ! mass before subroutine.
      real(r8) :: dt_local

      !output:
      real(r8) :: coltend (pcols, gas_pcnst)      ! ??
      real(r8) :: coltend_dummy(pcols, gas_pcnst) ! ??
      real(r8) :: nuclrate_pbl(pcols,pver)        ! [kg/kg] tracer lost
      real(r8) :: nuclrate(pcols,pver)            ! [kg/kg] tracer lost
      real(r8) :: formrate_pbl(pcols,pver)        ! [kg/kg] tracer lost
      real(r8) :: formrate(pcols,pver)            ! [kg/kg] tracer lost
      real(r8) :: h2so4nucl(pcols,pver)           ! h2so4 in nucleation code
      real(r8) :: grh2so4(pcols,pver)             ! growth rate h2so4
      real(r8) :: coagnucl(pcols,pver)            ! coagulation in nucleation
      real(r8) :: nuclnum(pcols,pver)             ! [#/m3/s] nucleation number rate from RM's parameterization

      real(r8), allocatable :: numconc_new(:,:,:) ! [#/m3] number concentration new
      real(r8), allocatable :: leaveSec(:,:,:)    ! [kg/kg] tracer lost
      real(r8), allocatable :: leaveSec_dummy(:,:,:) ![kg/kg] tracer lost
      real(r8), pointer     :: tmp_num(:,:)

      logical  :: notDone ! if not done, continues
      logical  :: split_dt ! whether timestep is split or not
      integer  :: nr_dt, cnt,i,j,k  !number of runs, counter, counter, counter
      integer  :: ibin ! index for bin

! TODO: what to do with this?
      real(r8) :: tend(pcols, pver, gas_pcnst) ! tendency for each tracer

      ! local variables for output
      logical                        :: history_aerosol
!      character(128)                 :: long_name
!      character(8)                   :: unit
!      character(18) :: fieldname_receiver

      !-----------------------------------------------------------------------------------

      allocate(numconc_new(ncol, pver, aero_props%nbins()))
      allocate(leaveSec_dummy(ncol, pver, aero_props%nbins()))
      allocate(leaveSec(ncol, pver, aero_props%nbins()))
      !initialization
      numconc_new = 0.0_r8

      q_t0(:ncol,:,:)= q(:,:,:) ! in case timestep needs to be decreased.
      coltend(:,:)   = 0.0_r8
      coltend_dummy(:,:) = 0.0_r8
      nuclrate_pbl(:,:) = 0.0_r8
      nuclrate(:,:)  = 0.0_r8
      formrate_pbl(:,:) = 0.0_r8
      formrate(:,:)  = 0.0_r8
      h2so4nucl(:,:) = 0.0_r8
      grh2so4(:,:)   = 0.0_r8
      coagnucl(:,:)  = 0.0_r8
      dt_local       = dt / 2.0_r8 ! always half timestep
      notDone        = .TRUE.
      split_dt       = .FALSE.
      nr_dt          = 2
      cnt            = 1

      ! get bin number concentrations
      tmp_num => null()
      do ibin = 1, aero_props%nbins()
          call aero_state%get_ambient_num(ibin, tmp_num)
          numconc_new(:,:,ibin) = tmp_num(:,:)
      end do

      ! run until no need to split time any longer
      do while (notDone )

           call condtend_sub(lchnk, q, cond_vap_gasprod, temperature, &
                           nuclrate,nuclrate_pbl,formrate, formrate_pbl, coagnucl, &
                           h2so4nucl, grh2so4, &
                           coltend_dummy, split_dt, &
                           leaveSec_dummy,&
                           pmid, pdel, dt_local, ncol, pblh, zm, qh20, aero_props, aero_state, numconc_new, nuclnum)
           coltend = coltend + coltend_dummy*dt_local ! divides by timestep at end
           leaveSec = leaveSec + leaveSec_dummy*dt_local

           if (split_dt) then ! If split timestep: split timestep
               dt_local = dt_local/2.0_r8
               cnt      = 1
               nr_dt    = nr_dt*2
               q(:ncol,:,:) = q_t0(:ncol,:,:)
               numconc_new(:ncol, :, :) = aero_state%bin_numconc(:ncol,:,:)    ! is reset if timestep is decreased.
               coltend(:,:) = 0.0_r8
               coltend_dummy(:,:) = 0.0_r8
               nuclrate_pbl(:,:) = 0.0_r8
               nuclrate(:,:) = 0.0_r8
               formrate_pbl(:,:) = 0.0_r8
               formrate(:,:)  = 0.0_r8
               h2so4nucl(:,:) = 0.0_r8
               grh2so4(:,:)   = 0.0_r8
               coagnucl(:,:)  = 0.0_r8
               notDone        = .TRUE.

           else if (nr_dt .eq. cnt) then ! if not, check if count is eq to number of splits
               !if (nr_dt .eq. cnt) then
               notDone=.FALSE.
               ! update bin with numconc from nucleation
! TODO: add nucleation
! TODO: add ranges
               do i = 1, ncol
                   do k = 1, pver
                       call aero_state%update_bin(1, i, k, 0._r8, -nuclnum(i,k)*dt, 0, dt, tend)
                   end do
                end do
           else ! if not done and no need to split timestep again, add one to count.
               cnt=cnt+1
           end if
      end do

      ! divide by timestep:
      leaveSec =   leaveSec/dt
      coltend  =   coltend/dt
      nuclrate =   nuclrate/dt
      nuclrate_pbl = nuclrate_pbl/dt
      formrate =   formrate/dt
      formrate_pbl = formrate_pbl/dt
      h2so4nucl =  h2so4nucl/dt
      grh2so4  =   grh2so4/dt
      coagnucl =  coagnucl/dt

      ! TODO: write output

      call phys_getopts(history_aerosol_out = history_aerosol)

 !     if(history_aerosol)then
! TODO: more nice output

      deallocate(numconc_new)

   end subroutine condtend_sub_super

!==========================================================================================
! Condensation helper routine
! sub method.
! calculate the sulphate nucleation rate, and condensation rate of
! aerosols used for parameterising the transfer of externally mixed
! aitken mode particles into an internal mixture.
! note the parameterisation for conversion of externally mixed particles
!  used the h2so4 lifetime onto the particles, and not a given
! increase in particle radius. will be improved in future versions of the model
! added input for h2so4 and soa nucleation: soa_lv_gasprod, soa_sv_gasprod, pblh,zi,qh20 (cka)
!==========================================================================================

   subroutine condtend_sub(lchnk,  q, cond_vap_gasprod, temperature,            &
           !smb++sectional
                nuclrate, nuclrate_pbl_o, formrate, formrate_pbl_o, coagnucl_o,  &
                h2so4nucl_o, grh2so4_o,                     &
                coltend_o, split_dt,                                            &
                leaveSec,                                                       &
           !smb--sectional
               pmid, pdel, dt, ncol, pblh,zm,qh20,                              &
               aero_props, aero_state, numconc_new, nuclnum)

      use cam_history,     only: outfld,fieldname_len
      use constituents,    only: pcnst  ! h2so4 and soa nucleation (cka)

      ! dummy arguments
      type(sectional_aerosol_properties), intent(in) :: aero_props
      type(sectional_aerosol_state), intent(in) :: aero_state

      real(r8), intent(inout) :: numconc_new(:,:,:) ![#/m3] number concentration new

       !++smb sectional
      real(r8), intent(inout)  :: nuclrate (:,:)            ! Nucleation rate output
      real(r8), intent(inout)  :: nuclrate_pbl_o (:,:)      ! Nucleation rate pbl output
      real(r8), intent(inout)  :: formrate(:,:)             ! Formation rate output
      real(r8), intent(inout)  :: formrate_pbl_o(:,:)        ! Formation rate pbl output

      real(r8), intent(inout)  :: coagnucl_o(:,:)           ! Coagulation sink for npf output
      real(r8), intent(inout)  :: h2so4nucl_o(:,:)          ! H2SO4 for nucleation output
      real(r8), intent(inout)  :: grh2so4_o(:,:)            ! GR from H2SO4 output

      real(r8), intent(out)    :: coltend_o(:,:)       ! column tendency output
      logical,  intent(out)    :: split_dt                          ! if true, time step needs to be split

      real(r8), intent(out)    :: leaveSec(:,:,:)                   ![kg/kg] tracer lost
      real(r8), intent(out)    :: nuclnum(:,:)                      ![#/m3/s] nucleation number rate from RM's parameterization
       !--smb sectional

      ! arguments
      integer,  intent(in) :: lchnk                      ! chunk identifier
      integer,  intent(in) :: ncol                       ! number of columns
      real(r8), intent(in) :: temperature(:,:)    ! Temperature (K)
      real(r8), intent(in) :: pmid(:,:)           ! [Pa] pressure at mid point
      real(r8), intent(in) :: pdel(:,:)           ! [Pa] difference in grid cell
      real(r8), intent(inout) :: q(:,:,:) ! TMR [kg/kg] including moisture

      real(r8), intent(in) :: cond_vap_gasprod(:,:) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in) :: dt                         ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(pcols)               ! pbl height (m)
      real(r8), intent(in)    :: zm(:,:)           ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(:,:)          ! specific humidity (kg/kg)

      ! local
      character(len=fieldname_len+3) :: fieldname
      integer :: i,k
      real(r8) :: sumCondensationSink(pcols,pver)       ![1/s] sum of condensation sink
      real(r8) :: totalLoss(pcols,pver,gas_pcnst)       ![kg/kg] tracer lost
      real(r8) :: coltend(pcols, gas_pcnst)
      real(r8) :: tracer_coltend(pcols)
      real(r8), allocatable :: condensationsink_sec(:)  ![1/s] loss rate per mode (mixture)
      real(r8), allocatable :: condensationsinkfraction_sec(:,:,:) ! [frc]
      real(r8), allocatable :: numberconcentration_sec(:,:,:) ![#/m3] number concentration
      real(r8), allocatable :: tend(:,:,:)

      real(r8)       :: intermediateConcentration(pcols,pver)
      real(r8)       :: rhoAir(pcols,pver)              ![kg/m3] density of air
      ! Volume of added  material from condensate;  surface area of core particle;
      real(r8)       :: volume_shell, area_core,vol_monolayer
      real(r8)       :: frac_transfer                   ! Fraction of hydrophobic material converted to an internally mixed mode

      ! needed for h2so4 and soa nucleation treatment
       integer  :: modeIndexReceiverCoag              ! Index of modes receiving coagulate
       integer  :: iCoagReceiver                      ! counter for species receiving coagulate
       real(r8) :: coagulationSink(pcols,pver)        ! [1/s] coaglation loss for SO4_n and soa_n
       real(r8), parameter :: lvocfrac=0.5            ! Fraction of organic oxidation products with low enough
                                                      ! volatility to enter nucleation mode particles (1-24 nm)
       real(r8), pointer :: tmp_num(:,:)
       real(r8) :: soa_lv_forNucleation(pcols,pver)   ! [kg/kg] soa gas available for nucleation
       real(r8) :: gasLost(pcols,pver)                ! [kg/kg] budget terms on H2SO4 (gas)
       real(r8) :: fracNucl(pcols,pver)               ! [frc] fraction of gas nucleated
       real(r8) :: firstOrderLossRateNucl(pcols,pver) ! [1/s] first order loss rate due to nucleation
       real(r8) :: nuclso4(pcols,pver)                ! [kg/kg/s] Nucleated so4 mass tendency from RM's parameterization
       real(r8) :: nuclsoa(pcols,pver)                ! [kg/kg/s] Nucleated soa mass tendency from RM's parameterization
       real(r8), allocatable :: rho_part(:, :, :)
       !smb++ sectional
       real(r8) :: dummy  !
       integer  :: ibin, irange ! indices

       ! local variables for output
       logical        :: history_aerosol
!       character(128) :: long_name                              ![-] needed for diagnostics

      !-----------------------------------------------------------------------------------

       allocate(condensationsink_sec(aero_props%nbins()))
       allocate(numberconcentration_sec(ncol,pver,aero_props%nbins()))
       allocate(condensationsinkfraction_sec(ncol,pver,aero_props%nbins()))
       allocate(rho_part(ncol,pver,aero_props%nbins()))
       allocate(tend(ncol, pver, gas_pcnst))
       !Initialize h2so4 and soa nucl variables
       coagulationSink = 0.0_r8
       condensationsinkfraction_sec = 0.0_r8
       numberconcentration_sec = 0.0_r8
       tmp_num => null()
       coltend_o = 0.0_r8
       leaveSec = 0.0_r8
       split_dt       = .FALSE.

       do k = 1, pver
           do i = 1, ncol
                !Air density
                rhoAir(i,k) = pmid(i,k)/rair/temperature(i,k)
           end do
        end do

       do ibin = 1, aero_props%nbins()
        ! No looping through species, mmr is added afterwards
!          call aero_state%get_ambient_num(ibin, tmp_num)
          numberconcentration_sec(:ncol,:,ibin) = numconc_new(:ncol, :, ibin) / rhoAir(:ncol,:)  ![#/m3] number concentration
          rho_part(:ncol,:,ibin) = aero_state%bin_dry_density(ibin, ncol)
       enddo

       do k=1,pver
           do i=1,ncol

                !smb++ sectional
                ! initialize condensation sink for sectional
                condensationSink_sec = 0.0_r8  !Sink to the coming "receiver" of any vapour
                !smb-- sectional

                !NB: The following is duplicated code, coordinate with koagsub!!
                !Initialize number concentration for this receiver


                !Go though all bins receiving condensation

                ! condensation sink to sectional bin:
                   do ibin = 1, aero_props%nbins()

                      !This is the loss rate a gas molecule will see due to aerosol surface area
                      condensationSink_sec(ibin)   = cond_sink_norm(ibin)  & ![m3/#/s]  per particle
                                                          * numberConcentration_sec(i,k,ibin)             ![#/m3]
                                                          !==> [1/s]
                   end do !Loop over receivers
                !smb-- sectional

                !Find concentration after condensation of all
                !condenseable vapours
                ! smb++ sectional edited:
                ! condensation sink fraction that goes to each bin in sectional scheme
                ! set to zero first:
                ! smb-- sectional

                    !sum of cond. sink for this vapour [1/s]
                    !smb++ sectional
                    ! Need to add condensation sink to sectional scheme for particles in sectional scheme
                    ! However, not all tracers may contribute:
                    ! assumes same order of gasses in sectional and other (1: H2SO4,2: SOA_LV, 3: SOA_SV)
                        sumCondensationSink(i,k) = sum(condensationSink_sec( :))
                        ! Keeps track of the fraction of the condensate to the sectional bins for each tracer
                        condensationSinkFraction_sec(i,k,:) = condensationSink_sec( :) &
                                /(sumCondensationSink(i,k)+1.e-30_r8)![frc]
                    !smb-- sectional


                !Solve the intermediate (end of timestep) concentration using
                !euler backward solution C_{old} + P *dt - L*C_{new}*dt = C_{new} ==>
                !Cnew -Cold = prod - loss ==>
                intermediateConcentration(i,k) = &
                                     ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                                     / (1.0_r8 + sumCondensationSink(i,k)*dt)


                !Assume only a fraction of ORG_LV left can contribute to nucleation

                modeIndexReceiverCoag = 0
                !Sum coagulation sink for nucleated so4 and soa particles over all receivers of coagulate. Needed for RM's nucleation code
                !OBS - looks like RM's coagulation sink is multiplied by 10^-12??


           end do !index i
       end do !index k

       !Calculate nucleated masses of so4 and soa (nuclso4, nuclsoa)
       !following RM's parameterization (cka)
! TODO: check, how is intermediateConcentration h2so4pc, and how is soa_lv_forNucleation = coagnuc?
       call aeronucl(lchnk, ncol, temperature, pmid, qh20, &
                   intermediateConcentration(:,:), soa_lv_forNucleation, &
                   coagulationSink, nuclnum, nuclso4, nuclsoa, zm, pblh, &
                   nuclrate, nuclrate_pbl_o, formrate, formrate_pbl_o, &
                   h2so4nucl_o, grh2so4_o, dt, &
                   bin_centers(1), aero_props &
               )

       !smb++ sectional
       coagnucl_o(:,:) = coagnucl_o(:,:) + coagulationSink(:,:)*dt
       !smb-- sectional
       firstOrderLossRateNucl = 0.0_r8
       do k=1,pver
          do i=1,ncol

             !First order loss rate (1/s) for nucleation
              !smb++ added check for 0
              if (intermediateConcentration(i,k) .eq. 0) then
                  firstOrderLossRateNucl(i,k)=0._r8
              else
                  firstOrderLossRateNucl(i,k) = nuclSo4(i,k)/intermediateConcentration(i,k)
              end if
! organics here
              !smb-- added check for 0
             !First order loss rate (1/s) for nucleation

                !Solve implicitly (again)
                !C_new - C_old =  PROD_{gas} - CS*C_new*dt - LR_{nucl}*C_new =>
                intermediateConcentration(i,k) = &
                               ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                               / (1.0_r8 + sumCondensationSink(i,k)*dt + firstOrderLossRateNucl(i,k)*dt)

                !fraction nucleated
                fracNucl(i,k) = firstOrderLossRateNucl(i,k) &
                                     /(firstOrderLossRateNucl(i,k) + sumCondensationSink(i,k))
               !From budget, we get: lost = prod -cnew + cold
                gasLost(i,k) = cond_vap_gasprod(i,k)*dt   & !Produced
                                     + q(i,k,l_h2so4_chem)            & !cold
                                     - intermediateConcentration(i,k)    !cnew

             !H2SO4 condensate
             do ibin=1, aero_props%nbins()
                ! bin_numconc (#/kg)
                ! gasLost: [kg/kg] mass lost from gas phase
                ! condensationSinkFraction_sec: [frc] fraction to the particular bin
                ! fracNucl: [frc] fraction nucleated
                ! get number from gasLost
! TODO: fix update_bin -> mmr to correct range
! TODO: condensation does not change the number!! -> add mass to ranges
!                numconc(i, k, ibin) = numconc(i, k, ibin) + gasLost(i,k) * (1._r8-fracNucl(i,k)) * condensationSinkFraction_sec(i,k, ibin) / aero_props%density(sulfate_specprop_ndx) / aero_props%particle_volume(ibin)
!                call aero_state%update_bin(ibin, i, k, 0._r8, -1._r8*gasLost(i,k) * (1._r8-fracNucl(i,k)) * condensationSinkFraction_sec(i,k, ibin) / aero_props%density(sulfate_specprop_ndx) / aero_props%particle_volume(ibin), 0, dt, tend)
                irange = aero_props%bins2ranges(ibin)

!                ispec = ??
                ! add mass to the range:
!                aero_state%aero_range_state(irange)%mmr(i,k,ispec) = aero_state%aero_range_state(irange)%mmr(i,k,ispec) &
!                     + gasLost(i,k)*(1.0_r8-fracNucl(i,k)) &
!                     *condensationSinkFraction_sec(i,k, ibin)

             end do

             !condenseable vapours
             q(i,k,l_h2so4_chem)  = intermediateConcentration(i,k)

             !smb++sectional grow particles in sectional scheme:
    ! removed median radius because not used
             q_mmr = q(i,k,:) * aero_props%particle_volume(:) * aero_state%bin_dry_density(:)    ! [kg/kg] mass in each tracer

             call sec_moveMass(q_mmr(i,k,:), numberConcentration_sec(i,k,:), leaveSec(i,k,:), &
                            rhoAir(i,k), split_dt, aero_props, aero_state, rho_part(i, k,:), numconc_new(i, k,:))

             ! Add nucleated number to smallest bin after redistribution (#/kg/s * s = #/kg)
             numconc_new(i, k, 1) = numconc_new(i, k, 1) + nuclnum(i,k)*dt
             ! Add nucleated mass to first bin of sectional scheme:
 !            q(i,k,chemistryIndex(secConstIndex(1,1))) =  q(i,k,chemistryIndex(secConstIndex(1,1)))       &
 !                        + gasLost(i,k)*fracNucl(i,k)

            ! Add mass from sectional scheme to so4_na and soa_na:
!             q(i,k,chemistryIndex(l_so4_na)) = q(i,k,chemistryIndex(l_so4_na))         &
!                                    +leaveSec(i,k, 1)

             !smb--sectional

          end do !physical index k
       end do    !physical index i

       !Output for diagnostics
!       call phys_getopts(history_aerosol_out = history_aerosol)

!       if(history_aerosol)then
!          coltend(:ncol,:) = 0.0_r8

          !smb++ sectional
          ! Remove so4_n ---> directly into so4_na
!          coltend(:ncol,chemistryIndex(secConstIndex(1,1))) = coltend(:ncol,chemistryIndex(secConstIndex(1,1))) + &
          !smb-- sectional
!                                                 sum(                                         &
!                                                    gasLost(:ncol,:)           &
!                                                    *fracNucl(:ncol,:)*pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt

          ! Remove so4_n ---> directly into so4_na
          !smb++ put in how much leaves sectional scheme:
!          coltend(:ncol,chemistryIndex(l_so4_na)) = coltend(:ncol,chemistryIndex(l_so4_na)) + &
!                                                 sum(                                         &
!                                                    leaveSec(:ncol,:)           &
!                                                    *pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt
          !smb--sectional
          !Take into account H2SO4 (gas) condensed in budget
!          coltend(:ncol,chemistryIndex(l_so4_a1)) = coltend(:ncol,chemistryIndex(l_so4_a1)) + &
!                                                 sum(                                         &
!                                                    gasLost(:ncol,:)           &
                                                    !smb++sectional subtract fraction to sectional scheme
!                                                    *(1-sum(condensationSinkFraction_sec(:ncol,:,:),3)) &
                                                    !smb--sectional
!                                                    *(1.0_r8 - fracNucl(:ncol,:))*pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt

          !Take into account soa_lv (gas) nucleated in budget
          !smb++ sectional
! organics here
          !smb++ sectional: putt in how condenses on sectional:
!          do ibin=1,aero_props%nbins()
!               coltend(:ncol,chemistryIndex(secConstIndex(1,ibin))) = coltend(:ncol, chemistryIndex(secConstIndex(1,ibin)))+ &
!                                                    sum(                    &
!                                                    gasLost(:ncol, :)    &
!                                                    *(condensationSinkFraction_sec(:ncol,:,ibin)) &
!                                                    *(1.0_r8 - fracNucl(:ncol, :))*pdel(:ncol,:),2) &
!                                                    /gravit/dt
!          end do
! organics here

!          coltend_o(:,:)=coltend(:,:)

!       endif


   end subroutine condtend_sub

!==========================================================================================
! Helper routine to move mass between bins from growth
!==========================================================================================

   subroutine sec_moveMass(massDistrib, numberConc_old, leave_sec, rhoAir, decrease_dt, aero_props, aero_state, rho_aer, numconc_new)

! TODO: add move mass for species
! TODO: add move number for bins
    ! Moves tracer mass from on bin to the other based on condensational/coagulation growth.
    ! Based on Jacobson Fundamentals of Atmospheric Modeling, second edition (2005),
    ! Chapter   13.5
    class(sectional_aerosol_properties), intent(in) :: aero_props
    class(sectional_aerosol_state), intent(in) :: aero_state

    real(r8), intent(in)    :: numberConc_old(:)    ! numbr concentration before growth
    real(r8), intent(in)    :: rho_aer(:)           ! Density of particle)
    real(r8), intent(out)   :: numconc_new(:)
    real(r8), intent(inout) :: massDistrib(:)       ! mass in each tracer
    real(r8), intent(out)   :: leave_sec(:)         ! the mass that leaves sectional scheme
    logical,  intent(out)   :: decrease_dt          ! if set to True, time step is divided
                                                    ! and the procedure is re run

    real(r8), allocatable   :: volume(:)            ! volume of particle in bin dimension(secNrBins)
    real(r8), allocatable   :: volume_new(:)        ! volume after growth dimension(secNrBins)
    real(r8), allocatable   :: volfrac(:)           ! dimension(secNrSpec,secNrBins)
    real(r8)                :: xfrac                ! fraction to stay in bin

    real(r8), parameter     :: pi = 3.141592654_r8
    real(r8)                :: rhoAir               ! Density of air
    integer                 :: ibin
    integer, allocatable    :: chem_ndx(:)
    !-----------------------------------------------------------------------------------

    allocate(volume(aero_props%nbins()))
    allocate(volume_new(aero_props%nbins()))
    allocate(volfrac(aero_props%nbins()))
    allocate(chem_ndx(aero_props%nbins()))

    decrease_dt=.FALSE.
    !compute volume in each bin with condensation (by mass) and by
    !numberconcentration
    do ibin = 1, aero_props%nbins()
       volume_new(ibin) = 0.0_r8
       volfrac(ibin) = 0.0_r8

       ! find chem index for num_ibin
       chem_ndx(ibin) = aero_state%num_transport_ndx(ibin) - imozart + 1

       if (numberConc_old(ibin)<1.e-30_r8) then
          volume_new(ibin)=0.0_r8
       else
! TODO: this is a bit nonsensical right now, since massDistrib is a number concentration, not a mass. Need to fix this.
          volume_new(ibin) = volume_new(ibin) + massDistrib(chem_ndx(ibin)) / &
                           rho_aer(ibin) * rhoAir / numberConc_old(ibin)

       end if

       volfrac(ibin) = massDistrib(chem_ndx(ibin))/&
                            rho_aer(ibin)*rhoAir
                            !kg/kg(air)*[kg(air)/m3(air)][kg/m3]--> m3/m3(air)
   !    volfrac(ibin)=volfrac(ibin)/(sum(volfrac(ibin))+1.E-50_r8)

       if (volfrac(ibin)<1.e-50_r8) then
          volfrac(ibin)=0.0_r8
       end if

       ! calculate volume in each bin by mass/density! m3
       volume(ibin) =  pi * bin_centers(ibin)**3/6._r8
            ! calculate volume in each bin by numberconcentration (volume from before condenstion)
    end do

    numconc_new(:) = 0._r8

    do ibin =  1, aero_props%nbins()-1
       ! fraction to stay in bin
       xfrac=(volume(ibin+1)-volume_new(ibin)) &
                            /(volume(ibin+1)-volume(ibin))
          if (numberConc_old(ibin)<1.e-30) then
             xfrac=1.0_r8
          end if

          if (xfrac .le. 0._r8) then        ! if the fraction to stay is equal to
                                            ! less than zero, then the
                                            ! aerosols have grown too large
                                            ! for the next bin and we will
                                            ! want to decrease the time step
                                            ! to avoid this.
             decrease_dt=.TRUE.
          end if

          if (xfrac .le. 0._r8) then
             decrease_dt=.TRUE.
          end if

          xfrac=max(0._r8, min(1._r8,xfrac))
          numconc_new(ibin) = numconc_new(ibin) + &
                                    xfrac*numberConc_old(ibin) &
                                    *volfrac(ibin)
          numconc_new(ibin+1) = numconc_new(ibin+1) + &
                                    (1-xfrac)*numberConc_old(ibin)   &
                                    *volfrac(ibin)
    end do

    xfrac = ( (bin_centers(aero_props%nbins())*2._r8)**3 * pi/6.0_r8 &
            - volume_new(aero_props%nbins())) &
            / ( (bin_centers(aero_props%nbins())*2._r8)**3 * pi/6.0_r8 &
            - volume(aero_props%nbins()))

    ! if less than or 0 % stays in bin, we must decrease timestep
    if (xfrac .le. 0._r8) then
            decrease_dt=.TRUE.
    end if

    xfrac=max(0._r8, min(1._r8,xfrac))

   ! do indSpec=1, secNrSpec
            numconc_new(aero_props%nbins()) = numconc_new(aero_props%nbins()) + &
                                            xfrac * numberConc_old(aero_props%nbins()) &
                                            * volfrac(aero_props%nbins())
! TODO: treat "leave_sec" differently, e.g. largest bin doesn't grow?
! rho_aer of last bin correct here?
            leave_sec = pi * ( bin_centers(aero_props%nbins())*2._r8 )**3 / 6.0_r8 * rho_aer(aero_props%nbins())/rhoAir &   ! [m3_aer/#]*[kg_aer/m3_aer]/[kg_air/m3_air]--> [kg_aer/kg_air/#][m3_air]
                                            * (1-xfrac) * numberConc_old(aero_props%nbins()) &              ! *[#/m3_air] --> kg_aer/kg_air
                                            * volfrac(aero_props%nbins())
   ! end do
    do ibin=1,aero_props%nbins()
        massDistrib(chem_ndx(ibin)) = rho_aer(ibin)/rhoAir &
                            * numconc_new(ibin) * pi * bin_centers(ibin)**3/6.0_r8
    end do

   end subroutine sec_moveMass

!==========================================================================================
! Helper routine for aerosol nucleation - Vehkamäki et al. (2002)
!==========================================================================================

   subroutine aeronucl(lchnk, ncol, t, pmid, h2ommr, h2so4pc, oxidorg, coagnuc, nuclnum, nuclso4, nuclorg, zm, pblht, &
                nuclrate, nuclrate_pbl_o, formrate, formrate_pbl_o,  &
                h2so4nucl_o, grh2so4_o, dt, &
                radius, aero_props)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use wv_saturation,  only: qsat_water
    use physconst,      only: avogad, rair
    use ppgrid,         only: pcols, pver, pverp
    use cam_history,    only: outfld
    use phys_control,   only: phys_getopts
    use chem_mods,      only: adv_mass
    use m_spc_id,       only : id_H2SO4
 !   use const,          only : volumeToNumber
    use shr_const_mod,  only: shr_const_rgas

    !-- Arguments
    class(sectional_aerosol_properties), intent(in) :: aero_props

    real(r8), intent(in)  :: dt                ! timestep (output is weighted by this)
    integer,  intent(in)  :: lchnk             ! chunk identifier
    integer,  intent(in)  :: ncol              ! number of atmospheric column
    real(r8), intent(in)  :: pmid(pcols,pver)  ! layer pressure (Pa)
    real(r8), intent(in)  :: h2ommr(:,:)       ! layer specific humidity
    real(r8), intent(in)  :: t(:,:)            ! Temperature (K)
    real(r8), intent(in)  :: h2so4pc(:,:)      ! Sulphuric acid concentration (kg kg-1)
    real(r8), intent(in)  :: oxidorg(:,:)      ! Organic vapour concentration (kg kg-1)
    real(r8), intent(in)  :: coagnuc(:,:)      ! Coagulation sink for nucleating particles [1/s]

    real(r8), intent(in)  :: zm(:,:)           ! Height at layer midpoints (m)
    real(r8), intent(in)  :: pblht(:)          ! Planetary boundary layer height (m)
    !smb++sectional
    real(r8), intent(in)  :: radius            ! Particle size at calculated formation rate [m]
    ! because the timestep may be divided in two, the output needs to be averaged over all timesteps
    ! therefore we track this in these variables
    real(r8), intent(inout)  :: nuclrate(:,:)       ! Nucleation rate output
    real(r8), intent(inout)  :: nuclrate_pbl_o(:,:) ! Nucleation in pbl rate output
    real(r8), intent(inout)  :: formrate(:,:)       ! formation rate output
    real(r8), intent(inout)  :: formrate_pbl_o(:,:) ! formation rate in pbl output

    real(r8), intent(inout)  :: h2so4nucl_o(:,:)    ! conc h2so4 for output
    real(r8), intent(inout)  :: grh2so4_o(:,:)      ! GR from H2SO4 for output
    !smb-- sectional
    real(r8), intent(out) :: nuclorg(:,:)      ! Nucleated mass (ORG)
    real(r8), intent(out) :: nuclso4(:,:)      ! Nucleated mass (H2SO4)
    real(r8), intent(out) :: nuclnum(:,:)    ! [#/m3/s] Nucleated mass (SO4)

! TODO: implement oxidorg, orgnucl_o, grsoa_o, nuclorg

    !-- Local variables
! TODO: are these in physics somewhere?
    real(r8), parameter   :: pi=3.141592654_r8
    !cka+
    real(r8), parameter   :: h2so4_dens=1841._r8       ! h2so4 density [kg m-3]
    real(r8), parameter   :: org_dens=2000._r8         ! density of organics [kg m-3], based on RM assumptions
    !cka -

    integer               :: i,k
    real(r8)              :: qs(pcols,pver)            ! Saturation specific humidity
    real(r8)              :: relhum(pcols,pver)        ! Relative humidity
    real(r8)              :: h2so4(pcols,pver)         ! Sulphuric acid concentration [#/cm3]
    real(r8)              :: nuclvolume(pcols,pver)    ! [m3/m3/s] Nucleated mass (SO4)
    real(r8)              :: rhoair(pcols,pver)        ! density of air [kg/m3] !cka
    real(r8)              :: pblht_lim(pcols)          ! Planetary boundary layer height (m) (500m<pblht_lim<7000m) (cka)

    real(r8)              :: nuclrate_bin(pcols,pver) ! Binary nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_bin(pcols,pver) ! Binary formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_bin(pcols,pver) ! Binary nucleation critical cluster size (m)
    real(r8)              :: nuclrate_pbl(pcols,pver) ! Boundary layer nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_pbl(pcols,pver) ! Boundary layer formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_pbl(pcols,pver) ! Boundary layer nucleation formation size (m)

    real(r8)              :: orgforgrowth(pcols,pver) ! Organic vapour mass available for growth
    real(r8)              :: gr(pcols,pver), grh2so4(pcols,pver), grorg(pcols,pver) !growth rates
    real(r8)              :: vmolh2so4, vmolorg       ! [m/s] molecular speed of condenseable gases
    real(r8)              :: frach2so4
    real(r8)              :: dummy

    integer               :: atm_nucleation           ! Nucleation parameterization for the whole atmosphere
    integer               :: pbl_nucleation           ! Nucleation parameterization for the boundary layer
    real(r8)              :: molmass_h2so4            ! molecular mass of h2so4 [g/mol]
    real(r8)              :: molmass_soa              ! molecular mass of soa [g/mol]
!TODO: fix this l_so4_na -> use aerosol tracers from state
 !   integer               :: l_so4_na
   ! Variables for binary nucleation parameterization
    real(r8)              :: zrhoa, zrh, zt, zt2, zt3, zlogrh, zlogrh2, zlogrh3, zlogrhoa, zlogrhoa2, zlogrhoa3, x, zxmole, zix
    real(r8)              :: zjnuc, zntot, zrc, zrxc

    !-----------------------------------------------------------------------------------

   !cka: OBS    call phys_getopts(pbl_nucleation_out=pbl_nucleation, atm_nucleation_out=atm_nucleation)
    !cka: testing by setting these flags:
    pbl_nucleation = 2 ! smb++ 3 use Riccobono 2014 for nucleation. -> cvb++ use 2, because we don't currently have organics
    atm_nucleation = 1

    nuclso4(:,:)=0._r8
    nuclorg(:,:)=0._r8

    ! set all organic stuff to 0
    grorg = 0._r8
    orgforgrowth = 0._r8
    vmolorg = 0._r8

    !-- The highest level in planetary boundary layer
    do i=1,ncol
        pblht_lim(i)=MIN(MAX(pblht(i),500._r8),7000._r8)
    end do

    !-- Get molecular mass of h2so4 and soa_lv (cka)
    molmass_h2so4=adv_mass(id_H2SO4)

! TODO: make this something sensible:
  !  l_so4_na = 1

    !-- Formation diameters (m). Nucleated particles are inserted to SO4(n), same size used for soa  (cka)

    !-- Conversion of H2SO4 from kg/kg to #/cm3
    !-- and calculation of relative humidity (needed by binary nucleation parameterization)
    do k=1,pver
        do i=1,ncol
            ! air density
            rhoair(i,k)=pmid(i,k)/(t(i,k)*rair)
            !avogad*1.e-3_r8 to get molec/mol instead of molec/kmol for h2so4 gas
            h2so4(i,k)=(1.e-6_r8*h2so4pc(i,k)*avogad*1.e-3_r8*rhoair(i,k))/(molmass_h2so4*1.E-3_r8)
!            orgforgrowth(i,k)=(1.e-6_r8*oxidorg(i,k)*avogad*1.e-3_r8*rhoair(i,k))/(molmass_soa*1.E-3_r8)
 !           orgforgrowth(i,k)=MAX(MIN(orgforgrowth(i,k),1.E10_r8),0._r8)

            call qsat_water(t(i,k),pmid(i,k),dummy,qs(i,k))

            relhum(i,k) = h2ommr(i,k)/qs(i,k)
            relhum(i,k) = max(relhum(i,k),0.0_r8)
            relhum(i,k) = min(relhum(i,k),1.0_r8)
        end do !ncol
    end do     !layers


    !-- Binary sulphuric acid-water nucleation rate
    if(atm_nucleation .EQ. 1) then
        do k=1,pver
            do i=1,ncol

                ! Calculate nucleation only for valid thermodynamic conditions:
                zrhoa = max(h2so4(i,k),1.E+4_r8)
                zrhoa = min(zrhoa,1.E11_r8)

                zrh   = max(relhum(i,k),1.E-4_r8)
                zrh   = min(zrh,1.0_r8)

                zt    = max(t(i,k),190.15_r8)
                zt    = min(zt,300.15_r8)

                zt2 = zt*zt
                zt3 = zt2*zt

                ! Equation (11) - molefraction of H2SO4 in the critical cluster

                zlogrh  = LOG(zrh)
                zlogrh2 = zlogrh*zlogrh
                zlogrh3 = zlogrh2*zlogrh

                zlogrhoa  = LOG(zrhoa)
                zlogrhoa2 = zlogrhoa*zlogrhoa
                zlogrhoa3 = zlogrhoa2*zlogrhoa

                x=0.7409967177282139_r8 - 0.002663785665140117_r8*zt   &
                + 0.002010478847383187_r8*zlogrh    &
                - 0.0001832894131464668_r8*zt*zlogrh    &
                + 0.001574072538464286_r8*zlogrh2        &
                - 0.00001790589121766952_r8*zt*zlogrh2    &
                + 0.0001844027436573778_r8*zlogrh3     &
                -  1.503452308794887e-6_r8*zt*zlogrh3    &
                - 0.003499978417957668_r8*zlogrhoa   &
                + 0.0000504021689382576_r8*zt*zlogrhoa

                zxmole=x

                zix = 1.0_r8/x

                ! Equation (12) - nucleation rate in 1/cm3s

                zjnuc=0.1430901615568665_r8 + 2.219563673425199_r8*zt -   &
                  0.02739106114964264_r8*zt2 +     &
                  0.00007228107239317088_r8*zt3 + 5.91822263375044_r8*zix +     &
                  0.1174886643003278_r8*zlogrh + 0.4625315047693772_r8*zt*zlogrh -     &
                  0.01180591129059253_r8*zt2*zlogrh +     &
                  0.0000404196487152575_r8*zt3*zlogrh +    &
                  (15.79628615047088_r8*zlogrh)*zix -     &
                  0.215553951893509_r8*zlogrh2 -    &
                  0.0810269192332194_r8*zt*zlogrh2 +     &
                  0.001435808434184642_r8*zt2*zlogrh2 -    &
                  4.775796947178588e-6_r8*zt3*zlogrh2 -     &
                  (2.912974063702185_r8*zlogrh2)*zix -   &
                  3.588557942822751_r8*zlogrh3 +     &
                  0.04950795302831703_r8*zt*zlogrh3 -     &
                  0.0002138195118737068_r8*zt2*zlogrh3 +    &
                  3.108005107949533e-7_r8*zt3*zlogrh3 -     &
                  (0.02933332747098296_r8*zlogrh3)*zix +     &
                  1.145983818561277_r8*zlogrhoa -    &
                  0.6007956227856778_r8*zt*zlogrhoa +    &
                  0.00864244733283759_r8*zt2*zlogrhoa -    &
                  0.00002289467254710888_r8*zt3*zlogrhoa -    &
                  (8.44984513869014_r8*zlogrhoa)*zix +    &
                  2.158548369286559_r8*zlogrh*zlogrhoa +   &
                  0.0808121412840917_r8*zt*zlogrh*zlogrhoa -    &
                  0.0004073815255395214_r8*zt2*zlogrh*zlogrhoa -   &
                  4.019572560156515e-7_r8*zt3*zlogrh*zlogrhoa +    &
                  (0.7213255852557236_r8*zlogrh*zlogrhoa)*zix +    &
                  1.62409850488771_r8*zlogrh2*zlogrhoa -    &
                  0.01601062035325362_r8*zt*zlogrh2*zlogrhoa +   &
                  0.00003771238979714162_r8*zt2*zlogrh2*zlogrhoa +    &
                  3.217942606371182e-8_r8*zt3*zlogrh2*zlogrhoa -    &
                  (0.01132550810022116_r8*zlogrh2*zlogrhoa)*zix +    &
                  9.71681713056504_r8*zlogrhoa2 -    &
                  0.1150478558347306_r8*zt*zlogrhoa2 +    &
                  0.0001570982486038294_r8*zt2*zlogrhoa2 +    &
                  4.009144680125015e-7_r8*zt3*zlogrhoa2 +    &
                  (0.7118597859976135_r8*zlogrhoa2)*zix -    &
                  1.056105824379897_r8*zlogrh*zlogrhoa2 +    &
                  0.00903377584628419_r8*zt*zlogrh*zlogrhoa2 -    &
                  0.00001984167387090606_r8*zt2*zlogrh*zlogrhoa2 +    &
                  2.460478196482179e-8_r8*zt3*zlogrh*zlogrhoa2 -    &
                  (0.05790872906645181_r8*zlogrh*zlogrhoa2)*zix -    &
                  0.1487119673397459_r8*zlogrhoa3 +    &
                  0.002835082097822667_r8*zt*zlogrhoa3 -    &
                  9.24618825471694e-6_r8*zt2*zlogrhoa3 +    &
                  5.004267665960894e-9_r8*zt3*zlogrhoa3 -    &
                  (0.01270805101481648_r8*zlogrhoa3)*zix

                zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]

                ! Equation (13) - total number of molecules in the critical cluster

                zntot=-0.002954125078716302_r8 - 0.0976834264241286_r8*zt +   &
                  0.001024847927067835_r8*zt2 - 2.186459697726116e-6_r8*zt3 -    &
                  0.1017165718716887_r8*zix - 0.002050640345231486_r8*zlogrh -   &
                  0.007585041382707174_r8*zt*zlogrh +    &
                  0.0001926539658089536_r8*zt2*zlogrh -   &
                  6.70429719683894e-7_r8*zt3*zlogrh -    &
                  (0.2557744774673163_r8*zlogrh)*zix +   &
                  0.003223076552477191_r8*zlogrh2 +   &
                  0.000852636632240633_r8*zt*zlogrh2 -    &
                  0.00001547571354871789_r8*zt2*zlogrh2 +   &
                  5.666608424980593e-8_r8*zt3*zlogrh2 +    &
                  (0.03384437400744206_r8*zlogrh2)*zix +   &
                  0.04743226764572505_r8*zlogrh3 -    &
                  0.0006251042204583412_r8*zt*zlogrh3 +   &
                  2.650663328519478e-6_r8*zt2*zlogrh3 -    &
                  3.674710848763778e-9_r8*zt3*zlogrh3 -   &
                  (0.0002672510825259393_r8*zlogrh3)*zix -    &
                  0.01252108546759328_r8*zlogrhoa +   &
                  0.005806550506277202_r8*zt*zlogrhoa -    &
                  0.0001016735312443444_r8*zt2*zlogrhoa +   &
                  2.881946187214505e-7_r8*zt3*zlogrhoa +    &
                  (0.0942243379396279_r8*zlogrhoa)*zix -   &
                  0.0385459592773097_r8*zlogrh*zlogrhoa -   &
                  0.0006723156277391984_r8*zt*zlogrh*zlogrhoa +   &
                  2.602884877659698e-6_r8*zt2*zlogrh*zlogrhoa +    &
                  1.194163699688297e-8_r8*zt3*zlogrh*zlogrhoa -   &
                  (0.00851515345806281_r8*zlogrh*zlogrhoa)*zix -    &
                  0.01837488495738111_r8*zlogrh2*zlogrhoa +   &
                  0.0001720723574407498_r8*zt*zlogrh2*zlogrhoa -   &
                  3.717657974086814e-7_r8*zt2*zlogrh2*zlogrhoa -    &
                  5.148746022615196e-10_r8*zt3*zlogrh2*zlogrhoa +    &
                  (0.0002686602132926594_r8*zlogrh2*zlogrhoa)*zix -   &
                  0.06199739728812199_r8*zlogrhoa2 +    &
                  0.000906958053583576_r8*zt*zlogrhoa2 -   &
                  9.11727926129757e-7_r8*zt2*zlogrhoa2 -    &
                  5.367963396508457e-9_r8*zt3*zlogrhoa2 -   &
                  (0.007742343393937707_r8*zlogrhoa2)*zix +    &
                  0.0121827103101659_r8*zlogrh*zlogrhoa2 -   &
                  0.0001066499571188091_r8*zt*zlogrh*zlogrhoa2 +    &
                  2.534598655067518e-7_r8*zt2*zlogrh*zlogrhoa2 -    &
                  3.635186504599571e-10_r8*zt3*zlogrh*zlogrhoa2 +    &
                  (0.0006100650851863252_r8*zlogrh*zlogrhoa2)*zix +   &
                  0.0003201836700403512_r8*zlogrhoa3 -    &
                  0.0000174761713262546_r8*zt*zlogrhoa3 +   &
                  6.065037668052182e-8_r8*zt2*zlogrhoa3 -    &
                  1.421771723004557e-11_r8*zt3*zlogrhoa3 +   &
                  (0.0001357509859501723_r8*zlogrhoa3)*zix

                zntot=EXP(zntot)  !  add. Eq. (13)

                ! Equation (14) - radius of the critical cluster in nm

                zrc=EXP(-1.6524245_r8+0.42316402_r8*x+0.33466487_r8*LOG(zntot))    ! [nm]

                !----1.2) Limiter

                IF(zjnuc<1.e-7_r8 .OR. zntot<4.0_r8) zjnuc=0.0_r8

                ! limitation to 1E+10 [1/cm3/s]

                nuclrate_bin(i,k)=MAX(MIN(zjnuc,1.E10_r8),0._r8)
                nuclsize_bin(i,k)=MAX(MIN(zrc,1.E2_r8),0.01_r8)

            end do
        end do
    else   !No atmospheric nucleation
        nuclrate_bin(:,:)=0._r8
        nuclsize_bin(:,:)=1._r8
    end if

    !-- Boundary layer nucleation
    do k=1,pver
        do i=1,ncol

            !-- Nucleation rate #/cm3/s
            if(pblht_lim(i)>zm(i,k) .AND. pbl_nucleation>0) then

                if(pbl_nucleation .EQ. 1) then

                    !-- Paasonen et al. (2010), eqn 10, Table 4
                    nuclrate_pbl(i,k)=(1.7E-6_r8)*h2so4(i,k)

                else if(pbl_nucleation .EQ. 2) then
                    !smb++ sectional : updated nucleation parameterization
                    !-- Paasonen et al. (2010)
                    !values from Table 3 in Paasonen et al (2010), modified version of eqn 14
                    nuclrate_pbl(i,k)=(6.1E-7_r8)*h2so4(i,k)+(0.39E-7_r8)*orgforgrowth(i,k) !(18)

                    !nuclrate_pbl(i,k)=(1.1E-14_r8)*h2so4(i,k)**2+(3.2E-14_r8)*h2so4(i,k)*orgforgrowth(i,k) !(19)
                    !nuclrate_pbl(i,k)=(1.4E-14_r8)*h2so4(i,k)**2+(2.6E-14_r8)*h2so4(i,k)*orgforgrowth(i,k) + (0.037E-14_r8)*orgforgrowth(i,k)**2 ! (20)
               else if(pbl_nucleation .EQ. 3) then
                    ! Riccobono 2014:
                  nuclrate_pbl(i,k)=3.27E-21_r8*h2so4(i,k)**2*orgforgrowth(i,k)

               end if

               nuclrate_pbl(i,k)=MAX(MIN(nuclrate_pbl(i,k),1.E10_r8),0._r8)

            else !Not using PBL-nucleation
               nuclrate_pbl(i,k)=0._r8
            end if
            !Size [nm] of particles in PBL
            nuclsize_pbl(i,k)=2._r8

         end do !horizontal points
      end do     !levels


      !-- Calculate total nucleated mass
      do k=1,pver
         do i=1,ncol

            !   Molecular speed and growth rate: H2SO4. Eq. 21 in Kerminen and Kulmala 2002
            vmolh2so4=SQRT(8._r8*shr_const_rgas*t(i,k)/(pi*molmass_h2so4*1.E-3_r8))
            grh2so4(i,k)=(3.E-9_r8/h2so4_dens)*(vmolh2so4*molmass_h2so4*h2so4(i,k))
            grh2so4(i,k)=MAX(MIN(grh2so4(i,k),10000._r8),1.E-10_r8)

            ! TODO: Molecular speed and growth rate: ORG. Eq. 21 in Kerminen and Kulmala 2002

            ! Combined growth rate (cka)
            gr(i,k)=grh2so4(i,k) !+grorg(i,k)

            !-- Lehtinen 2007 parameterization for apparent formation rate
            !   diameters in nm, growth rate in nm h-1, coagulation in s-1
            ! get formation rate: rate of formation of particles with smallest model bin radius #/cm3/s
            call appformrate(nuclsize_bin(i,k), radius*2.E9_r8, nuclrate_bin(i,k), formrate_bin(i,k), coagnuc(i,k), gr(i,k))
            call appformrate(nuclsize_pbl(i,k), radius*2.E9_r8, nuclrate_pbl(i,k), formrate_pbl(i,k), coagnuc(i,k), gr(i,k))

            formrate_bin(i,k)=MAX(MIN(formrate_bin(i,k),1.E3_r8),0._r8)
            formrate_pbl(i,k)=MAX(MIN(formrate_pbl(i,k),1.E3_r8),0._r8)

            !   Number of mol nucleated per g air per second.
            nuclvolume(i,k) = (formrate_bin(i,k) + formrate_pbl(i,k)) & ! [particles/cm3/s]
                            *1.0e6_r8                                 & !==> [particles / m3 /]
                            *2._r8*radius**3*pi/6._r8                 & !==> [m3_{aer} / m3_{air} / sec]
                            / rhoair(i,k)                               !==> m3_{aer} / kg_{air} /sec


            !Estimate how much is organic based on growth-rate
            ! TODO: implement organics here
       !     if(gr(i,k)>1.E-10_r8) then
       !       frach2so4=grh2so4(i,k)/gr(i,k)
       !     else
            frach2so4=1._r8
       !     end if

            ! Nucleated so4 and soa mass mixing ratio per second [kg kg-1 s-1]
            ! used density of particle phase, not of condensing gas
            nuclso4(i,k)=aero_props%density(sulfate_specprop_ndx)*nuclvolume(i,k)*frach2so4
            nuclnum(i,k) = (formrate_bin(i,k) + formrate_pbl(i,k)) * 1.0e6_r8 / rhoair(i,k) ! [#/kg/s]

         end do
      end do

      !-- Diagnostic output
      nuclrate(:,:)=nuclrate(:,:)+(nuclrate_pbl(:,:)+nuclrate_bin(:,:))*dt
      nuclrate_pbl_o(:,:)= nuclrate_pbl_o(:,:)+ nuclrate_pbl(:,:)*dt
      formrate(:,:)=formrate(:,:)+(formrate_pbl(:,:)+formrate_bin(:,:))*dt
      formrate_pbl_o(:,:)=formrate_pbl_o(:,:)+formrate_pbl(:,:)*dt
      grh2so4_o(:,:)=grh2so4_o(:,:)+grh2so4(:,:)*dt
      h2so4nucl_o(:,:)=h2so4nucl_o(:,:)+h2so4pc(:,:)*dt
! TODO: some nice outflds here
   end subroutine aeronucl

!==========================================================================================
! Helper routine to calculate the apparent formation rate of particles of size dx
! from nucleation rate of particles of size d1
!==========================================================================================

   subroutine appformrate(d1, dx, j1, jx, CoagS_dx, gr)
      !-- appformrate calculates the formation rate jx of dx sized particles from the nucleation rate j1 (d1 sized particles)
      !-- Formation rate is parameterized according to Lehtinen et al. (2007), JAS 38:988-994
      !-- Parameterization takes into account the loss of particles due to coagulation
      !-- Growth by self-coagulation is not accounted for
      !-- Typically, 1% of 1 nm nuclei make it to 12 nm
      !-- Written by Risto Makkonen
      ! First estimate: 99% of particles are lost during growth from 1 nm to 12 nm

      !-- Arguments
      real(r8), intent(in)  :: d1                  ! Size of nucleation-sized particles (nm)
      real(r8), intent(in)  :: dx                  ! Size of calculated apparent formation rate (nm)
      real(r8), intent(in)  :: j1                  ! Nucleation rate of d1 sized particles (# cm-3 s-1)
      real(r8), intent(in)  :: CoagS_dx            ! Coagulation term for nucleating particles (s-1)
      real(r8), intent(in)  :: gr                  ! Particle growth rate (nm h-1)
      real(r8), intent(out) :: jx                  ! Formation rate of dx sized particles (# cm-3 s-1)

      !-- Local variables
      real(r8)              :: m
      real(r8)              :: gamma
      real(r8)              :: CoagS_d1            ! Coagulation term for nucleating particles, calculated from CoagS_dx

      !-----------------------------------------------------------------------------------

      ! In Hyytiala, typically 80% of the nuclei are scavenged onto larger background particles while they grow from 1 to 3 nm

      !-- (Eq. 6) Exponent m, depends on background distribution
      ! m=log(CoagS_dx/CoagS_d1)/log(dx/d1)
      ! Or, if we dont want to calculate CoagS_d1, lets assume a typical value for m (-1.5 -- -1.9) and calculate CoagS_d1 from Eq.5
      m=-1.6_r8
      CoagS_d1=CoagS_dx*(d1/dx)**m
      CoagS_d1=MAX(MIN(CoagS_d1,1.E2_r8),1.E-10_r8)

      gamma=(1._r8/(m+1._r8))*((dx/d1)**(m+1._r8)-1._r8)
      gamma=MAX(MIN(gamma,1.E2_r8),1.E-10_r8)

      !-- (Eq. 7) CoagS_d1 is multiplied with 3600 to get units h-1
      jx=j1*exp(-gamma*d1*CoagS_d1*3600._r8/gr)

   end subroutine appformrate

end module condtend