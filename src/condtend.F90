module condtend

!
! script adapted from OsloAeroSec condtend.F90 and aeronucl.F90 by smb
!
! TODO: make copy of state for substepping
! TODO: get rid of: use aero_sectional,     only: secConstIndex index for aerosol in q array
! TODO: get rid of chemistryindex -> find physics index, then - imozart ( chemistryindex = phys_index - imozart + 1 )
! TODO: re-implement organics.. somehow :P
! TODO: make sure unit of bin_mmr is correct, in kg/kg, not number concentration
! TODO: make diagnostic output, addfld and outfld for nucleation and condensation tendencies
! TODO: verify all the calculations, this is just a quick cleanup of the osloaerosec version, and I have not checked the physics yet.
! TODO: check for unused variables, and remove them
! TODO: make an attribute in aerosol properties for the precursor gas, and de CAM-ify the code, so that it is not hardcoded to H2SO4
  ! CAM/NorESM specific stuff

   use phys_control, only: phys_getopts
   use chem_mods,    only: gas_pcnst
   use mo_tracname,  only: solsym
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld, addfld, add_default, fieldname_len, horiz_only
   use physconst,    only: rair, gravit, pi
   use chem_mods,    only: adv_mass             ! molecular weights from mozart
   use ppgrid,       only: pcols, pver, pverp
   use constituents, only: cnst_get_ind
   use spmd_utils,   only: masterproc

   ! Sectional aerosol specific stuff
   use aerosol_properties_mod,       only: aerosol_properties
   use fanci_aerosol_properties_mod, only: fanci_aerosol_properties
   use aerosol_state_mod,            only: aerosol_state, ptr2d_t
   use fanci_aerosol_state_mod,      only: fanci_aerosol_state

   implicit none

   ! public subroutines, used in aero_model
   public :: condensation_init
   public :: advance_aero_condtend

   private

   ! module variables
   real(r8), allocatable :: condsink_per_particle(:) ! [m3/#/s] condensation sink per particle in bin i
   real(r8), allocatable :: bin_centers(:)           ! [m] bin centers
   integer  :: l_h2so4_chem                          ! index for H2SO4 in chemical mechanism
   integer  :: imozart                               ! index shift for MOZART chemistry species
   integer  :: l_h2so4                               ! index for H2SO4 in q array in CAM
   integer  :: sulfate_specprop_ndx                  ! index for sulfate species in aero_props

contains

!==========================================================================================
! Initialize Condensation
!==========================================================================================

   subroutine condensation_init(aero_props)
      ! Find the index for the sulfate species (l_h2so4) in the chemistry
      ! Find sulfate properties from the aerosol properties object
      ! calculate molecular properties (vapor diffusion coefficient, mean free path, thermal velocity)
      ! calculate diffusion coefficient for each size bin, including correction for non-continuum effects
      ! Stores per-particle condensation sink coefficient for each bin (condsink_per_particle)
      !    Needed for vapor loss rate
      ! condensation coefficients:
      ! Theory: Poling et al, "The properties of gases and liquids"
      ! 5th edition, eqn 11-4-4

      ! dummy arguments
      type(fanci_aerosol_properties), intent(in) :: aero_props

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
      real(r8), parameter :: one_third  = 1.0_r8/3.0_r8
      real(r8), parameter :: cm2_to_m2 = 1.e-4_r8       ! convert from cm2 ==> m2

      ! local variables
      real(r8), allocatable :: diffusion_coeff_bin(:)          ! [m2/s] Diffusion coefficient sectional

      real(r8) :: mean_free_path  ! [m] mean free path
      real(r8) :: diffusion_coeff_vap  ! [m2/s] diffusion coefficient for cond. vap (vapor?)
      real(r8) :: Mdual           ! [molec/amu] 1/M_1 + 1/M_2
      real(r8) :: r_molecule      ! [m] radius molecule
      real(r8) :: thermal_vel     ! [m/s] thermal velocity
      integer  :: ibin, ispecprop ! indexes for loops
      character(len=16) :: type   ! type of species in aero_props

      ! local variables for output
      logical  :: history_aerosol
      character(128)                 :: fldname

      !-----------------------------------------------------------------------------------

      if (.not. aero_props%is_active('sulfate')) then
          return
      end if
      ! Initialize only if sulfate is present
      sulfate_specprop_ndx = 0

      do ispecprop = 1, aero_props%nspecies_tot() !TODO: make nicer to have several species
        type = aero_props%spectype(ispecprop)
        if (trim(type) == 'sulfate') then
            sulfate_specprop_ndx = ispecprop
            exit
        end if
      end do

      ! Allocate stuff
      allocate(condsink_per_particle(aero_props%nbins())) ! [m3/#/s]
      allocate(diffusion_coeff_bin(aero_props%nbins()))
      allocate(bin_centers(aero_props%nbins()))

      ! Find indices for precursor gas
      ! Couple the condenseable vapours to chemical species for properties and indexes
      ! add dimension for several species
! TODO: get this from properties?
      call cnst_get_ind('H2SO4',l_h2so4, abort=.true.)           ! CAM tracer for H2SO4 (gas)
      call cnst_get_ind(trim(solsym(1)), imozart, abort=.true.)  ! find first Mozart tracer imozart
      l_h2so4_chem = l_h2so4 - imozart + 1                       ! Find the chemistry index for gas H2SO4

      ! bin center radii [m] from aerosol properties
      bin_centers = aero_props%bin_centers(aero_props%nbins())

      ! pick up densities and weights from aerosol properties
      r_molecule = ( 3.0_r8*aero_props%molecular_weight(sulfate_specprop_ndx) * aunit &
                   / ( 4.0_r8*pi*aero_props%density(sulfate_specprop_ndx) ) )**one_third     ! Radius of molecule (straight forward assuming spherical)
      Mdual  = 2.0_r8/(1.0_r8/Mair+1.0_r8/aero_props%molecular_weight(sulfate_specprop_ndx)) ! factor of [1/m_1 + 1_m2]

      ! thermal velocity for H2SO4 in air (m/s)
      ! https://en.wikipedia.org/wiki/Thermal_velocity
      thermal_vel = sqrt(8.0_r8*boltz*t0/ ( pi*aero_props%molecular_weight(sulfate_specprop_ndx) * aunit ) )

      ! calculating microphysical parameters from equations in Ch. 8 of Seinfeld & Pandis (1998):
      ! mean free path for molec in air (m)
      mean_free_path = 1.0_r8 / ( pi*sqrt( 1.0_r8+aero_props%molecular_weight(sulfate_specprop_ndx)/Mair ) &
                    * (radair+r_molecule)**2 * p0/(boltz*t0) )

      ! Solve eqn 11-4.4 in Poling et al
      ! (A bit hard to follow units here, but result in the book is in cm2/s)..
      ! so scale by "cm2_to_m2" to get m2/sec

      diffusion_coeff_vap = cm2_to_m2   &
         * 0.00143_r8 * t0**1.75_r8     &
         / ((p0/1.0e5_r8) * sqrt(Mdual)   &
         * (((Vad)**one_third+(Vadair)**one_third)**2))

      do ibin = 1, aero_props%nbins()         !all bins receive condensation
          ! Correct for non-continuum effects, formula is from
          ! Chuang and Penner, Tellus, 1995, sticking coeffient from
          ! Vignati et al, JGR, 2004
          diffusion_coeff_bin(ibin) = diffusion_coeff_vap  &    !original diffusion coefficient
               / ( bin_centers(ibin) / (bin_centers(ibin) + mean_free_path )  &  ! non-continuum correction factor
               + 4.0_r8 * diffusion_coeff_vap / ( 1._r8 * thermal_vel * bin_centers(ibin) ) )
      enddo

      !Find sink per particle in bin "ibin"
      !Eqn 13 in Kulmala et al, Tellus 53B, 2001, pp 479
      !http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.530411.x/abstract

      condsink_per_particle = 0.0_r8

      do ibin = 1, aero_props%nbins()
         ! Since we do not sum over bins, it is slightly different than above (normnk=1, no summing)
         ! Condensation sink per particle in bin "ibin" [m3/s] per particle or [m3/#/s]
         condsink_per_particle(ibin) = 4.0_r8*pi * diffusion_coeff_bin(ibin) & ! [m2/s] diffusion coefficient
                                       * bin_centers(ibin)                     ! [m] radius of bin
      end do

      ! Initialize output
      call phys_getopts(history_aerosol_out = history_aerosol)

! condensation tendency for SO4_R1, maybe num_1, ...
    !Need to add so4_a1, soa_na, so4_na, soa_a1 also (which are not parts of the donor-receiver stuff)
    !aero_state%state%q(:ncol,:pver, aero_state%num_transport_ndx(ibin)) = &
    !aero_state%state%q(:ncol,:pver, aero_state%aero_range_state(irange)%transport_ndx(:)) = &

    ! TODO: for future maybe only add fields for the smallest bin/range
    do ibin = 1, aero_props%nbins()
        fldname = 'num'//trim(int2str(ibin))//'_condTend'
        call addfld( trim(fldname), horiz_only, 'A', "1/kg/s", "condensation/nucleation tendency of number in bin "//trim(int2str(ibin)))
        if(history_aerosol)then
            call add_default( trim(fldname), 1, ' ' )
        end if
    end do

    do irange = 1, aero_props%nranges()
        ! Currently only for sulfate!
        fldname = 'SO4_R'//trim(int2str(irange))
        call addfld( trim(fldname)//'_condTend', horiz_only, 'A', "kg/kg/s" , "condensation/nucleation tendency of SO4 mmr in range R"//trim(int2str(irange)) )
        if(history_aerosol)then
           call add_default( trim(fldname)//'_condTend', 1, ' ' )
        end if
    end do

   end subroutine condensation_init

!==========================================================================================
! Condensation + Nucleation called from aero_model
! manages time sub-stepping by calling condtend_substep
!==========================================================================================

   subroutine advance_aero_condtend(lchnk, q, cond_vap_gasprod, temperature, &
               pmid, pdel, dt, ncol, pblh, zm, qh20, aero_props, aero_state)
      ! Outer driver for aerosol condensation and nucleation processes
      !    Calls condtend_substep, retries if substep is too large
      ! Calculates nucleation rate and condensation rate of aerosols
      ! Accumulates tendencies and updates aerosol_state
! TODO: bin_numconc_new is modified but update_bin passes delmmr_sum=0 and tnd_ndx = 0
      ! This method calls the condtend method. If the timestep needs to be split in two,
      ! this will be done here and the condtend method will be called several times.
      ! This method also writes output once condend is done.

      use cam_abortutils, only: endrun
      use cam_history,    only: outfld,fieldname_len
      use constituents,   only: pcnst  ! h2so4 and soa nucleation (cka)

      ! dummy arguments
      type(fanci_aerosol_properties), intent(in) :: aero_props
      type(fanci_aerosol_state), intent(inout)   :: aero_state

      integer,  intent(in)    :: lchnk, ncol           ! chunk identifier, nr of columns
      real(r8), intent(in)    :: temperature(:,:)      ! Temperature (K)
      real(r8), intent(in)    :: pmid(:,:)             ! [Pa] pressure at mid point
      real(r8), intent(in)    :: pdel(:,:)             ! [Pa] difference in grid cell
      real(r8), intent(in)    :: cond_vap_gasprod(:,:) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in)    :: dt                    ! Time step
      real(r8), intent(in)    :: pblh(:)               ! pbl height (m)
      real(r8), intent(in)    :: zm(:,:)               ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(:,:)             ! specific humidity (kg/kg)
      real(r8), intent(inout) :: q(:,:,:)              ! TMR [kg/kg] including moisture

      integer, parameter :: max_substeps = 1024 ! TODO: makes sense?

      real(r8) :: rho_air(pcols,pver)                  ! [kg/m3] density of air
      real(r8) :: dt_substep
      real(r8) :: substep_output(pcols, gas_pcnst)

      real(r8), allocatable :: q_initial(:,:,:)      ! mass before subroutine.
      real(r8), allocatable :: bin_numconc_initial(:,:,:)
      real(r8), allocatable :: bin_numconc_new(:,:,:)  ! [#/m3] number concentration new
      real(r8), pointer     :: tmp_num(:,:)

      !output:
      real(r8) :: coltend (pcols, gas_pcnst)      ! ??

      real(r8) :: nuclrate_accum(pcols,pver)      ! [kg/kg] tracer lost
      real(r8) :: nuclrate_accum_pbl(pcols,pver)  ! [kg/kg] tracer lost
      real(r8) :: formrate_accum(pcols,pver)      ! [kg/kg] tracer lost
      real(r8) :: formrate_accum_pbl(pcols,pver)  ! [kg/kg] tracer lost
      real(r8) :: h2so4nucl(pcols,pver)           ! h2so4 in nucleation code
      real(r8) :: growth_rate_h2so4(pcols,pver)   ! growth rate h2so4
      real(r8) :: coag_sink_nucl(pcols,pver)            ! coagulation in nucleation
      real(r8) :: nuclnum_tend(pcols,pver)        ! [#/m3/s] nucleation number rate from RM's parameterization

      integer  :: ibin, irange ! index for bin and range
      logical  :: split_dt ! whether timestep is split or not
      integer  :: nsubsteps, isubstep

      ! aero props for sulfate
      real(r8), allocatable :: nucleated_numconc_added(:,:,:)
      real(r8), allocatable :: range_mmr_tend(:,:,:)
      real(r8) :: sulfate_density
      integer :: sulfate_range_species_ndx
      integer :: nrange_species
      integer :: ispec

      ! local variables for output
      logical                        :: history_aerosol
!      character(128)                 :: long_name
!      character(8)                   :: unit
!      character(18) :: fieldname_receiver

      !-----------------------------------------------------------------------------------
      if (.not. aero_props%is_active('sulfate')) then
          return
      end if

      ! TODO: an attempt to make this slightly less cam dependent
      allocate(q_initial(ncol, size(q,2), size(q,3)))
      allocate(bin_numconc_initial(ncol, pver, aero_props%nbins()))
      allocate(bin_numconc_new(ncol, pver, aero_props%nbins()))
      allocate(nucleated_numconc_added(ncol, pver, aero_props%nbins()))

      nucleated_numconc_added = 0._r8

      rho_air(:ncol, :) = pmid(:ncol,:) / (rair * temperature(:ncol,:))

      q_initial(:ncol,:,:) = q(:ncol,:,:)

      sulfate_density = aero_props%density(sulfate_specprop_ndx)

      nullify(tmp_num)
      ! get bin number concentrations from number mixing ratios
      do ibin = 1, aero_props%nbins()
         call aero_state%get_ambient_num(ibin, tmp_num)
         bin_numconc_initial(:,:,ibin) = tmp_num(:,:) * rho_air(:,:)
         nullify(tmp_num)
      end do

      nsubsteps = 2

      do
         if (nsubsteps > max_substeps) then
            call endrun('advance_aero_condtend: exceeded maximum number of substeps')
         end if

         q(:ncol,:,:) = q_initial(:ncol,:,:)
         bin_numconc_new = bin_numconc_initial

         nuclrate_accum_pbl(:,:) = 0.0_r8
         nuclrate_accum(:,:)     = 0.0_r8
         formrate_accum_pbl(:,:) = 0.0_r8
         formrate_accum(:,:)     = 0.0_r8
         h2so4nucl(:,:)          = 0.0_r8
         growth_rate_h2so4(:,:)  = 0.0_r8
         coag_sink_nucl(:,:)     = 0.0_r8

         dt_substep = dt / real(nsubsteps, r8)
         split_dt = .false.

         do isubstep = 1, nsubsteps
            substep_output(:,:) = 0._r8

           ! run until no need to split time any longer
           call condtend_substep(lchnk, q, cond_vap_gasprod, temperature, &
                           nuclrate_accum, nuclrate_accum_pbl, &
                           formrate_accum, formrate_accum_pbl, coag_sink_nucl, &
                           h2so4nucl, growth_rate_h2so4, &
                           substep_output, split_dt, &
                           pmid, pdel, dt_substep, ncol, pblh, zm, qh20, &
                           aero_props, aero_state, bin_numconc_new, nuclnum_tend)

           if (split_dt) then
           ! we need to increase the number of substeps

              exit
           end if

         end do

         if (.not. split_dt) then
            ! substep completed without needing to split further
            exit
         end if

         if (nsubsteps > max_substeps / 2) then
            call endrun('advance_aero_condtend: exceeded maximum number of substeps')
         end if

         nsubsteps = 2 * nsubsteps

      end do

     ! update aerosol_state
     ! numconc #/m3-air back to num mixing ratio #/kg-air.
     do ibin = 1, aero_props%nbins()
        aero_state%state%q(:ncol,:pver,aero_state%num_transport_ndx(ibin)) = &
                   bin_numconc_new(:ncol,:,ibin) / rho_air(:ncol,:)
     end do

     do irange = 1, aero_props%nranges()

   ! TODO: use update_bin:
                   !do i = 1, ncol
                   !do k = 1, pver
                   !    call aero_state%update_bin(1, i, k, 0._r8, -nuclnum_tend(i,k)*dt, 0, dt, tend)
                   !end do
                !end do
     end do

     ! Accumulated values for diagnostics * substep duration.
     ! Divide by the full timestep -> avg for dt
     nuclrate_accum     = nuclrate_accum     / dt
     nuclrate_accum_pbl = nuclrate_accum_pbl / dt
     formrate_accum     = formrate_accum     / dt
     formrate_accum_pbl = formrate_accum_pbl / dt
     h2so4nucl          = h2so4nucl          / dt
     growth_rate_h2so4  = growth_rate_h2so4  / dt
     coag_sink_nucl     = coag_sink_nucl     / dt

               ! update bin with numconc from nucleation
! TODO: add nucleation
! TODO: add ranges

      call phys_getopts(history_aerosol_out = history_aerosol)

      if(history_aerosol)then
      ! Write output
         do ibin = 1, aero_props%nbins()
            fieldname = 'num_'//trim(aero_props%bin_name(ibin))//'_condTend'
            call outfld(fieldname, coltend(:ncol, aero_state%num_transport_ndx(ibin)-imozart), pcols, lchnk)

    !aero_state%state%q(:ncol,:pver, aero_state%num_transport_ndx(ibin)) = &
    !aero_state%state%q(:ncol,:pver, aero_state%aero_range_state(irange)%transport_ndx(:)) = &
         end do

         do irange = 1, aero_props%nranges()
            fldname = 'SO4_R'//trim(aero_props%range_name(irange))//'_condTend'
            call outfld(fldname, coltend(:ncol, aero_state%aero_range_state(irange)%transport_ndx(sulfate_specprop_ndx)-imozart), pcols, lchnk)
         end do
      end if

      deallocate(bin_numconc_initial, bin_numconc)

   end subroutine advance_aero_condtend

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
! 1. calculate air density from pressure and temperature
! 2. Calculate each bin's particle concentration
! 3. Multiply the per-particle sink by concentration, then normalize to get fraction of condensation going to each bin
! 4. Calculates an implicit updated H₂SO₄ concentration after gas production and condensation.
! 5. Calls aerosol_nucleation to calculate nucleation and formation rates.
! 6. Accounts for H₂SO₄ consumed by nucleation and condensation.
! 7. Calls move_mass_and_number to move particles as their size grows.
! 8. Adds newly formed particles to the smallest bin.
!==========================================================================================

   subroutine condtend_substep(lchnk,  q, cond_vap_gasprod, temperature,            &
                nuclrate_accum, nuclrate_accum_pbl_o, formrate_accum, formrate_accum_pbl_o, coagnucl_o,  &
                h2so4nucl_o, growth_rate_h2so4_o,                     &
                coltend_o, split_dt,                                            &
                pmid, pdel, dt, ncol, pblh,zm,qh20,                              &
                aero_props, aero_state, bin_numconc_new, nuclnum_tend)

      use cam_history,     only: outfld,fieldname_len
      use constituents,    only: pcnst

      ! dummy arguments
      type(fanci_aerosol_properties), intent(in) :: aero_props
      type(fanci_aerosol_state), intent(in) :: aero_state

      real(r8), intent(inout) :: bin_numconc_new(:,:,:)      ! [#/m3] number concentration new

       !++smb sectional
      real(r8), intent(inout) :: nuclrate_accum (:,:)       ! Nucleation rate output
      real(r8), intent(inout) :: nuclrate_accum_pbl_o (:,:) ! Nucleation rate pbl output
      real(r8), intent(inout) :: formrate_accum(:,:)        ! Formation rate output
      real(r8), intent(inout) :: formrate_accum_pbl_o(:,:)  ! Formation rate pbl output

      real(r8), intent(inout) :: coagnucl_o(:,:)            ! Coagulation sink for npf output
      real(r8), intent(inout) :: h2so4nucl_o(:,:)           ! H2SO4 for nucleation output
      real(r8), intent(inout) :: growth_rate_h2so4_o(:,:)   ! GR from H2SO4 output

      real(r8), intent(out)   :: coltend_o(:,:)             ! column tendency output
      logical,  intent(out)   :: split_dt                   ! if true, time step needs to be split

      real(r8), intent(out)   :: nuclnum_tend(:,:)          ! [#/m3/s] nucleation number rate from RM's parameterization
       !--smb sectional

      ! arguments
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: ncol                  ! number of columns
      real(r8), intent(in) :: temperature(:,:)      ! Temperature (K)
      real(r8), intent(in) :: pmid(:,:)             ! [Pa] pressure at mid point
      real(r8), intent(in) :: pdel(:,:)             ! [Pa] difference in grid cell
      real(r8), intent(inout) :: q(:,:,:)           ! TMR [kg/kg] including moisture

      real(r8), intent(in) :: cond_vap_gasprod(:,:) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in) :: dt                    ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in) :: pblh(pcols)        ! pbl height (m)
      real(r8), intent(in) :: zm(:,:)            ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in) :: qh20(:,:)          ! specific humidity (kg/kg)

      ! local
      character(len=fieldname_len+3) :: fieldname
      integer :: i,k
      real(r8) :: condensation_sink_tot(pcols,pver)     ! [1/s] sum of condensation sink
      real(r8) :: totalLoss(pcols,pver,gas_pcnst)       ! [kg/kg] tracer lost
      real(r8) :: coltend(pcols, gas_pcnst)
      real(r8) :: tracer_coltend(pcols)
      real(r8), allocatable :: bin_condensation_sink(:) ! [1/s] loss rate per mode (mixture)
      real(r8), allocatable :: bin_condensation_fraction(:,:,:) ! [frc]
      real(r8), allocatable :: bin_numconc_old(:,:,:)   ! [#/m3] number concentration
      real(r8), allocatable :: tend(:,:,:)

      real(r8)       :: h2so4_mmr_tmp(pcols,pver)       ! intermediate concentration for h2so4
      real(r8)       :: rho_air(pcols,pver)             ! [kg/m3] density of air
      ! Volume of added  material from condensate;  surface area of core particle;
      ! TODO: add real(r8)       :: volume_shell, area_core,vol_monolayer

      ! needed for h2so4 and soa nucleation treatment
       integer  :: modeIndexReceiverCoag                ! Index of modes receiving coagulate
       integer  :: iCoagReceiver                        ! counter for species receiving coagulate
       real(r8) :: coagulationSink(pcols,pver)          ! [1/s] coaglation loss for SO4_n and soa_n
       real(r8), pointer :: tmp_num(:,:)                ! gets aerosol number mixing ratio through aero_state [#/kg_air]
       real(r8) :: h2so4_mmr_loss(pcols,pver)           ! [kg/kg] budget terms on H2SO4 (gas)
       real(r8) :: frac_loss_to_nucl(pcols,pver)        ! [frc] fraction of gas nucleated
       real(r8) :: firstOrderLossRateNucl(pcols,pver)   ! [1/s] first order loss rate due to nucleation
       real(r8) :: nuclso4_tend(pcols,pver)             ! [kg/kg/s] Nucleated so4 mass tendency from RM's parameterization
       real(r8), allocatable :: rho_aerosol(:, :, :)
       !smb++ sectional
       real(r8) :: dummy  !
       integer  :: ibin, irange ! indices
       real(r8) :: q_mmr(pcols,pver,aero_props%nbins()) ! [kg/kg] mass in each tracer -> rename to bin_mmr

       ! local variables for output
       logical        :: history_aerosol
!       character(128) :: long_name                     ! [-] needed for diagnostics

      !-----------------------------------------------------------------------------------

       allocate(bin_condensation_sink(aero_props%nbins()))
       allocate(bin_numconc(ncol,pver,aero_props%nbins()))
       allocate(bin_condensation_fraction(ncol,pver,aero_props%nbins()))
       allocate(rho_aerosol(ncol,pver,aero_props%nbins()))
       allocate(tend(ncol, pver, gas_pcnst))
       !Initialize h2so4 and soa nucl variables
       coagulationSink = 0.0_r8
       bin_condensation_fraction = 0.0_r8
       bin_numconc_old = 0.0_r8
       tmp_num => null()
       coltend_o = 0.0_r8
       split_dt       = .FALSE.

       do k = 1, pver
           do i = 1, ncol
                !Air density
                rho_air(i,k) = pmid(i,k)/rair/temperature(i,k)
           end do
        end do

       do ibin = 1, aero_props%nbins()
        ! No looping through species, mmr is added afterwards
!          call aero_state%get_ambient_num(ibin, tmp_num)
          bin_numconc_old(:ncol,:,ibin) = bin_numconc_new(:ncol, :, ibin) !* rho_air(:ncol,:)  ![#/m3] number concentration
          rho_aerosol(:ncol,:,ibin) = aero_state%bin_dry_density(ibin, ncol)
       enddo

       do k=1,pver
           do i=1,ncol

                !smb++ sectional
                ! initialize condensation sink for sectional
                bin_condensation_sink = 0.0_r8  !Sink to the coming "receiver" of any vapour
                !smb-- sectional

                !NB: The following is duplicated code, coordinate with koagsub!!
                !Initialize number concentration for this receiver


                !Go though all bins receiving condensation

                ! condensation sink to sectional bin:
                   do ibin = 1, aero_props%nbins()

                      !This is the loss rate a gas molecule will see due to aerosol surface area
                      bin_condensation_sink(ibin)   = condsink_per_particle(ibin)  & ![m3/#/s]  per particle
                                                          * bin_numconc_old(i,k,ibin)             ![#/m3]
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
                        condensation_sink_tot(i,k) = sum(bin_condensation_sink( :))
                        ! Keeps track of the fraction of the condensate to the sectional bins for each tracer
                        bin_condensation_fraction(i,k,:) = bin_condensation_sink( :) &
                                /(condensation_sink_tot(i,k)+1.e-30_r8)![frc]
                    !smb-- sectional


                !Solve the intermediate (end of timestep) concentration using
                !euler backward solution C_{old} + P *dt - L*C_{new}*dt = C_{new} ==>
                !Cnew -Cold = prod - loss ==>
                h2so4_mmr_tmp(i,k) = &
                                     ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                                     / (1.0_r8 + condensation_sink_tot(i,k)*dt)


                ! TODO: add organics contribution to nucleation ORG_LV

                modeIndexReceiverCoag = 0
                !Sum coagulation sink for nucleated so4 and soa particles over all receivers of coagulate. Needed for RM's nucleation code
                !OBS - looks like RM's coagulation sink is multiplied by 10^-12??


           end do !index i
       end do !index k

       !Calculate nucleated masses of so4 and soa (nuclso4_tend)
       !following RM's parameterization (cka)
! TODO: check, how is h2so4_mmr_tmp h2so4_mmr
       call aerosol_nucleation(lchnk, ncol, temperature, pmid, qh20, h2so4_mmr_tmp(:,:), &
                   coagulationSink, nuclnum_tend, nuclso4_tend, zm, pblh, &
                   nuclrate_accum, nuclrate_accum_pbl_o, formrate_accum, formrate_accum_pbl_o, &
                   h2so4nucl_o, growth_rate_h2so4_o, dt, &
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
              if (h2so4_mmr_tmp(i,k) .eq. 0) then
                  firstOrderLossRateNucl(i,k)=0._r8
              else
                  firstOrderLossRateNucl(i,k) = nuclso4_tend(i,k)/h2so4_mmr_tmp(i,k)
              end if
! organics here
              !smb-- added check for 0
             !First order loss rate (1/s) for nucleation

                !Solve implicitly (again)
                !C_new - C_old =  PROD_{gas} - CS*C_new*dt - LR_{nucl}*C_new =>
                h2so4_mmr_tmp(i,k) = &
                               ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                               / (1.0_r8 + condensation_sink_tot(i,k)*dt + firstOrderLossRateNucl(i,k)*dt)

                !fraction nucleated
                frac_loss_to_nucl(i,k) = firstOrderLossRateNucl(i,k) &
                                     /(firstOrderLossRateNucl(i,k) + condensation_sink_tot(i,k))
               !From budget, we get: lost = prod -cnew + cold
                h2so4_mmr_loss(i,k) = cond_vap_gasprod(i,k)*dt   & !Produced
                                     + q(i,k,l_h2so4_chem)            & !cold
                                     - h2so4_mmr_tmp(i,k)    !cnew

             !H2SO4 condensate
             do ibin=1, aero_props%nbins()
                ! bin_numconc_old (#/m3)
                ! h2so4_mmr_loss: [kg/kg] mass lost from gas phase
                ! bin_condensation_fraction: [frc] fraction to the particular bin
                ! frac_loss_to_nucl: [frc] fraction nucleated
                ! get number from h2so4_mmr_loss
! TODO: fix update_bin -> mmr to correct range
! TODO: condensation does not change the number!! -> add mass to ranges
!                bin_numconc_old(i, k, ibin) = bin_numconc_old(i, k, ibin) + h2so4_mmr_loss(i,k) * (1._r8-frac_loss_to_nucl(i,k)) * bin_condensation_fraction(i,k, ibin) / aero_props%density(sulfate_specprop_ndx) / aero_props%particle_volume(ibin)
!                call aero_state%update_bin(ibin, i, k, 0._r8, -1._r8*h2so4_mmr_loss(i,k) * (1._r8-frac_loss_to_nucl(i,k)) * bin_condensation_fraction(i,k, ibin) / aero_props%density(sulfate_specprop_ndx) / aero_props%particle_volume(ibin), 0, dt, tend)
                irange = aero_props%bins2ranges(ibin)

!                ispec = ??
                ! add mass to the range:
!                aero_state%aero_range_state(irange)%mmr(i,k,ispec) = aero_state%aero_range_state(irange)%mmr(i,k,ispec) &
!                     + h2so4_mmr_loss(i,k)*(1.0_r8-frac_loss_to_nucl(i,k)) &
!                     *bin_condensation_fraction(i,k, ibin)

             end do

             !condenseable vapours
             q(i,k,l_h2so4_chem)  = h2so4_mmr_tmp(i,k)

             !smb++sectional grow particles in sectional scheme:
    ! removed median radius because not used
             do ibin=1, aero_props%nbins()
                 ! get mmr for sulfate per bin?
! TODO: use chem_ndx or whatever
                 q_mmr(i,k,ibin) = q(i,k,XXX) * aero_props%particle_volume(ibin) * rho_aerosol(i, k, ibin)    ! [kg/kg] mass in each tracer
             end do

             call move_mass_and_number(q_mmr(i,k,:), bin_numconc_old(i,k,:), &
                            rho_air(i,k), split_dt, aero_props, rho_aerosol(i, k,:), q_mmr_new(i,k,:), bin_numconc_new(i, k,:))
             ! Add nucleated number to smallest bin after redistribution (#/kg/s * s = #/kg)
             bin_numconc_new(i, k, 1) = bin_numconc_new(i, k, 1) + nuclnum_tend(i,k) * dt
             ! Add nucleated mass to first bin of sectional scheme:
 !            q(i,k,chemistryIndex(secConstIndex(1,1))) =  q(i,k,chemistryIndex(secConstIndex(1,1)))       &
 !                        + h2so4_mmr_loss(i,k)*frac_loss_to_nucl(i,k)

            ! Add mass from sectional scheme to so4_na and soa_na:
!             q(i,k,chemistryIndex(l_so4_na)) = q(i,k,chemistryIndex(l_so4_na))

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
!                                                    h2so4_mmr_loss(:ncol,:)           &
!                                                    *frac_loss_to_nucl(:ncol,:)*pdel(:ncol,:) , 2 &
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
!                                                    h2so4_mmr_loss(:ncol,:)           &
                                                    !smb++sectional subtract fraction to sectional scheme
!                                                    *(1-sum(bin_condensation_fraction(:ncol,:,:),3)) &
                                                    !smb--sectional
!                                                    *(1.0_r8 - frac_loss_to_nucl(:ncol,:))*pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt

          !Take into account soa_lv (gas) nucleated in budget
          !smb++ sectional
! organics here
          !smb++ sectional: putt in how condenses on sectional:
!          do ibin=1,aero_props%nbins()
!               coltend(:ncol,chemistryIndex(secConstIndex(1,ibin))) = coltend(:ncol, chemistryIndex(secConstIndex(1,ibin)))+ &
!                                                    sum(                    &
!                                                    h2so4_mmr_loss(:ncol, :)    &
!                                                    *(bin_condensation_fraction(:ncol,:,ibin)) &
!                                                    *(1.0_r8 - frac_loss_to_nucl(:ncol, :))*pdel(:ncol,:),2) &
!                                                    /gravit/dt
!          end do
! organics here

!          coltend_o(:,:)=coltend(:,:)

!       endif


   end subroutine condtend_substep

!==========================================================================================
! Helper routine to move mass between bins from growth
!==========================================================================================

   subroutine move_mass_and_number(bin_mmr_old, bin_numconc_old, rho_air, decrease_dt, aero_props, &
                                     rho_aerosol, bin_mmr_new, bin_numconc_new)

    ! bin_mmr: mass in each tracer (kg_aer/kg_air) after growth
    ! Moves tracer mass from on bin to the other based on condensational/coagulation growth.
    ! Based on Jacobson Fundamentals of Atmospheric Modeling, second edition (2005),
    ! Chapter   13.5
    class(fanci_aerosol_properties), intent(in) :: aero_props

    ! inputs
    real(r8), intent(in)    :: bin_mmr_old(:)        ! mass in each tracer after growth, per bin (kg_aerosol/kg_air)
    real(r8), intent(in)    :: bin_numconc_old(:)    ! number concentration, after growth, before redistribution (particles/m^3)
    real(r8), intent(in)    :: rho_aerosol(:)        ! Density of particle (kg_aerosol/m3_aerosol)
    real(r8), intent(in)    :: rho_air               ! Density of air (kg_air/m3_air)

    ! outputs
    real(r8), intent(out)   :: bin_mmr_new(:)        ! mass in each tracer after growth, per bin (kg/kg_air)
    real(r8), intent(out)   :: bin_numconc_new(:)    ! output particle concentrations, after growth
    logical,  intent(out)   :: decrease_dt           ! if set to True, time step is divided
                                                     ! and the procedure is re run
    ! local variables
    real(r8)                :: particle_volume_grown ! volume after growth dimension(secNrBins)
    real(r8)                :: fraction_remaining    ! fraction to stay in bin
    integer                 :: ibin
    !-----------------------------------------------------------------------------------

    ! initialize
    bin_mmr_new(:) = 0._r8
    bin_numconc_new(:) = 0._r8
    decrease_dt=.FALSE.

    !compute volume in each bin with condensation (by mass) and by number concentration
    do ibin = 1, aero_props%nbins() - 1
        particle_volume_grown = 0._r8

        if (bin_numconc_old(ibin)<1.e-30_r8) then
            fraction_remaining = 1.0_r8
        else
            ! volume of single particle after growth
            particle_volume_grown = bin_mmr_old(ibin) * rho_air / rho_aerosol(ibin) &
                         / bin_numconc_old(ibin)

            ! Check that the particles don't grow beyond the next bin, decrease time step
            if (particle_volume_grown >= aero_props%particle_volume(ibin + 1)) then
                decrease_dt = .TRUE.
                return
            end if

            ! fraction to stay in bin
            fraction_remaining = (aero_props%particle_volume(ibin + 1) - particle_volume_grown) &
                  / (aero_props%particle_volume(ibin + 1) - aero_props%particle_volume(ibin))


            fraction_remaining=max(0._r8, min(1._r8,fraction_remaining))
        end if

        ! Calculate new number concentrations and mass mixing ratios
        bin_numconc_new(ibin)     = bin_numconc_new(ibin)     + fraction_remaining * bin_numconc_old(ibin)
        bin_numconc_new(ibin + 1) = bin_numconc_new(ibin + 1) + (1._r8 - fraction_remaining) * bin_numconc_old(ibin)

        bin_mmr_new(ibin)     = bin_mmr_new(ibin)     + fraction_remaining * bin_mmr_old(ibin)
        bin_mmr_new(ibin + 1) = bin_mmr_new(ibin + 1) + (1._r8 - fraction_remaining) * bin_mmr_old(ibin)
    end do

    ! Aerosol cannot leave the largest bin -> keep it in there
    ibin = aero_props%nbins()
    bin_numconc_new(ibin) = bin_numconc_new(ibin) + bin_numconc_old(ibin)
    bin_mmr_new(ibin) = bin_mmr_new(ibin) + bin_mmr_old(ibin)

   end subroutine move_mass_and_number

!==========================================================================================
! Helper routine for aerosol nucleation - Vehkamäki et al. (2002)
!==========================================================================================

   subroutine aerosol_nucleation(lchnk, ncol, t, pmid, h2ommr, h2so4_mmr, &
                coag_sink_nucl, nuclnum_tend, nuclso4_tend, zm, pblht, &
                nuclrate_accum, nuclrate_accum_pbl_o, formrate_accum, formrate_accum_pbl_o,  &
                h2so4nucl_o, growth_rate_h2so4_o, dt, &
                target_particle_radius, aero_props)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use wv_saturation,  only: qsat_water
    use physconst,      only: avogad, rair
    use ppgrid,         only: pcols, pver, pverp
    use cam_history,    only: outfld
    use phys_control,   only: phys_getopts
    use chem_mods,      only: adv_mass
    use m_spc_id,       only: id_H2SO4
    use shr_const_mod,  only: shr_const_rgas

    !-- Arguments
    class(fanci_aerosol_properties), intent(in) :: aero_props

    real(r8), intent(in)    :: dt                        ! timestep (output is weighted by this)
    integer,  intent(in)    :: lchnk                     ! chunk identifier
    integer,  intent(in)    :: ncol                      ! number of atmospheric column
    real(r8), intent(in)    :: pmid(pcols,pver)          ! layer pressure (Pa)
    real(r8), intent(in)    :: h2ommr(:,:)               ! layer specific humidity
    real(r8), intent(in)    :: t(:,:)                    ! Temperature (K)
    real(r8), intent(in)    :: h2so4_mmr(:,:)            ! Sulphuric acid (gas) mass mixing ratio (kg kg-1)
    real(r8), intent(in)    :: coag_sink_nucl(:,:)       ! Coagulation sink for nucleating particles [1/s]

    real(r8), intent(in)    :: zm(:,:)                   ! Height at layer midpoints (m)
    real(r8), intent(in)    :: pblht(:)                  ! Planetary boundary layer height (m)
    real(r8), intent(in)    :: target_particle_radius    ! Particle size at calculated formation rate [m]

    ! because the timestep may be divided in two, the output needs to be averaged over all timesteps
    ! therefore we track this in these variables
    real(r8), intent(inout) :: nuclrate_accum(:,:)       ! Nucleation rate output
    real(r8), intent(inout) :: nuclrate_accum_pbl_o(:,:) ! Nucleation in pbl rate output
    real(r8), intent(inout) :: formrate_accum(:,:)       ! formation rate output
    real(r8), intent(inout) :: formrate_accum_pbl_o(:,:) ! formation rate in pbl output

    real(r8), intent(inout) :: h2so4nucl_o(:,:)          ! conc h2so4 for output
    real(r8), intent(inout) :: growth_rate_h2so4_o(:,:)  ! GR from H2SO4 for output
    !smb-- sectional
    real(r8), intent(out)   :: nuclso4_tend(:,:)         ! Nucleated mass (H2SO4)
    real(r8), intent(out)   :: nuclnum_tend(:,:)         ! [#/m3/s] Nucleated mass (SO4)

! TODO: implement oxidorg, orgnucl_o, grsoa_o, nuclorg when adding organics

    !-- Local variables
    ! TODO: is there a way to not hard code the H2SO4 density?
    real(r8), parameter   :: h2so4_dens=1841._r8            ! h2so4 density [kg m-3]

    integer               :: i,k
    real(r8)              :: qs(pcols,pver)                 ! Saturation specific humidity
    real(r8)              :: relhum(pcols,pver)             ! Relative humidity
    real(r8)              :: h2so4_num_cm3(pcols,pver)      ! Sulphuric acid concentration [#/cm3]
    real(r8)              :: nucleated_volume(pcols,pver)   ! [m3/m3/s] Nucleated mass (SO4)
    real(r8)              :: rho_air(pcols,pver)            ! density of air [kg/m3] !cka
    real(r8)              :: pblht_lim(pcols)               ! Planetary boundary layer height (m) (500m<pblht_lim<7000m) (cka)

    real(r8)              :: nuclrate_accum_bin(pcols,pver) ! Binary nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_accum_bin(pcols,pver) ! Binary formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_bin(pcols,pver)       ! Binary nucleation critical cluster size (m)
    real(r8)              :: nuclrate_accum_pbl(pcols,pver) ! Boundary layer nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_accum_pbl(pcols,pver) ! Boundary layer formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_pbl(pcols,pver)       ! Boundary layer nucleation formation size (m)

! TODO: introduce orgforgrowth again, when implementing organics
    real(r8)              :: gr(pcols,pver), growth_rate_h2so4(pcols,pver) !growth rates
    real(r8)              :: vmolh2so4                      ! [m/s] molecular speed of condenseable gases
    real(r8)              :: frach2so4
    real(r8)              :: dummy

    integer               :: atm_nucleation                 ! Nucleation parameterization for the whole atmosphere
    integer               :: pbl_nucleation                 ! Nucleation parameterization for the boundary layer
    real(r8)              :: molmass_h2so4                  ! molecular mass of h2so4 [g/mol]
!TODO: fix this l_so4_na -> use aerosol tracers from state
 !   integer               :: l_so4_na
   ! Variables for binary nucleation parameterization
    real(r8)              :: zrhoa, zrh, zt, zt2, zt3, zlogrh, zlogrh2, zlogrh3, zlogrhoa, zlogrhoa2, zlogrhoa3, x, zxmole, zix
    real(r8)              :: zjnuc, zntot, zrc, zrxc

    !-----------------------------------------------------------------------------------

   !cka: OBS    call phys_getopts(pbl_nucleation_out=pbl_nucleation, atm_nucleation_out=atm_nucleation)
    !cka: testing by setting these flags:
    pbl_nucleation = 1 ! smb++ 3 use Riccobono 2014 for nucleation. -> cvb++ use 2, because we don't currently have organics
    atm_nucleation = 1

    nuclso4_tend(:,:) = 0._r8

    !-- The highest level in planetary boundary layer
    do i=1,ncol
        pblht_lim(i) = MIN( MAX( pblht(i), 500._r8), 7000._r8 )
    end do

    !-- Get molecular mass of h2so4 and soa_lv (cka)
    molmass_h2so4 = adv_mass(id_H2SO4)

! TODO: make this something sensible:
  !  l_so4_na = 1

    !-- Formation diameters (m). Nucleated particles are inserted to SO4(n), same size used for soa  (cka)

    !-- Conversion of H2SO4 from kg/kg to #/cm3
    !-- and calculation of relative humidity (needed by binary nucleation parameterization)
    do k=1,pver
        do i=1,ncol
            ! air density
            rho_air(i,k) = pmid(i,k) / ( t(i,k) * rair )
            !avogad*1.e-3_r8 to get molec/mol instead of molec/kmol for h2so4 gas
            h2so4_num_cm3(i,k) = ( 1.e-6_r8 * h2so4_mmr(i,k) * avogad * 1.e-3_r8 * rho_air(i,k) ) / ( molmass_h2so4 * 1.E-3_r8 )

            call qsat_water(t(i,k), pmid(i,k), dummy, qs(i,k))

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
                zrhoa = max(h2so4_num_cm3(i,k),1.E+4_r8)
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

                nuclrate_accum_bin(i,k)=MAX(MIN(zjnuc,1.E10_r8),0._r8)
                nuclsize_bin(i,k)=MAX(MIN(zrc,1.E2_r8),0.01_r8)

            end do
        end do
    else   !No atmospheric nucleation
        nuclrate_accum_bin(:,:)=0._r8
        nuclsize_bin(:,:)=1._r8
    end if

    !-- Boundary layer nucleation
    do k=1,pver
        do i=1,ncol

            !-- Nucleation rate #/cm3/s
            if(pblht_lim(i)>zm(i,k) .AND. pbl_nucleation>0) then

                if(pbl_nucleation .EQ. 1) then

                    !-- Paasonen et al. (2010), eqn 10, Table 4
                    nuclrate_accum_pbl(i,k)=(1.7E-6_r8)*h2so4_num_cm3(i,k)
                ! TODO: introduce other schemes when organics are available
                end if

               nuclrate_accum_pbl(i,k)=MAX(MIN(nuclrate_accum_pbl(i,k),1.E10_r8),0._r8)

            else !Not using PBL-nucleation
               nuclrate_accum_pbl(i,k)=0._r8
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
            growth_rate_h2so4(i,k)=(3.E-9_r8/h2so4_dens)*(vmolh2so4*molmass_h2so4*h2so4_num_cm3(i,k))
            growth_rate_h2so4(i,k)=MAX(MIN(growth_rate_h2so4(i,k),10000._r8),1.E-10_r8)

            ! TODO: Molecular speed and growth rate: ORG. Eq. 21 in Kerminen and Kulmala 2002

            ! Combined growth rate (cka) -> add growth_rate_organic, when available
            growth_rate_tot(i,k)=growth_rate_h2so4(i,k)

            !-- Lehtinen 2007 parameterization for apparent formation rate
            !   diameters in nm, growth rate in nm h-1, coagulation in s-1
            ! get formation rate: rate of formation of particles with smallest model bin radius #/cm3/s
            call apparent_formation_rate(nuclsize_bin(i,k), target_particle_radius*2.E9_r8, nuclrate_accum_bin(i,k), formrate_accum_bin(i,k), coag_sink_nucl(i,k), growth_rate_tot(i,k))
            call apparent_formation_rate(nuclsize_pbl(i,k), target_particle_radius*2.E9_r8, nuclrate_accum_pbl(i,k), formrate_accum_pbl(i,k), coag_sink_nucl(i,k), growth_rate_tot(i,k))

            formrate_accum_bin(i,k)=MAX(MIN(formrate_accum_bin(i,k),1.E3_r8),0._r8)
            formrate_accum_pbl(i,k)=MAX(MIN(formrate_accum_pbl(i,k),1.E3_r8),0._r8)

            !   Number of mol nucleated per g air per second.
            nucleated_volume(i,k) = (formrate_accum_bin(i,k) + formrate_accum_pbl(i,k)) & ! [particles/cm3/s]
                            *1.0e6_r8                                 & !==> [particles / m3 /]
                            *2._r8*target_particle_radius**3*pi/6._r8                 & !==> [m3_{aer} / m3_{air} / sec]
                            / rho_air(i,k)                               !==> m3_{aer} / kg_{air} /sec


            !Estimate how much is organic based on growth_rate_totowth-rate
            ! TODO: implement organics here -> fraction of growth due to H2SO4
            frach2so4=1._r8

            ! Nucleated so4 and soa mass mixing ratio per second [kg kg-1 s-1]
            ! used density of particle phase, not of condensing gas
            nuclso4_tend(i,k) = aero_props%density(sulfate_specprop_ndx) * nucleated_volume(i,k) * frach2so4
            nuclnum_tend(i,k) = (formrate_accum_bin(i,k) + formrate_accum_pbl(i,k)) * 1.e6_r8 ! ! [#/kg/s]
         end do
      end do

      !-- Diagnostic output
      nuclrate_accum(:,:)       = nuclrate_accum(:,:)       + (nuclrate_accum_pbl(:,:) + nuclrate_accum_bin(:,:)) * dt
      nuclrate_accum_pbl_o(:,:) = nuclrate_accum_pbl_o(:,:) + nuclrate_accum_pbl(:,:)                             * dt
      formrate_accum(:,:)       = formrate_accum(:,:)       + (formrate_accum_pbl(:,:) + formrate_accum_bin(:,:)) * dt
      formrate_accum_pbl_o(:,:) = formrate_accum_pbl_o(:,:) + formrate_accum_pbl(:,:)                             * dt
      growth_rate_h2so4_o(:,:)  = growth_rate_h2so4_o(:,:)  + growth_rate_h2so4(:,:)                              * dt
      h2so4nucl_o(:,:)          = h2so4nucl_o(:,:)          + h2so4_mmr(:,:)                                      * dt
! TODO: some nice outflds here
   end subroutine aerosol_nucleation

!==========================================================================================
! Helper routine to calculate the apparent formation rate of particles of size dx
! from nucleation rate of particles of size d1
!==========================================================================================

   subroutine apparent_formation_rate(d1, dx, j1, jx, CoagS_dx, growth_rate_tot)
      !-- apparent_formation_rate calculates the formation rate jx of dx sized particles from the nucleation rate j1 (d1 sized particles)
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
      real(r8), intent(in)  :: growth_rate_tot     ! Particle growth rate (nm h-1)
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
      m = -1.6_r8
      CoagS_d1 = CoagS_dx * (d1/dx)**m
      CoagS_d1 = MAX( MIN( CoagS_d1, 1.E2_r8 ), 1.E-10_r8 )

      gamma = ( 1._r8 / (m+1._r8) ) * ( (dx/d1)**(m+1._r8) - 1._r8 )
      gamma = MAX( MIN( gamma, 1.E2_r8 ), 1.E-10_r8 )

      !-- (Eq. 7) CoagS_d1 is multiplied with 3600 to get units h-1
      jx = j1 * exp(-gamma*d1*CoagS_d1*3600._r8 / growth_rate_tot)

   end subroutine apparent_formation_rate

end module condtend