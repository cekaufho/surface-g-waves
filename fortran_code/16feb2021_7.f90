!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! Feb 16, 2021 - Test 7
!=========================================================================================================

module config_module
 use main_module
 real*8 :: fac=2., mix=1e-6, L0, H0, B0
 real*8, allocatable :: T0(:,:,:), u0(:,:,:)
 real *8 :: N_0 = 2* pi /10.
 real *8 :: OM0 = 1./(1.5*10), x0, y0
end module config_module

!=========================================================================================================
! SUBROUTINE SET PARAMETERS
!=========================================================================================================

subroutine set_parameter
 use main_module   
 use config_module   
 use diagnostics_module  
 use tke_module   
 implicit none
  
  nx=int(60*fac)
  nz= 2 !int(2*fac)
  ny=int(60*fac)

  dt_tracer = 20*0.025/fac
  dt_mom    = 20*0.025/fac
  
  enable_free_surface    = .true.
  enable_conserve_energy = .true.
  coord_degree           = .false.
  enable_cyclic_x        = .false.
  enable_hydrostatic     = .true.
  enable_cyclic_y        = .true.
  eq_of_state_type       =  1

  congr_epsilon = 1e-12
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-9
  congr_max_itts_non_hydro = 5000

  enable_bottom_friction = .false.
  r_bot = 0.3

  ! Vertical friction
  enable_explicit_vert_friction = .true.
  kappam_0 = mix/fac**2
  
  ! Horizontal friction
  enable_hor_friction = .true.
  a_h = mix/fac**2
  
  enable_superbee_advection = .true.
  enable_tempsalt_sources =  .true.
  enable_momentum_sources =  .true.
  
  enable_hor_diffusion = .true.
  kappah_0 = mix/fac**2
  k_h = mix/fac**2

  runlen = 150 
  enable_diag_ts_monitor = .true.; ts_monint = 0.5 !dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 5.0 !1
  
end subroutine set_parameter

!=========================================================================================================
! SUBROUTINE SET GRID
!=========================================================================================================

subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 
 dxt(:)=0.25/fac
 dyt(:)=0.25/fac
 dzt(:)=0.25/fac
 
 L0 = dxt(is_pe)*nx  ! Length of the tank
 H0 = dzt(2)*nz      ! Height of the tank 
 B0 = dyt(js_pe)*ny  ! Width of the tank
 
end subroutine set_grid

!=========================================================================================================
! SUBROUTINE CORIOLIS
!=========================================================================================================

subroutine set_coriolis
use main_module   
implicit none

coriolis_t = 0

end subroutine set_coriolis

!=========================================================================================================
! SUBROUTINE INITIAL CONDITIONS
!=========================================================================================================

subroutine set_initial_conditions

 use main_module   
 use config_module   
 implicit none
 integer :: i,j,k
 allocate( T0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); T0 = 0.0
 allocate( u0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u0 = 0.0
 
 T0 = 12. 
 temp(:, :, :, taum1) = T0
 temp(:, :, :, tau)   = T0
 
  x0 = xu(1)
  
  do j = js_pe, je_pe
   do i = is_pe, ie_pe
    u0(i, j, 2) = maskU(i, j, 2) * 1./(100 * 60. * dt_tracer) * exp(-(xu(i)-x0)**2/(dxt(is_pe)*1)**2)
   enddo
  enddo

end subroutine set_initial_conditions

!=========================================================================================================
! SUBROUTINE FORCING
!=========================================================================================================

subroutine set_forcing

 use main_module   
 use config_module   
 implicit none
 integer :: i,j,k
 
 ! Wavemaker 
 u_source = u0 * sin( 2 * pi * OM0 * itt * dt_tracer )

end subroutine set_forcing

!=========================================================================================================
! SUBROUTINE TOPOGRAPHY
!=========================================================================================================

subroutine set_topography

 ! Define variables in advance
 use main_module   
 use config_module   
 implicit none
 integer :: i, j

 kbot = 1

!---------------------------------------------------------------------------------------------------------
! Double slits

do i = is_pe, ie_pe
    do j = js_pe, je_pe
        if (i == 40 .AND. j .LE. 18 ) kbot(i,j) =0
        if (i== 40 .AND. j .GE. 21 .AND. j .LE. 39) kbot(i,j) =0
        if (i ==40 .AND. j .GE. 42) kbot(i,j)=0
    enddo
enddo

!---------------------------------------------------------------------------------------------------------
! Single slits

!do i = is_pe, ie_pe
!    do j = js_pe, je_pe
!        if (i == 40 .AND. j .LE. 27 ) kbot(i,j) =0
!                if (i == 40 .AND. j .GE. 33 ) kbot(i,j) =0
!    enddo
!enddo

end subroutine set_topography

!=========================================================================================================
! END
!=========================================================================================================

subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles