!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! Feb 16 - 9
!=========================================================================================================

module config_module
 use main_module
 real*8 :: fac=1.0, p=2, mix=1e-03, L0, H0, B0
 real*8, allocatable :: T0(:,:,:),u0(:,:,:)
 real *8 :: N_0 = 2* pi /10.
 real *8 :: OM0 = 1./(1.5*10),x0,y0 !Note: Defined here x0, y0 for wavemaker location
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
  
  nx=int(60*p)
  nz=int(2*p)
  ny=int(60*p)

  dt_tracer=0.25/p
  dt_mom   =0.25/p
  
  enable_free_surface = .true.
  enable_conserve_energy = .false.
  coord_degree           = .false.
  enable_cyclic_x        = .false.
  enable_cyclic_y		 = .true.
  enable_hydrostatic     = .true.
  eq_of_state_type       =  1

  congr_epsilon = 1e-9
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-9
  congr_max_itts_non_hydro = 5000

  enable_bottom_friction = .true.
  r_bot = 0.2

  ! Vertical friction
  enable_explicit_vert_friction = .true.
  kappam_0 = mix/p**2
  
  ! Horizontal friction
  enable_hor_friction = .true.
  a_h = mix/p**2
  
  enable_superbee_advection = .true.
  enable_tempsalt_sources   =  .true.
  enable_momentum_sources   =  .true.
  
  !enable_hor_diffusion = .true.
  !kappah_0 = mix/p**2
  !k_h = mix/p**2

  runlen = 150
  enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 1
  
end subroutine set_parameter

!=========================================================================================================
! SUBROUTINE SET GRID
!=========================================================================================================

subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 
 dxt(:) = 0.1/p
 dyt(:) = 0.1/p
 dzt(:) = 0.1/p
 
 L0 = dxt(is_pe)*nx  ! Length of the tank
 H0 = dzt(2)*nz      ! Height of the tank 
 B0 = dyt(js_pe)*ny  ! Width of the tank
 
end subroutine set_grid

!=========================================================================================================
! SUBROUTINE CORIOLISs
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
 integer :: i, j
 allocate( T0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); T0 = 0.0
 allocate( u0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u0 = 0
  
 T0 = 12. 
 temp(:,:,:,taum1) = T0
 temp(:,:,:,tau)   = T0

 ! Note: here u0 is now allocated, meaning memory is made available to the dimensions of the tank 
 ! Specify here where you want the wavemaker

  x0 = xu(1)

  do j = js_pe, je_pe
   do i = is_pe, ie_pe
    u0(i, j, 2) = maskU(i, j, 2) * 1./(100 * 60. * dt_tracer) * exp( -(xu(i)-x0)**2/(dxt(is_pe)*1)**2 )
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
 
 !Wavemaker
 u_source=u0 * sin (2 * pi * OM0 * itt * dt_tracer )

end subroutine set_forcing

!=========================================================================================================
! SUBROUTINE TOPOGRAPHY
!=========================================================================================================

subroutine set_topography

 use main_module   
 use config_module   
 implicit none
 integer :: i, j

kbot = 1

!---------------------------------------------------------------------------------------------------------
! Double slits

do i = is_pe, ie_pe
    do j = js_pe, je_pe
        if (i == 40*p .AND. j .LE. 18*p ) kbot(i, j) =0
        if (i== 40*p .AND. j .GE. 21*p .AND. j .LE. 39*p) kbot(i, j) =0
        if (i ==40*p .AND. j .GE. 42*p) kbot(i, j)=0
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