!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! Feb 28 - 06 - Final Tank Scale
!=========================================================================================================

module config_module
 use main_module
 real*8 :: p=1, mix=1e-06
 real*8, allocatable :: T0(:,:,:), u0(:,:,:)
 real *8 :: N_0 = 2.* pi/10.
 real *8 :: OM0 = 1./(1.5*10.), x0 !Note: Defined here x0, y0 for wavemaker location
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
  
  nx=int(100)
  ny=int(100)
  nz=int(5)

  dt_tracer = 0.01
  dt_mom    = 0.01
  
  enable_free_surface = .true.
  enable_conserve_energy = .true.
  coord_degree           = .false.
  enable_cyclic_x        = .false.
  enable_cyclic_y		 = .true.
  enable_hydrostatic     = .true.
  eq_of_state_type       =  1

  congr_epsilon = 1e-8
  congr_max_iterations = 5000
  congr_epsilon_non_hydro=   1e-8
  congr_max_itts_non_hydro = 5000

  enable_bottom_friction = .true.
  r_bot = 0.3

  ! Vertical friction
  enable_explicit_vert_friction = .true.
  kappam_0 = mix/p**2
  
  ! Horizontal friction
  enable_hor_friction = .true.
  a_h = mix/p**2
  
  enable_superbee_advection = .true.
  enable_tempsalt_sources   =  .true.
  enable_momentum_sources   =  .true.
  
  enable_hor_diffusion = .true.
  kappah_0 = mix/p**2
  k_h = mix/p**2

  runlen = 2
  enable_diag_ts_monitor = .true.; ts_monint = dt_tracer
  enable_diag_snapshots  = .true.; snapint  = dt_tracer
  
end subroutine set_parameter

!=========================================================================================================
! SUBROUTINE SET GRID
!=========================================================================================================

subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 
 dxt(:) = 0.006
 dyt(:) = 0.006
 dzt(:) = 0.004
 
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
    u0(i, j, 2) = maskU(i, j, 2) * 1./(500*60.0 * dt_tracer) * exp( -(xu(i)-x0)**2/(dxt(is_pe)*1)**2 )
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
 u_source = u0 * sin (100 * pi * OM0 * itt *dt_tracer)

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

!do i = is_pe, ie_pe
!    do j = js_pe, je_pe
!        if (i == 78 .AND. j .LE. 29 ) kbot(i, j) = 0
!        if (i== 78 .AND. j .GE. 29+6 .AND. j .LE. 29+6+30) kbot(i, j) = 0
!        if (i == 78 .AND. j .GE. 29+6+30+6) kbot(i, j) = 0
!    enddo
!enddo

!---------------------------------------------------------------------------------------------------------
! Single slits

do i = is_pe, ie_pe
    do j = js_pe, je_pe
        if (i == 67 .AND. j .LE. 30 ) kbot(i,j) = 0
                if (i == 67 .AND. j .GE. 30+40 ) kbot(i,j) = 0
    enddo
enddo

end subroutine set_topography

!=========================================================================================================
! END
!=========================================================================================================

subroutine set_diagnostics
end subroutine set_diagnostics

subroutine set_particles
end subroutine set_particles
