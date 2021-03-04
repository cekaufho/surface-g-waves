!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! December 18, 2020 - Test 2
!=========================================================================================================

module config_module
 use main_module
 real*8 :: fac=1.0,mix=2e-06,L0,H0,B0
 real*8, allocatable :: T0(:,:,:),u0(:,:,:),usource(:,:,:)
 real *8 :: N_0 = 2* pi /10.
 real *8 :: OM0 = 1./(1.5*10)
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
 nz=int(20*fac)
 ny=int(60*fac)

 dt_tracer=0.25/fac
 dt_mom   =0.25/fac
 runlen = 2000  !86400
 enable_diag_ts_monitor = .true.; ts_monint =dt_tracer
 enable_diag_snapshots  = .true.; snapint  = 1
    
 enable_free_surface = .true.
 enable_conserve_energy = .false.
 coord_degree           = .false.
 enable_cyclic_x        = .true.
 enable_hydrostatic     = .false.  
 
 eq_of_state_type       =  1
 congr_epsilon = 1e-8
 congr_max_iterations = 5000 ! Changed from 5000
 congr_epsilon_non_hydro=   1e-8
 congr_max_itts_non_hydro = 5000 ! Changed from 5000
 enable_bottom_friction = .false.
 r_bot = 0.7

 ! Vertical friction
 enable_explicit_vert_friction = .true.
 kappam_0 = mix/fac**2
  
 ! Horizontal friction
 enable_hor_friction = .true.
 a_h = mix/fac**2
 enable_superbee_advection = .true.
 enable_tempsalt_sources =  .true. ! restoring zone
 enable_momentum_sources =  .true. ! restoring zone
  
end subroutine set_parameter

!=========================================================================================================
! SUBROUTINE SET GRID
!=========================================================================================================

subroutine set_grid
 use main_module   
 use config_module   
 implicit none
 
 dxt(:)=0.01/fac
 dyt(:)=0.01/fac
 dzt(:)=0.01/fac
 
 L0 = dxt(is_pe)*nx ! Length of the tank
 H0 = dzt(20)*nz     ! Height of the tank
 B0 = dyt(js_pe)*ny ! Width of the tank
 
end subroutine set_grid

!=========================================================================================================
! SUBROUTINE CORIOLIS
!=========================================================================================================

subroutine set_coriolis
 use main_module   
 implicit none

 coriolis_t = 0.
 coriolis_h = 0.

end subroutine set_coriolis

!=========================================================================================================
! SUBROUTINE INITIAL CONDITIONS
!=========================================================================================================

subroutine set_initial_conditions
 use main_module   
 implicit none
 integer :: i,j,k
 
 do k = 1, nz
  do j = js_pe-onx, je_pe+onx
   do i= is_pe-onx, ie_pe+onx
    u(i,j,k,:) = sin (2* 3.14 * itt * dt_tracer )
   enddo
  enddo
 enddo

end subroutine set_initial_conditions

!=========================================================================================================
! SUBROUTINE FORCING
!=========================================================================================================

subroutine set_forcing
 use main_module   
 implicit none
 integer :: i,j,k

 !do j= 1, ny
 ! do i= 1, nx
 !  usource(i,j,:) = maskU(i,j,:)*(5 * u(i,j,:,tau))
 ! enddo
 !enddo
 
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
   if (i == 30 .AND. j .LE. 18 ) kbot(i,j) = 0
   if (i== 30 .AND. j .GE. 22 .AND. j .LE. 40) kbot(i,j) = 0
   if (i ==30 .AND. j .GE. 44) kbot(i,j)= 0
 enddo
 enddo

!---------------------------------------------------------------------------------------------------------
! Single slits

!do i = is_pe, ie_pe
!    do j = js_pe, je_pe
!        if (i == 30 .AND. j .LE. 28 ) kbot(i,j) =0
!                if (i == 30 .AND. j .GE. 33 ) kbot(i,j) =0
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

! December 18, 2020 - Test 2
