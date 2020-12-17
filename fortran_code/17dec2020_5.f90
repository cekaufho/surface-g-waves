!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! December 17, 2020 - Test 5
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
  
  enable_free_surface = .true.
  enable_conserve_energy = .true.
  coord_degree           = .false.
  enable_cyclic_x        = .true.
  enable_hydrostatic     = .true.  
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
  
  !enable_hor_diffusion = .true
  !kappah_0 = mix/fac**2
  !k_h = mix/fac**2

  runlen = 2000  !86400
  enable_diag_ts_monitor = .true.; ts_monint =dt_tracer
  enable_diag_snapshots  = .true.; snapint  = 1
  
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
 H0 = dzt(1)*nz     ! Height of the tank
 B0 = dyt(js_pe)*ny ! Width of the tank
 
end subroutine set_grid

!=========================================================================================================
! SUBROUTINE CORIOLIS
!=========================================================================================================

subroutine set_coriolis
use main_module   
implicit none

!---------------------------------------------------------------------------------------------------------
! Beta Plane

! real*8 :: phi0 , betaloc
! phi0 = 10.0 /180. *pi
! betaloc = 2*omega*cos(phi0)/radius
! do j=js_pe-onx,je_pe+onx
!    coriolis_t(:,j) = 2*omega*sin(phi0) +betaloc*yt(j)
! enddo
!tank omega = 2 pi * rounds/60sec 
!coriolis_t = 0.1047 !1rpm
!coriolis_t = 0.2094 !2rpm
!coriolis_t = 0.8378 !8rpm
!coriolis_t =  0.4189 !4rpm
coriolis_t = 0

end subroutine set_coriolis

!=========================================================================================================
! SUBROUTINE INITIAL CONDITIONS
!=========================================================================================================

subroutine set_initial_conditions
 use main_module   
 use config_module   
 implicit none
 integer :: i,j,k,slope,ic_0(2),jc_0(2), gg
 real*8 :: temp_init,shelf
 allocate( T0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); T0 = 0.0
 allocate( u0(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz) ); u0 = 0
 
 T0 = 12 
 temp(:,:,:,taum1) = T0
 temp(:,:,:,tau)   = T0
surf_press(:,:,1) = 0.2 + 0.01

!  temp = 20.
!  salt = 0.
!  do i=is_pe-onx,ie_pe+onx
!     if (xt(i) <=25./80.*L0) salt(i,:,:,:)=2.
!  enddo
! u0 = ...
! u(:,:,:,tau) = 
! u(:,:,:,taum1) = 

end subroutine set_initial_conditions

!=========================================================================================================
! SUBROUTINE FORCING
!=========================================================================================================

subroutine set_forcing
 use main_module   
 use config_module   
 implicit none
 !real*8 :: T_rest, t_star
 integer :: i,j,k
 
 !t_star = 5 !Temperature to restore to
 !T_rest = 15.*dt_mom ! Restoring time scale
 !temp_source(1:5,:,nz) = maskT(1:5,:,nz) *1./T_rest*(t_star-temp(1:5,:,nz,tau))
 
!allocate ( u0 ( is_pe - onx : ie_pe + onx , js_pe - onx : je_pe + onx , nz ) ); u0 = 0
 
 !real*8 :: rho0=1024 , tau0 = -0.05, y
 !surface_flux
 !surface_tau
 !temp_source
 u_source=u0 * sin (2* 3.14 * OM0 * itt * dt_tracer ) ! changed pi to 3.14
 !     do j=1,ny
 !      do i=1,nx
 !       y = (yu(j)-yu(1))/ (yu(ny)-yu(1))
 !       surface_taux(i,j) = tau0/rho_0 *cos(2*pi*y)*maskU(i,j,nz-1)
 !       surface_tauy(i,j) = 0.0*maskV(i,j,nz-1)
 !      enddo
 !     enddo
 !if (itt*dt_mom < 5) then
 !  u_source = 0.005
 !else
 !  u_source = 0.
 !endif
 
 
  ! implement effect of background state

 ! update density, etc of last time step
! temp(:,:,:,tau) = temp(:,:,:,tau) + t0
! call calc_eq_of_state(tau)
! temp(:,:,:,tau) = temp(:,:,:,tau) - t0

 ! advection of background temperature
! call advect_tracer(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,nz,t0,dt0(:,:,:,tau))
! temp_source = (1.5+ab_eps)*dt0(:,:,:,tau) - ( 0.5+ab_eps)*dt0(:,:,:,taum1)

 ! wavemaker
! u_source= u0*sin(2*pi*OM0*itt*dt_tracer)
 

end subroutine set_forcing

!=========================================================================================================
! SUBROUTINE TOPOGRAPHY
!=========================================================================================================

subroutine set_topography

 ! Define variables in advance
 use main_module   
 use config_module   
 implicit none
 integer :: i, k, j, slope
 real*8 :: fxa, shelf, a, b

! kbot(i,j): e.g. kbot=0 means land 
! kbot=1 everywhere means no topographic features

kbot =1

!---------------------------------------------------------------------------------------------------------
! Double slits

do i = is_pe, ie_pe
    do j = js_pe, je_pe
        if (i == 30 .AND. j .LE. 18 ) kbot(i,j) =0
        if (i== 30 .AND. j .GE. 22 .AND. j .LE. 40) kbot(i,j) =0
        if (i ==30 .AND. j .GE. 44) kbot(i,j)=0
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

! December 17, 2020 - Test 5