!=========================================================================================================
! SURFACE WAVES EXPERIMENT
! By Johannes and Christine 
! Fortran Code for Calculating Surface Gravity Waves
! Jan 29, 2021 - Test 1
! Changed mixing term to lower
!=========================================================================================================

module config_module
 use main_module
 real*8 :: fac=1.0,p=4,mix=1e-07,L0,H0,B0
 real*8, allocatable :: T0(:,:,:),u0(:,:,:)!,u_source(:,:,:)
! note: u_source is already defined and allocated in the main module, so no
! need to do it here
! whereas u0 is a newly defined variable that needs allocation
 real *8 :: N_0 = 2* pi /10.
 real *8 :: OM0 = 1./(1.5*10),x0,y0 ! note: defined here x0,y0 for wavemaker location
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
  enable_cyclic_x        = .false. !note: this was set to true, meaning that you have open
                                  ! boundaries in x direction, so everthing that goes out
                                  ! on the right side comes back in on the left or vice versa,
                                  ! so set it to false for your closed basin
  enable_hydrostatic     = .true. ! Should enable hydrostatic as horizontal>vertical (increase performance?)  
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
  kappah_0 = mix/fac**2
  k_h = mix/fac**2

  runlen = 100  !86400
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
 
 dxt(:)=0.1/p
 dyt(:)=0.1/p
 dzt(:)=0.1/p
 
 L0 = dxt(is_pe)*nx! Length of the tank
 H0 = dzt(2)*nz    ! Height of the tank 
 B0 = dyt(js_pe)*ny ! Width of the tank
 
end subroutine set_grid

!=========================================================================================================
! SUBROUTINE CORIOLISs
!=========================================================================================================

subroutine set_coriolis
use main_module   
implicit none

! Tank Scale Rotation
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
! note: here u0 is now allocated, meaning memory is made available to the dimensions of the tank
 
 T0 = 12. 
 temp(:,:,:,taum1) = T0
 temp(:,:,:,tau)   = T0
! note: nothing happened because u0=0.., so specify here where you want the wavemaker
! and u0 is the amplitude
  x0=xu(1)
  y0=yu(30)
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    u0(i,j,2)= maskU(i,j,2)*1./(100*60.*dt_tracer)*exp( -(xu(i)-x0)**2/(dxt(is_pe)*1)**2 ) !&
                                                 ! *exp( -(yu(j)-y0)**2/(dyt(1)*1)**2 )
   enddo
  enddo

!surf_press(:,:,1) = 0.2 + 0.01 

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
 integer :: i,j,k
 
 ! Wavemaker 
! allocate(u_source) 
! note: u_source already allocated, and to allocated you need to specify the
! dimensions you want to allocate to, i.e. how much space
 u_source=u0 * sin (2* 3.14 * OM0 * itt * dt_tracer ) ! changed pi to 3.14

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

kbot = 1

!---------------------------------------------------------------------------------------------------------
! Double slits

do i = is_pe, ie_pe
    do j = js_pe, je_pe
        if (i == 40*p .AND. j .LE. 18*p ) kbot(i,j) =0
        if (i== 40*p .AND. j .GE. 21*p .AND. j .LE. 39*p) kbot(i,j) =0
        if (i ==40*p .AND. j .GE. 42*p) kbot(i,j)=0
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

! Jan 29, 2021 - Test 1