module constants_module
    
    use kind_module, only: wp8 => SHR_KIND_R8, wp4 => SHR_KIND_R4

    implicit none
    save
    public
    
    real(wp4), parameter:: prestompa=1.0e-6
    
    real(wp4), parameter:: Pi = 3.1415926,    & !Pi-constant
                       pip180 = Pi/180.0
    
    real(wp8), parameter:: dPi = 3.14159265358979d0,   &  !Pi-constant(double precision)
                       dpip180 = dPi/180.0d0
    
    real(wp8), parameter:: lat_extr= 89.99999d0
    
    real(wp4), parameter:: RadEarth = 6371000.0,    & !Earth radius[m]
                        EarthAngVel = 7.2921159e-5, & !earth angular velocity[rad/sec]
                       HeatCapWater = 4000.0,       & !heat capacity of water[J/kg/�C)]for 35%.
                             RefDen = 1025.0,       & !reference density[kg/m**3]
                        FreeFallAcc = 9.8,          & !free fall acceleration[m/s**2]
                           DenFresh = 1000.0,       & !Fresh water density[kg/m^3]
                            rgas    = 287.04,       & !J/kg/K Gas constant
                            cp_air  = 1000.5,       & !J/kg/K Air heat capacity
                            lambda_f= 3.337e5,      & !J/kg   Heat of snow/ice fusion
                            lambda_v= 2.5e6,        & !J/kg Heat of evaporation
                            lambda_s= 2.839e6         !J/kg Heat of sublimation
    
    real(wp4), parameter:: q1_o=0.98*640380.0,       &  !kg/m^3   ! saturated vapour parameters
                           q2_o=-5107.4                 !K        ! over ocean
    
    real(wp4), parameter:: q1_i=11637800.0,          &  !kg/m^3   ! saturated 
                           q2_i=-5897.8                 !K        ! over ice/snow
    
    ! Albedo
    real(wp4), parameter::   AlbOpWater=0.066        ! open water albedo
    
    real(wp4), parameter:: StephBoltz = 5.670367e-8  ! stephan-boltzman (W/m**2/K^4)
    
    real(wp4), parameter:: EmissWater=1.0            !water emissitivity
    real(wp4), parameter:: Emissice = 0.95           !water emissitivity
    
    real(wp4), parameter:: vonkarman=0.4
    
    !turbulent model parameters
    
    real(wp4), parameter:: dimdepth=3.0    !upper layer depth for setting high diffusivity
    
    real(wp4), parameter:: prandtl_lat=0.1
    
    real(wp4), parameter:: epsrho=1.0e-3,angle_max=0.01
    
    real(wp4), parameter ::   A1_t=0.92,   &
                              A2_t=0.74,   &
                              B1_t=16.6,   &
                              B2_t=10.1,   &
                              C1_t=0.08,   &
                              E1_t=1.80,   &
                              E2_t=1.33
    
    real(wp4), parameter:: tur_factor_nu=0.41,      &
                           tur_var_min=1.0e-8
    
    real(wp4), parameter:: Pice_cr = 2.75e+4,     &  !Ice pressure parameter(Pa)
                             Cstar = 20.0            !Ice pressure exponent scale
    
    real(wp4), parameter:: def_rate_min=1.0e-11,  &  !Minimum of ice deformation rate (1/s)
                           ice_mass_damp=615.0       !Typical ice mass for damping
    
    real(wp4), parameter:: den_ice  = 917.0,       & !ice density  (kg/m^3)
                           den_snow = 330.0,       & !snow density (kg/m^3)
                           Stenton  = 0.006,       & !Stenton number
                           z0star   = 0.005,       & !drag scale parameter (m)
                           h0       = 3.0,         & !Ice height parameter for drag (m) 
                           deltaz   = 1.0,         & !Wall proximity scale in ice drag (m)
                           extr2    = 4.0            !square of rheology excentricity
    
    real(wp4), parameter:: ice_mass_min=100.0        !Minimum of ice mass (kg/m^2)
    real(wp4), parameter:: tmelt=0.0                 !Ice/Snow melting temperature
    real(wp4), parameter:: kice=2.03, ksnow=0.3      !Heat conductivity of fresh ice and snow (W/m/�C)
    real(wp4), parameter:: accur_flux=0.01           !Accuracy for ice heat balance (W/m^2)
    real(wp4), parameter:: accur_tem =0.00001        !Accuracy for ice/snow temperature (K)
    
    real(wp4), parameter:: aimin=0.001,          &   !Minimum of ice compactness(0-1)
                           aimax=0.999,          &   !Maximum of ice compactness(0-1)
                           himin=0.001,          &   !Minimum of ice  height(m)
                           hsmin=0.001               !Minimum of snow height(m)
    
    real(wp4), parameter:: aimin_dyn=0.1
    
    real(wp4), parameter:: hsmax=0.3,            &   !maximum show height (if more, snow becomes ice)
                           gamma_stoice=1.0e-7       !time scale of snow to ice gravity transformation
    
    real(wp4), parameter:: cd_n10_i    = 1.63e-3,           &
                           cd_n10_i_rt = sqrt(cd_n10_i),    &
                           ce_n10_i    = 1.63e-3,           &
                           ch_n10_i    = 1.63e-3    
    
    real(wp4), parameter:: h_swi = 1.4,       &   !depth of solar radiation penetration in ice (m),
                           i0=0.2                 !part of solar radiation penetrated in ice
    
    real(wp4), parameter:: hmax_sugar=50.0,     &  !Maximum ocean depth to form a new ice at open water
                           href_min = 0.01,     &  !new ice height low  limit (m)
                           href_max = 0.5,      &  !new ice height high limit (m)
                           tref_min = -20.0,    &  !air temperature low limit (C)
                           tref_max = -10.0,    &  !air temperature high limit (C)
                           cmelt= 0.5              !lateral melting coefficient (0-1)
    
    real(wp4), parameter:: cda_ref=2.75e-3,      &
                           cdw_ref=5.50e-3,      &
                           airden_ref = 1.3   
    
    real(wp4), parameter:: freeze_frac_sal=0.25    !fraction of sea ice salinity relative to water formed it
    real(wp4), parameter:: rpart = 0.58, efold1 = 0.35, efold2 = 23.0   !SW-radiation penetration parameters
    
    !bottom friction parameters
    
    !  Type_fric - flag:
    integer, parameter:: type_fric = 2  ! 0 - no bottom fric
                                        ! 1 - linear
                                        ! 2 - nonlinear
    
    !Cb_l  - coef of    linear bottom friction (typical value 5e-4 m/s)
    !Cb_nl - coef of nonlinear bottom friction (typical value 1.0e-3 )
    !Ebottom - bottom turbulent kinetic Energy (typical value 0: (5d-2 m/s)**2)
    real(wp4), parameter::    Cb_l    = 5e-4,        & 
                              Cb_nl   = 2.5e-3,      &   !Popov S.K.
    !                         Cb_nl   = 1.0e-3,      &   !FRAM
                              Ebottom = 25.0e-4          
    
    
end module constants_module