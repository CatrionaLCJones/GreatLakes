k.dieke.mel = function(wnd.z, Kd, lat, lake.area, atm.press, datetime, z.aml, airt, wnd, rh, sw, lwnet,
                       wvht,fetch, Ts){
  
  # define constants used in function
  Kelvin <- 273.15 # temp mod for deg K
  emiss <- 0.972 # emissivity;
  S_B <- 5.67E-8 # Stefan-Boltzman constant (?K is used)
  vonK <- 0.41 # von Karman  constant
  dT <- 0.5   # change in temp for mixed layer depth
  C1 <- 114.278 # from Soloviev et al. 2007
  nu <- 0.29 # proportionality constant from Zappa et al. 2007, lower bounds
  KeCrit <- 0.18     # constant for wave age = 20 (Soloviev et al. 2007)
  albedo_SW <- 0.07
  swRat <- 0.46 # percentage of SW radiation that penetrates the water column
  g <- 9.81 # gravity
  C_w <- 4186 # J kg-1 ?C-1 (Lenters et al. 2005)
  mnWnd <- 0.2 # minimum wind speed
  
  # impose limit on wind speed
  rpcI <- wnd < mnWnd
  wnd[rpcI] <- mnWnd
  
  # calculate sensible and latent heat fluxes
 mm <- calc.zeng(datetime,Ts,airt,wnd,rh,atm.press,wnd.z)
  
  C_D <- mm$C_D # drag coefficient for momentum
  E <- mm$alh # latent heat flux
  H <- mm$ash # sensible heat flux
  

  # calculate total heat flux
  dUdt <- sw*0.93 - E - H + lwnet
  Qo <- sw*(1-albedo_SW)*swRat #PAR
  
  # calculate the effective heat flux
  q1 <- 2-2*exp(z.aml*-Kd)
  q2 <- z.aml*Kd
  q3 <- exp(z.aml*-Kd)
  H_star <- dUdt-Qo*(q1/q2-q3) #Effective surface heat flux Kim 1976
  
  # calculate water density
  rho_w <- water.density(Ts)
  
  # calculate u*
  if (wnd.z != 10) {
    e1 <- sqrt(C_D)
    wnd <- wnd/(1-e1/vonK*log(10/wnd.z))
  }
  rhoAir <- 1.2 #  air density
  tau <- C_D*wnd^2*rhoAir
  uSt <- sqrt(tau/rho_w)
  uSta <- sqrt(tau/rhoAir)  #friction velocity in air
  
  # calculate the thermal expansion coefficient
  thermalExpFromTemp <- function(Ts) {
    V <- 1
    dT <- 0.001
    T1 <- water.density(Ts)
    T2 <- water.density(Ts+dT)
    V2 <- T1/T2
    dv_dT <- (V2-V)/dT
    alpha <- dv_dT
    return (alpha)
  }
  tExp <- thermalExpFromTemp(Ts)
  
  # calculate buoyancy flux and w*
  B1 <- H_star*tExp*g #Hstar * coefficient of thermal expansion * gravity
  B2 <- rho_w*C_w
  Bflx <- B1/B2
  Bflx[Bflx>0] = 0
  
  wSt <- (-Bflx*z.aml)^1/3
  
  # calculate kinematic viscosiy
  getKinematicVis <- function(Ts) {
    # from Mays 2005, Water Resources Engineering
    tempTable <- seq(0,100,by=5)
    # table in m2/s E-6
    visTable <- c(1.792,1.519,1.308,1.141,1.007,0.897,
                  0.804,0.727,0.661,0.605,0.556,0.513,0.477,0.444,
                  0.415,0.39,0.367,0.347,0.328,0.311,0.296)
    v <- approx(tempTable,visTable,xout = Ts)$y
    v <- v*1e-6
    return(v)
  }
  kinV <- getKinematicVis(Ts)
  kinVa <- getKinematicVis(airt)
  
  KeDe <- (kinV*g)
  KeNm <- uSt^3
  Ke <- KeNm/KeDe
  tau <- tau    # estimate of total tau (includes wave stress)
  euPw <- (1+Ke/KeCrit)  # tau_shear = tau/(1+Ke/Kecr) Ke is the Keulegan number
  # Could calculate KeCrit (critical Keulegan number) from wave age
  #KeCrit <- (kinVa/kinV)*((rhoAir/rho_w)^1.5)*(Rbcr/Aw) # Eq1.16-Soloviev et al(2007)
  
  tau_t <- tau/euPw      # tau_t = tangential wind stress, tau = total wind stress
  uTanS <- tau_t/rho_w
  uTanS <- uTanS^0.5
  
  # calculate viscous sublayer
  Sv <- C1*kinV/uTanS  # effective thickness of the aqueous viscous sublayer
  eu_N <- uTanS^3      # e_u(0) = (tau_t/rho)^1.5/(vonK*Sv)
  eu_D <- vonK*Sv      # denominator
  eu_0 <- eu_N/eu_D    # in m2/s3
  ec_0 <- -1.0*Bflx       # buoyancy flux, but only when outward
  
  #ewave_0 turbulence due to wave breaking
 # Fetch <- Calculated.Fetch(latitude, longitude, wnd.dir) # fetch as function of wind direction
    Hs1 <- wvht #use actual significant wave height data from buoys. If significant wave height data 
    #unavailable, use Calculated.Fetch function followed by:
    Hs2 <- 0.0163*(fetch^0.5)*wnd
    
  Aw <- (1/(2*pi))*(( (g*Hs1*rhoAir)/(0.062*rho_w*
                                       uSt^2))^(2/3)) # wave age - eqn 1.11 Soloviev et al. (2007)
  
  W <- 3.8e-6*wnd^3.4 # simplified whitecap fraction (Fariall et al. 2000)
  
  
  Ap <- 2.45*W*((1/(W^0.25))-1)
  alphaW <- 100 # p. 185 - Soloviev et al. (2007)
  B <- 16.6 # p. 185 - Soloviev et al. (2007)
  Sq <- 0.2 # p. 185 - Soloviev et al. (2007)
  cT <- 0.6 # p. 188 - Soloviev et al. (2007)
  ewave_0 <- ((Ap^4)*alphaW)*((3/(B*Sq))^0.5) *
    (((Ke/KeCrit)^1.5)/((1+Ke/KeCrit)^1.5)) *
    (uSt*g*kinV)/(0.062*vonK*cT*((2*pi*Aw)^1.5)) *
    (rhoAir/rho_w)
  
  #------------------------------------
  e_0 <- ec_0+eu_0+ewave_0    # e(0) from Soloviev (w/o wave effects)
  Kc <- ec_0*kinV*100^4*3600^4      # convective component now in cm4/hr4  (Total)
  Ku <- eu_0*kinV*100^4*3600^4 # shear component now in cm4/hr4  (Total)
  Kwave <- ewave_0*kinV*100^4*3600^4 # wave component now in cm4/hr4  (Total)
  Kall <- e_0*kinV*100^4*3600^4       # turbulent kinetic energy now in cm4/hr4  (Total)
  
  #Schmidt number could be calculated as temperature dependent
  #Sc <- 1568+(-86.04*Ts)+(2.142*Ts^2)+(-0.0216*Ts^3)
  k600org <- nu*600^(-0.5)*(Kc+Ku)^(1/4)   # in cm/hr (Total)
  k600org <- k600org*24/100 #now in units in m d-1
  
  k600 <- nu*600^(-0.5)*Kall^(1/4)   # in cm/hr (Total)
  k600 <- k600*24/100 #now in units in m d-1
  
  # ---Breaking Wave Component, Author: Dieke and Melville
  
  # AB = dimensional fitting coefficient. Per DM18 = 10^-5 m2/s2),
  AB <- 10^(-5)
  # R = ideal gas constant - j/Kmol (check!). 
  R <- 8.314
  # a = O2 gas solubility coefficient in freshwater = ~1.22 × 10−3 mol dm−3
  # at 25C and 1 atm of pressure
  a <- 1.22 * (10^(-3))
  # Ts = surface temperature of water, from data
  Ts <- Ts
  # Os = Ostwald gas soluability as a function of R, a, and Ts
  Os <- R*a*Ts
  ABOs <- AB/Os
  # u*atm = wind friction - from Bye et al 2014: u* = sqrt(K10)*U10 where
  # K10 = 10m drag coefficient (see Edson et al 2011)
  U10 <- wnd.z #from data
  K10 <- if(U10 > 8) {0.035
  } else {
    0.03
  }
  # Split u* by wind speed - if greater than 8m/s use 0.035, if less than 8m/s use 0.03
  Ustar <- K10*U10
  U <- Ustar^(5/3)
  # g = gravitational acceleration in m/s2
  g <- 9.80665
  # Hs = significant wave height, from data
  Hs <- wvht
  gHs <- (g*Hs1)^(2/3)

  # kB = full bubble coefficient
  kb <- ABOs*U*gHs
  
  kb <- kb*24/100 #now in units in m d-1
  #----------------------------------------------------------------
  
  # Sc = Schmidt's number for oxygen in freshwater (from Wanninkof)
  Sc <- 530
  SC <- (Sc/600)^(-(1/2))
  
  k600b = (k600+kb)*SC
  allks = data.frame(Ku,Kc,Kwave,kb,k600org,k600,k600b)
  colnames(allks) = c("shear","convective","wave","bubble",
                      "k600org","k600",'k600b')
  return(as.numeric(k600b))
}