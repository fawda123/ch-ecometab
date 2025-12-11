#' Calculate gas transfer velocity for Wanninkhof 2009 equation, which has 0 wind speed intercept coefficient of 3
#'
#' @param Temp numeric for water temperature (°C)
#' @param Sal numeric for salinity (ppt)
#' @param WSpd numeric for wind speed (m/s)
#'
#' @return numeric vector of Kw in m/d
#' @export
#'
#' @details
#' Output is Kw vector that is alternative to calculating KL using Thiebault et al. 2008 (\code{\link{f_calcKL}}).
#' Interpreted as oxygen mass transfer coefficient.
#'
#' @references
#' Wanninkhof 2009. *Advances in quantifying air–sea gas exchange and environmental forcing.* Annu. Rev. Mar. Sci. 1:213–244.
#' Wanninkhof, R. 2014. *Relationship between wind speed and gas exchange over the ocean revisited.*
#' Limnology and Oceanography: Methods, 12(6):351–362. https://doi.org/10.4319/lom.2014.12.351
#'
#' @examples
#' data(SAPDC)
#' f_calcWan09(SAPDC$Temp, SAPDC$Sal, SAPDC$WSpd)

f_calcWan09 <- function(Temp, Sal, WSpd) {
  
  # Calculate from Wanninkhof (2009) equation see figure 6 equation 36
  # windmag is actually k660 (cm hr-1) vs U10 (m s-1)
  # windmag <- 3 + (0.1*u10 + 0.064*u10^2 + 0.011*u10^3)
# Yields Kw to ecometab in units of m/d
  
  windmag <- 3 + (0.1 * WSpd + 0.064 * WSpd^2 + 0.011 * WSpd^3)
  
  # Schmidt number for oxygen
  sc <- oxySchmidt(Temp, Sal)
  a <- 3.0  # DO NOT CHANGE — from Wanninkhof 2009 
  
  # k660 = 0.24 * U^2
  kw <- windmag * ((sc / 660) ^ (-0.5))  # cm/hr
  
  # Convert to m/day
  kw <- kw / (100 * 60 * 60)  # cm/hr → m/s
  kw <- kw * 60 * 60 * 24     # m/s → m/day
  
  return(kw)
}
