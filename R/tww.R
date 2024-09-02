#'
#' TWW Growth Model
#'
#' @description Calculates the 3-, 4-, and 5-parameter TWW Growth model estimates. For those who
#'     use the cycle number and fluorescence intensity to analyze real-time, or quantitative polymerase
#'     chain reaction (qPCR), this function will calculate the TWW cycle threshold (\eqn{C_{TWW}}).
#' @usage tww(x, y, start = list(alpha,theta,beta,delta = NULL,phi = NULL), ...)
#' @param x A numeric vector that must be same length as \code{y}
#' @param y A numeric vector that must be same length as \code{x}
#' @param start A numeric list.
#'              The supplied list of numbers are designated as starting parameters, or initial conditions,
#'              inserted into the \pkg{nls} function as \eqn{\alpha}, \eqn{\theta}, \eqn{\beta}, \eqn{\delta},
#'              and \eqn{\phi}, respectively. The length of the list determines which model will be used.
#'              List length should be between 3 and 5. See Details for more information.
#' @param ... Additional optional arguments passed to the \pkg{nls} function.
#' @details
#' The initialized parameters are inserted as a list in \code{start} and are passed to the \pkg{nls} function using the Gauss-Newton algorithm.
#' If you intend to use a 3-parameter model, insert values for \eqn{\alpha}, \eqn{\theta}, and \eqn{\beta} only. If you plan to use the
#' 4-parameter model, you must insert values for \eqn{\delta} in addition to \eqn{\alpha}, \eqn{\theta}, and \eqn{\beta}.
#' If you intend to use the 5-parameter model, you need to insert initial values for all five parameters. The parameters always follows the
#' order \eqn{\alpha}, \eqn{\theta}, \eqn{\beta}, \eqn{\delta}, and \eqn{\phi}. The number of items in the list
#' determines your choice of model.
#' The 3-parameter growth model has the form
#'   \deqn{F(x)=\alpha\ e^{-ArcSinh\left(\theta e^{-\beta x}\right)}}
#' while the 4-parameter growth model follows the equation
#'   \deqn{F(x)=\alpha\ e^{-ArcSinh\left(\theta e^{-\beta x}\right)}+\delta}
#' and the 5-parameter growth model is given by
#'   \deqn{F(x)=\alpha\ e^{-\phi ArcSinh\left(\theta e^{-\beta x}\right)}+\delta}
#' In each of these models, \eqn{\theta} > 0. In the 5-parameter model, \eqn{\phi} > 0.
#' \eqn{C_{TWW}} is only applicable to qPCR data and should not be considered in other cases.
#' @return This function is designed to calculate the parameter estimates, standard errors, and p-values
#'         for the TWW Growth (Decay) Model as well as estimating \eqn{C_{TWW}}, inflection point (poi) coordinates,
#'         sum of squares error (SSE), Akaike information criterion (AIC), and Bayesian information criterion (BIC).
#' @seealso
#'   \code{\link{nls}} to determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model.
#' @examples
#' #Data source: Guescini, M et al. BMC Bioinformatics (2008) Vol 9 Pg 326
#' fluorescence <- c(-0.094311625, -0.022077977, -0.018940959, -0.013167045,
#'                   0.007782761,  0.046403221,  0.112927418,  0.236954113,
#'                  0.479738750,  0.938835708,  1.821600610,  3.451747880,
#'                  6.381471101, 11.318606976, 18.669664284, 27.684433343,
#'                  36.269197588, 42.479513622, 46.054327283, 47.977882896,
#'                  49.141536806, 49.828324910, 50.280629676, 50.552338600,
#'                  50.731472869, 50.833299572, 50.869115345, 50.895051731,
#'                  50.904097158, 50.890804989, 50.895911798, 50.904685027,
#'                  50.899942221, 50.876866864, 50.878926417, 50.876938783,
#'                  50.857835844, 50.858580957, 50.854100495, 50.847128383,
#'                  50.844847982, 50.851447716, 50.841698121, 50.840564351,
#'                  50.826118614, 50.828983069, 50.827490974, 50.820366077,
#'                  50.823743224, 50.857581865)
#'
#' cycle_number <- 1:50
#'
#' #3-parameter model
#' tww(x = cycle_number, y = fluorescence, start = list(40,15.5,0.05))
#'
#' #4-parameter model
#' tww(x = cycle_number, y = fluorescence, start = list(40,15.5,0.05,0),
#'     algorithm = "port")$c_tww
#'
#' #5-parameter model
#' summary(tww(x = cycle_number, y = fluorescence, start = list(40,15.5,0.05,0,1),
#'             algorithm = "port",
#'             control = nls.control(maxiter = 250)))
#'
#' @import stats
#' @export tww

tww <- function(x, y, start = list(alpha,theta,beta,delta = NULL,phi = NULL), ...) {

  if (length(as.vector(x)) != length(as.vector(y))){
    stop("x and y must be the same length.")
  }
  if (!is.numeric(as.vector(x)[1])) {
    stop("x must be a vector, data.frame, or matrix with numeric elements.")
  }
  if (!is.list(start)) {
    stop("start must be a list.")
  }
  if (start[[2]] <= 0) {
    stop("Theta must be greater than 0.")
  }
  if (length(start) == 5) {
    if (start[[5]] <= 0){
      stop("Phi must be greater than 0.")
    }
  }

  start <- start[!sapply(start,is.null)]
  names(start) <- c("alpha","theta","beta","delta","phi")[1:length(start)]

  if (!(length(start) %in% 3:5)) {
    stop("Please ensure the correct number of starting parameters are supplied.
         There must be between 3 and 5 numbers entered in the list.")
  }

  selection <- length(start)-2
  switch(selection,
         {f3 <- function(x,alpha,theta,beta) {alpha * exp(-asinh(theta * exp(-beta * x)))}
         model <- nls(y ~ f3(x,alpha,theta,beta), start = start, ...)
         par = summary(model)$coefficients
         alpha = par["alpha","Estimate"]
         theta = par["theta","Estimate"]
         beta = par["beta","Estimate"]
         C_TWW <- (-sqrt(2*(3+sqrt(5)))+log((1+sqrt(5))/2)+2*log(theta))/(2*beta)
         SSE <- sum((y-f3(x,alpha,theta,beta))^2)
         poi_x <- log(theta**2*(1+sqrt(5))/2)/(2*beta)
         poi_y <- alpha*exp(-asinh(sqrt(2)/sqrt(1+sqrt(5))))},

         {f4 <- function(x,alpha,theta,beta,delta) {alpha * exp(-asinh(theta * exp(-beta * x))) + delta}
         model <- nls(y ~ f4(x,alpha,theta,beta,delta), start = start, ...)
         par = summary(model)$coefficients
         alpha = par["alpha","Estimate"]
         theta = par["theta","Estimate"]
         beta = par["beta","Estimate"]
         delta = par["delta","Estimate"]
         C_TWW <- (-sqrt((3+sqrt(5))/2)*(alpha+delta*exp(asinh(sqrt(2/(1+sqrt(5))))))+alpha*log((1+sqrt(5))*theta/2))/(alpha*beta)
         SSE <- sum((y-f4(x,alpha,theta,beta,delta))^2)
         poi_x <- log(theta**2*(1+sqrt(5))/2)/(2*beta)
         poi_y <- alpha*exp(-asinh(sqrt(2))/sqrt(1+sqrt(5)))+delta},

         {f5 <- function(x,alpha,theta,beta,delta,phi) {alpha * exp(-phi * asinh(theta * exp(-beta * x))) + delta}
         model <- nls(y ~ f5(x,alpha,theta,beta,delta,phi), start = start, ...)
         par = summary(model)$coefficients
         alpha = par["alpha","Estimate"]
         theta = par["theta","Estimate"]
         beta = par["beta","Estimate"]
         delta = par["delta","Estimate"]
         phi = par["phi","Estimate"]
         C_TWW <- -(sqrt(2+phi**2+phi*sqrt(4+phi**2)))*
           (alpha+delta*exp(phi*asinh(sqrt(2)/sqrt(phi*(phi+sqrt(4+phi**2))))))+
           log(sqrt(phi*theta**2*(phi+sqrt(4+phi**2)))/sqrt(2))/beta
         SSE <- sum((y-f5(x,alpha,theta,beta,delta,phi))^2)
         poi_x <- log(sqrt(phi**2*theta**2+phi*sqrt(4+phi**2)*theta**2)/sqrt(2))/beta
         poi_y <- alpha*exp(-phi*asinh(sqrt(2)/sqrt(phi*(phi+sqrt(4+phi**2)))))}
  )

  AIC <- length(x)*log(SSE/length(x)) + 2*length(start)

  BIC <- length(x)*log(SSE/length(x)) + log(length(x))*length(start)

  out <- structure(c(model, list(call = model$call, formula = formula, terms = rownames(par),
                                 coefficients = par),
                     "c_tww" = C_TWW,
                     "sse" = SSE,
                     "aic" = AIC,
                     "bic" = BIC),
                   class = c("nls"))
  return(out)
}

