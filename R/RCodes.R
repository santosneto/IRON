#'@name IRON
#'
#'@aliases diron
#'
#'@title The IRON distribution.
#'
#'@description Density function, distribution function, quantile function and random generation for the IRON distribution.
#'
#'@param x,q vetor of quantiles.
#'@param p vector of probabilities.
#'@param n number of observations.
#'@param dist kernel used.
#'@param tau quantile value.
#'@param shape shape parameter.
#'@param scale scale parameter.
#'@param df degrees of freedom.
#'@param kappa shape parameter.
#'@param log,log.p logical; if TRUE, the log-density or log(p) is used.  
#'@param ... additional arguments to be passed.
#'
#'@details The density function of the IRON distribution is
#'
#'\deqn{g(t|\lambda,\beta,\alpha) = \alpha f(a_t)[F(a_t)]^{\alpha-1} \vdot \frac{t^{-3/2}(t+\beta)}{2\lambda\sqrt{\beta}}\qc t,\lambda,\beta,\alpha,}
#'where \eqn{f(.)} and \eqn{F(.)} are, respectively, the density function and cumulative distribution of a symmetric distribution.
#'
#'@return diron gives the density, piron gives the distribution function, qiron gives the quantile function and riron generate pseudo-random numbers.
#' 
#'@author Manoel Santos-Neto \email{manoel.ferreira at professor.ufcg.edu.br} and Diego I. Gallardo \email{diego.gallardo.mateluna at gmail.com}
#'
#' 
#'@examples 
#' diron(x = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
#' piron(q = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
#' qiron(p = 0.5, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
#' riron(n = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
#' @export
#' 
#' @importFrom stats pgamma 
#'@import gamlss

diron <- function(x, dist = "norm", tau, shape, scale, kappa = NULL, log = FALSE, ...){
  g <- get(paste("d", dist, sep = ""), mode = "function")
  G <- get(paste("p", dist, sep = ""), mode = "function")
  
  a   <- -log(tau)/log(2)
  a.x <- (1.0/shape) * (sqrt(x/scale) - sqrt(scale/x))
  #A.x <- x^(-3/2) * (x + scale)/(2 * shape * sqrt(scale))
  
  if (dist != 'PE2') {
    log.f <- log(a) + g(x = a.x, log = TRUE, ...) + (a - 1)*G(q = a.x, log.p = TRUE, ...) - (3/2) * log(x) + log(x + scale) - log(2) - log(shape) - (1/2) * log(scale) 
  } else{
    term1 <- 0.5
    term2 <- (sign(a.x)/2)*(1/gamma(1/kappa))*pgamma(q = a.x, shape = 1/kappa)
    log.f <- log(a) + log(kappa) - log(2) - lgamma(1/kappa) - exp(kappa*log(abs(a.x))) + (a - 1)*(log(term1 + term2)) - (3/2) * log(x) + log(x + scale) - log(2) - log(shape) - (1/2) * log(scale) 
  }
  if (!log) log.f <- exp(log.f)
  
  log.f
}

#'@name IRON
#'
#'@aliases qiron
#'
#'@export
qiron <- function(p, shape, scale, tau, dist = "norm", ...){ 
  Q <- get(paste("q", dist, sep = ""), mode = "function")
  a   <- -log(tau)/log(2)
  p.new <- p^(1/a) 
  z.p <- Q(p = p.new, ...)
  term1 <- (shape/2) * z.p
  f_q <- scale * ((term1 + sqrt( (term1^2) + 1))^2)
  
  f_q
}

#'@name IRON
#'
#'@aliases piron
#'
#'@export
piron <- function(q, shape, scale, tau, dist = "norm", log.p = FALSE, ...){ 
  CDF <- get(paste("p", dist, sep = ""), mode = "function")
  a   <- -log(tau)/log(2)
  a.q <- (1.0/shape) * (sqrt(q/scale) - sqrt(scale/q))
  cdf <- CDF(a.q, log.p = log.p, ...)^a
  
  cdf 
}


#'@name IRON
#'
#'@aliases riron
#'
#'
#'@export
#'
#'@importFrom stats runif
#'@import gamlss

riron <- function(n, shape, scale, tau, dist = "normal", df = 4, kappa = NULL, ...)
{
  
  unif_number <- runif(n)
  
  if (dist == 'normal') {
    qexpbs_norm <- qiron
    prn <- qexpbs_norm(p = unif_number, shape = shape, scale = scale, tau = tau)
  }else if (dist == 't') {
    qexpbs_t <- function(p, ...) qiron(p, dist = 't', df = df, ...)
    prn <- qexpbs_t(p = unif_number, shape = shape, scale = scale, tau = tau)
  }else if (dist == 'pe') {
    qexpbs_pe <- function(p, ...) qiron(p, dist = 'PE2', ...)
    prn <- qexpbs_pe(p = unif_number, shape = shape, scale = scale, nu = kappa, tau = tau)
  }else if (dist == 'cauchy') {
    qexpbs_ca <- function(p, ...) qiron(p, dist = 't', df = 1, ...)
    prn <- qexpbs_ca(p = unif_number, shape = shape, scale = scale, tau = tau)
  } else{
    qexpbs_logis <- function(p, ...) qiron(p, dist = 'logis', ...)
    prn <- qexpbs_logis(p = unif_number, shape = shape, scale = scale, tau = tau)
  }
  
  prn
}

#'@name fit
#'
#'@title Fitting Quantile Regression Models
#'
#'@description quant_reg is used to fit quantile regression models, specified by giving a symbolic description of the linear predictor and a description of the kernel.
#'
#'@aliases quant_reg
#'
#'@param formula an object of class "formula".
#'@param link link function to be used in the model.
#'@param family a description of the kernel.
#'@param tau quantile.
#'@param data an optional data frame containing the variables of the model.
#'@param out.list logical.
#'@param estimate.df logical.
#'@param df degrees of freedom.
#'@param model object of class iron.
#'@param npoints number of points to be identified.
#'@param plot logical. 
#'@param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.
#'@param k number of replications for envelope construction. 
#'@param color a specification for the default envelope color.
#'@param xlabel a label for the x axis.
#'@param ylabel a label for the y axis.
#'@param font the name of a font family for x and y axis.
#'@param cex.axis The magnification to be used for axis annotation relative to the current setting of cex.
#'@param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#'@param ... additional arguments to be passed.
#'
#'
#'@examples
#'
#'set.seed(2)
#'y <- VGAM::rbisa(100, shape = 4, scale = 1)
#'quant_reg(y ~ 1, link = 'log', family = 'normal', tau = 0.9)
#'quant_reg(y ~ 1, link = 'log', family = 'pe', tau = 0.5)
#'quant_reg(y ~ 1, link = 'log', family = 't', tau = 0.5)
#'quant_reg(y ~ 1, link = 'log', family = 'logis', tau = 0.5)
#'quant_reg(y ~ 1, link = 'log', family = 'pe', tau = 0.2)
#'quant_reg(y ~ 1, link = 'log', family = 't', tau = 0.8)
#'
#'@export
#'
#'@importFrom stats AIC
#'@importFrom stats Gamma
#'@importFrom stats glm
#'@importFrom stats make.link
#'@importFrom stats model.matrix
#'@importFrom stats uniroot
#'@importFrom stats var
#'@importFrom maxLik maxLik
#'@import  gamlss

quant_reg <- function(formula, link = 'identity', family = 'normal', tau = 0.5, data = NULL, out.list = FALSE, estimate.df = FALSE, df = 4, ...){
  
  linkstr   <- link
  linkobj   <- make.link(linkstr)
  linkinv   <- linkobj$linkinv
  mu.eta    <- linkobj$mu.eta
  model_aux <- glm(formula, family = Gamma(link = link), data = data)  
  Y <- model_aux$y
  X <- model.matrix(model_aux)
  
  if (length(Y) < 1) 
    stop("empty model")
  if ((min(Y) < 0)) 
    stop("invalid dependent variable, all observations must be positive") 
  
  n <- length(Y)
  p <- ncol(X)
  xb <- mean(Y)
  xh <- 1/mean(1/Y)
  shape_start <- sqrt(2*(sqrt(xb/xh) - 1 ))    
  beta_start <- model_aux$coefficients[1:p]
  
  
  if (family == 'normal') {
    
    pdf_expbs_norm <- diron
    
    log_lik <- function(par, y, x, tau, link){
      linkstr   <- link
      linkobj   <- make.link(linkstr)   
      linkinv   <- linkobj$linkinv   
      shape     <- par[1]
      coef_beta <- par[-1]
      eta       <- as.vector(x %*% coef_beta)
      scale     <- linkinv(eta)
      
      pdf_expbs_norm(x = y, tau = tau, shape = shape, scale = scale, log = TRUE)
    }
    
    start_par <- c(shape = shape_start, beta_start)
    
    
    options(warn = -1)   # Disable warning messages globally  
    fit <- maxLik(logLik = log_lik, start = start_par, method = "BHHH", y = Y, x = X, tau = tau, link = link, control = list(iterlim = 600))
    
    
  } else if (family == 't') {
    
    if (estimate.df == FALSE) {
      df.t <- df
      pdf_expbs_t <- function(x, ...) diron(x = x, dist = 't', df = df.t, ...) 
      
      log_lik <- function(par, y, x, tau, link){
        linkstr   <- link
        linkobj   <- make.link(linkstr)   
        linkinv   <- linkobj$linkinv   
        shape     <- par[1]
        coef_beta <- par[-1]
        eta       <- as.vector(x %*% coef_beta)
        scale     <- linkinv(eta)
        
        pdf_expbs_t(x = y, tau = tau, shape = shape, scale = scale, log = TRUE)
      }
      start_par <- c(shape = shape_start, beta_start)
      
    } else{
      pdf_expbs_t <- function(x, ...) diron(x = x, dist = 't', ...)  
      
      log_lik <- function(par, y, x, tau, link){
        linkstr   <- link
        linkobj   <- make.link(linkstr)   
        linkinv   <- linkobj$linkinv   
        shape     <- par[1]
        df <- par[2] 
        coef_beta <- par[-(1:2)]
        eta       <- as.vector(x %*% coef_beta)
        scale     <- linkinv(eta)
        
        pdf_expbs_t(x = y, tau = tau, shape = shape, scale = scale, df = df, log = TRUE)
      }
      start_par <- c(shape = shape_start, df = 4, beta_start)
      
    }
    
    options(warn = -1)   # Disable warning messages globally  
    fit <- maxLik(logLik = log_lik, start = start_par, method = "BHHH", y = Y, x = X, tau = tau, link = link, control = list(iterlim = 600))
    
    
  } else if (family == 'cauchy') {
    
    pdf_expbs_ca <- function(x, ...) diron(x = x, dist = 't', df = 1, ...)
    
    log_lik <- function(par, y, x, tau, link){
      linkstr   <- link
      linkobj   <- make.link(linkstr)   
      linkinv   <- linkobj$linkinv   
      shape     <- par[1]
      coef_beta <- par[-1]
      eta       <- as.vector(x %*% coef_beta)
      scale     <- linkinv(eta)
      
      pdf_expbs_ca(x = y, tau = tau, shape = shape, scale = scale, log = TRUE)
    }
    
    start_par <- c(shape = shape_start, beta_start)
    
    
    options(warn = -1)   # Disable warning messages globally  
    fit <- maxLik(logLik = log_lik, start = start_par, method = "BHHH", y = Y, x = X, tau = tau, link = link, control = list(iterlim = 600))
    
    
  } else if (family == 'pe') {
    
    pdf_expbs_pe <- function(x, ...) diron(x = x, dist = 'PE2', ...)
    
    log_lik <- function(par, y, x, tau, link){
      linkstr   <- link
      linkobj   <- make.link(linkstr)   
      linkinv   <- linkobj$linkinv   
      shape     <- par[1]
      kappa     <- par[2]
      coef_beta <- par[-(1:2)]
      eta       <- as.vector(x %*% coef_beta)
      scale     <- linkinv(eta)
      
      pdf <- pdf_expbs_pe(x = y, tau = tau, shape = shape, scale = scale, kappa = kappa, log = TRUE)
      
      pdf
    }
    
    
    var.y <- var(Y)
    fnToFindRoot = function(kappa,vx){
      return(gamma(3/kappa)/gamma(1/kappa) - vx)
    }
    kappa_est <- uniroot(fnToFindRoot, c(0.01,1000), tol = 0.0001, vx = var.y)$root 
    
    fit0 <- quant_reg(formula, family = 'normal', link = link, tau = tau, data = data, out.list = TRUE)$coef
    start_par <- c(fit0[1], kappa = kappa_est, fit0[-1])
    
    
    options(warn = -1)   # Disable warning messages globally  
    fit <- maxLik(logLik = log_lik, start = start_par, method = "BHHH", y = Y, x = X, tau = tau, link = link, control = list(iterlim = 600))
    
    
  }  else{
    
    
    pdf_expbs_logis <- function(x, ...) diron(x = x, dist = 'logis', ...)
    
    log_lik <- function(par, y, x, tau, link){
      linkstr   <- link
      linkobj   <- make.link(linkstr)   
      linkinv   <- linkobj$linkinv   
      shape     <- par[1]
      coef_beta <- par[-1]
      eta       <- as.vector(x %*% coef_beta)
      scale     <- linkinv(eta)
      
      pdf_expbs_logis(x = y, tau = tau, shape = shape, scale = scale, log = TRUE)
    }
    
    start_par <- c(shape = shape_start, beta_start)
    
    
    options(warn = -1)   # Disable warning messages globally  
    fit <- maxLik(logLik = log_lik, start = start_par, method = "BHHH", y = Y, x = X, tau = tau, link = link, control = list(iterlim = 600))
    
  }
  
  if (out.list == TRUE) {
    coef <- fit$estimate
    npar <- length(coef)
    par.df <- npar - p 
    eta <- drop(X %*% coef[(par.df + 1):npar])
    mu <- linkinv(eta)
    residuals <- (Y - mu)/mu.eta(eta)
    conv <- fit$code
    aic <- AIC(fit)
    if (estimate.df == FALSE) out.df <- df else out.df <- fit$estimate[2]
    
    list(coefficients = coef, residuals = residuals, fitted.values = mu, family = family, 
         linear.predictors = eta, y = Y, converged = conv, hess = fit$hessian, vcov = solve(-fit$hessian), 
         AIC = aic, tau = tau, link = link, formula = formula, est.df = estimate.df, df = out.df)
  } else return(summary(fit))
  
}

#'@name fit
#'
#'@aliases res_quant
#'
#'
#'@export
#'
#'@importFrom stats qnorm
res_quant <- function(model){
  
  y <- model$y
  tau <- model$tau
  estimates <- model$coefficients
  scale_est <- model$fitted.values
  
  if (model$family == 'pe') {
    shape_est <- estimates[1]
    kappa_est <- estimates[2]
  } else if (model$family == 't') {
    if (model$est.df == TRUE) {
      shape_est <- estimates[1]
      df_est <- estimates[2]
    } else{
      df_value <- model$df  
      shape_est <- estimates[1]
    }
  } else{
    shape_est <- estimates[1]
  }
  
  if (model$family == 'normal') {
    pexpbs_norm <- piron
    prob <- pexpbs_norm(q = y, shape = shape_est, scale = scale_est, tau = tau)
  } else if (model$family == 't') {
    if (model$est.df == FALSE) {
      pexpbs_t <- function(q, ...) piron(q, dist = 't', df = df_value, ...)
      prob <- pexpbs_t(q = y, shape = shape_est, scale = scale_est, tau = tau)
    } else{
      pexpbs_t <- function(q, ...) piron(q, dist = 't', ...)
      prob <- pexpbs_t(q = y, shape = shape_est, scale = scale_est, tau = tau, df = df_est)
    }
  } else if (model$family == 'pe') {
    pexpbs_pe <- function(q, ...) piron(q, dist = 'PE2', ...)
    prob <- pexpbs_pe(q = y, shape = shape_est, scale = scale_est, nu = kappa_est, tau = tau)
  } else if (model$family == 'cauchy') {
    pexpbs_ca <- function(q, ...) piron(q, dist = 't', df = 1, ...)
    prob <- pexpbs_ca(q = y, shape = shape_est, scale = scale_est, tau = tau)
  } else{
    pexpbs_logis <- function(q, ...) piron(q, dist = 'logis', ...)
    prob <- pexpbs_logis(q = y, shape = shape_est, scale = scale_est, tau = tau)
  }
  
  
  res <- qnorm(prob)
  
  res
  
}

#'@name fit
#'
#'@aliases cooks_dist
#'
#'
#'@export
#'
#'@importFrom graphics par
#'@importFrom grDevices dev.new
#'@importFrom graphics identify

cooks_dist <- function(model, npoints = 0, plot = FALSE, data, cex=1, ...)
{
  
  y <- model$y
  theta <- model$coefficients
  inv.lpp <- model$vcov
  n <- length(y)    
  values <- vector()
  npar <- length(theta)
  
  if (model$family == 't')
  {
    if (model$est.df == FALSE) {
      for (i in 1:n) {
        datai <- as.data.frame(data[-i,])
        fit <- quant_reg(model$formula, link = model$link, family = model$family, tau = model$tau, out.list = TRUE, data = datai) 
        thetai <- fit$coefficients
        values[i] <- (1/(npar))*t(theta - thetai) %*% inv.lpp %*% (theta - thetai)
      }
    } else{
      for (i in 1:n) {
        datai <- as.data.frame(data[-i,])
        fit <- quant_reg(model$formula, link = model$link, family = model$family, tau = model$tau, estimate.df = TRUE, out.list = TRUE, data = datai)  
        thetai <- fit$coefficients
        values[i] <- (1/(npar))*t(theta - thetai) %*% inv.lpp %*% (theta - thetai)
      }
    }
  } else{
    for (i in 1:n)
    {
      datai <- as.data.frame(data[-i,])
      fit <- quant_reg(model$formula, link = model$link, family = model$family, tau = model$tau, out.list = TRUE, data = datai) 
      thetai <- fit$coefficients
      values[i] <- (1/(npar))*t(theta - thetai) %*% inv.lpp %*% (theta - thetai)
    }
    
  }
  if (plot == TRUE)
  {
    dev.new()
    par(mar = c(4.0,4.0,0.1,0.1))   
    plot(values, ...)
    if (npoints != 0) identify(values,n = npoints, cex = cex)
  } else return(values)
}




#'@name fit
#'
#'@aliases envelope
#'
#'@export
#'
#'@importFrom ggplot2 ggplot
#'@importFrom stats qqnorm
#'@importFrom stats rnorm

envelope <- function(model, k = 100, color = "grey50", xlabel = "Theorical Quantile", ylabel = "Empirical Quantile", font = "serif", cex.axis=1, cex.lab=1)
{
  par(mar = c(4,4.5,0.1,0.1))
  
  n   <- length(model$y)
  td  <- res_quant(model)
  re  <- matrix(0,n,k)
  
  for (i in 1:k)
  {
    y1 <- rnorm(n)
    re[, i] <- sort(y1)
  }
  
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)
  
  for (l in 1:n)
  {
    eo <- sort(re[l,])
    e10[l] <- eo[ceiling(k*0.01)]
    e20[l] <- eo[ceiling(k*(1 - 0.01))]
    e11[l] <- eo[ceiling(k*0.05)]
    e21[l] <- eo[ceiling(k*(1 - 0.05))]
    e12[l] <- eo[ceiling(k*0.1)]
    e22[l] <- eo[ceiling(k*(1 - 0.1))]
  }
  
  a   <- qqnorm(e10, plot.it = FALSE)$x
  r   <- qqnorm(td, plot.it = FALSE)$x
  xb  <- apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x
  
  df  <- data.frame(r = r,xab = a,emin = cbind(e10,e11,e12),emax = cbind(e20,e21,e22),xb = xb,td = td,rxb = rxb)
  ggplot(df, aes(r,td)) + geom_ribbon(aes(x = xab, ymin = emin.e10, ymax = emax.e20),fill = color,alpha = 0.5)  + geom_ribbon(aes(x = xab, ymin = emin.e11, ymax = emax.e21),fill = color,alpha = 0.5) + geom_ribbon(aes(x = xab, ymin = emin.e12, ymax = emax.e22),fill = color,alpha = 0.5) + scale_fill_gradient(low = "grey25", high = "grey75") + geom_point() + geom_line(aes(rxb,xb),lty = 2) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())  + theme(text = element_text(size = 10,family = font))
}

