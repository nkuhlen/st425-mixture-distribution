# Part II #################################################################

# Clean workspace.
rm(list=ls())

# install.packages("ggplot2")
library("ggplot2")

# Set random seed for replicability.
set.seed(425)

# Q.1 =====================================================================

# Define functions --------------------------------------------------------

mixture_pdf <- function(x, mu) {
  # Calculate the value of the probability density function (pdf) at x given
  # parameter mu.
  # Arguments:
  #     x : double
  #       Point that density is calculated at.
  #     mu : double
  #       Shape parameter of the density.

  1 / (4 * (1 - exp(-2))) * x * exp(-x ^ 2 / 2) * ((0 < x) & (x < 2)) + 
    3 / (4 * sqrt(2 * pi)) * exp(-(x - mu) ^ 2 / 2)
}

mixture_cdf <- function(x, mu) {
  # Calculate the value of the cumulative distribution function (cdf) at x given
  # parameter mu.
  # Arguments:
  #     x : double
  #       Point that cdf is calculated at.
  #     mu : double
  #       Shape parameter of the cdf.

  if (x <= 0) {
    3 / 4 * pnorm(x, mean=mu, sd=1)
  }
  else if (x >= 2) {
    1 / 4 + 3 / 4 * pnorm(x, mean=mu, sd=1)
  }
  else {
    (1 - exp(-x ^ 2 / 2)) / (4 * (1 - exp(-2))) + 3 / 4 * pnorm(x, mean=mu, sd=1)
  }
}

mixture_qf <- function(p, mu) {
  # Calculate the value of the distribution function for a given quantile p and
  # a given parameter mu.
  # Arguments:
  #     p : double
  #       Probability value for calculation of quantile. p has to lie in the
  #       interval [0, 1].
  #     mu : double
  #       Shape parameter of the quantile function.

  if (p < 0 || p > 1) {
    stop("Quantile function is only defined for values of p between 0 and 1.")
  }
  else if (p == 0) {
    return(-Inf)
  }
  else if (p == 1) {
    return(Inf)
  }
  else {
    # Solve for the x for which mixture_cdf(x,mu) - p = 0 to find p-th quantile. 
    # Define helper function for given mu that captures this condition.
    helper_function <- function(x) {
      mixture_cdf(x, mu) - p
    }

    # Find quantile as the root to helper_function.
    quantile <- uniroot(helper_function, lower=-10, upper=10, tol=1e-10)$root
    
    return(quantile)
  }
}

# Vectorise functions defined above.
mixture_pdf_v <- Vectorize(mixture_pdf)  
mixture_cdf_v <- Vectorize(mixture_cdf)  
mixture_qf_v <- Vectorize(mixture_qf)  

# Generate random numbers -------------------------------------------------

rmixture <- function(n, mu) {
  # Generate a random sample from the mixture distribution specified by mu using
  # the inverse transform method.
  # Arguments:
  #     n : double
  #       Number of random numbers to be generated.
  #     mu : double
  #       Shape parameter of the quantile function.

  # Draw uniform sample.
  quantiles_sample <- runif(n)

  # Recover sample from original distribution via quantile function.
  x <- mixture_qf_v(quantiles_sample, mu)
  
  return(x)
}

# Q.2 =====================================================================

# Define values for mu.
values_mu <- c(-2, 0, 2)

for (i in 1:length(values_mu)) {
  # Define temporary callable functions for given mu.
  pdf <- function(x) mixture_pdf_v(x, mu=values_mu[i])
  cdf <- function(x) mixture_cdf_v(x, mu=values_mu[i])
  qf <- function(x) mixture_qf_v(x, mu=values_mu[i])
  
  # Plot the density for given mu.
  print(
    ggplot(data.frame(x=c(-6, 6)), aes(x=x)) + ylim(0, 0.5) +
      stat_function(fun=pdf, geom="line") + xlab("x") + ylab("f(x)") + 
      labs(title=paste("PDF for mu =", values_mu[i])) +
      theme(plot.title = element_text(size=6)) + theme_bw()
  )
  
  # Plot distribution for given mu.
  print(
    ggplot(data.frame(x=c(-6, 6)), aes(x=x)) + ylim(0, 1) +
      stat_function(fun=cdf, geom="line") + xlab("x") + ylab("F(x)") + 
      labs(title=paste("CDF for mu =", values_mu[i])) +
      theme(plot.title = element_text(size=6)) + theme_bw()
  )
  
  # Plot quantile function for given mu.
  print(
    ggplot(data.frame(x=c(0, 1)), aes(x=x)) +
      stat_function(fun=qf, geom="line") + xlab("p") + ylab("Q(p)") + 
      labs(title=paste("QF for mu =", values_mu[i])) +
      theme(plot.title = element_text(size=6)) + theme_bw()
  )
}

# Q.3 =====================================================================

# Simulation --------------------------------------------------------------

# Set parameters for Monte Carlo simulation.
par_mu <- -1
n <- 10000

# Create empty list to store samples.
x <- list()

# Draw samples for three random variables following the mixture distribution.
for (i in 1:3) {
  x[[i]] <- rmixture(n, par_mu)
}

# Calculate sample for new random variable z.
z <- (x[[1]] - x[[2]]) / (1 + abs(x[[3]]))

# Monte Carlo estimate for the expectation from exercise is the mean of z.
expec <- mean(z)

# Calculate standard error of estimates.
std_err <- sd(z) / sqrt(n)

# Print results.
cat("Monte Carlo Estimate for Expectation:", expec, "\nStandard Error:", std_err)


# Visualisation -----------------------------------------------------------

# Create empty list to store estimates.
mc <- list()

for (i in 1:length(z)) {
  # Calculate mean of the first i values of z.
  mc[i] <- mean(z[1:i])
}

# Plot the estimated means of z for the corresponding number of draws.
plot(1:length(mc), mc, type = "l", ylab=expression(bar(Xn)), xlab="n")

# Q.4 =====================================================================


# Expectation Approach ----------------------------------------------------

mu_hat <- function(data){
  # Estimate shape parameter mu of mixture distribution based on given data 
  # sample.
  # Arguments:
  #     data : double
  #       Data sample containing values used for the estimation.

  # Calculate integral for the first part of the distribution.
  integral <- function(x) {
    x ^ 2 * 1 / (4 * (1 - exp(-2))) * exp(-x ^ 2 / 2)
  }
  integral_value <- integrate(integral, lower=0, upper=2)$value

  # Compute estimate for mu.
  mu <- 4 / 3 * (mean(data) - integral_value)

  return(mu)
}

# Open interactive window to choose and import data sample.
dat <- read.table(file.choose())

# Apply function to estimate mu.
mu_hat(dat$V1)


# Maximum Likelihood Estimation -------------------------------------------

mu_hat_mle <- function(
  data, lower=ceiling(min(dat_num)/10)*10-10, upper=ceiling(max(dat_num)/10)*10,
  step_size=0.1
){
  # Estimate shape parameter mu of mixture distribution based on given data 
  # sample by Maximum Likelihood Estimation (MLE).
  # Arguments:
  #     data : double
  #       Data sample containing values used for the estimation.
  #     lower : double
  #       Lower bound of interval used to find maximum of likelihood function.
  #     upper : double
  #       Upper bound of interval used to find maximum of likelihood function.
  #     step_size : double
  #       Stepsize for the grid used to find optimal value of mu.

  # Define grid to search for maximum value of likelihood function.
  mu <- seq(lower, upper, step_size)

  # Create vector of same length as grid to store values of likelihood function.
  likelihood_f <- 1:length(mu)
  
  # Create vector of same length as data set.
  f <- 1:length(dat_num)

  # Loop over different values of mu.
  for(j in 1:length(mu)) {
    likelihood_f[j] <- 1
      
      # Loop over observations to calculate value of likelihood function given
      # mu from outer loop.
      for(i in 1:length(dat_num)) {
        f[i] <- mixture_pdf(dat_num[i], mu[j])
        likelihood_f[j] <- likelihood_f[j] * f[i]
      }
  }

  # Find mu that maximises likelihood function.
  opt_mu <- mu[which.max(likelihood_f)]

  return(opt_mu)
}

# Convert data frame to numeric matrix.
dat_num <- data.matrix(dat) 

# Apply function to estimate mu.
mu_hat_mle(dat_num)

# Plots -------------------------------------------------------------------

# Store estimated mu from expectation approach.
estimated_mu <- mu_hat(dat$V1)

# Define callable function.
pdf <- function(x) mixture_pdf_v(x, mu=estimated_mu)

# Plot the density for given mu.
ggplot(dat, aes(x=V1)) + geom_histogram(aes(y=..density..), bins=40) +
  stat_function(fun=pdf, geom="line", colour="red") + xlab("x") +
  theme_classic() 

# Plot kernel density estimate for given data sample.
ggplot(dat, aes(x=V1)) + geom_histogram(aes(y=..density..), bins=40) + 
  geom_density(colour="red") + xlab("x") + theme_classic() 
