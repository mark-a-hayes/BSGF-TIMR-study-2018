###################################################################################
##
## BSGF TIMR FATALITY ESTIMATES 
##
## Using the "Carcass" package FOR R and the Korner, Huso, and Erickson estimators to
## estimate fatalities using control and treatment approaches at the BSGF site.
## Please contact Mark Hayes with comments and questions at the email below.
## 
## Mark A. Hayes
## Senior Bat Ecologist
## Normandeau Associates, Inc. 
## Date: 3/12/2018 
## Email: mhayes@normandeau.com
##
## Posted on MAH's GitHub repository: https://github.com/mark-a-hayes 
##
###################################################################################

## Preliminaries

install.packages("carcass")
library("carcass", lib.loc="~/R/win-library/3.4")

remove()
rm(list=ls())

###################################################################################

## Creating a custom estimateN function. 
# Note on plots: The estimateN plot function uses xlim=c(0,50). 
# So, J. Champion created the following custom 'estimateN_custom' function using the estimateN source code

estimateN_custom <- function (count, p = NA, p.lower = NA, p.upper = NA, f = NA, 
                              f.lower = NA, f.upper = NA, s = NA, s.lower = NA, s.upper = NA, 
                              arrival = "discrete", a = 1, a.lower = 1, a.upper = 1, pform = "korner", 
                              d = 1, n = NA, J = NA, maxn = 1000, nsim = 1000, plot = TRUE, 
                              postdist = FALSE, k = 1, x = c(1:10), 
                              plot_lims=c(0,50)) #new parameter for plot limits
{
  if (is.finite(p) & is.finite(f)) {
    warning("f, f.lower, f.upper, s, s.lower, s.upper, and d will be ignored")
  }
  if (a.lower < a.upper) {
    a.a <- shapeparameter(a, a.lower, a.upper)$a
    a.b <- shapeparameter(a, a.lower, a.upper)$b
  }
  if (is.finite(f)[1]) {
    if (f.lower >= f.upper) 
      stop("Something is wrong with the CI of f. If f is known without uncertainty, use posteriorN instead of estimateN")
    if (s.lower >= s.upper) 
      stop("Something is wrong with the CI of s. If s is known without uncertainty, use posteriorN instead of estimateN")
    if (a.lower > a.upper) 
      stop("Something is wrong with the CI of a.")
    if (length(s) == 1 & s == 1) 
      stop("s=1: no removal. This is very unlikely. The formulas implemented here are not made for this case.")
    f.a <- shapeparameter(f, f.lower, f.upper)$a
    f.b <- shapeparameter(f, f.lower, f.upper)$b
    s.a <- shapeparameter(s, s.lower, s.upper)$a
    s.b <- shapeparameter(s, s.lower, s.upper)$b
    Npostdist <- numeric(maxn + 1)
    if (length(s) == 1 & arrival == "discrete") {
      cat("You assume that all carcasses arrive simultaneously \nonce a day right after the time the search takes place (if there is one) \nand that the probability to persist from death to the first day \n(the time the searches normally are performed) is equal to the probability\nto persist from one day to the next. To change that, \nuse a different persistence probability for the first day \nor choose arrival='uniform' to allow carcasses to arrive continuously.\n")
    }
    if (arrival == "uniform" & length(s) > 1) 
      cat("You assume constant (continuous) carcass arrival. To change that, use arrival='discrete'")
    if (arrival == "uniform" & length(s) == 1) 
      cat("You assume that persistence probability is independent of age and that carcasses arrive constantly (continuously).\nTo change that, use different persistence probabilities for every day after \ndeath and/or arrival='discrete'\n")
    if (arrival == "discrete" & length(s) > 1) {
      cat("You assume that all carcasses arrive simultaneously and that the first \npersistence probability is the probability to persist from death to \nthe first day at the time when searches take place (if there is a search at that day).\n")
    }
    for (i in 1:nsim) {
      fr <- rbeta(length(f), f.a, f.b)
      if (length(s) > 1) {
        sr <- numeric(length(s))
        sr[1] <- rbeta(1, s.a[1], s.b[1])
        psr <- pbeta(sr[1], s.a[1], s.b[1])
        srelevant <- s > 0.001
        sr[srelevant] <- qbeta(psr, s.a[srelevant], s.b[srelevant])
        sr[!srelevant] <- 0
      }
      if (length(s) == 1) 
        sr <- rbeta(1, s.a, s.b)
      sr[sr > 0.9999] <- 0.9999
      if (arrival == "uniform") 
        sr <- integrate.persistence(sr, n = n, d = d)
      ar <- ifelse(a.lower < a.upper, rbeta(1, a.a, a.b), 
                   a)
      if (pform == "korner") 
        pr <- pkorner(s = sr, f = fr, d = d, n = n, k = k, 
                      search.efficiency.constant = ifelse(k == 1, 
                                                          TRUE, FALSE))
      if (pform == "huso") 
        pr <- phuso(s = sr, f = fr, d = d)
      if (pform == "erickson") 
        pr <- perickson(t.bar = -1/log(sr), f = fr, d = d)
      if (pform == "etterson") {
        if (length(J) == 1) 
          if (is.na(J)) 
            stop("Please, provide J")
        if (length(sr) == 1) 
          sr <- rep(sr, sum(J))
        if (length(fr) == 1) 
          fr <- rep(fr, length(J))
        pr <- ettersonEq14v2(s = sr, f = fr, J = J)
      }
      postNtemp <- posteriorN(nf = count, p = pr * ar, 
                              maxN = maxn, plot = FALSE, dist = TRUE)
      Npostdist <- Npostdist + postNtemp$pN
    }
    if (pform == "korner") 
      pm <- pkorner(s = s, f = f, d = d, n = n, k = k, 
                    search.efficiency.constant = ifelse(k == 1, TRUE, 
                                                        FALSE))
    if (pform == "huso") 
      pm <- phuso(s = s, f = f, d = d)
    if (pform == "erickson") 
      pm <- perickson(t.bar = -1/log(s), f = f, d = d)
    if (pform == "etterson") {
      if (length(J) == 1) 
        if (is.na(J)) 
          stop("Please, provide J")
      if (length(s) == 1) 
        s <- rep(s, sum(J))
      if (length(f) == 1) 
        f <- rep(f, length(J))
      pm <- ettersonEq14v2(s = s, f = f, J = J)
    }
    HT.estimate <- count/(pm * a)
  }
  if (is.finite(p)) {
    if (p.lower >= p.upper) 
      stop("Something is wrong with the CI of p. If p is known without error, use posteriorN instead of estimateN")
    p.a <- shapeparameter(p, p.lower, p.upper)$a
    p.b <- shapeparameter(p, p.lower, p.upper)$b
    HT.estimate <- count/(p * a)
    Npostdist <- numeric(maxn + 1)
    for (i in 1:nsim) {
      pr <- rbeta(1, p.a, p.b)
      ar <- ifelse(a.lower < a.upper, rbeta(1, a.a, a.b), 
                   a)
      postNtemp <- posteriorN(nf = count, p = pr * ar, 
                              maxN = maxn, plot = FALSE, dist = TRUE)
      Npostdist <- Npostdist + postNtemp$pN
    }
  }
  Npostdist.sc <- Npostdist/nsim
  indexLower <- cumsum(Npostdist.sc) < 0.025
  indexMedian <- cumsum(Npostdist.sc) < 0.5
  indexUpper <- cumsum(Npostdist.sc) < 0.975
  lower <- min(c(0:maxn)[!indexLower])
  estimate.median <- min(c(0:maxn)[!indexMedian])
  upper <- min(c(0:maxn)[!indexUpper])
  plarger <- 1 - cumsum(Npostdist.sc)[is.element(c(0:maxn), 
                                                 x)]
  names(plarger) <- paste0("x=", x)
  if (plot) 
    plot(0:maxn, Npostdist.sc, type = "h", lwd = 5, lend = "butt", 
         xlab = "Number of fatalities", ylab = "Posterior density", 
         xlim = plot_lims) #use new plot limits parameter
  if (!postdist) {
    result <- list(estimate = estimate.median, lower = lower, 
                   upper = upper, HT.estimate = HT.estimate, P.true.larger.x = plarger)
  }
  if (postdist) {
    result <- list(estimate = estimate.median, lower = lower, 
                   upper = upper, HT.estimate = HT.estimate, postdist = Npostdist.sc, 
                   P.true.larger.x = plarger)
  }
  return(result)
}


## Test that the new function works for a new xlim

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
          s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000,
          plot=TRUE)

estimateN_custom(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
                 s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000,
                 plot=TRUE, plot_lims=c(0,500))

###################################################################################

## Using the BSGF carcass results for counts (count), searcher efficiency (f), and carcass
## persistence probability (s) estimates to estimate total fatalities (detected and not detected).

# Fatality estimator parameters are for 95% confidence intervals, and are estimated 
# using the results reported by Gruver (2016). M. Hayes took these estimates, which
# were for 90% CI and converted to 95% CI's (MAH, calculations on 11/27/2017).

# f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890, s.upper=0.950

# NOTE: the estimateN function uses xlim=c(0,50), so here I use a costum version of the function
# so that the limits of the x axis can be expanded and changed as needed.  


## Pooled (all species) for control turbine fatalities, with count = 189:
  
# The Erickson estimator, pooled then per turbine
  
  
estimateN_custom(count=189, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
            s.upper=0.950, d=1, pform="erickson", n=78, maxn=600, nsim=10000,
            plot=TRUE, plot_lims=c(0,600))

dev.print(tiff, "pooled_c_e.tiff", height=4, width=6.5, units='in', res=300)


# The Huso estimator

estimateN_custom(count=189, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=600, nsim=10000,
          plot=TRUE, plot_lims=c(0,600))

dev.print(tiff, "pooled_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator, with count as total number for the 10 control turbines.

estimateN_custom(count=189, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=600, nsim=10000,
          plot=TRUE, plot_lims=c(0,600))

dev.print(tiff, "p_c_K.tiff", height=4, width=6.5, units='in', res=300)

# Pooled treatment fatalities, with count = 31:

# The Erickson estimator

estimateN_custom(count=31, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500))

dev.print(tiff, "p_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=31, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500))

dev.print(tiff, "p_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator, with count as total number for the 10 treatment turbines.

estimateN_custom(count=31, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500))

dev.print(tiff, "p_t_k.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################

## LACI control fatalities, with count = 48:

# The Erickson estimator

estimateN_custom(count=48, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=48, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=48, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LACI treatment fatalities, with count = 11:

# The Erickson estimator

estimateN_custom(count=11, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=11, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_h.tiff", height=4, width=6.5, units='in', res=300)


# The Korner estimor

estimateN_custom(count=11, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_k.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################

## LABO control fatalities, with count = 39:

# The Erickson estimator

estimateN_custom(count=39, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_e.tiff", height=4, width=6.5, units='in', res=300)


# The Huso estimator

estimateN_custom(count=39, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=39, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LABO treatment fatalities, with count = 6:

# The Erickson estimator

estimateN_custom(count=6, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=6, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=6, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_k.tiff", height=4, width=6.5, units='in', res=300)

###################################################################################

## LANO control fatalities, with count = 44:

# The Erickson estimator

estimateN_custom(count=44, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=44, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=44, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LANO treatment fatalities, with count = 1:

# The Erickson estimator

estimateN_custom(count=1, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_e.tiff", height=4, width=6.5, units='in', res=300)


# The Huso estimator

estimateN_custom(count=1, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=1, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_k.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################

## EPFU control fatalities, with count = 27:

# The Erickson estimator

estimateN_custom(count=27, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=27, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=27, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_k.tiff", height=4, width=6.5, units='in', res=300)

## EPFU treatment fatalities, with count = 7:

# The Erickson estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_k.tiff", height=4, width=6.5, units='in', res=300)

###################################################################################

## MYLU control fatalities, with count = 29:

# The Erickson estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_k.tiff", height=4, width=6.5, units='in', res=300)

## MYLU treatment fatalities, with count = 3:

# The Erickson estimator

estimateN_custom(count=3, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=3, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=3, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_t_k.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################
##
## End
##
###################################################################################