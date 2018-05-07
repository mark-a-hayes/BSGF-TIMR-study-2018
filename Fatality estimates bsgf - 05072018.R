###################################################################################
##
## BSGF TIMR FATALITY ESTIMATES AND ANALYSIS 
##
## Using the "Carcass" package FOR R and the Korner, Huso, and Erickson estimators to
## estimate fatalities using control and treatment approaches at the BSGF site.
## Please contact Mark Hayes with comments and questions at the email below.
## 
## Mark A. Hayes
## Senior Bat Ecologist
## Normandeau Associates, Inc. 
## Date: 5/7/2018 
## Email: mhayes@normandeau.com
##
## Posted on MAH's GitHub repository: https://github.com/mark-a-hayes 
##
###################################################################################

## Preliminaries

install.packages("carcass")
library("carcass")

remove()
rm(list=ls())

###################################################################################

## Some example code for EstimateN, pp 8 of the carcass pdf:

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
          s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000,
          plot=TRUE)


estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
          s.upper=0.94, d=2, pform="huso", maxn=500, nsim=1000, plot=TRUE)


## Creating a custom estimateN function. 
# Note on plots: The estimateN plot function uses xlim=c(0,50). 
# J. Champion created the following custom 'estimateN_custom' function using the estimateN source code
# MAH added the "plot_limsy" variable to allow for adjustment of the  y-axis limits when plotting the posterior densities.

estimateN_custom <- function (count, p = NA, p.lower = NA, p.upper = NA, f = NA, 
                              f.lower = NA, f.upper = NA, s = NA, s.lower = NA, s.upper = NA, 
                              arrival = "discrete", a = 1, a.lower = 1, a.upper = 1, pform = "korner", 
                              d = 1, n = NA, J = NA, maxn = 1000, nsim = 1000, plot = TRUE, 
                              postdist = FALSE, k = 1, x = c(1:10), 
                              plot_lims=c(0,50), plot_limsy=c(0,1)) #new parameter for plot limits
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
         xlim = plot_lims, ylim = plot_limsy) #use new plot limits parameter
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


## Test that the new function works for a new xlim, and ylim.

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
          s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000,
          plot=TRUE)

estimateN_custom(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64,
                 s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000,
                 plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,1))

################
##
## Checking carcass distances from turbines in control and treatment groups
## to look for evidence of a difference in the underlying distributions.
## Using proportion of carcasses in each 5 m increment from turbine
## and Wilcoxon signed rank test for paired vectors.
##
################

distances = read.table(file.choose(), header=TRUE, sep=",")
attach(distances)
y1 <- proportion_control
y2 <- proportion_treatment


# dependent 2-group Wilcoxon Signed Rank Test 

wilcox.test(y1,y2,paired=TRUE) # where y1 and y2 are numeric

# Conclusion: there is not evidence for differences in our distance vectors.

# Checking result with almost identical vectors

v1 <- c(1:20)
v2 <- c(1,2,4,3,5,6,7,8,10,9,11,12,13,14,15,17,16,18,19,20)

wilcox.test(v1,v2,paired=TRUE) # where v1 and v2 are numeric

# ...and very different vectors

v3 <- c(41:60)

wilcox.test(v1,v3,paired=TRUE) # where v1 and v3 are numeric

# Conclusion: there is not evidence for differences in our distance vectors.

# Some plotting:

par(mfrow=c(1,2))

d1 <- number_control 
d2 <- number_treatment
d1
d2

plot(range,y1, ylim = c(0,0.22))
plot(range,y2, ylim = c(0,0.22))
hist(d1)

par(op)

################
##
## Wilcoxon tests
##
################

daily = read.table(file.choose(), header=TRUE, sep=",")
attach(daily)
summary(daily)

labo_t

wilcox.test(pooled_c, pooled_t, paired=TRUE)
wilcox.test(labo_c, labo_t, paired=TRUE)
wilcox.test(laci_c, laci_t, paired=TRUE)
wilcox.test(lano_c, lano_t, paired=TRUE)
wilcox.test(epfu_c, epfu_t, paired=TRUE)
wilcox.test(mylu_c, mylu_t, paired=TRUE)

###################################################################################

## Using the BSGF carcass results for counts (count), searcher efficiency (f), and carcass
## persistence probability (s) estimates to estimate total fatalities (detected and not detected).

# Fatality estimator parameters are for 95% confidence intervals, and are estimated 
# using the results reported by Gruver (2016). M. Hayes took these estimates, which
# were for 90% CI and converted to 95% CI's (MAH, calculations on 11/27/2017).

# f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890, s.upper=0.950

# NOTE: the estimateN function uses xlim=c(0,50), so here I use a costume version of the function
# so that the limits of the x-axis can be expanded and changed as needed. Same with y-axis. 


## Pooled (all species) for control turbine fatalities, with total count = 187, new count = 138:
  
# The Erickson estimator, pooled for 10 turbines

par(mfrow=c(2,3))  
  
estimateN_custom(count=187, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
            s.upper=0.950, d=1, pform="erickson", n = 78, maxn=600, nsim=10000,
            plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "pooled_c_e.tiff", height=4, width=6.5, units='in', res=300)


# The Huso estimator

estimateN_custom(count=138, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n = 78, maxn=600, nsim=10000,
          plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "pooled_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=187, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=600, nsim=10000,
          plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "p_c_K.tiff", height=4, width=6.5, units='in', res=300)


# Pooled treatment fatalities, with total count = 29, new count = 22:

# The Erickson estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "p_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=22, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "p_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,500), plot_limsy = c(0,0.05))

dev.print(tiff, "fig2.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################

## LACI control fatalities, with total count = 48, new count = 40:

# The Erickson estimator

estimateN_custom(count=48, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=40, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=48, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LACI treatment fatalities, with total count = 9, new count = 7:

# The Erickson estimator

estimateN_custom(count=9, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimor

estimateN_custom(count=9, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "laci_t_k.tiff", height=4, width=6.5, units='in', res=300)


###################################################################################

## LABO control fatalities, with total count = 37, new count = 21:

# The Erickson estimator

estimateN_custom(count=37, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=21, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=37, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LABO treatment fatalities, with total count = 6, new count = 5:

# The Erickson estimator

estimateN_custom(count=6, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=5, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=6, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "labo_t_k.tiff", height=4, width=6.5, units='in', res=300)

###################################################################################

## LANO control fatalities, with total count = 45, new count = 36:

# The Erickson estimator

estimateN_custom(count=45, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=36, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=45, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_c_k.tiff", height=4, width=6.5, units='in', res=300)

## LANO treatment fatalities, with total count = 4, new count = 3:

# The Erickson estimator

estimateN_custom(count=4, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=3, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=4, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "lano_t_k.tiff", height=4, width=6.5, units='in', res=300)

###################################################################################

## EPFU control fatalities, with total count = 27, new count = 18:

# The Erickson estimator

estimateN_custom(count=27, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=18, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=27, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_c_k.tiff", height=4, width=6.5, units='in', res=300)

## EPFU treatment fatalities, with total count = 7, new count = 5:

# The Erickson estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=5, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=7, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "epfu_t_k.tiff", height=4, width=6.5, units='in', res=300)

###################################################################################

## MYLU control fatalities, with total count = 29, new count = 23:

# The Erickson estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=23, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="huso", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_h.tiff", height=4, width=6.5, units='in', res=300)

# The Korner estimator

estimateN_custom(count=29, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="korner", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_c_k.tiff", height=4, width=6.5, units='in', res=300)

## MYLU treatment fatalities, with total count = 3, new count = 2:

# The Erickson estimator

estimateN_custom(count=3, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
          s.upper=0.950, d=1, pform="erickson", n=78, maxn=500, nsim=10000,
          plot=TRUE, plot_lims=c(0,150))

dev.print(tiff, "mylu_t_e.tiff", height=4, width=6.5, units='in', res=300)

# The Huso estimator

estimateN_custom(count=2, f=0.600, f.lower=0.421, f.upper=0.779, s=0.920, s.lower=0.890,
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
