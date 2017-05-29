
##################################################################################
# an R script to fit the parameters of a structured SEIR model to norovirus data
# data by minimizing the Least Squares statistic
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Feb 6, 2017
#
# Copyright Sherry Towers, 2017
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
#
# https://www.ncbi.nlm.nih.gov/pubmed/20508526
#
##################################################################################

##################################################################################
# par(pch=20) sets a solid round dot style for the plots
# The chron package contains utilities to calculate dates
##################################################################################
rm(list = ls(all = TRUE))  # resets R to fresh
require("sfsmisc")
require("deSolve")
par(pch=20)  
set.seed(362552)
ltwo_cruise = 1
lfixed_gamma_kappa = 0 # if true, then don't sweep gamma and kappa
lpoisson = 0 # if 1 then do Poisson likelihood fit.  Otherwise do Negative Binomial likelihood fit
source("norovirus_utils.R")

##################################################################################
##################################################################################
##################################################################################
# read in the data, and select cruises #1 and #2
##################################################################################
adat = read.table("cruise_norovirus_outbreak_data.csv",header=T,as.is=T,sep=",")

if (ltwo_cruise){
  bdat = subset(adat,cruise_number<=2)
}else{
  bdat = subset(adat,cruise_number<=1)
}
incidence_observed1 = bdat$npass
incidence_observed2 = bdat$ncrew
times_of_observed = bdat$day

time_binning = min(diff(times_of_observed))

##################################################################################
##################################################################################
# now, set up the iterations of the Monte Carlo method
# At each iteration, we will randomly sample a hypothesis for the
# reproduction number, R0, and the time-of-introduction of the virus to the
# population, t0.
#
# With these hypotheses, we will solve for the model predicted incidence
# and calculate the negative log likelihood statistic comparing this to the observed
# incidence.  We store the hypotheses for R0 and t0 in the vectors vR0_direct and vt0,
# and the resulting goodness of fit statistic in the vector vnegloglike
#
# At each iteration, we'll check if the predicted incidence is the best fit
# so far, and if so, we'll store that in vbest_negloglike_fit_incidence_prediction
#
# best_negloglike_so_far keeps track of the best-fit GoF statistic so far obtained
# in the iterations.
##################################################################################
vkappa = numeric(0)
vgamma = numeric(0)
valpha = numeric(0)
vr     = numeric(0)
vq     = numeric(0)
vxi    = numeric(0)
vR0_direct    = numeric(0)
vR0_env   = numeric(0)
vR0_tot   = numeric(0)

vwkappa = numeric(0)
vwgamma = numeric(0)
vwalpha = numeric(0)
vwr     = numeric(0)
vwq     = numeric(0)
vwR0_direct    = numeric(0)
vwR0_env   = numeric(0)
vwxi    = numeric(0)


veta = numeric(0)
vnegloglike = numeric(0) 
vfrac_conf1 = numeric(0)
vfrac_conf2 = numeric(0)
vfcrewcrew = numeric(0)
vratio = numeric(0)
vt0 = numeric(0)

vB11 = numeric(0)
vB12 = numeric(0)
vB21 = numeric(0)
vB22 = numeric(0)

best_negloglike_so_far = 1e10 
vbest_negloglike_fit_incidence_prediction = rep(0,length(incidence_observed1))

R0_direct_best = 6.7
q_best = 0.96
r_best = 0.07
alpha_best = 0.36

niter = 1000000
for (iter in 1:niter){

   ###############################################################################
   # This process is computationally intensive, so once in a while during the
   # iterations it is nice to inform user the script is doing something, 
   # and not hung
   ###############################################################################
   if (iter%%100==0){
     cat("Doing iteration ",iter," out of ",niter,"\n")  
   }

   ###############################################################################
   ################################################################################
   # set up the model parameters
   # npop is approximately the population of IL IN MI MN OH WI (CDC region 5)
   # I_0 is the initial number infected
   # R_0 is the inital number recovered and immune
   # S_0 is the susceptibles
   #
   # 1/gamma is the average recovery period of the disease
   # R0_direct      is the reproduction number
   # t0      is the time-of-introduction of the disease to the population, measured
   #         in days from Jan 1, 2007
   #
   # For the SIR model, R0=eta/gamma, thus given our hypotheses for gamma and R0,
   # we calculate eta as eta=R0*gamma (note that eta and gamma appear in the model
   # equations, not gamma and R0, which is why we need to calculate eta given R0
   # and gamma).
   ###############################################################################
   npass = 2318
   ncrew = 999
   npop = npass+ncrew
   f = c(npass,ncrew)/npop
   N = f*npop
   nsub_population = length(f)

   ###############################################################################
   # https://www.ncbi.nlm.nih.gov/pubmed/20508526 
   #     zelner2010infections
   #     infectious period 1.17 [1,1.88] days 95% CI
   #     shape parameter for infectious period was found to be 1 with 95% CI [1-2]
   #     based on likelihood analysis of household outbreak data
   #     mean incubation period estimated to be 1.7 days, with SE 0.048
   #     incubation period used a Gamma distribution, with shape parameter 3.73
   #
   # https://www.researchgate.net/profile/Rowena_Bull/publication/5434332_Norovirus_Excretion_in_an_Aged-Care_Setting/links/004635268419bbd3f3000000.pdf 
   #     tu2009norovirus
   #     half life of viral shedding in the elderly after onset of symptoms was T1/2=2.5 days, meaning 1/gamma = 3.6 days
   #     assumed Exponential distribution for viral shedding profile
   #     These were elderly people.
   #
   #  https://academic.oup.com/cid/article/33/5/622/465934/Clinical-Spectrum-and-Transmission-Characteristics
   #     outbreak in childcare setting... looked at incubation time, and serial interval
   #     figure 3, incubation time is clearly Gamma distributed
   #     v = c(12,24,36,48,60,72)
   #     w = c(3,13,58,32,18,2) 
   #     mean = 41.2 hours with se 1.27 hours
   #     sd of the distribution is 14.25 hours
   #     thus the parameters of the Gamma distribution are theta = 4.93, k=8.4
   #     In Figure 4, the serial interval also appears to be Gamma distributed
   #     v = c(12,24,36,48,60,72,84,96)
   #     w = c(1,4,1,7,7,4,4,3)
   #     mean 58.45 hours with se 4.42 hours
   #     sd of the distribution is 24.65 hours
   #     thus the parameters of the Gamma distribution are theta = 10.4, k=5.6
   #
   #  http://www2.math.su.se/matstat/reports/seriea/2005/rep14/report.pdf
   #     svennson2007note
   #     if X is incubation period with either exponential or uniform distribution with mean_x
   #     and Y is infectious period with either exponential or uniform distribution with mean_y
   #         then the serial interval is mean_x+mean_y (eqn 5.10)
   #     on the other hand, if Y is Gamma distributed with mean_x and shape parameter alpha, then
   #         the serial interval is mean_x+mean_y (alpha+1)/(2*alpha)  (eqn 5.9)
   ###############################################################################
   gamma_cen = 1.17
   egammapl = (1.88-gamma_cen)/1.96
   egammami = (gamma_cen-1)/1.96

   kappa_cen = 1.15
   ekappapl = 0.10/1.96/2
   ekappami = 0.10/1.96/2


   if (!lfixed_gamma_kappa){
     kmin = max(0,kappa_cen-4*ekappami)
     kmax = kappa_cen+4*ekappapl
     gmin = max(0,gamma_cen-4*egammami)
     gmax = gamma_cen+4*egammapl
     #kappa = 1/runif(1,kmin,kmax)
     #gamma = 1/runif(1,gmin,gmax)

     kappa = 1/rnorm(1,kappa_cen,ekappapl)
     wkappa = dnorm(1/kappa,kappa_cen,ekappapl)

     atest = runif(1,0,1)
     mu = 1.2
     if (atest>0.5){
       sigma =0.35
       gamma = mu+abs(rnorm(1,0,sigma))
       wgamma = 0.5*dnorm(gamma-mu,0,sigma)
       gamma = 1/gamma
     }else{
       sigma = 0.15
       gamma = mu-abs(rnorm(1,0,sigma))
       while(gamma<0) gamma = mu-abs(rnorm(1,0,sigma))
       wgamma = 0.5*dnorm(mu-gamma,0,sigma)
       gamma = 1/gamma
     }

   }else{
     kappa = 1/kappa_cen
     gamma = 1/gamma_cen
   }

   ###############################################################################
   # randomly sample R0_direct and t0 uniformly
   ###############################################################################
   t0 = min(times_of_observed)-1

   atest = runif(1,0,1)
   mu = 0.97
   if (atest>0.5){
     sigma = 0.01
     q = mu+abs(rnorm(1,0,sigma))
     while(q>1) q = mu+abs(rnorm(1,0,sigma))
     wq = 0.5*dnorm(q-mu,0,sigma)
   }else{
     sigma = 0.030
     q = mu-abs(rnorm(1,0,sigma))
     while(q<0) q = mu-abs(rnorm(1,0,sigma))
     wq = 0.5*dnorm(mu-q,0,sigma)
   }

   #R0_env = 0
   #R0_env = runif(1,0.0,3.0)
   #wR0_env = 1

   atest = runif(1,0,1)
   mu = 0.90
   if (atest>0.5){
     sigma = 0.55
     R0_env = mu+abs(rnorm(1,0,sigma))
     wR0_env = 0.5*dnorm(R0_env-mu,0,sigma)
   }else{
     sigma = mu/3
     R0_env = mu-abs(rnorm(1,0,sigma))
     while(R0_env<0) R0_env = mu-abs(rnorm(1,0,sigma))
     wR0_env = 0.5*dnorm(mu-R0_env,0,sigma)
   }
   

   mu = 0.607
   sigma = 1.5*0.117
   xi = rnorm(1,mu,sigma)
   wxi = dnorm(xi,mu,sigma)
   #xi = runif(1,0,1.5)
   #wxi = 1

   atest = runif(1,0,1)
   mu = 11
   if (atest>0.5){
     sigma = 2.0
     R0_direct = mu+abs(rnorm(1,0,sigma))
     wR0_direct = 0.5*dnorm(R0_direct-mu,0,sigma)
   }else{
     sigma = 1.5
     R0_direct = mu-abs(rnorm(1,0,sigma))
     while(R0_direct<0) R0_direct = mu-abs(rnorm(1,0,sigma))
     wR0_direct = 0.5*dnorm(mu-R0_direct,0,sigma)
   }



   atest = runif(1,0,1)
   mu = 0.24
   if (atest>0.5){
     sigma = 0.2
     alpha = mu+abs(rnorm(1,0,sigma))
     walpha = 0.5*dnorm(alpha-mu,0,sigma)
   }else{
     sigma = 0.075
     alpha = mu-abs(rnorm(1,0,sigma))
     while(alpha<0) alpha = mu-abs(rnorm(1,0,sigma))
     walpha = 0.5*dnorm(mu-alpha,0,sigma)
   }

   lcalculate_transmission_probability = 1
   C = matrix(0,nrow=nsub_population,ncol=nsub_population)
   C[1,1] = q
   C[1,2] = 1-q
   C[2,1] = C[1,2]*f[1]/f[2]
   r = runif(1,0,0.20)
   if (!lpoisson) r = runif(1,0,0.35)
   r = abs(rnorm(1,0,0.15))
   wr = 0.5*dnorm(r,0,0.15)
   C[2,2] = r

   if (lcalculate_transmission_probability==1){
      M = C
      M[1,1] = C[1,1]*f[1]/f[1]
      M[1,2] = C[1,2]*f[1]/f[2]
      M[2,1] = C[2,1]*f[2]/f[1]
      M[2,2] = C[2,2]*f[2]/f[2]
      eig = eigen(M)
      # reverse engineer eta from the R0_direct and gamma
      eta = R0_direct*gamma/max(Re(eig$values))
      eta = eta
   }else{
      eta = 0.05
   }

   B = C
   B[1,1] = eta*C[1,1]
   B[1,2] = eta*C[1,2]
   B[2,1] = eta*C[2,1]
   B[2,2] = eta*C[2,2]

   
   

   ###############################################################################
   # Now solve the system of differential equations numerically with lsoda in the 
   # deSolve package.  Put the results in solved_model
   ###############################################################################
   tmin = max(bdat$day[bdat$cruise_number==1])
   vt = seq(0,tmin)
   I = c(1,0)
   E = c(0,0)
   R = c(0,0)
   S = npop*f-E-I-R
   newI = c(0,0)
   W = 0
   inits = c(S=S,E=E,I=I,R=R,newI=newI,W=W)
   vparameters = c(gamma=gamma,eta=eta,kappa=kappa,xi=xi,beta_w=R0_env*gamma,B=B,fquar=0)
   solved_model = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters)) 

   npass = 2400
   n = nrow(solved_model)
   I = c(0,solved_model$I2[n])
   E = c(0,solved_model$E2[n])
   R = c(0,solved_model$R2[n])
   S = c(npass,solved_model$S2[n])
   newI = c(0,solved_model$newI2[n])
   W = (solved_model$W[n])
   vt = seq(tmin,15)
   inits = c(S=S,E=E,I=I,R=R,newI=newI,W=W)
   vparameters = c(gamma=gamma,eta=eta,kappa=kappa,xi=xi,beta_w=R0_env*gamma,B=B,fquar=0)
   solved_modelb = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters)) 

   ###############################################################################
   # The incidence over some time step is the difference in S over that time
   # step (in a model with no births or deaths, immigration or emigration, or 
   # recurring susceptibility).
   #
   # The time step are derived from the dates at the end of the week of each
   # data point (vtime_data)
   #
   # solved_model$time%in%vtime_data returns the indices of the elements of simodel$time that
   # are also found in vtime_data
   ###############################################################################
   tmin_data = min(times_of_observed)-time_binning
   tmax_data = max(times_of_observed)
   vtime_data = seq(tmin_data,tmax_data,time_binning)
   
   susceptible_predicted1 = solved_model$newI1[solved_model$time%in%vtime_data]  
   incidence_predicted1 = +diff(susceptible_predicted1)
   susceptible_predicted2 = solved_model$newI2[solved_model$time%in%vtime_data]  
   incidence_predicted2 = +diff(susceptible_predicted2)

   if (ltwo_cruise){
     susceptible_predicted1 = solved_modelb$newI1[solved_modelb$time%in%vtime_data]  
     incidence_predicted1 = c(incidence_predicted1,+diff(susceptible_predicted1))
     susceptible_predicted2 = solved_modelb$newI2[solved_modelb$time%in%vtime_data]  
     incidence_predicted2 = c(incidence_predicted2,+diff(susceptible_predicted2))
   }
   
   ###############################################################################
   # from the model estimate of the incidence and the data, we can 
   # estimate the fraction of cases that were confirmed
   ###############################################################################
   frac_confirmed1 = sum(incidence_observed1)/sum(incidence_predicted1)
   frac_confirmed2 = sum(incidence_observed2)/sum(incidence_predicted2)
   frac_confirmed = sum(incidence_observed1+incidence_observed2)/sum(incidence_predicted1+incidence_predicted2)
   #incidence_predicted1 = incidence_predicted1*frac_confirmed
   #incidence_predicted2 = incidence_predicted2*frac_confirmed
   incidence_predicted1 = incidence_predicted1*frac_confirmed1
   incidence_predicted2 = incidence_predicted2*frac_confirmed2


   ###############################################################################
   # now calculate the negative log likelihood
   # statistic that compares the data to this model calculated
   # under a particular hypothesis of R0_direct and t0
   ###############################################################################

   if (length(incidence_predicted1)==length(incidence_observed1)
      &!is.na(sum(incidence_predicted1))){

       ########################################################################### 
       # calculate the goodness of fit statistic
       ########################################################################### 
       if (lpoisson){
         negloglike = sum(-(incidence_observed1*log(incidence_predicted1)-incidence_predicted1))
         negloglike = negloglike + sum(-(incidence_observed2*log(incidence_predicted2)-incidence_predicted2))
       }else{
         negloglike = sum(my_lnbinom(incidence_observed1,incidence_predicted1,alpha))
         negloglike = negloglike+sum(my_lnbinom(incidence_observed2,incidence_predicted2,alpha))
       }
       if (1/gamma>gamma_cen){
         negloglike = negloglike + 0.5*(1/gamma-gamma_cen)^2/(egammapl)^2
       }else{
         negloglike = negloglike + 0.5*(1/gamma-gamma_cen)^2/(egammami)^2
       }
       if (1/kappa>kappa_cen){
         negloglike = negloglike + 0.5*(1/kappa-kappa_cen)^2/(ekappapl)^2
       }else{
         negloglike = negloglike + 0.5*(1/kappa-kappa_cen)^2/(ekappami)^2
       }
       negloglike = negloglike + 0.5*(xi-0.607)^2/(0.117)^2

       #negloglike = negloglike - log(dgamma(1/kappa+1/gamma,shape=3.35,scale=1.09)) 

       vB11 = append(vB11,B[1,1])
       vB12 = append(vB12,B[1,2])
       vB21 = append(vB21,B[2,1])
       vB22 = append(vB22,B[2,2])
       beta_W = R0_env*gamma
       a = calculate_R0_matrix(B,gamma,beta_W,N)
       vR0_direct = append(vR0_direct,R0_direct)
       vR0_env = append(vR0_env,R0_env)
       vR0_tot = append(vR0_tot,a$R0_tot)
       vxi = append(vxi,xi)
       vt0 = append(vt0,t0)
       vr = append(vr,r)
       veta = append(veta,eta)
       vq = append(vq,q)
       vnegloglike = append(vnegloglike,negloglike)
       valpha = append(valpha,alpha)
       vgamma = append(vgamma,gamma)
       vkappa = append(vkappa,kappa)
       vfrac_conf1 = c(vfrac_conf1,frac_confirmed1)
       vfrac_conf2 = c(vfrac_conf2,frac_confirmed2)
       vfcrewcrew = c(vfcrewcrew,B[2,2]/(B[2,1]+B[2,2]))
       vratio = c(vratio,(B[2,1]+B[2,2])/(B[1,1]+B[1,2]))
 
       vwkappa = c(vwkappa,wkappa)
       vwgamma = c(vwgamma,wgamma)
       vwR0_direct    = c(vwR0_direct   ,wR0_direct   )
       vwalpha = c(vwalpha,walpha)
       vwr     = c(vwr    ,wr    )
       vwq     = c(vwq    ,wq    )
       vwR0_env   = c(vwR0_env  ,wR0_env  )
       vwxi    = c(vwxi   ,wxi   )

       if (negloglike<best_negloglike_so_far){
         best_negloglike_so_far = negloglike
         vbest_negloglike_fit_incidence_prediction1 = incidence_predicted1
         vbest_negloglike_fit_incidence_prediction2 = incidence_predicted2
         frac_confirmed1_best = frac_confirmed1
         frac_confirmed2_best = frac_confirmed2
         solved_model_best = solved_model
         solved_modelb_best = solved_modelb
         alpha_best = alpha
         R0_direct_best = R0_direct
         R0_env_best = R0_env
         xi_best = xi
         t0_best = t0
         r_best = r
         eta_best = eta
         q_best = q
         B_best = B
         gamma_best = gamma
         kappa_best = kappa
         fcrewcrew_best = B[2,2]/(B[2,1]+B[2,2])
         ratio_best = (B[2,1]+B[2,2])/(B[1,1]+B[1,2])
         cat("The best value of R0_direct so far is:",R0_direct_best,"\n")
         cat("The best value of q so far is:",q_best,"\n")
         cat("The best value of r so far is:",r_best,"\n")
         cat("The best value of eta   so far is:",eta_best,"\n")
         cat("The best value of frac_confirmed1 is:",frac_confirmed1,"\n")
         cat("The best value of frac_confirmed2 is:",frac_confirmed2,"\n")
         print(B_best)
       }

       ######################################################################## 
       # plot the best-fit results every once in a while
       # cex is the point size
       ######################################################################## 
       if (iter%%100==0){
         vdat = data.frame(R0_direct=vR0_direct,eta=veta,t0=vt0,r=vr,q=vq,like=vnegloglike,alpha=valpha,kappa=vkappa,gamma=vgamma,frac_conf1=vfrac_conf1,frac_conf2=vfrac_conf2,fcrewcrew=vfcrewcrew,ratio=vratio,xi=vxi,beta_W=vR0_env*vgamma,R0_tot=vR0_tot,R0_env=vR0_env,wkappa=vwkappa,wgamma=vwgamma,wR0_direct=vwR0_direct,walpha=vwalpha,wr=vwr,wq=vwq,wR0_env=vwR0_env,wxi=vwxi,B11=vB11,B12=vB12,B21=vB21,B22=vB22)
         write.table(vdat,"fitting_results_z.txt",row.names=F)
         weight = dnorm(sqrt(2*(vnegloglike-min(vnegloglike))))
         weight = weight/vwkappa/vwgamma/vwR0_direct/vwalpha/vwr/vwq/vwR0_env/vwxi

         R0_direct_best_weight = weighted.mean(vR0_direct,weight)
         R0_env_best_weight = weighted.mean(vR0_env,weight)
         q_best_weight = weighted.mean(vq,weight)
         r_best_weight = weighted.mean(vr,weight)
         alpha_best_weight = weighted.mean(valpha,weight)
         xi_best_weight = weighted.mean(vxi,weight)
         if (lpoisson){
           A = cov.wt(as.matrix(cbind(vR0_direct,vq,vr,vR0_env)),weight,cor=T)
         }else{
           A = cov.wt(as.matrix(cbind(vR0_direct,vq,vr,vR0_env,valpha)),weight,cor=T)
         }
         Vcov_weight_method = A$cov
         eR0_direct_best_weight = sqrt(Vcov_weight_method[1,1])
         eq_best_weight = sqrt(Vcov_weight_method[2,2])
         er_best_weight = sqrt(Vcov_weight_method[3,3])
         eR0_env_best_weight = sqrt(Vcov_weight_method[4,4])
         exi_best_weight = sqrt(Vcov_weight_method[5,5])
         n = ncol(Vcov_weight_method)
         if (!lpoisson) ealpha_best_weight = sqrt(Vcov_weight_method[n,n])

         text_main = paste("Confirmed norovirus cases:\n result of",iter,"Monte Carlo fit iterations")
         if (lpoisson){
           text_main = paste(text_main,"\n Poisson Likelihood Fit")
         }else{
           text_main = paste(text_main,"\n Negative Binomial Likelihood Fit")
         }

         mult.fig(12,main=text_main,oma=c(1,2,4,1),mar=c(3,3,4,1))

         num_points_to_show = 250
         adel = 0.5*3^2
         ymax = max(vnegloglike)
         if (length(vnegloglike)>num_points_to_show) ymax = sort(vnegloglike)[num_points_to_show]
         if (ymax<=(min(vnegloglike)+adel)) ymax = min(vnegloglike)+adel
         l = which(vnegloglike<=ymax)

         lmin = which.min(vnegloglike)

         lgoodb = which(vnegloglike<=(min(vnegloglike)+0.5))
         fminplushalf = min(vnegloglike[lgoodb])+0.5
         fminplus196 = min(vnegloglike[lgoodb])+1.96^2/2

         amin = min(vR0_direct[lgoodb])
         amax = max(vR0_direct[lgoodb])
         aminb = R0_direct_best_weight-1.96*eR0_direct_best_weight
         amaxb = R0_direct_best_weight+1.96*eR0_direct_best_weight
         del = (amax-amin)/2
         plot(vR0_direct[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="R0_direct hypothesis"
             ,main=paste("Best-fit R0_direct so far fmin+1/2: ",round(R0_direct_best,2),"+/-",round(del,2)," [",round(amin,2),",",round(amax,2),"], \n"
                        ,"Best-fit R0_direct so far weighted mean: ",round(R0_direct_best_weight,2),"+/-",round(eR0_direct_best_weight,2)," [",round(aminb,2),",",round(amaxb,2),"]",sep="")
             ,col.main=3,cex.main=0.75)
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)
         points(vR0_direct[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }

         amin = min(vR0_env[lgoodb])
         amax = max(vR0_env[lgoodb])
         aminb = R0_env_best_weight-1.96*eR0_env_best_weight
         amaxb = R0_env_best_weight+1.96*eR0_env_best_weight
         del = (amax-amin)/2
         plot(vR0_env[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="R0_env hypothesis"
             ,main=paste("Best-fit R0_env so far fmin+1/2: ",round(R0_env_best,2),"+/-",round(del,2)," [",round(amin,2),",",round(amax,2),"], \n"
                        ,"Best-fit R0_env so far weighted mean: ",round(R0_env_best_weight,2),"+/-",round(eR0_env_best_weight,2)," [",round(aminb,2),",",round(amaxb,2),"]",sep="")
             ,col.main=3,cex.main=0.75)
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)
         points(vR0_env[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }

       

         amin = min(vq[lgoodb])
         amax = max(vq[lgoodb])
         aminb = q_best_weight-1.96*eq_best_weight
         amaxb = q_best_weight+1.96*eq_best_weight
         del = (amax-amin)/2
         plot(vq[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="q hypothesis"
             ,main=paste("Best-fit q so far fmin+1/2: ",round(q_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n"
                        ,"Best-fit q so far weighted mean: ",round(q_best_weight,3),"+/-",round(eq_best_weight,3)," [",round(aminb,3),",",round(amaxb,3),"]",sep="")
             ,col.main=3,cex.main=0.75)
         points(vq[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)

         amin = min(vr[lgoodb])
         amax = max(vr[lgoodb])
         aminb = r_best_weight-1.96*er_best_weight
         amaxb = r_best_weight+1.96*er_best_weight
         del = (amax-amin)/2
         plot(vr[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="r hypothesis"
             ,main=paste("Best-fit r so far fmin+1/2: ",round(r_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n"
                        ,"Best-fit r so far weighted mean: ",round(r_best_weight,3),"+/-",round(er_best_weight,3)," [",round(aminb,3),",",round(amaxb,3),"]",sep="")
             ,col.main=3,cex.main=0.75)
         points(vr[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)

         amin = min(vratio[lgoodb])
         amax = max(vratio[lgoodb])
         del = (amax-amin)/2
         plot(vratio[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="ratio hypothesis"
             ,main=paste("Best-fit ratio so far fmin+1/2: ",round(ratio_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n" ,sep="")
             ,col.main=3,cex.main=0.75)
         points(vratio[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)

         amin = min(1/vkappa[lgoodb])
         amax = max(1/vkappa[lgoodb])
         del = (amax-amin)/2
         plot(1/vkappa[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="1/kappa hypothesis"
             ,main=paste("Best-fit 1/kappa so far fmin+1/2: ",round(1/kappa_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n" ,sep="")
             ,col.main=3,cex.main=0.75)
         points(1/vkappa[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)



         amin = min(1/vgamma[lgoodb])
         amax = max(1/vgamma[lgoodb])
         del = (amax-amin)/2
         plot(1/vgamma[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="1/gamma hypothesis"
             ,main=paste("Best-fit 1/gamma so far fmin+1/2: ",round(1/gamma_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n" ,sep="")
             ,col.main=3,cex.main=0.75)
         points(1/vgamma[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)

         amin = min(vxi[lgoodb])
         amax = max(vxi[lgoodb])
         aminb = xi_best_weight-1.96*exi_best_weight
         amaxb = xi_best_weight+1.96*exi_best_weight
         del = (amax-amin)/2
         plot(vxi[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="xi hypothesis"
             ,main=paste("Best-fit xi so far fmin+1/2: ",round(xi_best,2),"+/-",round(del,2)," [",round(amin,2),",",round(amax,2),"], \n"
                        ,"Best-fit xi so far weighted mean: ",round(xi_best_weight,2),"+/-",round(exi_best_weight,2)," [",round(aminb,2),",",round(amaxb,2),"]",sep="")
             ,col.main=3,cex.main=0.75)
         points(vxi[lmin],vnegloglike[lmin],col=3,cex=2)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)




         if (!lpoisson){
         amin = min(valpha[lgoodb])
         amax = max(valpha[lgoodb])
         del = (amax-amin)/2
         aminb = alpha_best_weight-1.96*ealpha_best_weight
         amaxb = alpha_best_weight+1.96*ealpha_best_weight
         plot(valpha[l]
             ,vnegloglike[l]
             ,ylab="negative log likelihood"
             ,xlab="alpha hypothesis"
             ,main=paste("Best-fit alpha so far fmin+1/2: ",round(alpha_best,3),"+/-",round(del,3)," [",round(amin,3),",",round(amax,3),"], \n"
                        ,"Best-fit alpha so far weighted mean: ",round(alpha_best_weight,3),"+/-",round(ealpha_best_weight,3)," [",round(aminb,3),",",round(amaxb,3),"]",sep="")
             ,col.main=3,cex.main=0.75)
         points(valpha[lmin],vnegloglike[lmin],col=3,cex=2)
         lines(c(-1e6,1e6),c(vnegloglike[lmin],vnegloglike[lmin]),col=4)
         lines(c(-1e6,1e6),c(fminplushalf,fminplushalf),col=4,lty=3)
         lines(c(-1e6,1e6),c(fminplus196,fminplus196),col=2,lty=3)
         if (amax>amin){
         arrows(amin,fminplushalf,amax,fminplushalf,col=3,length=0.2,lwd=2,code=3)
         }
         }

         ymax = max(c(incidence_observed1,vbest_negloglike_fit_incidence_prediction1))
         plot(times_of_observed
             ,incidence_observed1
             ,ylim=c(0,1.2*ymax)
             ,xlab="Time, in days"
             ,ylab="Incidence"
             ,cex=2,main="Identified cases in passengers")
         lines(times_of_observed
              ,vbest_negloglike_fit_incidence_prediction1
              ,col=2
              ,lwd=5) 

         ymax = max(c(incidence_observed2,vbest_negloglike_fit_incidence_prediction2))
         plot(times_of_observed
             ,incidence_observed2
             ,ylim=c(0,1.2*ymax)
             ,xlab="Time, in days"
             ,ylab="Incidence"
             ,cex=2,main="Identified cases in crew")
         lines(times_of_observed
              ,vbest_negloglike_fit_incidence_prediction2
              ,col=2
              ,lwd=5) 

       }
     
   } # end check that the predicted incidence vector is the same length as the observed and doesn't contain NA's
} # end loop over the Monte Carlo iterations
   

