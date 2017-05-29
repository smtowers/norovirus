##################################################################################
##################################################################################
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
rm(list = ls(all = TRUE))
require("sfsmisc")
require("deSolve")
par(pch=20)  
ltwo_cruise = 1
source("norovirus_utils.R")

##################################################################################
##################################################################################
##################################################################################
if (1){
  vdat  = read.table("fitting_results_y.txt",header=T,as.is=T)
  vdat = subset(vdat,like<=(min(vdat$like)+18))
  write.table(vdat,"fitting_results_within_6_sigma.txt",row.names=F)
}else{
  vdat  = read.table("fitting_results_within_6_sigma.txt",header=T,as.is=T)
}

##################################################################################
##################################################################################
##################################################################################
if (0){
  vdat  = read.table("fitting_results_d.txt",header=T,as.is=T)
  vdatb = read.table("fitting_results_e.txt",header=T,as.is=T)
  vdat = rbind(vdat,vdatb)
  vdat = subset(vdat,like<=(min(vdat$like)+18))
  write.table(vdat,"fitting_results_within_6_sigma.txt",row.names=F)
}else{
  vdat  = read.table("fitting_results_within_6_sigma.txt",header=T,as.is=T)
}

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
##################################################################################
vparameters=fill_vparameters(vdat)
a = solve_model(bdat,vparameters)
incidence_predicted_raw1 = a$incidence_predicted_raw1
incidence_predicted_raw2 = a$incidence_predicted_raw2


vrel1 = numeric(0)
vrel2 = numeric(0)
vreltot = numeric(0)
wscale = seq(0,1,0.05)
n = length(wscale)
vinterv = c(rep("fquar",n),rep("env",n),rep("contact",n))
vscale = rep(wscale,3)
f_symptomatic = 0.5
for (i in 1:length(vinterv)){
  vparametersb = as.list(vparameters)
  scale = vscale[i]
  if (vinterv[i]=="fquar") vparametersb$fquar = scale*f_symptomatic
  if (vinterv[i]=="env"){
    vparametersb$R0_env = scale*vparametersb$R0_env
    if (scale>0){
       vparametersb$xi = vparametersb$xi/scale
    }else{
       vparametersb$xi = 1e9
    }
  }
  if (vinterv[i]=="contact"){
    vparametersb$R0_direct = scale*vparametersb$R0_direct
    vparametersb$R0_env = scale*vparametersb$R0_env
  }

  vparametersb = unlist(vparametersb)
  b = solve_model(bdat,vparametersb)
  rel1 = sum(b$incidence_predicted_raw1)/sum(incidence_predicted_raw1)
  rel2 = sum(b$incidence_predicted_raw2)/sum(incidence_predicted_raw2)
  reltot = sum(b$incidence_predicted_raw1+b$incidence_predicted_raw2)/sum(incidence_predicted_raw1+incidence_predicted_raw2)
  vrel1 = c(vrel1,rel1)
  vrel2 = c(vrel2,rel2)
  vreltot = c(vreltot,reltot)
}

for (iter in 1:3){
  if (iter==1) jpeg('interventions.jpg')
  if (iter==2){
    setEPS()
    postscript("interventions.eps")
  }


mult.fig(4,main="Assessment of relative efficacy of potential intervention strategies\n (relative to baseline, with no interventions)",oma=c(1,1,4,1))

ymax = 1.05
plot(wscale,vreltot[1:n],ylim=c(0,ymax),xlab="",ylab="Relative reduction in outbreak final size",type="l",col=1,lwd=5,xaxs="i",yaxs="i",main="Quarantine")
mtext(side=1,"Fraction of symptomatic cases quarantined\n (assuming 50% of cases are symptomatic)",cex=0.8,line=2.5)
lines(wscale,vrel1[1:n],col=2,lwd=5)
lines(wscale,vrel2[1:n],col=4,lwd=5)
grid()
legend("bottomleft",legend=c("All","Passengers","Crew"),col=c(1,2,4),bty="n",lwd=5)


plot(1-wscale,vreltot[(n+1):(2*n)],ylim=c(0,ymax),xlab="",ylab="Relative reduction in outbreak final size",type="l",col=1,lwd=5,xaxs="i",yaxs="i",main="Environmental cleaning")
mtext(side=1,line=2.5,"Relative reduction in environmental transmission\n due to cleaning of environment",cex=0.8)
lines(1-wscale,vrel1[(n+1):(2*n)],col=2,lwd=5)
lines(1-wscale,vrel2[(n+1):(2*n)],col=4,lwd=5)
grid()
par(new=T)
plot(c(0,1),c(0,1),col=0,axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
text(0.99,0.01,"Complete cleaning",adj=c(0,1),srt=90,cex=0.8,font=2)
text(-0.02,0.01,"No extra cleaning",adj=c(0,1),srt=90,cex=0.8,font=2)

plot(1-wscale,vreltot[(2*n+1):(3*n)],ylim=c(0,ymax),xlab="",ylab="Relative reduction in outbreak final size",type="l",col=1,lwd=5,xaxs="i",yaxs="i",main="Handwashing")
mtext(side=1,line=2.5,"Relative reduction in transmission\n due to hand washing",cex=0.8)
lines(1-wscale,vrel1[(2*n+1):(3*n)],col=2,lwd=5)
lines(1-wscale,vrel2[(2*n+1):(3*n)],col=4,lwd=5)
grid()
par(new=T)
plot(c(0,1),c(0,1),col=0,axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
text(0.99,0.01,"Complete handwashing",adj=c(0,1),srt=90,cex=0.8,font=2)
text(-0.02,0.01,"No extra handwashing",adj=c(0,1),srt=90,cex=0.8,font=2)

  if (iter<=2) dev.off()
}


