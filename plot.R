library(ggplot2)
library(scales)
load("plot.RData")

sd_threshold=0.01
alphaLevel=0.05
##############################################################
##################### plots in the paper #####################
##############################################################

###### empirical cdf for each group (Fig. 4) ######
jpeg(paste("empirical_cdf_delta",sd_threshold,"_alpha",alphaLevel,".jpeg",sep=""), quality = 100,width = 1000, height = 800, res = 200)
plot(error4_plot,type='l',log="x", xaxt = "n",cex.axis=1.2,cex.lab=1.2,
     ylab="Empirical CDF",xlab = "Absolute error",col=5,mgp=c(2,1,0))
title(bquote(list(delta==.(sd_threshold), alpha==.(alphaLevel))), line = 1)
xticks <- seq(-13L,0L,2L)
for (i in 1:length(xticks)){
  axis(1, at=10^xticks[i], labels=bquote(10^.(xticks[i])),cex.axis=1.2)
}
lines(error3_plot,type='l',col=4)
lines(error2_plot,type='l',col=3)
lines(error1_plot,type='l',col=2)
abline(v=sd_threshold)
abline(h=(1-alphaLevel),col=6,lty=2)
points(xpoints[4],y=(par("usr")[3]),col=5,pch=20,xpd = TRUE)
points(xpoints[3],y=(par("usr")[3]),col=4,pch=20,xpd = TRUE)
points(xpoints[2],y=(par("usr")[3]),col=3,pch=20,xpd = TRUE)
points(xpoints[1],y=(par("usr")[3]),col=2,pch=20,xpd = TRUE)
legend(x=10^(-11),y=0.95, legend = c("group 1","group 2","group 3", "group 4",as.expression(bquote(delta == .(sd_threshold))),paste((1-alphaLevel)*100,"th percentile",sep="")),
       lty = c(1,1,1,1,1,2),col=c(2:5,1,6),bty="n",cex=1)
dev.off()


###### RMSE comparison (Fig. 5) ######
jpeg("RMSE_model_separately.jpeg", quality = 100,width = 1300, height = 600, res = 200)
ggplot(dat_sep, aes(xmin=as.numeric(as.factor(Group))-.4,xmax=as.numeric(as.factor(Group))+.4, x=Group, ymin=1E-7, ymax=RMSE_rho, fill=Design)) +
  geom_rect(position=position_dodge(.8))+
  geom_point(aes(x=Group, y=RMSE_Omega,group=Design,color=Design))+
  geom_line(aes(x=Group, y=RMSE_Omega,group=Design,color=Design))+ 
  #geom_text(aes(x=Group, y=RMSE_rho,label=n_train,group=Design,vjust=0),position = position_dodge(width= .8))+
  scale_y_log10("RMSE",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_classic()+theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14))+
  scale_fill_brewer(palette="Pastel1")+scale_color_brewer(palette="Set1")
dev.off()

###### Group 4 comparison (Fig. 6) ######
#6(a)
jpeg("boxplot.jpeg", quality = 100,width = 850, height = 1000, res = 220)
par(mgp=c(2,1,0))
boxplot(RMSE~design,data=RMSE4,col=c(2,3,4,7),outpch=20,xlab='Design',ylab=expression("RMSE"[rho]),cex.lab=1.2)
dev.off()
#6(b)
jpeg("scatter.jpeg", quality = 100,width = 1000, height = 1000, res = 300)
plot(NULL,xlim=range(energy4_RS1,energy4_RS2),ylim=range(energy4_RS1,energy4_RS2),
     xlab=expression("true"~Omega),ylab=expression("predicted"~hat(Omega)),mgp=c(2,1,0))
points(energy4_RS1,pch=15,col=3,cex=0.3)
points(energy4_RS2,pch=17,col=4,cex=0.3)
points(energy4_Dopt,pch=18,col=7,cex=0.4)
points(energy4_ALEC,pch=16,col=2,cex=0.3)
legend("topleft",legend = c("RS1", "RS2","D-opt", "ALEC"),col=c(3,4,7,2),pch=c(15,17,18,16),cex=0.8,bty="n")
points(select4[1,1],select4[1,2],pch=1,col=1,cex=0.6)
points(select4[2,1],select4[2,2],pch=1,col=1,cex=0.6)
dev.off()
#6(c)
jpeg("density1.jpeg", quality = 100,width = 750, height = 500, res = 150)
plot(seq(-4.5,4.5,0.01),density1[1,],type="l",xlab="",ylab=expression(rho),ylim=c(0,0.75),mgp=c(2,1,0),cex.lab=1.2)
points(seq(-4.5,4.5,0.01),density1[2,],pch=15,col=3,cex=0.8)
points(seq(-4.5,4.5,0.01),density1[3,],pch=17,col=4,cex=0.9)
points(seq(-4.5,4.5,0.01),density1[4,],pch=18,col=7,cex=0.8)
points(seq(-4.5,4.5,0.01),density1[5,],pch=16,col=2,cex=0.6)
lines(seq(-4.5,4.5,0.01),density1[1,],lwd=1.4)
legend("topright",legend = c("truth","RS1","RS2","D-opt","ALEC"),pch=c(NA,15,17,18,16),lty=c(1,NA,NA,NA,NA),
       lwd=c(1.4,NA,NA,NA,NA),col=c(1,3,4,7,2),bty="n")
dev.off()
#6(d)
jpeg("density2.jpeg", quality = 100,width = 750, height = 500, res = 150)
plot(seq(-4.5,4.5,0.01),density2[1,],type="l",xlab="",ylab=expression(rho),ylim=c(0,0.7),mgp=c(2,1,0),cex.lab=1.2)
points(seq(-4.5,4.5,0.01),density2[2,],pch=15,col=3,cex=0.8)
points(seq(-4.5,4.5,0.01),density2[3,],pch=17,col=4,cex=0.9)
points(seq(-4.5,4.5,0.01),density2[4,],pch=18,col=7,cex=0.8)
points(seq(-4.5,4.5,0.01),density2[5,],pch=16,col=2,cex=0.7)
lines(seq(-4.5,4.5,0.01),density2[1,],lwd=1.4)
# legend("topright",legend = c("truth","RS1","RS2","D-opt","ALEC"),pch=c(NA,15,17,18,16),lty=c(1,NA,NA,NA,NA),
#        lwd=c(1.4,NA,NA,NA,NA),col=c(1,3,4,7,2),bty="n")
dev.off()

###### All groups RMSE comparison (Fig. 7) ######

jpeg("RMSE_model_all_seq.jpeg", quality = 100,width = 1300, height = 600, res = 200)
ggplot(dat_all, aes(xmin=as.numeric(as.factor(Group))-.4,xmax=as.numeric(Group)+.4, x=Group, ymin=1E-7, ymax=RMSE_rho, fill=Design, alpha = factor(Group))) +
  geom_rect(position=position_dodge(.8))+
  geom_point(aes(x=Group, y=RMSE_Omega,group=Design,color=Design))+
  geom_line(aes(x=Group, y=RMSE_Omega,group=Design,color=Design),alpha=0.6)+
  scale_alpha_manual(values = c("all groups"=1,"group 1"=0.6, "group 2"=0.6,"group 3"=0.6,"group 4"=0.6), guide='none')+
  scale_y_log10("RMSE",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_classic()+theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14))+
  scale_fill_brewer(palette="Pastel1")+scale_color_brewer(palette="Set1")
dev.off()

###### time comparison (Fig. 8) ######
jpeg("time_separately.jpeg", quality = 100,width = 1000, height = 500, res = 250)
ggplot(time_sep,aes(x=Method,y=Time,fill=Solver))+
  geom_col()+
  facet_grid(~ Group)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Time (s)")+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),axis.title=element_text(size=12))
dev.off()

jpeg("time_all.jpeg", quality = 100,width = 510, height = 500, res = 250)
ggplot(time_all,aes(x=Method,y=Time,fill=Solver))+
  geom_col()+
  facet_grid(~ Group)+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Time (s)")+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=8),axis.title=element_text(size=12))
dev.off()

###### Mixed/new group (Fig. 9) ######
#9 left
jpeg(paste("empirical_cdf_mix23all_delta",sd_threshold,"_alpha",alphaLevel,".jpeg",sep=""), quality = 100,width = 1000, height = 800, res = 200)
plot(error5w4_plot,type='l',log="x", xaxt = "n",cex.lab=1.2,
     ylab="Empirical CDF",xlab = "Absolute error",col=5,mgp=c(2,1,0))
xticks <- seq(-13L,0L,2L)
for (i in 1:length(xticks)){
  axis(1, at=10^xticks[i], labels=bquote(10^.(xticks[i])),cex.lab=1.2)
}
lines(error5w3_plot,type='l',col=4)
lines(error5w2_plot,type='l',col=3)
lines(error5w1_plot,type='l',col=2)
lines(error5wall_plot,type='l',col=7)
abline(v=sd_threshold)
abline(h=(1-alphaLevel),col=6,lty=2)
points(xpoints_mix[4],y=(par("usr")[3]),col=5,pch=20,xpd = TRUE)
points(xpoints_mix[3],y=(par("usr")[3]),col=4,pch=20,xpd = TRUE)
points(xpoints_mix[2],y=(par("usr")[3]),col=3,pch=20,xpd = TRUE)
points(xpoints_mix[1],y=(par("usr")[3]),col=2,pch=20,xpd = TRUE)
points(xpoints_mix[5],y=(par("usr")[3]),col=7,pch=20,xpd = TRUE)
legend("left", legend = c("group 1","group 2","group 3", "group 4","all groups",as.expression(bquote(delta == .(sd_threshold))),paste((1-alphaLevel)*100,"th percentile",sep="")),lty = c(1,1,1,1,1,1,2),col=c(2:5,7,1,6),bty="n",cex=1)
dev.off()
#9 middle
jpeg("coordinate_coverage_all.jpg", quality = 100,width = 1000, height = 800, res = 200)
plot(error_col4,type="l",col=5,lwd=1.5,ylab="Percentage",xaxt = "n",xlab="Coordinate",cex.lab=1.2,mgp=c(2,1,0))
lines(error_col3,col=4,lwd=1.5)
lines(error_col2,col=3,lwd=1.5)
lines(error_col1,col=2,lwd=1.5)
lines(error_colall,col=7,lwd=1.5)
abline(h=alphaLevel,col=6,lty=2,lwd=1.5)
legend("top", legend = c("group 1","group 2","group 3", "group 4","all groups",expression(paste(alpha==0.05))),lty = c(1,1,1,1,1,2),col=c(2:5,7,6),lwd=1.5,bty="n",cex=1)
dev.off()
# 9 right
jpeg(paste("empirical_cdf_mix23all_max_delta",sd_threshold,"_alpha",alphaLevel,".jpeg",sep=""), quality = 100,width = 1000, height = 800, res = 200)
plot(error5w4max_plot,type='l',log="x", xaxt = "n",cex.lab=1.2,
     ylab="Empirical CDF",xlab = "Absolute error",col=5,mgp=c(2,1,0))
xticks <- seq(-13L,0L,2L)
for (i in 1:length(xticks)){
  axis(1, at=10^xticks[i], labels=bquote(10^.(xticks[i])),cex.lab=1.2)
}
lines(error5w3max_plot,type='l',col=4)
lines(error5w2max_plot,type='l',col=3)
lines(error5w1max_plot,type='l',col=2)
lines(error5wallmax_plot,type='l',col=7)
abline(v=sd_threshold)
abline(h=(1-alphaLevel),col=6,lty=2)
points(xpoints_mix_max[4],y=(par("usr")[3]),col=5,pch=20,xpd = TRUE)
points(xpoints_mix_max[3],y=(par("usr")[3]),col=4,pch=20,xpd = TRUE)
points(xpoints_mix_max[2],y=(par("usr")[3]),col=3,pch=20,xpd = TRUE)
points(xpoints_mix_max[1],y=(par("usr")[3]),col=2,pch=20,xpd = TRUE)
points(xpoints_mix_max[5],y=(par("usr")[3]),col=7,pch=20,xpd = TRUE)
legend("left", legend = c("group 1","group 2","group 3", "group 4","all groups",as.expression(bquote(delta == .(sd_threshold))),paste((1-alphaLevel)*100,"th percentile",sep="")),lty = c(1,1,1,1,1,1,2),col=c(2:5,7,1,6),bty="n",cex=1)
dev.off()

# jpeg("coordinate_coverage_max_all.jpg", quality = 100,width = 1000, height = 800, res = 200)
# plot(error_colallmax,type="l",col=7,lwd=1.5,ylab="Percentage",xaxt = "n",xlab="Coordinate",cex.lab=1.2,mgp=c(2,1,0))
# lines(error_col3max,col=4,lwd=1.5)
# lines(error_col2max,col=3,lwd=1.5)
# lines(error_col1max,col=2,lwd=1.5)
# lines(error_col4max,col=5,lwd=1.5)
# abline(h=alphaLevel,col=6,lty=2,lwd=1.5)
# legend("top", legend = c("group 1","group 2","group 3", "group 4","all groups",expression(paste(alpha==0.05))),lty = c(1,1,1,1,1,2),col=c(2:5,7,6),lwd=1.5,bty="n",cex=1)
# dev.off()

