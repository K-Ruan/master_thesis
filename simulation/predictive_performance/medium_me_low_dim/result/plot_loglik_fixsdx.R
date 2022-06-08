library(ggplot2)
library(tidyverse)
mos=log10(c(1.5,3,6,20,40,70,100,1e5,1e10,1e20)) #here it is: max(me)/sdX

cvloglik_gmus_result<-read.csv(file="cvloglik_gmus_fixsdX.csv")
cvloglik_gmus0_result<-read.csv(file="cvloglik_gmus0_fixsdX.csv")
cvloglik_gmus00_result<-read.csv(file="cvloglik_gmus00_fixsdX.csv")
cvloglik_gds_result<-read.csv(file="cvloglik_gds_fixsdX.csv")

cvloglik_gmus_result<-cvloglik_gmus_result$x #1cvloglik(5-fold-cv)200samples * B=50 times)* 10 mos
cvloglik_gmus0_result<-cvloglik_gmus0_result$x #1cvloglik(5-fold-cv)200samples * B=50 times)* 10 mos
cvloglik_gmus00_result<-cvloglik_gmus00_result$x #1cvloglik(5-fold-cv)200samples * B=50 times)* 10 mos
cvloglik_gds_result<-cvloglik_gds_result$x #1cvloglik(5-fold-cv)200samples *  B=50 times)* 10 mos

#60samples * B=50 *10 mos

split1<-seq(1,length(cvloglik_gmus_result),by=200) #1:10,0000
split2<-seq(200,length(cvloglik_gmus_result),by=200)# 1:10,0000


#######################


mean_cvloglik_gmus<-c()
mean_cvloglik_gmus0<-c()
mean_cvloglik_gmus00<-c()
mean_cvloglik_gds<-c()


for(i in 1:length(split1)){
  mean_cvloglik_gmus<-c(mean_cvloglik_gmus,mean(cvloglik_gmus_result[split1[i]:split2[i]],na.rm=T))
}

for(i in 1:length(split1)){
  mean_cvloglik_gmus0<-c(mean_cvloglik_gmus0,mean(cvloglik_gmus0_result[split1[i]:split2[i]],na.rm=T))
}

for(i in 1:length(split1)){
  mean_cvloglik_gmus00<-c(mean_cvloglik_gmus00,mean(cvloglik_gmus00_result[split1[i]:split2[i]],na.rm=T))
}

for(i in 1:length(split1)){
  mean_cvloglik_gds<-c(mean_cvloglik_gds,mean(cvloglik_gds_result[split1[i]:split2[i]],na.rm=T))
}

split3<-seq(1,500,by=50)
split4<-seq(50,500,by=50)
#splitgds3<-seq(1,15000,by=1500)
#splitgds4<-seq(1500,15000,by=1500)

mos=log10(c(1.5,3,6,20,40,70,100,1e5,1e10,1e20))
df_cvloglik_gmus<-as.data.frame(cbind(mean_cvloglik_gmus, rep(mos,each=50)))
colnames(df_cvloglik_gmus)<-c("cvloglik_gmus","mos")
df_cvloglik_gmus0<-as.data.frame(cbind(mean_cvloglik_gmus0, rep(mos,each=50)))
colnames(df_cvloglik_gmus0)<-c("cvloglik_gmus0","mos")

df_cvloglik_gmus00<-as.data.frame(cbind(mean_cvloglik_gmus00, rep(mos,each=50)))
colnames(df_cvloglik_gmus00)<-c("cvloglik_gmus00","mos")
df_cvloglik_gds<-as.data.frame(cbind(mean_cvloglik_gds,rep(mos,each=50)))
colnames(df_cvloglik_gds)<-c("cvloglik_gds","mos")


ggplot(df_cvloglik_gmus,aes(x=log(mos),y=cvloglik_gmus, fill=mos,group=mos)) +
  geom_boxplot() +
  theme(legend.position = "none")


ggplot(df_cvloglik_gmus0,aes(x=log(mos),y=cvloglik_gmus0, fill=mos,group=mos)) +
  geom_boxplot() +
  theme(legend.position = "none")

ggplot(df_cvloglik_gmus00,aes(x=log(mos),y=cvloglik_gmus00, fill=mos,group=mos)) +
  geom_boxplot() +
  theme(legend.position = "none")

ggplot(df_cvloglik_gds,aes(x=log(mos),y=cvloglik_gds, fill=mos,group=mos)) +
  geom_boxplot() +
  theme(legend.position = "none")

mean_cvllgmus_permos<-df_cvloglik_gmus %>%
  group_by(mos) %>%
  summarise(mean_cvll = mean(cvloglik_gmus, na.rm = TRUE))

mean_cvllgmus0_permos<-df_cvloglik_gmus0 %>%
  group_by(mos) %>%
  summarise(mean_cvll = mean(cvloglik_gmus0, na.rm = TRUE))

mean_cvllgmus00_permos<-df_cvloglik_gmus00 %>%
  group_by(mos) %>%
  summarise(mean_cvll = mean(cvloglik_gmus00, na.rm = TRUE))

mean_cvllgds_permos<-df_cvloglik_gds %>%
  group_by(mos) %>%
  summarise(mean_cvll = mean(cvloglik_gds, na.rm = TRUE))



plot(x=mean_cvllgmus_permos$mos,y=mean_cvllgmus_permos$mean_cvll,type="l",col="blue")
lines(x=mean_cvllgmus0_permos$mos,y=mean_cvllgmus0_permos$mean_cvll,col="red")
lines(x=mean_cvllgmus00_permos$mos,y=mean_cvllgmus00_permos$mean_cvll,col="orange")

xxx1<-cbind(mean_cvloglik_gmus, rep(mos,each=50)) #gmus B=50 * 10 mos
colnames(xxx1)<-c()
xxx2<-cbind(mean_cvloglik_gmus0, rep(mos,each=50)) #gmus0 B=50 * 10 mos
colnames(xxx2)<-c()
xxx3<-cbind(mean_cvloglik_gmus00, rep(mos,each=50)) #gmus00 B=50 * 10 mos
colnames(xxx3)<-c()
xxx4<-cbind(mean_cvloglik_gds,rep(mos,each=50)) #gds B=50 * 10 mos
colnames(xxx3)<-c()

df_combined<-as.data.frame(rbind(xxx1,xxx2,xxx3,xxx4))
colnames(df_combined)<-c("mean_loglik","mos")
df_combined$mos<-round(df_combined$mos,3)
df_combined$type<-rep(c("CGMUS","CGMUS0","GMUS","GDS"),each=500)
#"gmus","gmus0","gmus00","gds"

df_combined$mos<-as.factor(df_combined$mos)
df_combined$type<-as.factor(df_combined$type)


ggplot(df_combined,aes(x=mos,y=mean_loglik, fill=type)) +
  geom_boxplot()+
  labs(title = "ME setting 1, n=200", x="Constant C")
#theme(legend.position = "top")

############ paired ttest
lvl<-c("0.176","0.477","0.778","1.301","1.602","1.845","2","5" ,"10","20")
df_combined[df_combined$mos==lvl[1],]

# loglik_gmus - loglik_gmus0
#mos1_diff<-df_combined[df_combined$type=="gmus",]$mean_loglik[1:50]-df_combined[df_combined$type=="gmus0",]$mean_loglik[1:50]
##############

#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[1:50],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[1:50],paired=T)
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[51:100],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[51:100],paired=T)
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[101:150],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[101:150],paired=T)
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[151:200],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[151:200],paired=T)
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[201:250],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[201:250],paired=T)
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[251:300],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[251:300],paired=T)
#
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[1:400],
#       y=df_combined[df_combined$type=="gmus0",]$mean_loglik[1:400],paired=T,alternative = "greater")
#
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[1:400],
#       y=df_combined[df_combined$type=="gds",]$mean_loglik[1:400],paired=T,alternative = "greater")
#
#t.test(x=df_combined[df_combined$type=="gmus",]$mean_loglik[51:100],
#       y=df_combined[df_combined$type=="gds",]$mean_loglik[51:100],paired=T,alternative = "two.sided")

#df_combined$paired<-as.numeric(df_combined$mos)
#ggplot(df_combined,aes(x=mos,y=mean_loglik, fill=type)) +
#  geom_boxplot() + 
#  geom_line(aes(group=paired), position = position_dodge(0.2)) 