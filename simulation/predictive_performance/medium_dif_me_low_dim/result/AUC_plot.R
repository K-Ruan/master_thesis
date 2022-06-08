library(ggplot2)
AUC_mat<-read.csv(file="AUC_Gmus_Gmus0_Gmus00_Gds_fixsdX.csv",row.names = 1,header=T) #50rows * 40columns
mos=log10(c(1.5,3,6,20,40,70,100,1e5,1e10,1e20))

temp<-as.matrix(AUC_mat)
temp_df<-rbind(temp[,c(1,2,3,4)],    temp[,c(5,6,7,8)],    temp[,c(9,10,11,12)],  temp[,c(13,14,15,16)],  temp[,c(17,18,19,20)],
               temp[,c(21,22,23,24)],temp[,c(25,26,27,28)],temp[,c(29,30,31,32)], temp[,c(33,34,35,36)],  temp[,c(37,38,39,40)])

df<-cbind(temp_df,rep(round(mos,digits = 3),each=50))

long_df<-as.data.frame( cbind( rbind(df[,c(1,5)],df[,c(2,5)],df[,c(3,5)],df[,c(4,5)]), rep(c("CGMUS","CGMUS0","GMUS","GDS"),each=500) ) )
colnames(long_df)<-c("AUC","mos","type")



long_df$type<-as.factor(long_df$type)
long_df$mos<-as.factor(as.numeric(long_df$mos))
long_df$AUC<-as.numeric(long_df$AUC)

ggplot(long_df,aes(x=mos,y=AUC, fill=type)) +
  geom_boxplot() +
  labs(title="ME setting 2, n=200", x="Constant C")




ggplot(long_df, aes(x=mos, y=AUC, group=type,color=type)) + 
  geom_point()