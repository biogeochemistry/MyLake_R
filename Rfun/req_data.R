load("albedot1.bin")

obj<-readMat("varsimp_1.mat")
Wt<-obj$Wt
tt<-obj$tt
Bio_par_names<-obj[[3]]
Bio_par_range<-obj[[4]]
Bio_par<-obj[[5]]
Phys_par_names<-obj[[6]]
Phys_par_range<-obj[[7]]
Phys_par<-obj[[8]]
Inflw<-obj[[9]]
Ice0<-obj[[10]]
In_FIM<-obj[[11]]
In_Chlz_sed<-obj[[12]]
In_TPz_sed<-obj[[13]]
In_DOCz<-obj[[14]]
In_Chlz<-obj[[15]]
In_DOPz<-obj[[16]]
In_TPz<-obj[[17]]
In_Sz<-obj[[18]]
In_Cz<-obj[[19]]
In_Tz<-obj[[20]]
In_Az<-obj[[21]]
In_Z<-obj[[22]]
rm(obj)

