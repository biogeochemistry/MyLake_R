# === MyLake model, version 1.2, 15.03.05 ===
# by Tom Andersen & Tuomo Saloranta, NIVA 2004
# Convection module
       
convection_v12<-function(Tz_in,Cz_in,Sz_in,Pz_in,Chlz_in,PPz_in,
	DOPz_in,DOCz_in,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,
	zz,swa_b0,tracer_switch,springautumn){

# Inputs (with extension "_in") and Outputs:
#       Tz   : Temperature profile
#       Cz   : Tracer profile
#       Sz   : Suspended inorg. matter profile
#       Pz   : Dissolved inorg. P profile
#       Chlz : Chlorophyll a profile
#       PPz  : Phosphorus bound to inorganic particles profile
#       DOPz  : Dissolved organic phosphorus profile
#       DOCz  : Particulate inorganic phosphorus profile

# Inputs: 
#       Tprof_prev   : Temperature profile from previous timestep
#       etc.

# These variables are still global and not transferred by functions


Trhomax<-3.98; #temperature of maximum water density (deg C) 
Nz<-length(zz); #total number of layers in the water column
dz<-zz[2]-zz[1]; #model grid step

  # Convective mixing adjustment
  # Mix successive layers until stable density profile
   Rho <- rho(0, Tz_in, 0) #old: polyval(ies80,max(0,Tz_in(:)))+min(Tz_in(:),0);	# Density (kg/m3)
   d_rho<-c(diff(Rho), 1); #d_rho = how much a layer is lighter than layer below; last cell in "d_rho" is always positive (sediment)
   while (length(which(d_rho<0))){ 	#old:any(d_rho < 0),
      blnUnstb_layers<-1*(d_rho <= 0); #1=layer is heavier or equal than layer below, 0=layer is lighter than layer below
      A_Unstb<-which(diff(c(0,blnUnstb_layers))==1); #layer index(es) where unstable/neutral column(s) start(s)
      B_Unstb<-which(diff(c(0,blnUnstb_layers))==-1)-1;#layer index(es) where unstable/neutral column(s) end(s)
        
      for (n in 1:length(A_Unstb)){
       j <- c(A_Unstb[n]:(B_Unstb[n]+1))
       Tz_in[j] <- weighted.mean(Tz_in[j], Vz[j])

        if (tracer_switch==1) Cz_in[j]<- weighted.mean(Cz_in[j], Vz[j])
        
       Sz_in[j]<- weighted.mean(Sz_in[j], Vz[j])

       Pz_in[j]<- weighted.mean(Pz_in[j], Vz[j])
       
	 Chlz_in[j]<- weighted.mean(Chlz_in[j], Vz[j])
       
       PPz_in[j]<- weighted.mean(PPz_in[j], Vz[j])
      
       DOPz_in[j]<- weighted.mean(DOPz_in[j], Vz[j])
       
       DOCz_in[j]<- weighted.mean(DOCz_in[j], Vz[j])
      }#'for'
             
      Rho <- rho(0, Tz_in, 0) # rho = polyval(ies80,max(0,Tz_in(:))) + min(Tz_in(:),0);
      d_rho<-c(diff(Rho), 1); #d_rho=[diff(rho); 1];
   }# 'while'
      if (springautumn==1){    # 'if -springautumn-'         
        # Spring/autumn turnover
        # don't allow temperature jumps over temperature of maximum density
        jumpinx<-which( ((Tprof_prev>=Trhomax)&(Tz_in<Trhomax))|((Tprof_prev<=Trhomax)&(Tz_in>Trhomax)) );
	  if (length(jumpinx)>0){# 'if -1-'
        sellay<-jumpinx[1]:length(Tz_in)
         intSign<-sign(Trhomax-Tz_in[jumpinx[1]]); #plus in autumn turnover, minus in spring turnover   
         XE_turn<-cumsum((Tz_in[sellay]-Trhomax)*Vz[sellay]*Cw*intSign); #always starts negative 
         Dummy<-which(XE_turn>0);
           if(length(Dummy)==0){ # 'if -2-'
             Tz_in[sellay]<-Trhomax;
                if (intSign==1){# 'if -3-'
                 Tz_in[jumpinx[1]]<-Tz_in[jumpinx[1]]+intSign*XE_turn[length(XE_turn)]/(Vz[jumpinx[1]]*Cw) #put overshoot on top layer
                } else {
                 distrib_key<-(f_par * exp(c(0, -lambdaz_wtot_avg) * c(zz, zz[length(zz)]+dz) ) + 
                  (1-f_par) * exp(c(0, rep(-swa_b0, l=Nz)) * c(zz, max(zz)+dz)))
                  Tz_in[sellay]<-Tz_in[sellay] + (-diff( intSign*XE_turn[length(XE_turn)] * 
			  (distrib_key[jumpinx[1]:length(distrib_key)]/distrib_key[jumpinx[1]]) ) /(Vz[sellay]*Cw))
                 #distribute overshoot as shortwave energy}
                }# 'if -3-'
           } else {    
             Tz_in[jumpinx[1]:(jumpinx[1]+Dummy[1]-1-1)]<-Trhomax;
             Tz_in[jumpinx[1]+Dummy[1]-1]<-Trhomax+intSign*(XE_turn[Dummy[1]])/(Vz[jumpinx[1]+Dummy[1]-1]*Cw); 
           } # 'if -2-'
        } # 'if -1-'
      }# 'if -springautumn-'
   
        Tz<-Tz_in
        Cz<-Cz_in
        Sz<-Sz_in
        Pz<-Pz_in
        PPz<-PPz_in
        Chlz<-Chlz_in
        DOPz<-DOPz_in
        DOCz<-DOCz_in    

return(list(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz))
}