tfun<-function(x){

if(a>0) {alb <- alb_melt_snow
}else{if (b>0) { alb <- alb_melt_ice 
} else alb <- albedo_mod(Trnsmiss, alt, alb_table) } 

return(alb)
}
