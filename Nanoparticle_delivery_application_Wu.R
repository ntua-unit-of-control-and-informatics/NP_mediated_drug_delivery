# Working directory

setwd("C:/Users/user/Desktop/Nanoparticle_delivery_application")


# *** metrics ***

# The metric used for the optimization
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

mape <- function(observed, predicted){
  mean(abs(observed-predicted)*100/observed)
}

rmse <- function(observed, predicted){
  sqrt(mean((observed-predicted)^2)) 
}

AAFE <- function(observations, predictions, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

SODI <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    Et <- 0 #relative error with observations
    St <- 0  #relative error with simulations
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate the N_obs vector
    for (j in 1:N){
      # sum of relative squared errors (error = observed - predicted)
      Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
      St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
    }
    
    # root mean of the square of observed values
    RMEt <- sqrt(Et/N)
    # root mean of the square of simulated values
    RMSt <- sqrt( St/N)
    
    I[i] <- (RMEt + RMSt)/2
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    # Give weight to compartments with more observations (more information)
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}



#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function( dose){
  NP_AS<-dose; NP_CV<-0;NP_CI <- 0; NP_OT <- 0;NP_EX <- 0;
  
  return(c( "NP_AS" = NP_AS, "NP_CV"=NP_CV,"NP_CI" = NP_CI,
            "NP_OT"=NP_OT, "NP_EX" = NP_EX))
}
#==============
#3. ODEs System
#==============
Wu_model <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Description:
    # k_AStCV: Administrtion Site to target Cell Vicinity
    # k_CVtAS: Cell Vicinity to  Administration Site
    # k_CVtCI: Cell Vicinity to  Cell Interior
    # k_AStOT: Administration Site to Off-Target Sites
    # k_OTtAS: Off-Target sites to Administration Site
    # k_CItCV: Cell Interior to Cell Vicinity
    # k_CVtOT: Cell Vicinity to Off-Target sites
    # k_OTtCV: Off-Target sites to Cell Vicinity
    # CL_CItEX: Maximum clearance from Cell Interior to excreta
    # Km_CItEX: km constant from Cell Interior to excreta
    # CL_OTtEX: Maximum clearance from off targets to excreta
    # Km_OTtEX: km constant from off-targets to excreta
    
    # NP_AS: nanoparticles in Administration Site
    # NP_CV: nanoparticles in target Cell Vicinity
    # NP_CI: nanoparticles in target Cell Interior
    # NP_OT: nanoparticles in Off-Target sites
    
    # Units:
    # k_TTC, k_C, k_ITC, k_ATC ---> 1/
    # NP_AS, NP_CV, NP_CI, NP_OT ---> 1/
    CL_CItEX <- 0
    Km_CItEX <- 1
    
    
    dNP_AS <-  -k_AStCV*NP_AS + k_CVtAS*NP_CV - k_AStOT*NP_AS + k_OTtAS*NP_OT
    dNP_CV <-  k_AStCV*NP_AS - k_CVtAS*NP_CV- k_CVtCI*NP_CV + k_CItCV*NP_CI -
               k_CVtOT*NP_CV + k_OTtCV*NP_OT
    dNP_CI <-  k_CVtCI*NP_CV - k_CItCV*NP_CI -(CL_CItEX* NP_CI)/(Km_CItEX+NP_CI)
    dNP_OT <-  k_AStOT*NP_AS -  k_OTtAS*NP_OT+ k_CVtOT*NP_CV -  k_OTtCV*NP_OT-
      (CL_OTtEX* NP_CI)/(Km_OTtEX+NP_CI) * NP_OT
    dNP_EX <- (CL_CItEX* NP_CI)/(Km_CItEX+NP_CI) +  (CL_OTtEX* NP_CI)/(Km_OTtEX+NP_CI) * NP_OT
    NP_tot <- NP_AS + NP_CV + NP_OT + NP_CI
    selectivity <- NP_CI/NP_tot
    
    return(list(c("dNP_AS" = dNP_AS,   "dNP_CV" = dNP_CV,
                  "dNP_CI" = dNP_CI, "dNP_OT" = dNP_OT, "dNP_EX" = dNP_EX), 
                "selectivity" = selectivity))
  })
}



obj_func <- function(x, dose, df, metric = "AAFE"){
  
  BodyBurden <- c(df$administration_site, df$cell_vicinity, df$cell_interior, df$off_target,
                  df$excreta)

  parms <- c("k_AStCV" = (x[1]),  "k_CVtAS" = (x[2]),
             "k_CVtCI" =  (x[3]),
             "k_AStOT" =  (x[4]) , "k_CItCV" =  (x[5]),
             "k_OTtAS" = (x[6]) , "k_CVtOT" =  (x[7]),
             "k_OTtCV"=  (x[8]),  "CL_OTtEX"=  (x[9]), 
             "Km_OTtEX"=  (x[10]) )
  
  sol_times <- c(seq(0,1, 0.001),seq(1.1,5, 0.1),  seq(6,28*24, 1))
  inits <- create.inits(unname(dose))
  solution <- data.frame(deSolve::ode(times = sol_times,  func = Wu_model,
                                      y = inits, parms = parms, method="lsodes",
                                      rtol = 1e-7, atol = 1e-7))
  if(sum(solution$time %in% df$time) == length( df$time)){
    results <- c(solution[solution$time %in% df$time, "NP_AS"],
                 solution[solution$time %in% df$time, "NP_CV"],
                 solution[solution$time %in% df$time, "NP_CI"],
                 solution[solution$time %in% df$time, "NP_OT"],
                 solution[solution$time %in% df$time, "NP_EX"])
  }else{
    results <-   BodyBurden *100
    
  }
  
  
  # Find the position of the current PFAS in the PFAS_names vector
  if(metric == "AAFE"){
    score <- AAFE(BodyBurden, results) 
  }else if (metric =="rmse"){
    score <- rmse(BodyBurden, results)
  }else if(metric == "SODI"){
    score <- SODI(list(BodyBurden), list(results))
  }       
  return(score)
}

##############################################################################
################################################################################
# Load the data for PFAS concentration
case <- openxlsx::read.xlsx ('case2.xlsx')
dose_per_weight <- 0.7 #mg/kg
dose <- dose_per_weight*229

opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",NLOPT_LN_SBPLX , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-7, 
              "ftol_rel" = 1e-7,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 5000,
              "print_level" = 1)

N_pars <-10
# Define initial values of fitted parameters to provide to the optimization routine
x0 <-  rep(10,N_pars)
set.seed(1221312)
optimization<- nloptr::nloptr(x0 = x0,
                              eval_f = obj_func,
                              lb	= rep(0.00001,N_pars),
                              ub =   rep(10000,N_pars),
                              opts = opts,
                              dose = dose,
                              df = case,
                              metric = "AAFE")

parms <- c("k_AStCV" =(optimization$solution[1]), 
           "k_CVtAS" =(optimization$solution[2]),
           "k_CVtCI" = (optimization$solution[3]),
           "k_AStOT" = (optimization$solution[4]) , 
           "k_CItCV" = (optimization$solution[5]),
           "k_OTtAS" = (optimization$solution[6]) ,
           "k_CVtOT" = (optimization$solution[7]),
           "k_OTtCV"= (optimization$solution[8]),
           "CL_OTtEX"=  (optimization$solution[9]), 
           "Km_OTtEX"=   (optimization$solution[10]))




sol_times <- seq(0,28*24, 1 )
inits <- create.inits(unname(dose))
solution <- data.frame(deSolve::ode(times = sol_times,  func = Wu_model,
                                    y = inits, parms = parms, method="lsodes",
                                    rtol = 1e-7, atol = 1e-7))



library(ggplot2)
# administration site plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_AS ), size=1.7)+
  geom_point(data = case, aes(x=time  , y=administration_site ), size=5)+
  
  labs(title = "NP mass in administration site", y = "Micrograms", x = "Time (hours)")+
  
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22)) + 
  
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))

########
#cell vicinity plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_CV ), size=1.7)+
  geom_point(data = case, aes(x=time  , y=cell_vicinity ), size=5)+
  
  labs(title = "NP mass in target cell vicinity", y = "Micrograms", x = "Time (hours)")+
  
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22)) + 
  
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))

########
# Cell interior plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_CI ), size=1.7)+
  geom_point(data = case, aes(x=time  , y=cell_interior ), size=5)+
  
  labs(title = "NP mass in target cell interior", y = "Micrograms", x = "Time (hours)")+
  
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22)) + 
  
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))


########
# Off_Target_sites plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_OT ), size=1.7)+
  geom_point(data = case, aes(x=time  , y= off_target ), size=5)+
  
  labs(title = "NP mass in off-target sites", y = "Micrograms", x = "Time (hours)")+
  
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22)) + 
  
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))

########
# Excreta plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_EX ), size=1.7)+
  geom_point(data = case, aes(x=time  , y= excreta ), size=5)+
  
  labs(title = "NP mass in excreta", y = "Micrograms", x = "Time (hours)")+
  
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22)) + 
  
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))



