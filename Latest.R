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

#####################################
### Function to create Parameters ###
#####################################
create.params <<- function(weight){
  # Physiological parameters units
  # V_blood, V_ven, V_art (ml): Volume of total blood, venous blood and arterial blood
  # w_i (g):                    mass of tissue or organ "i"
  # V_tis_i (ml):                volume of tissue or organ "i"
  # V_cap_i (ml):                volume of capillary blood in tissue "i"
  # Q_i, Q_total (ml/h):        regional blood flow of tissue or organ "i"
  
  # List with names of  compartments
  compartments <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", 
                       "Brain"="Brain", "Spleen"="Spleen",
                       "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", 
                       "Bone"="Bone", "Adipose"="Adipose", "Skin"="Skin",
                       "Muscles"="Muscles",
                       "GIT"="GIT") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*weight^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*weight+0.77 # Total blood volume (ml)
  
  #The following equation gives the  adipose % of body weight 
  fr_ad <- 0.0199*weight + 1.644 # w in g,  Brown et al.1997 p.420. 
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat_physiological_parameters.xlsx",
                                   sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  rownames(fractions) <- compartments
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  
  W_tis <- rep(0,length(compartments))
  V_tis <- rep(0,length(compartments))
  V_cap <- rep(0,length(compartments))
  
  # Retrieval of tissue weight fractions
  fw_heart <- Tissue_fractions[2]
  fw_kidneys <- Tissue_fractions[3]
  fw_brain <- Tissue_fractions[4]
  fw_spleen <-Tissue_fractions[5]
  fw_lungs <- Tissue_fractions[6]
  fw_liver <- Tissue_fractions[7]
  fw_skeleton <- Tissue_fractions[8]
  fw_uterus <- Tissue_fractions[9]
  fw_adipose <- Tissue_fractions[10]
  fw_skin <- Tissue_fractions[11]
  fw_muscles <- Tissue_fractions[12]
  fw_git <- Tissue_fractions[13]
  
  ### Calculation of tissue weights  
  W_tis[2] <- fw_heart*weight
  W_tis[3] <- fw_kidneys*weight
  W_tis[4] <- fw_brain*weight
  W_tis[5] <- fw_spleen*weight
  W_tis[6] <- fw_lungs*weight
  W_tis[7] <- fw_liver*weight
  W_tis[8] <- fw_uterus*weight
  W_tis[9] <- fw_skeleton*weight
  W_tis[10] <- fw_adipose*weight
  W_tis[11] <- fw_skin*weight
  W_tis[12] <- fw_muscles*weight
  W_tis[13] <- fw_git*weight
  
  for (i in 1:length(compartments)) {
    ###Calculation of tissue volumes
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
  }
  
  
  ### Calculations for "Soft tissue" compartment
  W_tis[1] <- weight - sum(W_tis[2:length(W_tis)], na.rm = TRUE)-Total_Blood
  V_tis[1] <- W_tis[1]     
  
  Vven=0.64*Total_Blood
  Vart=0.15*Total_Blood
  
  V_lu_is <- V_tis[6] 
  V_lu_cap <- V_cap[6] 
  V_rob_is <- sum(V_tis)-V_tis[6]
  V_rob_cap <- sum(V_cap)-V_cap[6]
  
  return(c("Q_blood_total"=Q_total,  "V_ven"=Vven, "V_art"=Vart,
           "V_lu_is" = V_lu_is, "V_rob_is" = V_rob_is,
           "V_lu_cap" = V_lu_cap, "V_rob_cap" = V_rob_cap))
}


#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function( dose){
  M_lu_cap<-0; M_lu_is<-0;M_lu_cell <- 0; M_lu_pc <- 0; M_rob_cap<-0; M_rob_is<-0; 
  M_rob_cell <- 0; M_rob_pc <- 0;
  M_ven <- dose; M_art<-0; M_excreta<-0
  
  return(c( "M_lu_cap" = M_lu_cap, "M_lu_is"=M_lu_is,"M_lu_cell" = M_lu_cell,
            "M_lu_pc" = M_lu_pc,
            "M_rob_cap"=M_rob_cap, "M_rob_is"=M_rob_is,"M_rob_cell" = M_rob_cell,
            "M_rob_pc" = M_rob_pc,
            
            "M_ven" = M_ven, "M_art" = M_art,"M_excreta" = M_excreta))
}
#==============
#3. ODEs System
#==============
Rat_model <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    #Estimate the interpolated weights
    if (time<=7*24){
      interpolated_weight = round(229+ time*(243-229)/(7*24))
    }else{
      interpolated_weight =  round(243 + (time-7*24)*(288-243)/((28-7)*24))
    }
    #Find the position of the interpolated weight
    which_weight <- which(all_weights == interpolated_weight)
    # Obtain physiological parameters based on interpolated weight
    parameters <- physiological_pars[[which_weight]]
    Q_blood_total <- parameters["Q_blood_total"]
    V_ven <- parameters["V_ven"]
    V_art <- parameters["V_art"]
    V_lu_is <-  parameters["V_lu_is"]
    V_lu_cap <-  parameters["V_lu_cap"]
    V_rob_is <- parameters["V_rob_is"]
    V_rob_cap <- parameters["V_rob_cap"]
    
    Q_total <- Q_blood_total*(1-Hct)
    QL_lu <- Q_total/500
    QL_rob <- Q_total/500
    # Concentrations (micro grams of NPs)/(mL tissue)
    C_lu_is <- M_lu_is/V_lu_is
    C_lu_cap <- M_lu_cap/V_lu_cap
    C_rob_is <- M_rob_is/V_rob_is
    C_rob_cap <- M_rob_cap/V_rob_cap
    C_ven <- M_ven/V_ven
    C_art <- M_art/V_art
    
    #Estimate lung's reflection coefficient
    a_lu <- np_size/lung_pore_size
    Phi_lu = (1-a_lu)^2
    F_lu <- (((1-a_lu^2)^(3/2))*Phi_lu)/(1+0.2*(a_lu^2)*(1-a_lu^2)^16)
    G_lu <- ((1- (2*a_lu^2)/3 - 0.20217*a_lu^5 )/ (1-0.75851*a_lu^5)) - (0.0431*(1-(1-a_lu^10)))
    sigma_lu <- 1-(1-(1-Phi_lu)^2)*G_lu+2*a_lu^2*Phi_lu*F_lu
    
    #Estimate rob reflection coefficient
    a_rob <- np_size/rob_pore_size
    Phi_rob = (1-a_rob)^2
    F_rob <- (((1-a_rob^2)^(3/2))*Phi_rob)/(1+0.2*(a_rob^2)*(1-a_rob^2)^16)
    G_rob <- ((1- (2*a_rob^2)/3 - 0.20217*a_rob^5 )/ (1-0.75851*a_rob^5)) - (0.0431*(1-(1-a_rob^10)))
    sigma_rob <- 1-(1-(1-Phi_rob)^2)*G_rob+2*a_rob^2*Phi_rob*F_rob
    
    # Lungs
    dM_lu_cap <-  Q_total*C_ven - (Q_total-QL_lu)*C_lu_cap - (1-sigma_lu)*QL_lu*C_lu_cap 
    dM_lu_is <- (1-sigma_lu)*QL_lu*C_lu_cap  - k_lu_in *M_lu_is + k_lu_out *M_lu_cell 
    dM_lu_cell <- k_lu_in *M_lu_is - k_lu_out *M_lu_cell
    dM_lu_pc <-  0
    
    #Rest of the body
    dM_rob_cap <-  Q_total*C_art - (Q_total-QL_rob)*C_rob_cap - (1-sigma_rob)*QL_rob*C_rob_cap -M_rob_cap*pc_endo
    dM_rob_is <- (1-sigma_rob)*QL_rob*C_rob_cap- (CLE*M_rob_is)/(Km+M_rob_is)
    dM_rob_cell <- 0
    dM_rob_pc <-  M_rob_cap*pc_endo
    
    # Venous Blood
    dM_ven <-(Q_total-QL_rob)*C_rob_cap + QL_rob*C_rob_is + QL_lu*C_lu_is- Q_total*C_ven
    
    # Arterial Blood
    dM_art <-  (Q_total-QL_lu)*C_lu_cap - Q_total*C_art
    
    # Excreta
    dM_excreta <- (CLE*M_rob_is)/(Km+M_rob_is)
    
    Whole_blood <- M_art+M_ven
    Lungs <- M_lu_cap + M_lu_is + M_lu_pc+ M_lu_cell
    Rob <- M_rob_cap + M_rob_is + M_rob_cell+M_rob_pc
    
    list(c("dM_lu_cap" = dM_lu_cap, "dM_lu_is" = dM_lu_is, 
           "dM_lu_cell" = dM_lu_cell,"dM_lu_pc" = dM_lu_pc,
           "dM_rob_cap" = dM_rob_cap, "dM_rob_is" = dM_rob_is, 
           "dM_rob_cell" = dM_rob_cell, 
           "dM_rob_pc" = dM_rob_pc,
           "dM_ven" = dM_ven,"dM_art" = dM_art,"dM_excreta" = dM_excreta), 
         "Whole_blood"=Whole_blood,"Lungs" = Lungs, "Rob" = Rob)
  })
}

obj_func <- function(x, dose, df, pars, metric = "AAFE"){
  BodyBurden <- c(df$lungs, df$rob, df$excreta, df$blood)
  parms <- c("rob_pore_size" = exp(x[1]), "lung_pore_size" = exp(x[2]),
             "CLE" = exp(x[3]),"Km" = exp(x[4]),
             "pc_endo" =  exp(x[5]) ,
             #"pc_exo" =  exp(x[6]), "Q_pc_sinusoids" =  exp(x[7]),"Q_pc_interstitial" =  exp(x[8]),
             # "Q_pc_lu" =  exp(x[9]),
             
             "Hct" = 0.45,  "np_size" = 6.55,
             'k_lu_in' = 0.400,  'k_lu_out' =0.0598,
             "k_rob_in" = 0.051, "k_rob_out" = 0.063,
             "physiological_pars" = pars$physiological_pars, 
             "all_weights" = pars$all_weights)
  sol_times <- c(seq(0,1, 0.001),seq(1.1,5, 0.1),  seq(6,28*24, 1))
  inits <- create.inits(unname(dose))
  solution <- data.frame(deSolve::ode(times = sol_times,  func = Rat_model,
                                      y = inits, parms = parms, method="bdf",
                                      rtol = 1e-4, atol = 1e-4))
  if(sum(solution$time %in% df$time) == length( df$time)){
    results <- c(solution[solution$time %in% df$time, "Lungs"],
                 solution[solution$time %in% df$time, "Rob"],
                 solution[solution$time %in% df$time, "M_excreta"],
                 solution[solution$time %in% df$time, "Whole_blood"])
  }else{
    results <-   BodyBurden *1000
    
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
setwd("C:/Users/user/Desktop/Nanoparticle_delivery_application")
data_concentration <- t(openxlsx::read.xlsx ('PEG-AU-NPs.xlsx', 
                                             sheet = "Kozics_concentration", 
                                             colNames = TRUE, rowNames = TRUE))
data_mass <- t(openxlsx::read.xlsx ('PEG-AU-NPs.xlsx', 
                                    sheet = "Kozics_mass", 
                                    colNames = TRUE, rowNames = TRUE)[1:5,])


dose_per_weight <- 0.7 #mg/kg
dose <- dose_per_weight*229
#The percentage of gold recovery in the analyzed samples
#detected 1 h after i.v. injection accounted for approximately 68.5% of the total injected
#dose.
#dose <- rowSums(data_mass)[1]/0.685 
df <- data.frame(time = c(1, 4, 24, 7*24, 28*24), lungs = unname(data_mass[,"lungs"]), 
                 rob = unname(data_mass[,"liver"]+ data_mass[,"spleen"] + 
                                data_mass[,"kidneys"]),  excreta = rep(NA, 5),
                 blood = unname(data_mass[,"blood"]))
df$excreta <- dose - (df$lungs+df$rob+df$blood)

# Vector of all possible interpolated weights
all_weights <- 229:288
physiological_pars <- list()
for(i in 1:length(all_weights)){
  physiological_pars[[i]] <- create.params(all_weights[i])
}

opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",NLOPT_LN_SBPLX , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-7, 
              "ftol_rel" = 1e-7,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 1000,
              "print_level" = 1)


# Define initial values of fitted parameters to provide to the optimization routine
x0 <-  c(log(10),log(8),5,1,1)#,-17,1,1,-5)
set.seed(435)
optimization<- nloptr::nloptr(x0 = x0,
                              eval_f = obj_func,
                              lb	=  c(log(6.551),log(6.551),-10,-10,-5),#-30,-10,-10,-10),
                              ub =   c(log(2500),log(100),7,7,7),#-5,7,7,2),
                              opts = opts,
                              dose = dose,
                              df = df,
                              metric = "SODI",
                              pars = list("all_weights" = all_weights,
                                          "physiological_pars" = physiological_pars ))


parms <- c("rob_pore_size" = exp(optimization$solution[1]), 
           "lung_pore_size" = exp(optimization$solution[2]),
           "CLE" = exp(optimization$solution[3]),
           "Km" =  exp(optimization$solution[4]),
           
           "pc_endo" =  exp(optimization$solution[5]) ,
           # "pc_exo" =  exp(optimization$solution[6]),
           # "Q_pc_sinusoids" =  exp(optimization$solution[7]) ,
           # "Q_pc_interstitial" =  exp(optimization$solution[8]),
           # "Q_pc_lu" =  exp(optimization$solution[9]) ,
           
           "k_rob_in" =  0.051 , 
           "k_rob_out" = 0.063,
           "Hct" = 0.45,"np_size" = 6.5,
           'k_lu_in' = 0.400,  'k_lu_out' =0.0598,
           "physiological_pars" = physiological_pars, 
           "all_weights" = all_weights)


sol_times <- seq(0,28*24, 1 )
inits <- create.inits(unname(dose))
solution <- data.frame(deSolve::ode(times = sol_times,  func = Rat_model,
                                    y = inits, parms = parms, method="lsodes",
                                    rtol = 1e-7, atol = 1e-7))
library(ggplot2)
# Lungs plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=Lungs ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=lungs ), size=5)+
  
  labs(title = "NP mass in Lungs", y = "Micrograms", x = "Time (hours)")+
  theme_bw()+
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
# Rob plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=Rob ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=rob ), size=5)+
  
  labs(title = "NP mass in Rob", y = "Micrograms", x = "Time (hours)")+
  theme_bw()+
  
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
  geom_line(data = solution, aes(x=time, y=M_excreta ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=excreta ), size=5)+
  
  labs(title = "NP mass in Excreta", y = "Micrograms", x = "Time (hours)")+
  theme_bw()+
  
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
# Blood plot
ggplot()+
  geom_line(data = solution, aes(x=time, y=Whole_blood ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=blood ), size=5)+
  
  labs(title = "NP mass in blood", y = "Micrograms", x = "Time (hours)")+
  theme_bw()+
  
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

###################
# Generate Data   #
###################
# times <- seq(0,28*24, 0.01)
# solution <- data.frame(deSolve::ode(times = times,  func = Rat_model,
#                                     y = inits, parms = parms, method="lsodes",
#                                     rtol = 1e-7, atol = 1e-7))
# select_times <-c(0.1, 0.5, 1, 2, 5, 10,20, 30, 40, 50,100,250,672)
# solution <- solution[solution$time %in% select_times,]
# 
# # Case 1: cell vicinity  is the whole organ
# case1 <- data.frame(time = solution$time, administration_site = solution$M_ven+solution$M_art, 
#                     cell_vicinity = solution$Lungs-solution$M_lu_cell,
#                     cell_interior = solution$M_lu_cell, 
#                     off_target = solution$Rob,
#                     excreta = solution$M_excreta)
# openxlsx::write.xlsx(case1, "case1.xlsx")
# 
# # Case 2: cell vicinity  is interstitial space
# case2 <- data.frame(time = solution$time, administration_site = solution$M_ven+solution$M_art, 
#                     cell_vicinity = solution$M_lu_is,
#                     cell_interior = solution$M_lu_cell, 
#                     off_target = solution$Rob+solution$M_lu_cap+solution$M_lu_pc,
#                     excreta = solution$M_excreta)
# openxlsx::write.xlsx(case1, "case2.xlsx")






















