# Required packages
my.packs <- c(
  'tidyverse',   # For data manipulation/plotting
  'reshape2',    # For data manipulation
  'jagsUI',      # For Bayesian analysis
  'lme4',        # For Frequentist analysis
  'viridis',     # For plotting
  'cowplot'#,     # For plotting
  #'optimx'       # For fitting models
)

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}

lapply(my.packs, require, character.only = TRUE) # Load packages

rm(list=ls())

#------------------------------------------------------------------------------------
# Load data and prepare for analysis
#------------------------------------------------------------------------------------

all.summary = data.frame()

for (current.species in "ALFL"){#c("ALFL","AMRO","CHSP","GRAJ","PHVI","REVI","WETA")){ 
  
  # Load data
  dat = read.csv(paste0("./data/",current.species,".csv"))
  
  #-------------------------------------
  # Restructure data for stand-level analysis
  #-------------------------------------
  data.stand = expand.grid(  round = sort(unique(dat$Round)),
                             stand = sort(unique(dat$StandNumber)),
                             
                             year = sort(unique(dat$YEAR)))
  
  data.stand$nstation = NA # Will store number of stations that were visited on each visit (round) in each year
  data.stand$count = NA    # Will store summed count
  
  # Loop through the round/stand/year combinations and sum up the counts from all stations
  for (i in 1:nrow(data.stand)){
    
    i.year = data.stand$year[i]
    i.stand = data.stand$stand[i]
    i.round = data.stand$round[i]
    
    # Extract relevant data for this stand/year/round combination
    i.data = subset(dat, YEAR == i.year & StandNumber == i.stand & Round == i.round)
    
    if (nrow(i.data) == 0) next
    
    # Number of stations surveyed on this visit
    data.stand$nstation[i] = length(unique(i.data$Station))
    
    # Total count on this visit (sum of counts at all stations)
    data.stand$count[i] = sum(i.data$Abundance.x)
  }
  
  #-------------------------------------
  # Analyze stand-level data
  #-------------------------------------
  data.stand$year_abs = data.stand$year
  data.stand$year = data.stand$year-min(data.stand$year)
  data.stand$standyear = as.numeric(as.factor(paste0(data.stand$stand,data.stand$year)))
  data.stand$observation = 1:nrow(data.stand)
  
  # Remove stands where the species was never observed
  stand.max.count = aggregate(count ~ stand, data = data.stand, FUN = max)
  stand.observed = subset(stand.max.count, count > 0)
  glmer.data = subset(data.stand, stand %in% stand.observed$stand)
  glmer.data = subset(glmer.data, !is.na(count))
  
  # Fit the model using maximum likelihood
  intercept = trend = spatial.sd = temporal.sd = residual.sd = NA
  tryCatch({
  fit.glmer = glmer(count ~ year + (1|stand) + (1|year) + (1|observation) +
                      offset(log(nstation)), 
                    data = glmer.data, family = "poisson", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  
  # Extract estimates of intercept and trend
  intercept = as.numeric(fixef(fit.glmer)[1]) # intercept on log scale
  trend = as.numeric(fixef(fit.glmer)[2])     # trend on log scale
  
  # Extract estimates of spatial, temporal, and residual variance
  spatial.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`stand`))
  temporal.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`year`))
  residual.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`observation`))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # Extract confidence interval (precision of trend estimate)
  trend.ci = ci.width = NA
  tryCatch({
    
    trend.ci = confint(fit.glmer,"year")
    ci.width = trend.ci[2] - trend.ci[1]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # #------------------------------------------------------------------------------------
  # # Goodness-of-fit assessment (use a 'parametric boostrap')
  # #    See Kery and Royle's 2016 book, pp. 196-197.
  # #       Steps are:
  # #          1. Compute some useful discrepancy measure for actual dataset
  # #          2. Generate simulated dataset using fitted parameters
  # #          3. Refit the model to the new simulated data
  # #          4. Calculate discrepancy measure for simulated dataset
  # #          5. Repeat many times, evaluate where "actual" discrepancy measures falls
  # #------------------------------------------------------------------------------------
  # 
  # # Step 1. Compute some relevant discrepancy measures for actual dataset
  # observed = glmer.data$count
  # expected = predict(fit.glmer, type = "response")
  # sqerror.obs = (observed - expected)^2
  # chi2.obs = sum((observed-expected)^2 / (expected + 0.0001)) # fit statistic 1
  # max.obs  = max(observed)                                    # fit statistic 2
  # se.obs = summary(fit.glmer)$coefficients[2,2]
  # trend.obs = summary(fit.glmer)$coefficients[2,1]
  # 
  # # Steps 2-5. Generate simulated dataset using estimated parameters
  # simrep = 1
  # chi2.sim = max.sim = rep(NA,simrep)
  # for (i in min(which(is.na(chi2.sim))):simrep){
  #   set.seed(i)
  #   observed.sim = rpois(n = length(observed), lambda = expected)
  # 
  #   tryCatch({
  #     fit.glmer.sim = glmer(observed.sim ~ year + (1|stand) + (1|year) + (1|observation) +
  #                             offset(log(nstation)), data = glmer.data, family = "poisson", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  #   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # 
  #   expected.sim = predict(fit.glmer.sim, type = "response")
  #   chi2.sim[i] = sum((observed.sim-expected.sim)^2 / (expected.sim + 0.0001)) # fit statistic 1
  #   max.sim[i] = max(observed.sim) # fit statistic 2
  # 
  #   # Update progress and plot histogram of chi-square fit statistics every 50 iterations
  #   cat(paste0(100*round(i/simrep,3),"% complete \n"))
  #   if (i%%50 == 0){
  #     hist(chi2.sim, breaks = 50, col = "gray80", xlim = range(c(chi2.sim,chi2.obs), na.rm = TRUE),
  #          xlab = "Chi-square goodness of fit statistic", main = paste0("Parametric boostrap results for ",
  #                                                                       current.species,"\nc-hat =",round(chi2.obs / mean(chi2.sim,na.rm=TRUE),3)))
  #     abline(v = chi2.obs, col = "blue", lwd = 2)
  #   }
  # 
  # }
  # 
  # chat = chi2.obs / mean(chi2.sim,na.rm=TRUE)            # Mean overdispersion
  # chi2.perc = mean(chi2.obs > chi2.sim, na.rm = TRUE)    # Probability that observed test statistic is greater than simulated test statistic (ideally will equal 0.5)
  
  # Store relevant parameters for this species in a dataframe
  species.summary = data.frame (species = current.species,
                                
                                # Goodness-of-fit statistics
                                #chat = chat,
                                #chi2.perc = chi2.perc,
                                
                                intercept = intercept,
                                trend = trend,
                                spatial.sd = spatial.sd,
                                temporal.sd = temporal.sd,
                                residual.sd = residual.sd,
                                trend.lci = trend.ci[1],
                                trend.uci = trend.ci[2],
                                ci.width = ci.width
  )
  
  all.summary = rbind(all.summary,species.summary)
  
  write.csv(all.summary, file = "all_species_summary.csv", row.names = FALSE)
}

print(all.summary)

#---------------------------------------------------------
#---------------------------------------------------------
# Liard precision analysis
#---------------------------------------------------------
#---------------------------------------------------------

# Steps: 
# 1. Extract parameter estimates for a species
# 2. Simulate historical and future data at liard on particular sampling schedules
# 3. Estimate trend up to that point and record precision

summary_data = data.frame()
for (current.species in c("ALFL")){ 
  
  # Load data
  dat = read.csv(paste0("./data/",current.species,".csv"))
  
  #-------------------------------------
  # Restructure data for stand-level analysis
  #-------------------------------------
  data.stand = expand.grid(  round = sort(unique(dat$Round)),
                             stand = sort(unique(dat$StandNumber)),
                             
                             year = sort(unique(dat$YEAR)))
  
  data.stand$nstation = NA # Will store number of stations that were visited on each visit (round) in each year
  data.stand$count = NA    # Will store summed count
  
  # Loop through the round/stand/year combinations and sum up the counts from all stations
  for (i in 1:nrow(data.stand)){
    
    i.year = data.stand$year[i]
    i.stand = data.stand$stand[i]
    i.round = data.stand$round[i]
    
    # Extract relevant data for this stand/year/round combination
    i.data = subset(dat, YEAR == i.year & StandNumber == i.stand & Round == i.round)
    
    if (nrow(i.data) == 0) next
    
    # Number of stations surveyed on this visit
    data.stand$nstation[i] = length(unique(i.data$Station))
    
    # Total count on this visit (sum of counts at all stations)
    data.stand$count[i] = sum(i.data$Abundance.x)
  }
  
  #-------------------------------------
  # Analyze stand-level data
  #-------------------------------------
  data.stand$year_abs = data.stand$year
  data.stand$year = data.stand$year-min(data.stand$year)
  data.stand$standyear = as.numeric(as.factor(paste0(data.stand$stand,data.stand$year)))
  data.stand$observation = 1:nrow(data.stand)
  
  # Remove stands where the species was never observed
  stand.max.count = aggregate(count ~ stand, data = data.stand, FUN = max)
  stand.observed = subset(stand.max.count, count > 0)
  glmer.data = subset(data.stand, stand %in% stand.observed$stand)
  glmer.data = subset(glmer.data, !is.na(count))
  
  #***************************************************************************************************************************
  # Note: Could also just pull the parameters from the all.summary file created in part 1 here, rather than run lines 238-262
  #***************************************************************************************************************************
  # Fit the model using maximum likelihood
  intercept = trend = spatial.sd = temporal.sd = residual.sd = NA
  tryCatch({
    fit.glmer = glmer(count ~ year + (1|stand) + (1|year) + (1|observation) +
                        offset(log(nstation)), 
                      data = glmer.data, family = "poisson", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    
    # Extract estimates of intercept and trend
    intercept = as.numeric(fixef(fit.glmer)[1]) # intercept on log scale
    trend = as.numeric(fixef(fit.glmer)[2])     # trend on log scale
    
    # Extract estimates of spatial, temporal, and residual variance
    spatial.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`stand`))
    temporal.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`year`))
    residual.sd = sqrt(as.numeric(summary(fit.glmer)$varcor$`observation`))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # Extract confidence interval (precision of trend estimate)
  trend.ci = ci.width = NA
  tryCatch({
    
    trend.ci = confint(fit.glmer,"year")
    ci.width = trend.ci[2] - trend.ci[1]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #***************************************************************************************************************************
  # The code above could be commented out, and replaced by just reading from the all.summary file if it's already been created
  #***************************************************************************************************************************
  
  
  #--------------
  # Note that up to this point, the analysis is identical to the code in the first part of this script
  #--------------
  
  #---------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------
  # Conduct precision analysis on new simulated data
  #---------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------
  
  max.year = 50            # Final year of sampling
  interannual.schedule = 5 # Sampling every X years
  n.visit = 2              # Number of visits per season
  
  # Final year of sampling at liard up to this point
  final.year.at.liard = max(glmer.data$year)
  
  # Future years of sampling under the proposed sampling scheme
  additional.years = seq((final.year.at.liard+interannual.schedule),max.year,interannual.schedule)
  
  # Create new dataframe to store new simulated data
  newdata = expand.grid(round = 1:n.visit,
                        stand = unique(glmer.data$stand),
                        year = additional.years,
                        nstation = 3)
  
  newdata$count = NA
  newdata$year_abs = newdata$year + min(glmer.data$year_abs)
  newdata$standyear = paste0(newdata$stand,newdata$year)
  newdata$observation = NA
  newdata = rbind(glmer.data, newdata)
  newdata$observation = 1:nrow(newdata)
  
  n.years = max(newdata$year+1)
  n.stands = max(newdata$stand)
  n.obs = max(newdata$observation)
  year_sequence = rev(unique(newdata$year)[which(unique(newdata$year) >= final.year.at.liard)])
  
  #-------------------------------------------------------------
  # Repeatedly generate new data and analyze it
  #-------------------------------------------------------------
  nrep = 100 # Repeat the process XX times (probably want to do this 100 times at least)
  
  # Matrix to store results (width of confidence interval) and convergence code for each iteration
  convergence.matrix = ci.width.matrix = matrix(NA, nrow = nrep, ncol = length(year_sequence)) 
  for (reps in 1:nrep){
    
    year.effects = rnorm(n.years,0,temporal.sd)
    stand.effects = rnorm(n.stands, 0, spatial.sd)
    residual.effects = rnorm(n.obs,0,residual.sd)
    
    # Fill in count data with new simulated data
    newdata$lambda = exp(intercept + trend*newdata$year + stand.effects[newdata$stand] + year.effects[newdata$year+1] + residual.effects[newdata$observation])
    newdata$count = rpois(nrow(newdata),newdata$lambda)
    
    # Fit the statistical model to the simulated counts and estimate precision
    # Work backwards from the last simulated year (max.year) to the last year of liard sampling (final.year.at.liard)
    for (last.year in year_sequence){
      
      # Bit of code that ensures the script doesn't break if it encounters an error message
      ci.width.new = NA
      tryCatch({
        fit.glmer.new = glmer(count ~ year + (1|stand) + (1|year) + (1|observation)+
                                offset(log(nstation)), 
                              data = subset(newdata, year <= last.year), 
                              family = "poisson", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
        
        trend.ci.new = confint(fit.glmer.new, "year")
        ci.width.new = trend.ci.new[2] - trend.ci.new[1]
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      
      # Store precision estimates in matrix
      ci.width.matrix[reps, which(year_sequence == last.year)] = ci.width.new
      
    }# end year sequence loop
    
    print(paste0(round(100*reps/nrep,3),"% complete for ", current.species))
  } # end reps loop
  
  #Summary of precision analysis for this species
  summary_data = rbind(summary_data,data.frame(species = current.species,
                                               max.year = max.year,
                                               interannual.schedule = interannual.schedule,
                                               n.visit = n.visit,
                                               nreps = nrep,
                                               ci.width = apply(ci.width.matrix,2,mean, na.rm=TRUE), # Calculate mean precision across the 100 repetitions
                                               year.max = year_sequence,
                                               year.abs = year_sequence + min(glmer.data$year_abs),
                                               years.into.future = year_sequence - final.year.at.liard))
  
} # end species loop


#---------------------------------------------
# Plot resulting precision curves for this run of the precision analysis
#---------------------------------------------

ggplot(data = summary_data) +
  
  geom_line(aes(x = year.abs, y = ci.width, col = species))+
  geom_hline(yintercept = 0.035, linetype = 2)+ # Arbitrary precision threshold
  geom_text(data = subset(summary_data, year.max == max(summary_data$year.max)), aes(x = year.abs, y = ci.width, label = species, col = species))+
  
  theme_bw()

