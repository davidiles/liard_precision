rm(list=ls())

# Required packages
my.packs <- c(
  'tidyverse',   # For data manipulation/plotting
  'reshape2',    # For data manipulation
  'jagsUI',      # For Bayesian analysis
  'lme4',        # For Frequentist analysis
  'viridis',     # For plotting
  'cowplot',#,     # For plotting
  #'optimx'       # For fitting models
  'parallel',
  'doParallel'
)

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}

lapply(my.packs, require, character.only = TRUE) # Load packages


#---------------------------------------------------------
#---------------------------------------------------------
# Liard precision analysis

# Steps: 
# 1. Extract parameter estimates for a species
# 2. Simulate historical and future data at liard on particular sampling schedules
# 3. Estimate trend up to that point and record precision

#---------------------------------------------------------
#---------------------------------------------------------

species_to_run <- c("CHSP","ALFL","GRAJ")

#----------------------
# Simulation parameters
nrep = 10 # Repeat the process XX times (probably want to do this 100 times at least)
max.year = 30            # Final year of sampling
interannual.schedule = 5 # Sampling every X years
n.visit = 2              # Number of visits per season
#----------------------


#----------------------
# Parallel processing setup
numCores <- detectCores() # Detect number of cores on machine
numCores <- numCores - 2  # Reserve 2 cores for other tasks
registerDoParallel(numCores) # Register cores to use for parallel processing
#----------------------

all_results = foreach(yy = 1:length(species_to_run), .combine = rbind, .packages = my.packs) %dopar% {
  
  current.species <- species_to_run[yy]
  
  species_results <- data.frame()
  
  # Load data
  dat = read.csv(paste0("../data/",current.species,".csv"))
  
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
  
  #---------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------
  # Conduct precision analysis on new simulated data
  #---------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------
  
  # Final year of sampling at liard up to this point
  final.year.at.liard = max(data.stand$year)
  
  # Future years of sampling under the proposed sampling scheme
  additional.years = seq((final.year.at.liard+interannual.schedule),max.year,interannual.schedule)
  
  # Create new dataframe to store new simulated data
  newdata = expand.grid(round = 1:n.visit,
                        stand = unique(data.stand$stand),
                        year = additional.years,
                        nstation = 3)
  
  newdata$count = NA
  newdata$year_abs = newdata$year + min(data.stand$year_abs)
  newdata$standyear = paste0(newdata$stand,newdata$year)
  newdata$observation = NA
  newdata = rbind(data.stand, newdata)
  newdata$observation = 1:nrow(newdata)
  
  n.years = max(newdata$year+1)
  n.stands = max(newdata$stand)
  n.obs = max(newdata$observation)
  year_sequence = rev(unique(newdata$year)[which(unique(newdata$year) >= final.year.at.liard)])
  
  #-------------------------------------------------------------
  # Repeatedly generate new data and analyze it
  #-------------------------------------------------------------
  
  #------------
  # If results file already exists, load it and continue progress
  #------------
  # Check if species_results already exists
  if (file.exists(paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""))){
    species_results = read.csv(paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""))
    completed.reps = max(species_results$rep)
  } else{
    completed.reps <- 0
  }
  
  # If progress is incomplete, continue processing from where the script left off
  if (completed.reps < nrep){
    
    # Simulate new data
    for (rep in (completed.reps + 1):nrep){
      
      newdata$count = NA
      
      # Resample existing data in random order (preserving stand-level variation)
      year.assignments = sample(unique(data.stand$year),n.years, replace = TRUE)
      stand.assignments = sample(unique(data.stand$stand),n.stands, replace = TRUE)
      
      newdata_years = unique(newdata$year)
      newdata_stands = unique(newdata$stand)
      
      for (i in 1:length(newdata_years)){
        for (j in 1:length(unique(newdata$stand))){
          
          assigned.data = subset(data.stand, stand == stand.assignments[j] & year == year.assignments[i])
          
          for (v in 1:n.visit){
            
            newdata[which(newdata$stand == newdata_stands[j] & newdata$year == newdata_years[i] & newdata$round == v),c("nstation","count")] = subset(assigned.data,round == v)[,c("nstation","count")]
            
          }
          
        }
      }
      
      # Fit the statistical model to the simulated counts and estimate precision
      # Work backwards from the last simulated year (max.year) to the last year of liard sampling (final.year.at.liard)
      for (last.year in year_sequence){
        
        current.data = subset(newdata, year <= last.year)
        
        # Eliminate sites where species was never observed
        stands.observed = aggregate(count~stand, data = current.data, FUN = sum) %>% subset(., count > 0)
        
        current.data = subset(current.data, stand %in% stands.observed$stand)
        
        # Bit of code that ensures the script doesn't break if it encounters an error message
        intercept = trend = spatial.sd = temporal.sd = residual.sd = ci.width.new = NA
        tryCatch({
          fit.glmer.new = glmer(count ~ year + (1|stand) + (1|year) + (1|observation)+
                                  offset(log(nstation)), 
                                data = current.data, 
                                family = "poisson", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
          
          trend.ci.new = confint(fit.glmer.new, "year")
          ci.width.new = trend.ci.new[2] - trend.ci.new[1]
          
          intercept = as.numeric(fixef(fit.glmer.new)[1]) # intercept on log scale
          trend = as.numeric(fixef(fit.glmer.new)[2])     # trend on log scale
          
          # Extract estimates of spatial, temporal, and residual variance
          spatial.sd = sqrt(as.numeric(summary(fit.glmer.new)$varcor$`stand`))
          temporal.sd = sqrt(as.numeric(summary(fit.glmer.new)$varcor$`year`))
          residual.sd = sqrt(as.numeric(summary(fit.glmer.new)$varcor$`observation`))
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        # Check if species_results already exists
        if (file.exists(paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""))){
          species_results = read.csv(paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""))
        }
        
        
        species_results = rbind(species_results,data.frame(species = current.species,
                                                           max.year = max.year,
                                                           interannual.schedule = interannual.schedule,
                                                           n.visit = n.visit,
                                                           rep = rep,
                                                           n.stands.observed.at = nrow(stands.observed),
                                                           total.n = nrow(current.data),
                                                           year.max = last.year,
                                                           year.abs = last.year + min(data.stand$year_abs),
                                                           years.into.future = last.year - final.year.at.liard,
                                                           ci.width = ci.width.new, # Calculate mean precision across the 100 repetitions
                                                           intercept = intercept,
                                                           trend = trend,
                                                           spatial.sd = spatial.sd,
                                                           temporal.sd = temporal.sd,
                                                           residual.sd = residual.sd))
        
        
        # Save results
        write.csv(species_results, file = paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""), row.names = FALSE)
        
      }# end year sequence loop
      
      
    } # end reps loop
    
  }
  
} # end species (parallelized) loop



#------------------------------------------
# Once species runs are complete, load all data and merge into a single dataframe (called all_results)
#------------------------------------------
all_results <- data.frame()
for (current.species in species_to_run){
  species_results <- read.csv(paste0("../results/",current.species,"_results_",nrep,"_",max.year,"_",interannual.schedule,"_",n.visit,".csv",sep=""))
  all_results <- rbind(all_results, species_results)
}

# Remove rows with infinite and missing values
all_results[all_results == Inf] = NA
all_results <- na.omit(all_results)

#------------------------------------------
# Plot all results for all species
#------------------------------------------

result_plot = ggplot(data = all_results) +
  
  #Results for each simulated dataset
  geom_point(aes(x = year.abs, y = ci.width, col = factor(rep)), alpha = 0.5)+
  geom_line(aes(x = year.abs, y = ci.width, col = factor(rep)), alpha = 0.5)+
  
  #Average results
  geom_line(data = aggregate(ci.width ~ year.abs + species, data = all_results, FUN = mean, na.rm = TRUE), aes(x = year.abs, y = ci.width), col = "black")+
  
  geom_hline(yintercept = 0.035, linetype = 2)+ # Arbitrary precision threshold
  
  theme_bw()+
  scale_color_manual(values = rep("gray50",length(unique(all_results$rep))), guide = FALSE)+
  xlab("Year")+
  facet_wrap(species~.)
print(result_plot)


#------------------------------------------
# Extract and plot average confidence interval widths
#------------------------------------------
summary_aggregated <- aggregate(ci.width~year.abs + species, data = all_results, FUN = mean, na.rm = TRUE)
aggregated_plot <- ggplot(data = summary_aggregated) +
  geom_line(aes(x = year.abs, y = ci.width))+
  geom_hline(yintercept = 0.035, linetype = 2)+ # Arbitrary precision threshold
  theme_bw()+
  xlab("Year")+
  facet_wrap(species~.)

print(aggregated_plot)
