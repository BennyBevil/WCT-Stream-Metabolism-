library(tidyverse)
library(FSA)
library(lubridate)

setwd("C:/Users/benny/OneDrive/Desktop/Masters.working")

#load data
CTTdata <- read.csv('CTT.Data.2025.F.csv')

CTTdata <- CTTdata %>%
  mutate(
    Date = dmy(Date),                # Convert String Date to Date Object
    Tag = as.character(Tag),         # Treat Tag as ID, not number
    Location = as.factor(Location)
  )

# Fit the log-log length-weight model
# Only uses rows where W and FL are both present and > 0
lw_mod <- lm(log(W) ~ log(FL),
             data = CTTdata,
             subset = !is.na(W) & !is.na(FL) & W > 0 & FL > 0)

# Predict weights for all rows
# predict() yields log(W), so we exponentiate to return to grams
CTTdata$Est.W <- exp(predict(lw_mod, newdata = CTTdata))

# Create completed weight column:
# Use observed W if it exists, otherwise use the estimated weight
CTTdata$C.W <- ifelse(!is.na(CTTdata$W), CTTdata$W, CTTdata$Est.W)


CTTdata <- CTTdata %>% 
  filter(!Location %in% c("Jake Canyon", "Alkali"))

# total length ~ fork length

lm_tl_fl <- lm(TL ~ FL,data=CTTdata)
summary(lm_tl_fl)


CTTdata <- CTTdata %>% mutate(logL = log(FL), logW = log(C.W))

library(lme4)
lmeLW <- lmer(logW ~ logL + (1|Location),data=CTTdata)
summary(lmeLW)

ranef(lmeLW)

CTTdata %>% 
  ggplot(aes(x=logL,y=logW)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_minimal()


### Length Frequency Histograms

CTTdata %>%
  ggplot(aes(x=FL)) +
  facet_wrap(~Location) +
  geom_histogram() +
  theme_minimal()

#### Growth

growth_CTT <- CTTdata %>%
  filter(!is.na(Tag)) %>%
  group_by(Tag) %>%
  filter(n() > 1) %>%
  arrange(Tag, Date) %>%
  mutate(
    Date_Diff = as.numeric(Date - lag(Date)),
    Growth_FL = FL - lag(FL),
    Growth_W  = C.W - lag(C.W),
    # Specific Growth Rate (SGR) = (ln(W2) - ln(W1)) / t * 100
    SGR = (log(C.W) - log(lag(C.W))) / Date_Diff * 100
  ) %>%
  filter(!is.na(Date_Diff)) %>%
  ungroup()

# Summary of Growth by Location
growth_summary <- growth_CTT %>%
  group_by(Location) %>%
  summarise(
    N_Recaps = n(),
    Avg_Growth_FL_mm = mean(Growth_FL, na.rm=TRUE),
    Avg_SGR_Percent = mean(SGR, na.rm=TRUE)
  )

print("--- Growth Summary by Location ---")
print(growth_summary)


# 5. Abundance & Biomass Estimation (Depletion Method) -------------------------
# We calculate Population (N) for each Location and Date using removal method

estimates <- CTTdata %>%
  group_by(Location, Date) %>%
  summarise(
    Passes = n_distinct(Pass),
    Total_Catch = n(),
    # Create a vector of catches per pass for the estimator
    Catch_Vector = list(tapply(Count, Pass, sum)),
    Mean_W_g = mean(C.W, na.rm = TRUE),
    .groups = "keep"
  ) %>%
  mutate(
    # Estimate N using FSA::removal (Carle-Strub or Zippin)
    # Only runs if Passes > 1, otherwise uses Total_Catch as minimum N
    Pop_Est = map(Catch_Vector, function(catch) {
      if(length(catch) > 1) {
        tryCatch({
          fit <- removal(catch, method="CarleStrub")
          return(fit$est["No"])
        }, error = function(e) return(sum(catch))) # Fallback to catch if error
      } else {
        return(sum(catch)) # If only 1 pass, N = Catch (conservative)
      }
    })
  ) %>%
  unnest(Pop_Est) %>%
  mutate(
    Biomass_g = Pop_Est * Mean_W_g
  )

# 6. Production Calculation ----------------------------------------------------
# Production (P) = Growth Rate (G) * Mean Biomass (B_bar)
# We assume 'Growth Rate' comes from the recaptured fish (Instantaneous G)
# If no recaptures, we can estimate G from change in mean weight of the population.

production_est <- estimates %>%
  group_by(Location) %>%
  arrange(Date) %>%
  filter(n() >= 2) %>% # Need at least 2 timepoints
  summarise(
    T1 = first(Date),
    T2 = last(Date),
    Interval_Days = as.numeric(last(Date) - first(Date)),
    B1 = first(Biomass_g),
    B2 = last(Biomass_g),
    Mean_B = (B1 + B2) / 2, # Linear approximation of mean biomass
    
    # Calculate G (Instantaneous Growth)
    # Option A: Use population mean weights (Cohort analysis)
    # G = (log(W_bar_2) - log(W_bar_1))
    Mean_W1 = first(Mean_W_g),
    Mean_W2 = last(Mean_W_g),
    G_Inst = log(Mean_W2) - log(Mean_W1),
    
    # Production
    Production_g = G_Inst * Mean_B,
    Production_g_per_day = Production_g / Interval_Days
  )

print("--- Production Estimates ---")
print(production_est)

# 7. Save Results --------------------------------------------------------------
write_csv(growth_summary, "Trout_Recapture_Growth.csv")
write_csv(production_est, "Trout_Production_Estimates.csv")

####Graph

# ==============================================================================
# Trout Growth and Production Visualization Script
# ==============================================================================

library(ggplot2)
library(patchwork) # Optional: for combining plots. Install if needed: install.packages("patchwork")

# Set a consistent theme
theme_set(theme_bw() + theme(text = element_text(size = 12)))


# ------------------------------------------------------------------------------
# 2. Recapture Growth Rates by Location (Boxplot)
# ------------------------------------------------------------------------------
# Visualizes the Specific Growth Rate (SGR) for recaptured fish
p2 <- ggplot(growth_CTT, aes(x = Location, y = SGR, fill = Location)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  labs(title = "Specific Growth Rates (SGR) by Location",
       subtitle = "Based on recaptured tagged fish",
       y = "SGR (% weight gain per day)",
       x = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# ------------------------------------------------------------------------------
# 3. Biomass Comparison (Start vs. End)
# ------------------------------------------------------------------------------
# Reshape production data to 'long' format for plotting before/after
biomass_long <- production_est %>%
  select(Location, B1, B2) %>%
  pivot_longer(cols = c(B1, B2), names_to = "Timepoint", values_to = "Biomass") %>%
  mutate(Timepoint = ifelse(Timepoint == "B1", "Start", "End"))

p3 <- ggplot(biomass_long, aes(x = Location, y = Biomass, fill = Timepoint)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Start" = "#A3A500", "End" = "#00BF7D")) +
  labs(title = "Estimated Population Biomass Change",
       y = "Total Biomass (g)",
       x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# ------------------------------------------------------------------------------
# 4. Total Production by Location
# ------------------------------------------------------------------------------
# Visualizes the total biomass produced over the interval
p4 <- ggplot(production_est, aes(x = reorder(Location, Production_g), y = Production_g)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(Production_g, 1)), hjust = -0.2, size = 3) +
  coord_flip() + # Horizontal bars for easier reading
  labs(title = "Total Trout Production Estimate",
       subtitle = "Production = G * Mean Biomass",
       y = "Production (g)",
       x = "Location") +
  ylim(min(production_est$Production_g) * 1.1, max(production_est$Production_g) * 1.15)

print(p4)

# ------------------------------------------------------------------------------
# Optional: Combine plots (requires 'patchwork' library)
# ------------------------------------------------------------------------------
# (p2 + p4) / (p3)

#### Depletion model single stream

library(R2jags)

sink("CTTDepletion_0_1.txt")
cat("
model {
    catch[1] ~ dbin(q,Np);
    
    Np ~ dpois(N);
    
    for(k in 2:Npass){
       catch[k] ~ dbin(q,(Np-sum(catch[1:(k-1)])));
    }
    
    q ~ dunif(0,1);
    N ~ dunif(minN,3000)
}
",fill = TRUE)
sink()

DepleteData <- CTTdata %>% 
  group_by(Location,Pass) %>% 
  summarise(Count = sum(Count),RL = median(RL),StreamW = median(StreamW))

DepleteData <- DepleteData %>%
  mutate(Location.Pass = paste(Location,Pass,sep='.')) %>%
  mutate(Count = ifelse(Location.Pass == 'Cow.2',Count+15,
                        ifelse(Location.Pass == 'Cow.3', 
                               Count - 15,Count)))

Stream_DepleteData <- DepleteData %>% filter(Location == 'CLT')

jags.params<-c('q','Np','N')
jags.data<-list(catch = Stream_DepleteData$Count,
                Npass = length(Stream_DepleteData$Pass),
                minN = sum(Stream_DepleteData$Count))

depletion1<- jags(jags.data, inits =NULL, parameters.to.save = jags.params,
                  model.file="CTTDepletion_0_1.txt",n.chains=3,
                  n.thin=10, n.iter=50000,n.burnin=25000,DIC=TRUE)

plot(depletion1)

depletion1$BUGSoutput$sims.list$Np %>% hist()
depletion1$BUGSoutput$sims.list$q %>% hist()


modelScript.name <- "DepletionMultiSite_Density_Biomass.txt"
jagsscript <- cat("
model {

    for(s in 1:Nstream){
      catch[s,1] ~ dbin(q[s],Np[s]);
      for(k in 2:Npass[s]){
        catch[s,k] ~ dbin(q[s],(Np[s]-sum(catch[s,(1:(k-1))])));
      }
    }
    
    lqm ~ dunif(-3,3);
    lqsd ~ dnorm(0,3) T(0,);
    for(s in 1:Nstream){
      Np[s] ~ dpois(N[s]);
      N[s] ~ dunif(minN[s],10000);
      lq[s] ~ dnorm(lqm,pow(lqsd,-2));
      q[s] <- exp(lq[s])/(1+exp(lq[s]));
      
      Biomass[s] <- Np[s] * AvgW[s];
      NDens[s] <- Np[s]/(SectionLength[s]*StreamWidth[s]);
      BDens[s] <- Biomass[s]/(SectionLength[s]*StreamWidth[s]);
    }
    
    
    
  
}  
", file = modelScript.name)

CatchMatrix <- DepleteData %>% xtabs(Count ~ Location + Pass,data=.) %>% as.matrix()
Npass_vec <- DepleteData %>% group_by(Location) %>% summarise(Npass = length(unique(Pass))) %>% pull(Npass)
Streams <- DepleteData %>% pull(Location) %>% unique()
Nstream <-  Streams %>% length()
minN_vec <- DepleteData %>% group_by(Location) %>% summarise(Count = sum(Count)) %>% pull(Count)
AvgW_vec <- CTTdata %>% group_by(Location) %>% summarise(AvgW = mean(W,na.rm=T)) %>% pull(AvgW)
SectionLengths <- DepleteData %>% group_by(Location) %>% summarise(RL = median(RL)) %>% pull(RL)
StreamWidths <- DepleteData %>% group_by(Location) %>% summarise(StreamW = median(StreamW)) %>% pull(StreamW)

jags.params<-c('q','Np','N','Biomass','NDens','BDens')
jags.data<-list(catch = CatchMatrix,
                Npass = Npass_vec,
                Nstream = Nstream,
                minN = minN_vec,
                AvgW = AvgW_vec,
                SectionLength = SectionLengths,
                StreamWidth = StreamWidths)

depletion2<- jags(jags.data, inits =NULL, parameters.to.save = jags.params,
                  model.file="DepletionMultiSite_Density_Biomass.txt",n.chains=3,
                  n.thin=10, n.iter=500000,n.burnin=250000,DIC=TRUE)

plot(depletion2)

## Plot
library(ggplot2)
library(dplyr)
library(tidyr)

# Extract posterior samples from JAGS object
mcmc_samples <- depletion2$BUGSoutput$sims.matrix
colnames(mcmc_samples)

# Identify columns corresponding to NDens
NDens_cols <- grep("^NDens\\[", colnames(mcmc_samples), value = TRUE)

# Convert to long format
NDens_df <- as.data.frame(mcmc_samples[, NDens_cols])
colnames(NDens_df) <- Streams#paste0("Stream_", 1:length(NDens_cols))

NDens_long <- NDens_df %>%
  mutate(iter = row_number()) %>%
  pivot_longer(cols=1:Nstream,names_to = "Stream",
               values_to = "NDens")

# Plot violin + boxplot for posterior densities
ggplot(NDens_long, aes(x = Stream, y = log(NDens))) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  theme_minimal() +
  labs(title = "Posterior Density (Fish/m²) Across Streams",
       x = "Stream",
       y = "Density (fish per m²)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####

BDens_cols <- grep("^BDens\\[", colnames(mcmc_samples), value = TRUE)

# Convert to long format
BDens_df <- as.data.frame(mcmc_samples[, BDens_cols])
colnames(BDens_df) <- Streams#paste0("Stream_", 1:length(NDens_cols))

BDens_df <- BDens_df[,order(apply(BDens_df,2,FUN=median))]

BDens_long <- BDens_df %>%
  mutate(iter = row_number()) %>%
  pivot_longer(cols=1:Nstream,names_to = "Stream",
               values_to = "BDens")

stream_order <- BDens_long %>%
  group_by(Stream) %>%
  summarise(mean_BDens = mean(BDens, na.rm = TRUE)) %>%
  arrange((mean_BDens)) %>%
  pull(Stream)

BDens_long <- BDens_long %>%
  mutate(Stream = factor(Stream, levels = stream_order))

# Plot violin + boxplot for posterior densities
ggplot(BDens_long, aes(x = Stream, y = log10(BDens))) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  theme_minimal() +
  labs(title = "Posterior Density (Fish/m²) Across Streams",
       x = "Stream",
       y = expression(Log[10]*"(Biomass Density, g/m"^2*")")
  )  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#write.csv(BDens_long, file = "FDens.csv")

# "Alkali"      "Buffalo"     "CCT"         "CLT"         "Cow"        
# # [6] "Dry"         "Jake Canyon" "Pintler"     "Plimpton"   
# TempC <- c(10.8,10.8,18.2,8.2,15.2,10.5,9.6,11.1,12.7)
# TempDF <- data.frame(Stream = Streams,Temp = TempC)
# 
# ggplot(BDens_long %>% left_join(TempDF,by='Stream'), aes(x = as.factor(Temp), y = log10(BDens))) +
#   geom_violin(fill = "skyblue", alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
#   theme_minimal() +
#   labs(title = "Posterior Density (Fish/m²) Across Streams",
#        x = "",
#        y = "Biomass Density (g per m²)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


########
###
# 
# modelScript.name <- "DepletionMultiSite_Density_Biomass.txt"
# jagsscript <- cat("
# model {
# 
#   # OBSERVATION MODEL
#   for(s in 1:Nstream){
#     N_remain[s,1] <- Np[s]
#     for(k in 2:maxNpass){
#       # If the stream had fewer passes, don't compute extra N_remain
#       N_remain[s,k] <- N_remain[s,k-1] - catch[s,k-1]
#     }
#     for(k in 1:Npass[s]){
#       catch[s,k] ~ dbin(q_sk[s,k], N_remain[s,k])
#     }
#   }
# 
#   # PASS-EFFECT MODEL: q[s,k] = logistic(logit_base[s] + delta[k])
#   for(s in 1:Nstream){
#     logit_base[s] ~ dnorm(lq_m, pow(lq_sd, -2))
#     
#     logit_q[s,1] <- logit_base[s];
#     q_sk[s,1] <- ilogit(logit_q[s,1])
#     
#     for(k in 2:maxNpass){
#       logit_q[s,k] <- logit_base[s] + delta[k-1]
#       q_sk[s,k] <- ilogit(logit_q[s,k])
#     }
#   }
# 
#   # Shared pass effect (e.g., decreasing efficiency)
#   for(k in 1:(maxNpass-1)){
#     delta[k] ~ dnorm(0, 0.001)  # centered prior for pass effect
#   }
# 
#   lq_m ~ dnorm(0, 1)
#   lq_sd ~ dnorm(0, 3) T(0,)
# 
#   # LATENT ABUNDANCE AND DERIVED QUANTITIES
#   for(s in 1:Nstream){
#     Np[s] ~ dpois(N[s])
#     N[s] ~ dunif(minN[s], 10000)
# 
#     Biomass[s] <- Np[s] * AvgW[s]
#     NDens[s] <- Np[s] / (SectionLength[s] * StreamWidth[s])
#     BDens[s] <- Biomass[s] / (SectionLength[s] * StreamWidth[s])
#   }
# }
# ", file = modelScript.name)
# 
# CatchMatrix <- DepleteData %>% xtabs(Count ~ Location + Pass,data=.) %>% as.matrix()
# Npass_vec <- DepleteData %>% group_by(Location) %>% summarise(Npass = length(unique(Pass))) %>% pull(Npass)
# Streams <- DepleteData %>% pull(Location) %>% unique()
# Nstream <-  Streams %>% length()
# minN_vec <- DepleteData %>% group_by(Location) %>% summarise(Count = sum(Count)) %>% pull(Count)
# AvgW_vec <- CTTdata %>% group_by(Location) %>% summarise(AvgW = median(W,na.rm=T)) %>% pull(AvgW)
# SectionLengths <- DepleteData %>% group_by(Location) %>% summarise(RL = median(RL)) %>% pull(RL)
# StreamWidths <- DepleteData %>% group_by(Location) %>% summarise(StreamW = median(StreamW)) %>% pull(StreamW)
# 
# jags.params<-c('q_sk','Np','N','Biomass','NDens','BDens','logit_base','delta')
# jags.data<-list(catch = CatchMatrix,
#                 Npass = Npass_vec,
#                 maxNpass = 3,
#                 Nstream = Nstream,
#                 minN = minN_vec,
#                 AvgW = AvgW_vec,
#                 SectionLength = SectionLengths,
#                 StreamWidth = StreamWidths)
# 
# depletion3<- jags(jags.data, inits =NULL, parameters.to.save = jags.params,
#                   model.file="DepletionMultiSite_Density_Biomass.txt",n.chains=3,
#                   n.thin=10, n.iter=500000,n.burnin=250000,DIC=TRUE)
# 
# plot(depletion3)
# 
# 
# ## Plot
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # Extract posterior samples from JAGS object
# mcmc_samples <- depletion3$BUGSoutput$sims.matrix
# colnames(mcmc_samples)
# 
# # Identify columns corresponding to NDens
# NDens_cols <- grep("^NDens\\[", colnames(mcmc_samples), value = TRUE)
# 
# # Convert to long format
# NDens_df <- as.data.frame(mcmc_samples[, NDens_cols])
# colnames(NDens_df) <- Streams#paste0("Stream_", 1:length(NDens_cols))
# 
# NDens_long <- NDens_df %>%
#   mutate(iter = row_number()) %>%
#   pivot_longer(cols=1:Nstream,names_to = "Stream",
#                values_to = "NDens")
# 
# # Plot violin + boxplot for posterior densities
# ggplot(NDens_long, aes(x = Stream, y = log(NDens))) +
#   geom_violin(fill = "skyblue", alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
#   theme_minimal() +
#   labs(title = "Posterior Density (Fish/m²) Across Streams",
#        x = "Stream",
#        y = "Density (fish per m²)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# ####
# 
# BDens_cols <- grep("^BDens\\[", colnames(mcmc_samples), value = TRUE)
# 
# # Convert to long format
# BDens_df <- as.data.frame(mcmc_samples[, BDens_cols])
# colnames(BDens_df) <- Streams#paste0("Stream_", 1:length(NDens_cols))
# 
# BDens_long <- BDens_df %>%
#   mutate(iter = row_number()) %>%
#   pivot_longer(cols=1:Nstream,names_to = "Stream",
#                values_to = "BDens")
# 
# # Plot violin + boxplot for posterior densities
# ggplot(BDens_long, aes(x = Stream, y = log(BDens))) +
#   geom_violin(fill = "skyblue", alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
#   theme_minimal() +
#   labs(title = "Posterior Density (Fish/m²) Across Streams",
#        x = "Stream",
#        y = "Biomass Density (g per m²)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
# 

modelScript.name <- "DepletionMultiSite_Density_Biomass.txt"
jagsscript <- cat("
model {
  
  # ---------------------------
  # Priors on detection model
  # ---------------------------
  alpha ~ dnorm(0, 0.1)         # Intercept (logit scale)
  beta ~ dnorm(0, 0.1)          # Length effect (logit scale)
  
  # ---------------------------
  # Stream-level priors
  # ---------------------------
  for(i in 1:Nstream){
    psi[i] ~ dbeta(1, 1)        # Inclusion probability for augmented individuals
  }
  
  # ---------------------------
  # Fish-level model
  # ---------------------------
  for(i in 1:Nstream){
    for(j in 1:M[i]){
      
      # Latent inclusion indicator (1 = real fish, 0 = pseudo-fish)
      z[i,j] ~ dbern(psi[i])
      
      # Length-based detection probability
      logit(p[i,j]) <- alpha + beta * Length[i,j]
      
      # Observation: only real fish can be detected
      y[i,j] ~ dbern(p[i,j] * z[i,j])
    }
    
    # Derived abundance
    N[i] <- sum(z[i,1:M[i]])
  }
  
}
", file = modelScript.name)

# Observed fish data (example)
obs_data <- CTTdata %>% arrange(Location,Pass)

# Set augmentation
max_augment <- 207 + 100  # or pick max(n_obs) + buffer
streams <- unique(obs_data$Location)
Nstream <- length(streams)

# Initialize data containers
M <- rep(max_augment, Nstream)
y_mat <- matrix(0, nrow = Nstream, ncol = max_augment)
Length_mat <- matrix(NA, nrow = Nstream, ncol = max_augment)

for (i in 1:Nstream) {
  fish_i <- obs_data %>% filter(Location == streams[i])
  n_obs <- sum(fish_i$Count)
  
  # Observed fish
  y_mat[i, 1:n_obs] <- 1
  Length_mat[i, 1:n_obs] <- fish_i$FL[1:n_obs]
  
  # Augmented: fill with sampled or mean length
  if (n_obs < max_augment) {
    Length_mat[i, (n_obs + 1):max_augment] <- rnorm(max_augment - n_obs,
                                                    mean = mean(log(fish_i$FL),na.rm=T),
                                                    sd = sd(log(fish_i$FL),na.rm=T)) %>% exp()
  }
}

jags_data <- list(
  y = y_mat,
  Length = Length_mat,
  M = M,
  Nstream = Nstream
)





