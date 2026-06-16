library(tidyverse)

#load data
CTTdata <- read.csv('CTT.Data.2025.csv')

# total length ~ fork length
lm_tl_fl <- lm(TL ~ FL,data=CTTdata)
summary(lm_tl_fl)


CTTdata <- CTTdata %>% mutate(logL = log(FL), logW = log(W))

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
       y = "Biomass Density (g per m²)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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





