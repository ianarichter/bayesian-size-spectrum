# Script to parameterize Bayesian size spectrum model 

#### Prep ####
library(tidyverse)
library(reshape)
library(runjags)
library(R2jags)
library(mcmcplots)
library(broom)
library(EnvStats)
library(FSA)
library(boot)
library(ggpubr)

seed <- 111
set.seed(seed)

#### Import data ####
sitedata <- read.csv("samplefishdata.csv")

# important model coefficients for size bias correction model (Richter et al. 2022)
b0.summary <- read.csv("b0.summary.csv")
beta.summary <- read.csv("beta.summary.csv")
b0.median <- b0.summary[2,2]
b.size <- beta.summary[5,2] 

#### Data ####
# raw data requires some manipulation to obtain the data necessary for model parameterization
# fish are assigned to log2 size class bins
intervals <- c(2^1, 2^1.5, 2^2, 2^2.5, 2^3, 2^3.5, 2^4, 2^4.5, 2^5, 2^5.5, 2^6,
               2^6.5, 2^7, 2^7.5, 2^8, 2^8.5, 2^9, 2^9.5, 2^10, 2^10.5, 2^11, 2^11.5, 2^12)
sitedata <- sitedata %>%
  dplyr::mutate(lencat = lencat(Weight, breaks = intervals),
                lencat_log2 = log2(lencat))

# calculate the total catches for each size bin along with the lower, upper, and midpoints of the size bins
sitedata_binned <- sitedata %>%
  group_by(lencat_log2, .drop = TRUE) %>%
  dplyr::summarise(catch = n()) %>%
  dplyr::mutate(lencat = 2^(lencat_log2), # lencat function from FSA package returns the lowest value for each size
                midpoint_log2 = lencat_log2 + 0.25,
                midpoint = 2^(midpoint_log2),
                upper_log2 = lencat_log2 + 0.5,
                upper = 2^(upper_log2),
                lencat_width = upper - lencat)

# calculate capture probability and corrected catches for each size bin
sitedata_binned <- sitedata_binned %>%
  dplyr::mutate(q_logit = b0.median + b.size * midpoint %>%
                  as.numeric(),
                q = inv.logit(q_logit))

# number of size bins
N_bins <- nrow(sitedata_binned) %>%
  as.numeric()

#### Model data ####
# include all of the necessary data in a list
df <- sitedata_binned
data <- list("catch" = df$catch, "N" = N_bins,
             "catch_total" = sum(df$catch), "q" = df$q,
             "lower" = df$lencat, "upper" = df$upper)

#### JAGS Model Parameters ####
# Set parameter values for the sampling algorithm
nc <- 3 # number of chains
ni <- 10000 # number of iterations/chain length
nb <- 1000 # burn in value
nt <- 20 # thinning rate

#### Size spectrum model file ####
# create a text file that will be used by jags
sink("sizespectrum.model.jags")
cat("model {
  
  # size bin follows a multinomial distribution
  catch[1:N] ~ dmulti(p[1:N], catch_total)
  
  # calculate probability of fish belonging to each size bin
  for (i in 1:N) {
  p_bin[i] <- (q[i]*(lower[i]^(-1 * alpha) - upper[i]^(-1 * alpha)))
  }
  
  # obtain total value to standardize the probabilities
  p_total <- sum(p_bin[])
  
  # standardization of probabilities
  for (j in 1:N) {
  p[j] <- (q[j]*(lower[j]^(-1 * alpha) - upper[j]^(-1 * alpha)))/p_total
  }
  
  alpha ~ dunif(0,50)

  }")
sink()

#### Run model ####
# identify the parameter (alpha) that the sampler will monitor across iterations
params <- c("alpha")
# run the model using the R2jags package
jagsfit <- R2jags::jags(data = data, parameters.to.save = params, n.chains = nc,
                        n.iter = ni, n.thin = nt, n.burnin = nb, 
                        model.file = "sizespectrum.model.jags")
# check model output
print(jagsfit)
output <- as.data.frame(jagsfit$BUGSoutput$summary)

#### Model output plot ####
# use the slope parameter (alpha) to create the size spectrum model
# create a dataframe including all size bins
dist <- seq(from = min(sitedata_binned$lencat), to = max(sitedata_binned$upper), by = 1) %>%
  as.data.frame() %>%
  `colnames<-`(c("Weight")) %>%
  dplyr::mutate(pdf = dpareto(Weight, shape = output["alpha",1], location = min(sitedata$Weight)),
                pdf_uci = dpareto(Weight, shape = output["alpha",7], location = min(sitedata$Weight)),
                pdf_lci = dpareto(Weight, shape = output["alpha",3], location = min(sitedata$Weight))) %>%
  dplyr::filter(pdf > 0) # omit any individuals with a 0 probability


(prob_density_plots <- ggplot(dist, aes(x = log2(Weight), y = log2(pdf))) +
    geom_line(size = 2) +
    labs(x = expression(Log[2]*" Weight"), y = expression(Log[2]*" Probability Density")) +
    stat_regline_equation(aes(label =  paste(after_stat(eq.label))), label.x.npc = 0.6, label.y.npc = 0.85, size = 6, show.legend = FALSE) +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(size = 16))
)

# Size spectrum model with 95% CI values for alpha
# sigma (location) parameter is fairly consistent
(prob_density_plots_faceted_ci <- ggplot(dist, aes(x = log2(Weight), y = log2(pdf))) +
    geom_line(size = 2) +
    geom_line(aes(y = log2(pdf_lci)), size = 2, linetype = 'dashed', colour = 'black') +
    geom_line(aes(y = log2(pdf_uci)), size = 2, linetype = 'dashed', colour = 'black') +
    labs(x = expression(Log[2]*" Weight"), y = expression(Log[2]*" Probability Density")) +
    stat_regline_equation(aes(label =  paste(..eq.label..)), label.x.npc = 0.2, label.y.npc = 0.4, size = 6, show.legend = FALSE) +
    scale_x_continuous(limits = c(min(log2(dist$Weight)), max(log2(dist$Weight)))) +
    scale_y_continuous(limits = c(min(log2(dist$pdf_uci)), max(log2(dist$pdf_uci)))) +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(size = 16),
          strip.background = element_rect(colour = 'black', fill = 'white'),
          strip.text = element_text(size = 14)))
