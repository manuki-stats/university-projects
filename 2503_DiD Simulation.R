############
#
#  Script:  2503_DiD Simulation.R
#
#  Authors: Manuel Soffici
# 
#  Purpose: Simulation for own paper (Research seminar on ML in Econometrics
#           2024/25): "Learning and Identification in DiD designs". Supervisor:
#           Prof. Daniel Wilhelm at LMU Munich.
#
#  Created: 2025-03-05
#
#  Version: 2025-03-16
# 
############


## Clean

rm(list = ls())
set.seed(153)

library(dplyr)
library(tidyr)
library(ggplot2)


## Learning from/about untreated outcomes ####

alpha <- 1200     # intercept for wages
lambda <- 100     # linear time trend
tau <- 700        # ATT
beta <- 0.5       # discount rate
rho <- 0.5        # learning rate
K1 <- 70000       # cost of first treatment in t=1 (-> force sharp design)
K01 <- 700        # cost of first treatment in t=2
K11 <- 700        # cost of second treatment
K0 <- 0           # cost of no treatment
N <- 10000        # population/sample size

mean_epsilon <- 1200
sd_epsilon <- 300
effort0 <- rnorm(N, mean = mean_epsilon, sd = sd_epsilon)     # individual initial efforts

sim_data <- data.frame(
  effort0 = effort0,
  Y0_untreated = alpha + effort0,
  Y0_treated = alpha + tau + effort0
)

summary(sim_data)     # check

bias_tauHat0 <- 0     # bias in treatment effect as predicted by workers
sd_tauHat0 <- 0       # variability in predictions
sim_data$tauHat0 <- rnorm(N, mean = tau + bias_tauHat0, sd = sd_tauHat0)

# Compute predicted outcomes for Y0

sim_data$Y0Hat_untreated <- alpha + sim_data$effort0
sim_data$Y0Hat_treated <- alpha + sim_data$tauHat0 + sim_data$effort0

summary(sim_data)        # check

# Compute predicted outcomes for Y1 (with information U0)

sim_data$Y1Hat_untreated_U0 <- sim_data$Y0Hat_untreated + lambda
sim_data$Y1Hat_treated_U0 <- sim_data$Y0Hat_treated + lambda

# Calculate utilities based on the given formulas

sim_data$utility_00_U0 <- sim_data$Y0Hat_untreated - K0 +
  beta * (sim_data$Y1Hat_untreated_U0 - K0)

sim_data$utility_01_U0 <- sim_data$Y0Hat_untreated - K0 +
  beta * (sim_data$Y1Hat_treated_U0 - K01)

sim_data$utility_10_U0 <- sim_data$Y0Hat_treated - K1 +
  beta * (sim_data$Y1Hat_untreated_U0 - K0)

sim_data$utility_11_U0 <- sim_data$Y0Hat_treated - K1 +
  beta * (sim_data$Y1Hat_treated_U0 - K11)

summary(sim_data)       # check

# Choose treatment and update predictions

# treatment0 = 0 if either utility_00_U0 or utility_01_U0 is larger than both utility_10_U0 and utility_11_U0, 1 otherwise
sim_data$treatment0 <- ifelse(
  (sim_data$utility_00_U0 > sim_data$utility_10_U0 & sim_data$utility_00_U0 > sim_data$utility_11_U0) |
    (sim_data$utility_01_U0 > sim_data$utility_10_U0 & sim_data$utility_01_U0 > sim_data$utility_11_U0),
  0, 1
)

# if treatment0==0, observed outcome is Y0_untreated; if treatment0==1, it's Y0_treated.
sim_data$Y0_observed <- ifelse(sim_data$treatment0 == 0,
                                      sim_data$Y0_untreated,
                                      sim_data$Y0_treated)

# tauHat1 is simply equal to tauHat0 in this scenario:
sim_data$tauHat1 <- sim_data$tauHat0

# adjust effort by adding rho*effort0; otherwise keep it like this
sim_data$effort1 <- (1 + rho) * sim_data$effort0

summary(sim_data)        # check

# Computing new potential outcomes (modified effort)

sim_data$Y1_untreated = alpha + lambda + sim_data$effort1
sim_data$Y1_treated = alpha + lambda + tau + sim_data$effort1

# Adding period 1 predictions

sim_data$Y1Hat_untreated_U1 <- alpha + lambda + sim_data$effort1
sim_data$Y1Hat_treated_U1 <- alpha + lambda + sim_data$tauHat1 + sim_data$effort1

# Compute utilities for period 1

sim_data$utility_X0_U1 <- sim_data$Y1Hat_untreated_U1 - K0

sim_data$utility_X1_U1 <- ifelse(
  sim_data$treatment0 == 1,
  sim_data$Y1Hat_treated_U1 - K11,
  sim_data$Y1Hat_treated_U1 - K01
)

summary(sim_data)        # check

# Adding treatment1 and observed Y1

# treatment1 is set to 0 if utility_X0_U1 exceeds utility_X1_U1; otherwise, it's 1.
sim_data$treatment1 <- ifelse(sim_data$utility_X0_U1 > sim_data$utility_X1_U1, 0, 1)

# if treatment1 equals 0, use Y1_untreated; if treatment1 equals 1, use Y1_treated.
sim_data$Y1_observed <- ifelse(sim_data$treatment1 == 0,
                                      sim_data$Y1_untreated,
                                      sim_data$Y1_treated)

summary(sim_data)        # check

# Estimating the ATT with the DiD estimator

att_DiD <- (mean(sim_data$Y1_observed[sim_data$treatment1 == 1]) -
                 mean(sim_data$Y0_observed[sim_data$treatment1 == 1])) -
  (mean(sim_data$Y1_observed[sim_data$treatment1 == 0]) -
     mean(sim_data$Y0_observed[sim_data$treatment1 == 0]))

cat("Direct DiD estimate of ATT:", att_DiD, "\n")

# Bootstrap CI

B <- 1000
boot_ATT <- numeric(B)

for(b in 1:B){
  boot_sample <- sim_data[sample(1:nrow(sim_data), replace = TRUE), ]
  boot_ATT[b] <- (mean(boot_sample$Y1_observed[boot_sample$treatment1 == 1]) -
                    mean(boot_sample$Y0_observed[boot_sample$treatment1 == 1])) -
    (mean(boot_sample$Y1_observed[boot_sample$treatment1 == 0]) -
       mean(boot_sample$Y0_observed[boot_sample$treatment1 == 0]))
}

CI <- quantile(boot_ATT, probs = c(0.025, 0.975))
cat("95% CI for ATT: ", CI, "\n")

# Restrict to individuals who were untreated in period 0
plot_data <- sim_data %>% 
  filter(treatment0 == 0) %>%
  mutate(group = ifelse(treatment1 == 1, "Switchers", "Never-treated"))

# Compute group-level averages for untreated potential outcomes in period 0 and period 1
plot_data_long <- plot_data %>%
  group_by(group) %>%
  summarise(
    Y0_avg = mean(Y0_untreated),
    Y1_avg = mean(Y1_untreated)
  ) %>%
  pivot_longer(cols = c(Y0_avg, Y1_avg), names_to = "period", values_to = "avg_outcome")

# Rename the period variable for better labeling
plot_data_long <- plot_data_long %>%
  mutate(period = ifelse(period == "Y0_avg", "t=0", "t=1"))

# Create the plot
ggplot(plot_data_long, aes(x = period, y = avg_outcome, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Untreated potential outcomes (group average over time)",
       x = "Time",
       y = "Untreated outcomes (average)",
       color = "Group") +
  theme_minimal()

###

# Compute group-level averages for observed outcomes in period 0 and period 1
plot_data_long <- plot_data %>%
  group_by(group) %>%
  summarise(
    Y0_avg = mean(Y0_observed),
    Y1_avg = mean(Y1_observed)
  ) %>%
  pivot_longer(cols = c(Y0_avg, Y1_avg), names_to = "period", values_to = "avg_outcome")

# Rename the period variable for better labeling
plot_data_long <- plot_data_long %>%
  mutate(period = ifelse(period == "Y0_avg", "t=0", "t=1"))

# Create the plot
ggplot(plot_data_long, aes(x = period, y = avg_outcome, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Observed outcomes (group average over time)",
       x = "Time",
       y = "Observed outcomes (average)",
       color = "Group") +
  theme_minimal()
