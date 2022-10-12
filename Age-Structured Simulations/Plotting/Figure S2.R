library(ggplot2)
library(grid)
library(gridExtra)

# load data
load("ABC_data.RData")

# priors for initial-year demographic parameters
a_r <- 1
b_r <- 1

a_s_juv <- 3
b_s_juv <- 3

a_s_ad <- 2
b_s_ad <- 0.4

a_b <- 1.25
b_b <- 1

accepted_ind <- which(keep_bin == 1)

df <- data.frame(r = r_mean[accepted_ind], s_juv = s_juv_mean[accepted_ind], s_ad = s_ad_mean[accepted_ind], b = b_mean[accepted_ind])
df_priors <- data.frame(r_prior = rbeta(1e5, a_r, b_r), s_juv_prior = rbeta(1e5, a_s_juv, b_s_juv), s_ad_prior = rbeta(1e5, a_s_ad, b_s_ad), b_prior = rbeta(1e5, a_b, b_b))

p1 <- ggplot(df, mapping = aes(x=r)) + 
  geom_histogram(mapping = aes(y = ..density..), colour = "black", fill = "#638475", bins = 40) +
  geom_density(df_priors, mapping = aes(x=r_prior, y = ..density..), colour = "#7B8CDE", lwd = 1.5) + ylab("density")

p2 <- ggplot(df, mapping = aes(x=s_juv)) + 
  geom_histogram(mapping = aes(y = ..density..), colour = "black", fill = "#FF5666", bins = 40) +
  geom_density(df_priors, mapping = aes(x=s_juv_prior, y = ..density..), colour = "#7B8CDE", lwd = 1.5) + ylab("")

p3 <- ggplot(df, mapping = aes(x=s_ad)) + 
  geom_histogram(mapping = aes(y = ..density..), colour = "black", fill = "#FFB86F", bins = 40) +
  geom_density(df_priors, mapping = aes(x=s_ad_prior, y = ..density..), colour = "#7B8CDE", lwd = 1.5) + ylab("density")

p4 <- ggplot(df, mapping = aes(x=b)) + 
  geom_histogram(mapping = aes(y = ..density..), colour = "black", fill = "#3E2F5B", bins = 40) +
  geom_density(df_priors, mapping = aes(x=s_juv_prior, y = ..density..), colour = "#7B8CDE", lwd = 1.5) + ylab("")

grid.arrange(p1, p2, p3, p4, nrow = 2, top = textGrob("Accepted Demographic Rates with Priors", gp = gpar(fontsize = 15)))
