

# for mapppd site CROZ
site = "CROZ"
load("CROZ_ADPE_obs.RData")

# number of years
years_obs <- site_ADPE_obs$year
year <- min(years_obs):max(years_obs)
num_years <- length(year)

# create count vector with NAs
count <- c()
for (yr in year){
  if (yr %in% years_obs){
    count <- c(count, site_ADPE_obs$count[site_ADPE_obs$year == yr])
  } else{
    count <- c(count, NA)
  }
}
Data <- as.data.frame(cbind(year, count))
colnames(Data) <- cbind("year", "count")
Data$count <- as.numeric(Data$count)

# define truncated normal function
rtruncnorm <- function(n, mean, sd) {
  min <- 0
  max <- 1
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  qnorm(u, mean, sd)
}

# number of runs
num_runs <- 300000 # need ~ 100,000 to get 1000 w/ acceptance rate of 0.01 (1%)

###### priors for demographic parameters (means)
a_r <- 1
b_r <- 1

a_s_juv <- 3
b_s_juv <- 3

a_s_ad <- 2
b_s_ad <- 0.4

a_b <- 1.25
b_b <- 1

# wiggle demographic parameters around selected value for each year
# truncated normal (0,1) with sigma value:
sigma_r <- 0.025
sigma_s_juv <- 0.025
sigma_s_ad <- 0.025
sigma_b <- 0.025

# initialize integers
num_skipped <- 0 ###################################

# initialize vectors
r_mean <- matrix(NA, nrow = num_runs)
s_juv_mean <- matrix(NA, nrow = num_runs)
s_ad_mean <- matrix(NA, nrow = num_runs)
b_mean <- matrix(NA, nrow = num_runs)

r_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
r_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
r_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
r_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

s_juv <- matrix(NA, nrow = num_runs, ncol = num_years)

s_2 <- matrix(NA, nrow = num_runs, ncol = num_years)
s_3 <- matrix(NA, nrow = num_runs, ncol = num_years)
s_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
s_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
s_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
s_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

b_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
b_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
b_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
b_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

# initialize matrices
C_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
C_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
C_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
C_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

N_1 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_2 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_3 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
N_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

B_4 <- matrix(NA, nrow = num_runs, ncol = num_years)
B_5 <- matrix(NA, nrow = num_runs, ncol = num_years)
B_6 <- matrix(NA, nrow = num_runs, ncol = num_years)
B_7p <- matrix(NA, nrow = num_runs, ncol = num_years)

B_total <- matrix(NA, nrow = num_runs, ncol = num_years) # total number of breeders for each run for each year (for plotting)


################### SIMULATE TIME SERIES

# start the clock
ptm <- proc.time()

for (i in 1:num_runs){
  #
  # on every run, draw values from priors for initial year
  #
  # reproductive success PER EGG
  r_mean[i] <- rbeta(1, a_r, b_r) ########################################
  r_4[i,1] <- r_mean[i] # make rbeta(1, 1, 1) if you want each r_i different
  r_5[i,1] <- r_mean[i]
  r_6[i,1] <- r_mean[i]
  r_7p[i,1] <- r_mean[i]
  # don't need initial year survival, but calculate the means to wiggle rest around
  # juvenile survival
  s_juv_mean[i] <- rbeta(1, a_s_juv, b_s_juv) ####################################
  # adult survival
  s_ad_mean[i] <- rbeta(1, a_s_ad, b_s_ad) ########################################
  # breeding probability
  b_mean[i] <- rbeta(1, a_b, b_b) ########################################
  b_4[i,1] <- b_mean[i] # make rbeta(1, 1, 1) if you want each b_i different
  b_5[i,1] <- b_mean[i]
  b_6[i,1] <- b_mean[i]
  b_7p[i,1] <- b_mean[i]
  
  #
  # calculate initial year populations (year 1)
  # 
  initial_count <- Data$count[1] # 987 (B_total)
  B_4[i,1] <- round(initial_count/4)
  B_5[i,1] <- round(initial_count/4)
  B_6[i,1] <- round(initial_count/4)
  B_7p[i,1] <- round(initial_count/4)
  B_total[i,1] <- sum(as.numeric(B_4[i,1]), as.numeric(B_5[i,1]), as.numeric(B_6[i,1]), as.numeric(B_7p[i,1])) # ~987
  #
  N_2[i,1] <- round((initial_count/4)/b_4[i,1])
  N_3[i,1] <- round((initial_count/4)/b_4[i,1])
  #
  N_4[i,1] <- round(B_4[i,1]/b_4[i,1])
  N_5[i,1] <- round(B_5[i,1]/b_5[i,1])
  N_6[i,1] <- round(B_6[i,1]/b_6[i,1])
  N_7p[i,1] <- round(B_7p[i,1]/b_7p[i,1])
  #
  # check if N_i is too large (close to int overflow), if it is, don't finish this simulation
  if (N_2[i,1] > 1e7){
    num_skipped <- num_skipped + 1
    next # move onto next simuation in outer for loop
  }
  
  #
  # egg 1 + egg 2
  C_4[i,1] <- sum(as.numeric(rbinom(2, size=B_4[i,1], prob=r_4[i,1])))
  C_5[i,1] <- sum(as.numeric(rbinom(2, size=B_5[i,1], prob=r_5[i,1])))
  C_6[i,1] <- sum(as.numeric(rbinom(2, size=B_6[i,1], prob=r_6[i,1])))
  C_7p[i,1] <- sum(as.numeric(rbinom(2, size=B_7p[i,1], prob=r_7p[i,1])))
  #
  N_1[i,1] <- sum(as.numeric(C_4[i,1]), as.numeric(C_5[i,1]), as.numeric(C_6[i,1]), as.numeric(C_7p[i,1]))
  
  
  # simulation for num_years, starting at year 2 (initial year is year 1)
  for (j in 2:num_years){
    #
    # every year, draw demographic parameter values from betas
    #
    # reproductive success PER EGG
    r <- rtruncnorm(1, r_mean[i], sigma_r) ##################################
    r_4[i,j] <- r # change if you want each r_i different
    r_5[i,j] <- r
    r_6[i,j] <- r
    r_7p[i,j] <- r
    # juvenile survival
    s_juv[i,j] <- rtruncnorm(1, s_juv_mean[i], sigma_s_juv) ###########################
    # adult survival
    s_ad <- rtruncnorm(1, s_ad_mean[i], sigma_s_ad) ###############################
    s_2[i,j] <- s_ad # change if you want each s_ad different
    s_3[i,j] <- s_ad
    s_4[i,j] <- s_ad
    s_5[i,j] <- s_ad
    s_6[i,j] <- s_ad
    s_7p[i,j] <- s_ad
    # breeding probability
    b <- rtruncnorm(1, b_mean[i], sigma_b) ##################################
    b_4[i,j] <- b # make change if you want each b_i different
    b_5[i,j] <- b
    b_6[i,j] <- b
    b_7p[i,j] <- b
    
    N_2[i,j] <- round(0.5*rbinom(1, size=N_1[i,j-1], prob=s_juv[i,j]))
    N_3[i,j] <- rbinom(1, size=N_2[i,j-1], prob=s_2[i,j])
    N_4[i,j] <- rbinom(1, size=N_3[i,j-1], prob=s_3[i,j])
    N_5[i,j] <- rbinom(1, size=N_4[i,j-1], prob=s_4[i,j])
    N_6[i,j] <- rbinom(1, size=N_5[i,j-1], prob=s_5[i,j])
    N_7p[i,j] <- sum(as.numeric(rbinom(1, size=N_6[i,j-1], prob=s_6[i,j])), as.numeric(rbinom(1, size=N_7p[i,j-1], prob=s_7p[i,j])))
    
    B_4[i,j] <- rbinom(1, size=N_4[i,j-1], prob=b_4[i,j])
    B_5[i,j] <- rbinom(1, size=N_5[i,j-1], prob=b_5[i,j])
    B_6[i,j] <- rbinom(1, size=N_6[i,j-1], prob=b_6[i,j])
    B_7p[i,j] <- rbinom(1, size=N_7p[i,j-1], prob=b_7p[i,j])
    
    B_total[i,j] <- sum(as.numeric(B_4[i,j]), as.numeric(B_5[i,j]), as.numeric(B_6[i,j]), as.numeric(B_7p[i,j]))
    
    # egg 1 + egg 2
    C_4[i,j] <- sum(as.numeric(rbinom(2, size=B_4[i,j-1], prob=r_4[i,j])))
    C_5[i,j] <- sum(as.numeric(rbinom(2, size=B_5[i,j-1], prob=r_5[i,j])))
    C_6[i,j] <- sum(as.numeric(rbinom(2, size=B_6[i,j-1], prob=r_6[i,j])))
    C_7p[i,j] <- sum(as.numeric(rbinom(2, size=B_7p[i,j-1], prob=r_7p[i,j])))
    N_1[i,j] <- sum(as.numeric(C_4[i,j]), as.numeric(C_5[i,j]), as.numeric(C_6[i,j]), as.numeric(C_7p[i,j]))
  }
  
}


################### ACCEPT / REJECT

# MAPE threshold (chosen to arrive at reasonable acceptance rate)
threshold <- 35

keep_bin <- rep(NA, num_runs)
kept <- 0

for (i in 1:num_runs){
  # after each run is finished, accept or reject
  #
  # mean absolute percent error
  MAPE_i <- mean(abs((Data$count-B_total[i,])/Data$count), na.rm=TRUE) * 100
  
  # we're already keeping all parameter values of all accepted OR rejected runs, so we just need to keep a binary yes/no if that numbered run was accepted (accepted = 1, rejected = 0)
  if (anyNA(B_total[i,])){ # if the simulation was stopped (any B_total are NA), reject it
    keep_bin[i] <- 0
  } else if (MAPE_i <= threshold) { # otherwise... accept if MAPE is within threshold
    keep_bin[i] <- 1
    kept <- kept + 1
  } else { # reject if MAPE is not within threshold
    keep_bin[i] <- 0
  }
}
# calculate the acceptance rate
acc_rate <- kept/num_runs # number of runs accepted / number of total runs so far
# want to shoot for an acceptance rate of 0.01




# calculate skipped rate
skipped_rate <- num_skipped/num_runs

# Stop the clock
proc.time() - ptm

# save data so you don't have to run this every time
save(B_total, keep_bin, acc_rate, skipped_rate, r_mean, s_juv_mean, s_ad_mean, b_mean, r_4, r_5, r_6, r_7p, s_juv, s_2, s_3, s_4, s_5, s_6, s_7p, b_4, b_5, b_6, b_7p, Data, file = "ABC_data.RData")

