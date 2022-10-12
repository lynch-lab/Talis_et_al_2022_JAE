
# load count data
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

num_runs_each <- 10


# load in the "known" mean for each demographic parameter (from the ABC)
load("ABC_parameter_med.RData")

r_mean <- r_med
s_juv_mean <- s_juv_med
s_ad_mean <- s_ad_med
b_mean <- b_med


# initialize variables/vectors/arrays
num_steps <- 15 # number of steps / values of sigma for each demographic parameter

# bound on sigma
sigma_steps <- seq(0.001, 0.1, length.out = num_steps+2)[2:(num_steps+1)] # values of standard deviation of each demographic parameter to be stepped through

year_forecast <- 1982:2060 # years of our simulation
num_years_sim <- length(year_forecast)

r_realized <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim)) # [h,i,j,k, n, year]
s_juv_realized <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim)) # NOTE: year 1 will always stay NA, since we don't need an initial year survival
s_ad_realized <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim)) # NOTE: year 1 will always stay NA, since we don't need an initial year survival
b_realized <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim))

B_total <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim))
N_adults_total <- array(NA, dim = c(num_steps, num_steps, num_steps, num_steps, num_runs_each, num_years_sim))

# define truncated normal function
rtruncnorm <- function(n, mean, sd) {
  min <- 0
  max <- 1
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  qnorm(u, mean, sd)
}


# start the clock
ptm <- proc.time()

for (h in 1:num_steps){ # r
  sigma_r <- sigma_steps[h]
  
  for (i in 1:num_steps){ # s_juv
    sigma_s_juv <- sigma_steps[i]
    
    for (j in 1:num_steps){ # s_ad
      sigma_s_ad <- sigma_steps[j]
      
      for (k in 1:num_steps){ # b
        sigma_b <- sigma_steps[k]
        
        for (n in 1:num_runs_each){ # n
          
          ###########################################################################
          # INNER LOOP
          
          # initialize vectors
          B_4 <- rep(NA, num_years_sim)
          B_5 <- rep(NA, num_years_sim)
          B_6 <- rep(NA, num_years_sim)
          B_7p <- rep(NA, num_years_sim)
          N_1 <- rep(NA, num_years_sim)
          N_2 <- rep(NA, num_years_sim)
          N_3 <- rep(NA, num_years_sim)
          N_4 <- rep(NA, num_years_sim)
          N_5 <- rep(NA, num_years_sim)
          N_6 <- rep(NA, num_years_sim)
          N_7p <- rep(NA, num_years_sim)
          C_4 <- rep(NA, num_years_sim)
          C_5 <- rep(NA, num_years_sim)
          C_6 <- rep(NA, num_years_sim)
          C_7p <- rep(NA, num_years_sim)
          
          #
          ###############################
          #
          # initial year populations (year 1)
          # initial year parameter draws (but only for r & b, we don't need survival in the initial year)
          r_realized[h,i,j,k,n,1] <- rtruncnorm(1, r_mean, sigma_r)
          r_4 <- r_realized[h,i,j,k,n,1] # for the age classes w/ diff param values
          r_5 <- r_realized[h,i,j,k,n,1] # change these
          r_6 <- r_realized[h,i,j,k,n,1]
          r_7p <- r_realized[h,i,j,k,n,1]
          # breeding probability
          b_realized[h,i,j,k,n,1] <- rtruncnorm(1, b_mean, sigma_b)
          b_4 <- b_realized[h,i,j,k,n,1] # change for diff param values
          b_5 <- b_realized[h,i,j,k,n,1] # change these
          b_6 <- b_realized[h,i,j,k,n,1]
          b_7p <- b_realized[h,i,j,k,n,1]
          
          # 
          initial_count <- Data$count[1]
          B_4[1] <- round(initial_count/4)
          B_5[1] <- round(initial_count/4)
          B_6[1] <- round(initial_count/4)
          B_7p[1] <- round(initial_count/4)
          #
          B_total[h,i,j,k,n,1] <- sum(as.numeric(B_4[1]), as.numeric(B_5[1]), as.numeric(B_6[1]), as.numeric(B_7p[1]))
          #
          N_2[1] <- round(initial_count/4)
          N_3[1] <- round(initial_count/4)
          #
          N_4[1] <- round(B_4[1]/b_4)
          N_5[1] <- round(B_5[1]/b_5)
          N_6[1] <- round(B_6[1]/b_6)
          N_7p[1] <- round(B_7p[1]/b_7p)
          #
          # check if N_i is too large (close to int overflow), if it is, don't finish this simulation
          if (N_2[1] > 1e7){
            num_skipped <- num_skipped + 1
            next # move onto next simulation in outer for loop
          }
          #
          N_adults_total[h,i,j,k,n,1] <- sum(as.numeric(N_4[1]), as.numeric(N_5[1]), as.numeric(N_6[1]), as.numeric(N_7p[1]))
          
          #
          # egg 1 + egg 2
          C_4[1] <- sum(as.numeric(rbinom(2, size=B_4[1], prob=r_4)))
          C_5[1] <- sum(as.numeric(rbinom(2, size=B_5[1], prob=r_5)))
          C_6[1] <- sum(as.numeric(rbinom(2, size=B_6[1], prob=r_6)))
          C_7p[1] <- sum(as.numeric(rbinom(2, size=B_7p[1], prob=r_7p)))
          #
          N_1[1] <- sum(as.numeric(C_4[1]), as.numeric(C_5[1]), as.numeric(C_6[1]), as.numeric(C_7p[1]))
          
          # simulation for num_years, starting at year 2 (initial year is year 1)
          for (yr in 2:num_years_sim){
            # draw parameter values (from beta with realized parameters)
            # reproductive success PER EGG
            r_realized[h,i,j,k, n,yr] <- rtruncnorm(1, r_mean, sigma_r)
            r_4 <- r_realized[h,i,j,k, n,yr] # for the age classes to have diff param values
            r_5 <- r_realized[h,i,j,k, n,yr] # change these each to rbeta(1, a_r_realized, b_r_realized)
            r_6 <- r_realized[h,i,j,k, n,yr]
            r_7p <- r_realized[h,i,j,k, n,yr]
            # juvenile survival
            s_juv_realized[h,i,j,k, n,yr] <- rtruncnorm(1, s_juv_mean, sigma_s_juv)
            s_juv <- s_juv_realized[h,i,j,k, n,yr]
            # adult survival
            s_ad_realized[h,i,j,k, n,yr] <- rtruncnorm(1, s_ad_mean, sigma_s_ad)
            s_2 <- s_ad_realized[h,i,j,k, n,yr] # change for diff param values
            s_3 <- s_ad_realized[h,i,j,k, n,yr]
            s_4 <- s_ad_realized[h,i,j,k, n,yr]
            s_5 <- s_ad_realized[h,i,j,k, n,yr]
            s_6 <- s_ad_realized[h,i,j,k, n,yr]
            s_7p <- s_ad_realized[h,i,j,k, n,yr]
            # breeding probability
            b_realized[h,i,j,k, n,yr] <- rtruncnorm(1, b_mean, sigma_b)
            b_4 <- b_realized[h,i,j,k, n,yr] # change for diff param values
            b_5 <- b_realized[h,i,j,k, n,yr]
            b_6 <- b_realized[h,i,j,k, n,yr]
            b_7p <- b_realized[h,i,j,k, n,yr]
            
            N_2[yr] <- round(0.5*rbinom(1, size=N_1[yr-1], prob=s_juv))
            N_3[yr] <- rbinom(1, size=N_2[yr-1], prob=s_2)
            N_4[yr] <- rbinom(1, size=N_3[yr-1], prob=s_3)
            N_5[yr] <- rbinom(1, size=N_4[yr-1], prob=s_4)
            N_6[yr] <- rbinom(1, size=N_5[yr-1], prob=s_5)
            N_7p[yr] <- sum(as.numeric(rbinom(1, size=N_6[yr-1], prob=s_6)), as.numeric(rbinom(1, size=N_7p[yr-1], prob=s_7p)))
            
            B_4[yr] <- rbinom(1, size=N_4[yr-1], prob=b_4)
            B_5[yr] <- rbinom(1, size=N_5[yr-1], prob=b_5)
            B_6[yr] <- rbinom(1, size=N_6[yr-1], prob=b_6)
            B_7p[yr] <- rbinom(1, size=N_7p[yr-1], prob=b_7p)
            
            B_total[h,i,j,k, n,yr] <- sum(as.numeric(B_4[yr]), as.numeric(B_5[yr]), as.numeric(B_6[yr]), as.numeric(B_7p[yr]))
            
            # egg 1 + egg 2
            C_4[yr] <- sum(as.numeric(rbinom(2, size=B_4[yr-1], prob=r_4)))
            C_5[yr] <- sum(as.numeric(rbinom(2, size=B_5[yr-1], prob=r_5)))
            C_6[yr] <- sum(as.numeric(rbinom(2, size=B_6[yr-1], prob=r_6)))
            C_7p[yr] <- sum(as.numeric(rbinom(2, size=B_7p[yr-1], prob=r_7p)))
            N_1[yr] <- sum(as.numeric(C_4[yr]), as.numeric(C_5[yr]), as.numeric(C_6[yr]), as.numeric(C_7p[yr]))
            
            N_adults_total[h,i,j,k,n,yr] <- sum(as.numeric(N_4[yr]), as.numeric(N_5[yr]), as.numeric(N_6[yr]), as.numeric(N_7p[yr]))
            
            
          } # close yr
        } # close n
      } # close k
    } # close j
  } # close i
} # close h

# Stop the clock
proc.time() - ptm


save(B_total, N_adults_total, r_mean, s_juv_mean, s_ad_mean, b_mean, r_realized, s_juv_realized, s_ad_realized, b_realized, file = "Agestructured_simulations_data.RData")


