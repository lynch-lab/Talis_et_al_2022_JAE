library(ggplot2)
library(reshape2)

# load data
load("ABC_data.RData")


# configure data frames
df_all <- melt(B_total, value.name = "B")
names(df_all)[names(df_all) == "Var1"] <- "iterate"
names(df_all)[names(df_all) == "Var2"] <- "year"
df_all$accept = rep(keep_bin, times = length(Data$year))
df_all$year_actual <- rep(Data$year, each = dim(B_total)[1])

df_subset_reject <- subset(df_all, df_all$iterate %in% sample(which(keep_bin == 0), 600))
df_accept <- subset(df_all, df_all$iterate %in% which(keep_bin == 1))

Data$iterate <- 1

# plot
p1 <- ggplot(df_subset_reject, mapping = aes(x=year_actual, y=B, group = iterate)) +
  geom_line(alpha = 0.2, colour = "#A1B0AB") + ylim(0,max(Data$count, na.rm = TRUE)+10000) +
  geom_line(df_accept, mapping = aes(x=year_actual, y=B, group = iterate), alpha = 0.3, colour = "#7B8CDE") +
  geom_line(Data, mapping = aes(x=year, y=count), colour = "#890620") +
  geom_point(Data, mapping = aes(x=year, y=count), colour = "#890620") +
  labs(y= "Breeding Abundance", x = "Year") 
p1 + ggtitle("Accepted and Rejected ABC-Simulated Time Series")