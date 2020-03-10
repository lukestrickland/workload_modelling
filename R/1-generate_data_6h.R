library(FIPS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

#Sim start and end (same for all)
simulation.start = lubridate::ymd_hms('2018-05-02 21:55:00', tz = "Australia/Perth")
simulation.end = lubridate::ymd_hms('2018-05-17 23:00:00', tz = "Australia/Perth")

#Set up common sleep schedule (last 3 days)
sleeptimes <- tibble::tibble(
  sleep.start = seq(
    from = lubridate::ymd_hms('2018-05-02 23:00:00', tz = "Australia/Perth"),
    to = lubridate::ymd_hms('2018-05-16 23:00:00', tz = "Australia/Perth"),
    by = '24 hours'),
  sleep.end = sleep.start + lubridate::dhours(7.96),
  sleep.id = rank(sleep.start))


simulation_df = FIPS::parse_sleeptimes(
  sleeptimes = sleeptimes ,
  series.start = simulation.start,
  series.end = simulation.end,
  sleep.start.col = "sleep.start",
  sleep.end.col = "sleep.end",
  sleep.id.col = "sleep.id",
  roundvalue = 5,
  grouping.factors = NULL)

simulation <- FIPS_simulation(
  FIPS_df = simulation_df , # The FIPS_df
  modeltype = "TPM",             # three process model
  pvec = TPM_make_pvec()      # paramater vector
)

sim2<- FIPS_simulation(
  FIPS_df = simulation_df , # The FIPS_df
  modeltype = "TPM",             # three process model
  pvec = TPM_make_pvec(g=-3.7)      # paramater vector
)

simulation$FIPS_df$sim2 <- sim2$FIPS_df$alertness

plot(simulation, plot_stat=c("alertness", "sim2"))

simulation <- simulation$FIPS_df

simulation <- simulation %>% mutate(sleep_group=sleep.id)%>% fill(sleep_group, .direction="down")

observed_times = 8:22

#observed_times = 1:24
#observed_times <- c(9,12,15,21)

#Add in data collection event at specified time
simulation <- simulation %>% mutate(observation=

                                            (hour(datetime) %in% observed_times &
                                               minute(datetime)==0 & sleep_group  <8)
)

simulation <- simulation %>% filter(wake_start==1|
                                            sleep_start==1|
                                            observation)

simulation <- simulation %>% mutate(time_since_previous =
                                            (as.numeric(datetime) - as.numeric(lag(datetime,1)))/3600,
                                          was_just_asleep = lag(sleep_start,1))

simulation$was_just_asleep[simulation$was_just_asleep==0] <- 2

simulation$event_number <- 1:length(simulation$datetime)
simulation


create_group_df <- function(sim_df, ns){
  for (i in ns) {
    tmp <- sim_df
    tmp$subject <- as.numeric(i)
    if (i==min(ns)) sim_df_ps <- tmp else
      sim_df_ps <- rbind(sim_df_ps, tmp)
  }
  sim_df_ps
}

full_df <- create_group_df(simulation, 1:40)

full_df$event <- factor(full_df$observation, levels=c("FALSE", "TRUE"),
                        labels=c("notobs", "observation"))


data_list_8h<- list(
  Nsubj=length(unique(full_df$subject)),
  Ntotal = length(full_df$datetime),
  subject = full_df$subject,
  event_number = full_df$event_number,
  previous_episode_type = full_df$was_just_asleep,
  time_since_previous = full_df$time_since_previous,
  timeofday = full_df$time,
  valid = which(full_df$observation),
  Nvalid= sum(full_df$observation)
)

data_list_8h

as.data.frame(data_list_8h[3:7])

save(full_df, data_list_8h, file="img/data_list_8h.RData")

table(full_df$time[full_df$observation],
      full_df$subject[full_df$observation])

table(full_df$datetime[full_df$wake_start==1],
      full_df$subject[full_df$wake_start==1])

table(full_df$time[full_df$sleep_start==1],
      full_df$subject[full_df$sleep_start==1])

plot(full_df$subject[full_df$sleep_start==1],
     full_df$time[full_df$sleep_start==1]
     )

plot(full_df$time_since_previous,
     full_df$was_just_asleep
)

plot(full_df$subject,
     full_df$taw
)

plot(full_df$subject,
     full_df$tas
)

table(full_df$was_just_asleep, full_df$subject)



