library(FIPS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

#Sim start and end (same for all)
simulation.start = lubridate::ymd_hms('2018-05-02 21:55:00', tz = "Australia/Perth")
simulation.end = lubridate::ymd_hms('2018-05-12 23:00:00', tz = "Australia/Perth")



#Set up common sleep schedule (last 3 days)
sleeptimes_common <- tibble::tibble(
  sleep.start = seq(
    from = lubridate::ymd_hms('2018-05-09 23:00:00', tz = "Australia/Perth"),
    to = lubridate::ymd_hms('2018-05-11 23:00:00', tz = "Australia/Perth"),
    by = '24 hours'),
  sleep.end = sleep.start + lubridate::dhours(7.95),
  sleep.id = rank(sleep.start)+7)


#wake_datetime <- lubridate::ymd_hms('2018-05-03 7:00:00', tz = "Australia/Perth")
# ndays=7
# sleep_hrs=3
# sim_start = simulation.start
# sim_end = simulation.end
# observed_times = c(9,12,15,21,0,2)


create_sleep_dep_df <- function(wake_datetime,
                                ndays,
                                sleep_hrs,
                                observed_times,
                                sim_start,
                                sim_end) {

  #Set up sleep schedule and build df for convenience
  sleeptimes <- rbind(tibble(
    sleep.start = seq(
      from = wake_datetime - hours(sleep_hrs),
      to = wake_datetime - hours(sleep_hrs) + days(ndays-1),
      by = '24 hours'),
    sleep.end = sleep.start + lubridate::dhours(sleep_hrs-0.08),
    sleep.id = rank(sleep.start)),
    sleeptimes_common)

  simulation_df = FIPS::parse_sleeptimes(
    sleeptimes = sleeptimes ,
    series.start = sim_start,
    series.end = sim_end,
    sleep.start.col = "sleep.start",
    sleep.end.col = "sleep.end",
    sleep.id.col = "sleep.id",
    roundvalue = 5,
    grouping.factors = NULL)
  
  simulation <- FIPS_simulation(
    FIPS_df = simulation_df , # The FIPS_df
    modeltype = "unified",             # three process model
    pvec = unified_make_pvec()      # paramater vector
  )
  
  plot(simulation)
  
  simulation <- simulation$FIPS_df

  simulation <- simulation %>% mutate(sleep_group=sleep.id)%>% fill(sleep_group, .direction="down")

  #Add in data collection event at specified time
  simulation <- simulation %>% mutate(observation=

                                              (hour(datetime) %in% observed_times &
                                                 minute(datetime)==0 & sleep_group  <8)|
                                              (hour(datetime) %in% c(9,12,15,21)&
                                                 minute(datetime)==0 & sleep_group  >7)
  )

  simulation <- simulation %>% filter(wake_start==1|
                                              sleep_start==1|
                                              observation)

  simulation <- simulation %>% mutate(time_since_previous =
                                              (as.numeric(datetime) - as.numeric(lag(datetime,1)))/3600,
                                            was_just_asleep = lag(sleep_start,1)

  )
  simulation$was_just_asleep[simulation$was_just_asleep==0] <- 2

  simulation$hour <- simulation$sim_hours - simulation$sim_hours[2]
  simulation <- simulation[-1,]
  simulation$time_since_previous[1] <- 0

  simulation$event_number <- 1:length(simulation$datetime)
  simulation


}

create_group_df <- function(sim_df, ns){
  for (i in ns) {
    tmp <- sim_df
    tmp$subject <- as.numeric(i)
    if (i==min(ns)) sim_df_ps <- tmp else
      sim_df_ps <- rbind(sim_df_ps, tmp)
  }
  sim_df_ps
}


first_wake <- lubridate::ymd_hms('2018-05-03 7:00:00', tz = "Australia/Perth")

simulation_3h <- create_sleep_dep_df( wake_datetime = first_wake,
                             ndays=7,
                             sleep_hrs=3,
                             sim_start = simulation.start,
                             sim_end = simulation.end,
                             observed_times = c(9,12,15,21,0,2))

simulation_3h_ps <- create_group_df(simulation_3h, 1:14)

simulation_5h <- create_sleep_dep_df( wake_datetime = first_wake,
                                      ndays=7,
                                      sleep_hrs=5,
                                      sim_start = simulation.start,
                                      sim_end = simulation.end,
                                      observed_times = c(9,12,15,21,0))

simulation_5h_ps <- create_group_df(simulation_5h, 15:26)

simulation_7h <- create_sleep_dep_df( wake_datetime = first_wake,
                                      ndays=7,
                                      sleep_hrs=7,
                                      sim_start = simulation.start,
                                      sim_end = simulation.end,
                                      observed_times = c(9,12,15,21))

simulation_7h_ps <- create_group_df(simulation_7h, 27:41)

simulation_9h <- create_sleep_dep_df( wake_datetime = first_wake,
                                      ndays=7,
                                      sleep_hrs=9,
                                      sim_start = simulation.start,
                                      sim_end = simulation.end,
                                      observed_times = c(9,12,15,21))

simulation_9h_ps <- create_group_df(simulation_9h, 42:57)


full_df <- rbind(simulation_3h_ps,
                   simulation_5h_ps,
                   simulation_7h_ps,
                   simulation_9h_ps)
full_df$event <- factor(full_df$observation, levels=c("FALSE", "TRUE"),
                        labels=c("notobs", "observation"))


data_list <- list(
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

data_list

as.data.frame(data_list[3:7])

save(full_df, data_list, file="img/data_list_sleepdep.RData")

quick_FIPS_plot <- function(simdf, ...) {
  simlist<- list(FIPS_df=simdf)
  class(simlist)<- c("FIPS_simulation", "list")
  plot(simlist, ...)
  
}

quick_FIPS_plot(simulation_3h, plot_stat="alertness")

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



