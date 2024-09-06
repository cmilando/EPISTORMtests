

library(outbreaks)
library(EpiEstim)
data("Flu2009")

# *****
daily_incidience_dates = Flu2009$incidence$dates
daily_incidence_values = as.integer(Flu2009$incidence$I)
daily_si_values = as.vector(Flu2009$si_distr)
# *****
source('00_rtcheck.R')
source('00_plotrtcheck.R')
source('00_getEpiEstim_R.R')

rt1 <- rtcheck(daily_incidience_dates, 
               daily_incidence_values, 
               daily_si_values,
               window_size = as.integer(7))
plot_rtcheck(rt1)