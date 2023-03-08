# Q1: Run the model for 5 weeks without intervention, repeat this 3 times.  Plot cumulative incidence vs. time for the three runs. 

# Parameters (no intervention): time in week
alpha = 1; # incubation period: 7 days = 1 week
gamma.h = 7/5; # from onset to hospitalization: 5 days
gamma.d = 7/9.6; # from onset to death: 9.6 days
gamma.i = 7/10; # from onset to end of infectiousness for survivors: 10 days
gamma.f = 7/2; # from death to traditional burial 2 days
gamma.ih = 7/(10-5); # from hospitalization to end of infectiousness for survivors
gamma.dh = 7/(9.6-5); # from hospitalization to death
theta1 = 0.67; # proportion infectious in the hospital
delta1 = 0.8; # CFR for non-hospitalized
delta2 = 0.8; # CFR for hospitalized
beta.I = 0.588; # transmission rate in the community
beta.H = 0.794; # transmission rate in the hospital
beta.F = 7.653; # transmission rate at funerals; 

# Initial Conditions
N = 2e5
E0 = H0 = FF0 = R0 = 0; 
I0 = 3; 
S0 = N-I0; 
tm = 0; # in week
S = S0; E = E0; I = I0; H = H0; FF = FF0; R = R0; cumI = I0;
res = c(tm, NA, S, E, I, H, FF, R, cumI); 
num_wk = 5;

# Simulate the epidemic for the first 20 week without intervention:
while((I>0 | E>0 | H>0 | FF>0) & S>0 & tm < num_wk){ # if all the conditions are true, execute the following:
  # step 1: compute the transition rates
  rate.se = S/N*(beta.I*I + beta.H*H + beta.F*FF);  # S->E
  rate.ei = alpha*E;  # E->I
  rate.ih = gamma.h*theta1*I; # I->H
  rate.hf = gamma.dh*delta2*H;  # H->F
  rate.fr = gamma.f*FF;  # F->R
  rate.ir = gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
  rate.if = delta1*(1-theta1)*gamma.d*I; # I->F
  rate.hr = gamma.ih*(1-delta2)*H; # H->R
  rate.tot = rate.se + rate.ei + rate.ih + rate.hf + rate.fr + rate.ir + rate.if + rate.hr  # M_t = hazard rate
  
  # step 2: draw a random number, u1, between 0 and 1 and 
  # calculate the time after which the next transition occurs
  u1 = runif(1); # draw a random number 
  tau_i = -log(u1)/rate.tot
  
  # step 3: compute the probability that each type of transition will occur based on the rates
  # use this to calculate the range in which a number drawn at random 
  # must lie for a given transition to occur
  p.se = rate.se/rate.tot;  # prob the next event is S->E
  p.ei = rate.ei/rate.tot;  # prob the next event is E->I
  p.ih = rate.ih/rate.tot;  # prob the next event is I->H
  p.hf = rate.hf/rate.tot;  # prob the next event is H->F
  p.fr = rate.fr/rate.tot;  # prob the next event is F->R
  p.ir = rate.ir/rate.tot;  # prob the next event is I->R
  p.if = rate.if/rate.tot;  # prob the next event is I->F
  p.hr = rate.hr/rate.tot;  # prob the next event is H->R
  
  # step 4: draw a random number, u2, to determine the transition event which occurs next.
  u2 = runif(1); # draw a random number 
  
  # step 5: use the result from step 4 to update the numbers in each compartment
  if(u2 < p.se) { # u2 lies in [0, p.se), the next event is S->E
    S = S-1; E = E+1;
  } else if (u2 < p.se + p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
    E = E-1; I = I+1;
    cumI = cumI+1; # record cumulative incidence
  } else if (u2 < p.se + p.ei + p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
    I = I-1; H = H+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
    H = H-1; FF = FF+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr){ # FF->R
    FF = FF-1; R = R+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir){ # I->R
    I = I-1; R = R+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir + p.if) { # I->FF
    I = I-1; FF = FF+1;
  } else { # H->R
    H = H-1; R = R+1;
  }
  tm = tm + tau_i; # update time
  # save the results
  res = rbind(res, c(tm, tau_i, S, E, I, H, FF, R, cumI))
} 
colnames(res) = c('time', 'tm_step', 'S', 'E', 'I', 'H', 'FF', 'R', 'cumI')

# Need to run 3 times, so save the "res" for each time:
res1 = res # run the above code 1st time
res2 = res # run the above code 2nd time
res3 = res # run the above code 3rd time

if(T){
  # aggregate the results to weekly interval
  num_wk = ceiling(res1[nrow(res1), 'time'])
  weekly.res1 = matrix(NA, num_wk, 9); colnames(weekly.res1) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
  num_events = NULL;
  for(wk in 1:num_wk){
    idx = tail(which(res1[, 'time'] < wk), 1);
    weekly.res1[wk, ] = res1[idx, ];
  }
  # compute the weekly incidence
  weekly_inci = weekly.res1[-1, 'cumI'] - weekly.res1[-nrow(weekly.res1), 'cumI']
  # plot it
  par(mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.3, 0), tck = -0.02)
  plot(weekly_inci, type = 'l', ylab = 'Weekly incidence', xlab = 'Week', main = 'Res1')
} 

# Chatgpt way to do this:
# aggregate the results to weekly interval for res1
num_wk = ceiling(res1[nrow(res1), 'time'])
weekly.res1 = matrix(NA, num_wk, 9); colnames(weekly.res1) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
num_events = NULL;
for(wk in 1:num_wk){
  idx = tail(which(res1[, 'time'] < wk), 1);
  weekly.res1[wk, ] = res1[idx, ];
}
# compute the weekly incidence for res1
weekly_inci1 = weekly.res1[-1, 'cumI'] - weekly.res1[-nrow(weekly.res1), 'cumI']

# aggregate the results to weekly interval for res2
num_wk = ceiling(res2[nrow(res2), 'time'])
weekly.res2 = matrix(NA, num_wk, 9); colnames(weekly.res2) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
num_events = NULL;
for(wk in 1:num_wk){
  idx = tail(which(res2[, 'time'] < wk), 1);
  weekly.res2[wk, ] = res2[idx, ];
}
# compute the weekly incidence for res2
weekly_inci2 = weekly.res2[-1, 'cumI'] - weekly.res2[-nrow(weekly.res2), 'cumI']

# aggregate the results to weekly interval for res3
num_wk = ceiling(res3[nrow(res3), 'time'])
weekly.res3 = matrix(NA, num_wk, 9); colnames(weekly.res3) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
num_events = NULL;
for(wk in 1:num_wk){
  idx = tail(which(res3[, 'time'] < wk), 1);
  weekly.res3[wk, ] = res3[idx, ];
}
# compute the weekly incidence for res3
weekly_inci3 = weekly.res3[-1, 'cumI'] - weekly.res3[-nrow(weekly.res3), 'cumI']

# compute the maximum weekly incidence across all data sets
max_weekly_inci = max(c(max(weekly_inci1), max(weekly_inci2), max(weekly_inci3)))

# plot the weekly incidence of all 3 results on the same plot, with y-axis limit set to the maximum weekly incidence
par(mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.3, 0), tck = -0.02)
plot(weekly_inci1, type = 'l', ylab = 'Weekly incidence', xlab = 'Week', ylim = c(0, max_weekly_inci))
lines(weekly_inci2, col = 'red')
lines(weekly_inci3, col = 'blue')
legend('topleft', legend=c('res1', 'res2', 'res3'), col=c('black', 'red', 'blue'), lty=1)


# Q2: How does the computing time change with number of weeks simulated? 
n_week = c(9, 11, 13, 15, 17)
run_time = c(0.02, 0.18, 0.39, 12, 30)

plot(x = n_week, y = run_time)


# Q3: Modify the code (no intervention included) to run the Ebola stochastic model with no intervention during the first 8 weeks and with intervention starting Week 9 (Day 63).

# Parameters (no intervention to Week 8): time in week
alpha = 1; # incubation period: 7 days = 1 week
gamma.h = 7/5; # from onset to hospitalization: 5 days
gamma.d = 7/9.6; # from onset to death: 9.6 days
gamma.i = 7/10; # from onset to end of infectiousness for survivors: 10 days
gamma.f = 7/2; # from death to traditional burial 2 days
gamma.ih = 7/(10-5); # from hospitalization to end of infectiousness for survivors
gamma.dh = 7/(9.6-5); # from hospitalization to death
theta1 = 0.67; # proportion infectious in the hospital
delta1 = 0.8; # CFR for non-hospitalized
delta2 = 0.8; # CFR for hospitalized
beta.I = 0.588; # transmission rate in the community
beta.H = 0.794; # transmission rate in the hospital
beta.F = 7.653; # transmission rate at funerals; 

# Parameters (with intervention starting from Week 9)
beta.I2 = beta.I*0.12; # transmission rate in the community after intervention
beta.H2 = 0; # transmission rate in the hospital after intervention
beta.F2 = 0; # transmission rate at funerals after intervention

# Intervention starting from Week 9 (Day 63)
z = 0.88; # effectiveness of intervention in the community
z.H = 1; # effectiveness of intervention in the hospital
z.F = 1; # effectiveness of intervention in safe burial

# Initial Conditions
N = 2e5
E0 = H0 = FF0 = R0 = 0; 
I0 = 3; 
S0 = N-I0; 
tm = 0; # in week
S = S0; E = E0; I = I0; H = H0; FF = FF0; R = R0; cumI = I0;
res = c(tm, NA, S, E, I, H, FF, R, cumI); 
num_wk = 20;

# Simulate the epidemic for the first 20 week without intervention:
while((I>0 | E>0 | H>0 | FF>0) & S>0){ # if all the conditions are true, execute the following:
  if(tm < 9) { # no intervention
    beta.I = 0.588; # transmission rate in the community
    beta.H = 0.794; # transmission rate in the hospital
    beta.F = 7.653; # transmission rate at funerals; 
  } else { # with intervention
    beta.I = beta.I2
    beta.H = beta.H2
    beta.F = beta.F2
  }
  # step 1: compute the transition rates
  rate.se = S/N*(beta.I*I + beta.H*H + beta.F*FF);  # S->E
  rate.ei = alpha*E;  # E->I
  rate.ih = gamma.h*theta1*I; # I->H
  rate.hf = gamma.dh*delta2*H;  # H->F
  rate.fr = gamma.f*FF;  # F->R
  rate.ir = gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
  rate.if = delta1*(1-theta1)*gamma.d*I; # I->F
  rate.hr = gamma.ih*(1-delta2)*H; # H->R
  rate.tot = rate.se + rate.ei + rate.ih + rate.hf + rate.fr + rate.ir + rate.if + rate.hr  # M_t = hazard rate
  
  # step 2: draw a random number, u1, between 0 and 1 and 
  # calculate the time after which the next transition occurs
  u1 = runif(1); # draw a random number 
  tau_i = -log(u1)/rate.tot
  
  # step 3: compute the probability that each type of transition will occur based on the rates
  # use this to calculate the range in which a number drawn at random 
  # must lie for a given transition to occur
  p.se = rate.se/rate.tot;  # prob the next event is S->E
  p.ei = rate.ei/rate.tot;  # prob the next event is E->I
  p.ih = rate.ih/rate.tot;  # prob the next event is I->H
  p.hf = rate.hf/rate.tot;  # prob the next event is H->F
  p.fr = rate.fr/rate.tot;  # prob the next event is F->R
  p.ir = rate.ir/rate.tot;  # prob the next event is I->R
  p.if = rate.if/rate.tot;  # prob the next event is I->F
  p.hr = rate.hr/rate.tot;  # prob the next event is H->R
  
  # step 4: draw a random number, u2, to determine the transition event which occurs next.
  u2 = runif(1); # draw a random number 
  
  # step 5: use the result from step 4 to update the numbers in each compartment
  if(u2 < p.se) { # u2 lies in [0, p.se), the next event is S->E
    S = S-1; E = E+1;
  } else if (u2 < p.se + p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
    E = E-1; I = I+1;
    cumI = cumI+1; # record cumulative incidence
  } else if (u2 < p.se + p.ei + p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
    I = I-1; H = H+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
    H = H-1; FF = FF+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr){ # FF->R
    FF = FF-1; R = R+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir){ # I->R
    I = I-1; R = R+1;
  } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir + p.if) { # I->FF
    I = I-1; FF = FF+1;
  } else { # H->R
    H = H-1; R = R+1;
  }
  tm = tm + tau_i; # update time
  # save the results
  res = rbind(res, c(tm, tau_i, S, E, I, H, FF, R, cumI))
} 
colnames(res) = c('time', 'tm_step', 'S', 'E', 'I', 'H', 'FF', 'R', 'cumI')

if(T){
  # aggregate the results to weekly interval
  num_wk = ceiling(res[nrow(res), 'time'])
  weekly.res = matrix(NA, num_wk, 9); colnames(weekly.res) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
  num_events = NULL;
  for(wk in 1:num_wk){
    idx = tail(which(res[, 'time'] < wk), 1);
    weekly.res[wk, ] = res[idx, ];
  }
  # compute the weekly incidence
  weekly_inci = weekly.res[-1, 'cumI'] - weekly.res[-nrow(weekly.res), 'cumI']
  # plot it
  par(mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.3, 0), tck = -0.02)
  plot(weekly_inci, type = 'l', ylab = 'Weekly incidence', xlab = 'Week')
} 


# Q5: Run the Ebola stochastic model with intervention at week 9 (code provided in the script). Run the model until the epidemic dies out 
# (i.e. E=I=H=FF=0, no source of transmission). Repeat for 1000 times. Plot the distribution of final epidemic size and epidemic peak.
set.seed(123)

num_run = 1000

out = matrix(0, num_run, 4)
colnames(out) = c('num_run', 'duration', 'peak', 'final_size')
out[ , 1] = 1:num_run

for(ir in 1:num_run){
  # Parameters (no intervention to Week 8): time in week
  alpha = 1; # incubation period: 7 days = 1 week
  gamma.h = 7/5; # from onset to hospitalization: 5 days
  gamma.d = 7/9.6; # from onset to death: 9.6 days
  gamma.i = 7/10; # from onset to end of infectiousness for survivors: 10 days
  gamma.f = 7/2; # from death to traditional burial 2 days
  gamma.ih = 7/(10-5); # from hospitalization to end of infectiousness for survivors
  gamma.dh = 7/(9.6-5); # from hospitalization to death
  theta1 = 0.67; # proportion infectious in the hospital
  delta1 = 0.8; # CFR for non-hospitalized
  delta2 = 0.8; # CFR for hospitalized
  beta.I = 0.588; # transmission rate in the community
  beta.H = 0.794; # transmission rate in the hospital
  beta.F = 7.653; # transmission rate at funerals; 
  
  # Parameters (with intervention starting from Week 9)
  beta.I2 = beta.I*0.12; # transmission rate in the community after intervention
  beta.H2 = 0; # transmission rate in the hospital after intervention
  beta.F2 = 0; # transmission rate at funerals after intervention
  
  # Intervention starting from Week 9 (Day 63)
  z = 0.88; # effectiveness of intervention in the community
  z.H = 1; # effectiveness of intervention in the hospital
  z.F = 1; # effectiveness of intervention in safe burial
  
  N = 2e5
  E0 = H0 = FF0 = R0 = 0; 
  I0 = 3; 
  S0 = N-I0; 
  tm = 0; # in week
  S = S0; E = E0; I = I0; H = H0; FF = FF0; R = R0; cumI = I0;
  res = c(tm, NA, S, E, I, H, FF, R, cumI); 
  num_wk = 20
  
  while((I>0 | E>0 | H>0 | FF>0) & S>0){ # if all the conditions are true, execute the following:
    if(tm < 9) { # no intervention
      beta.I = 0.588; # transmission rate in the community
      beta.H = 0.794; # transmission rate in the hospital
      beta.F = 7.653; # transmission rate at funerals; 
    } else { # with intervention
      beta.I = beta.I2
      beta.H = beta.H2
      beta.F = beta.F2
    }
    # step 1: compute the transition rates
    rate.se = S/N*(beta.I*I + beta.H*H + beta.F*FF);  # S->E
    rate.ei = alpha*E;  # E->I
    rate.ih = gamma.h*theta1*I; # I->H
    rate.hf = gamma.dh*delta2*H;  # H->F
    rate.fr = gamma.f*FF;  # F->R
    rate.ir = gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
    rate.if = delta1*(1-theta1)*gamma.d*I; # I->F
    rate.hr = gamma.ih*(1-delta2)*H; # H->R
    rate.tot = rate.se + rate.ei + rate.ih + rate.hf + rate.fr + rate.ir + rate.if + rate.hr  # M_t = hazard rate
    
    # step 2: draw a random number, u1, between 0 and 1 and 
    # calculate the time after which the next transition occurs
    u1 = runif(1); # draw a random number 
    tau_i = -log(u1)/rate.tot
    
    # step 3: compute the probability that each type of transition will occur based on the rates
    # use this to calculate the range in which a number drawn at random 
    # must lie for a given transition to occur
    p.se = rate.se/rate.tot;  # prob the next event is S->E
    p.ei = rate.ei/rate.tot;  # prob the next event is E->I
    p.ih = rate.ih/rate.tot;  # prob the next event is I->H
    p.hf = rate.hf/rate.tot;  # prob the next event is H->F
    p.fr = rate.fr/rate.tot;  # prob the next event is F->R
    p.ir = rate.ir/rate.tot;  # prob the next event is I->R
    p.if = rate.if/rate.tot;  # prob the next event is I->F
    p.hr = rate.hr/rate.tot;  # prob the next event is H->R
    
    # step 4: draw a random number, u2, to determine the transition event which occurs next.
    u2 = runif(1); # draw a random number 
    
    # step 5: use the result from step 4 to update the numbers in each compartment
    if(u2 < p.se) { # u2 lies in [0, p.se), the next event is S->E
      S = S-1; E = E+1;
    } else if (u2 < p.se + p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
      E = E-1; I = I+1;
      cumI = cumI+1; # record cumulative incidence
    } else if (u2 < p.se + p.ei + p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
      I = I-1; H = H+1;
    } else if (u2 < p.se + p.ei + p.ih + p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
      H = H-1; FF = FF+1;
    } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr){ # FF->R
      FF = FF-1; R = R+1;
    } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir){ # I->R
      I = I-1; R = R+1;
    } else if (u2 < p.se + p.ei + p.ih + p.hf + p.fr + p.ir + p.if) { # I->FF
      I = I-1; FF = FF+1;
    } else { # H->R
      H = H-1; R = R+1;
    }
    tm = tm + tau_i; # update time
    # save the results
    res = rbind(res, c(tm, tau_i, S, E, I, H, FF, R, cumI))
  } 
  colnames(res) = c('time', 'tm_step', 'S', 'E', 'I', 'H', 'FF', 'R', 'cumI')
  
  # compile the results to weekly interval
  num_wk = ceiling(res[nrow(res), 'time'])
  weekly.res = matrix(NA, num_wk, 9); 
  colnames(weekly.res) = c('time', 'tm_step', 'S','E','I','H','FF','R','cumI')
  num_events = NULL;
  for(wk in 1:num_wk){
    idx = tail(which(res[, 'time'] < wk), 1);
    weekly.res[wk, ] = res[idx, ];
  }
  # compute the weekly incidence
  weekly_inci = weekly.res[-1, 'cumI'] - weekly.res[-nrow(weekly.res), 'cumI']

  # save the result for this run
  out[ir, 'duration'] = num_wk
  out[ir, 'peak'] = max(weekly_inci)
  out[ir, 'final_size'] = weekly.res[nrow(weekly.res), 'cumI']
  print(c(paste('run #', ir), out[ir, 2:4]))
}

# plot histogram of epidemic peak:
hist(out[ , 'peak'], xlab = "Peak Numbers", main = "Distribution of Epidemic Peak")

# plot histogram of final epidemic size:
hist(out[ , 'final_size'], xlab = "Final Epidemic Size", main = "Distribution of Final Epidemic Size")
