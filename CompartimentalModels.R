library(rootSolve)
library(deSolve) 
library(phaseR)

# Parameters
beta <- 6*10^(-9) # infectious rate
gamma <- 0.17 # recovery rate
# gamma = 1/D where D = average duration of infection
sigma <- 0.07 # rate of exposed individuals becoming infectious
# 1/sigma = average duration of incubation

# Initial conditions
Si <- 70*10^6 # number of susceptible individuals
Ei <- 0 # number of exposed individuals
Ii <- 1000 # number of infectious individuals
Ri <- 0 # number of recovered individuals

############################## Utils ##################################

plot_simulation_SIR <- function(results){
  par(mar=c(6,6,3,0.75))
  matplot(results[ ,1], results[ ,2:4],
          type = "l", lwd = 2, col = c("Green","Red","Blue"), lty = 1,
          xlab = "Time", ylab = "Number of individuals", main = "SIR model",
          cex.main = 1.5, cex.lab = 1.2, cex.axis = 0.9,
          lab=c(10, 6, 2), las = 1, mgp=c(3.5, 1, 0))
  
  legend("right", c("Susceptible", "Infectious","Recovered"),
         col = c("Green","Red","Blue"), lty = 1, lwd = 2)
}

plot_simulation_SEIR <- function(results){
  par(mar=c(6,6,3,0.75))
  matplot(results[ ,1], results[ ,2:5],
          type = "l", lwd = 2, col = c("Green","Orange", "Red","Blue"), lty = 1,
          xlab = "Time", ylab = "Number of individuals", main = "SEIR model",
          cex.main = 1.5, cex.lab = 1.2, cex.axis = 0.9,
          lab=c(10, 6, 2), las = 1, mgp=c(3.5, 1, 0))
  
  legend("left", c("Susceptible", "Exposed", "Infectious","Recovered"),
         col = c("Green","Orange","Red","Blue"), lty = 1, lwd = 2)
}

#################################### Fonction SIR : ####################################
SIR <- function(t,y,parameters){ 
  ds <- -parameters[1]*y[1]*y[2]
  di <- parameters[1]*y[1]*y[2]-parameters[2]*y[2]
  dr <- parameters[2]*y[2]
  list(c(ds,di,dr))
}

initial_population_SIR <- c(Si, Ii, Ri)
parameters_SIR <- c(beta,gamma)

# Simulation time
time <- seq(0,30*5,by=0.01)
result_SIR <- ode(initial_population_SIR, time, SIR, parameters_SIR)

plot_simulation_SIR(result_SIR)


#################################### Fonction SEIR : ####################################
SEIR <- function(t,y,parameters){ 
  ds <- -parameters[1]*y[1]*y[3]
  de <- parameters[1]*y[1]*y[3] - parameters[2]*y[2]
  di <- parameters[2]*y[2] - parameters[3]*y[3]
  dr <- parameters[3]*y[3]
  list(c(ds,de,di,dr))
}

initial_population_SEIR <- c(Si, Ei, Ii, Ri)
parameters_SEIR <- c(beta, sigma, gamma)

# Simulation time
time <- seq(0,30*10,by=0.01)
result_SEIR <- ode(initial_population_SEIR, time, SEIR, parameters_SEIR)
plot_simulation_SEIR(result_SEIR)
