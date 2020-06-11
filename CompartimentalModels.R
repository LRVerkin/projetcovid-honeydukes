library(rootSolve)
library(deSolve) 
library(phaseR)

# Paramètres :
alpha <- 6*10^(-9)
rho <- 0.17

# Conditions initiales :
Si <- 70*10^6
Ii <- 1000
Ri <- 0

# Temps de simulation :
time <- seq(0,30*5,by=0.01)

#################################### Fonction SIR : ####################################
SIR <- function(t,y,parameters){ 
  ds <- -parameters[1]*y[1]*y[2]
  di <- parameters[1]*y[1]*y[2]-parameters[2]*y[2]
  dr <- parameters[2]*y[2]
  list(c(ds,di,dr))
}

initialisation1 <- c(Si, Ii, Ri)
parameters1 <- c(alpha,rho)

result1 <- ode(initialisation1, time, SIR, parameters1)

par(mar=c(6,6,3,0.75))
matplot(result1[ ,1], result1[ ,2:4],
        type = "l", lwd = 2, col = c("Green","Red","Blue"), lty = 1,
        xlab = "Temps", ylab = "Nombre d'individus", main = "Modèle SIR",
        cex.main = 1.5, cex.lab = 1.2, cex.axis = 0.9,
        lab=c(10, 6, 2), las = 1, mgp=c(3.5, 1, 0))

legend("right", c("Sains", "Infectés","resistants"),
       col = c("Green","Red","Blue"), lty = 1, lwd = 2)


#################################### Fonction SEIR : ####################################
SEIR <- function(t,y,parameters){ 
  ds <- -parameters[1]*y[1]*y[2]
  di <- parameters[1]*y[1]*y[2]-parameters[2]*y[2]
  dr <- parameters[2]*y[2]
  list(c(ds,di,dr))
}

initialisation <- c(Si, Ii, Ri)
parameters <- c(alpha,rho)

result <- ode(initialisation, time, SIR, parameters)
result
par(mar=c(6,6,3,0.75))
matplot(result[ ,1], result[ ,2:4],
        type = "l", lwd = 2, col = c("Green","Red","Blue"), lty = 1,
        xlab = "Temps", ylab = "Nombre d'individus", main = "Modèle SIR",
        cex.main = 1.5, cex.lab = 1.2, cex.axis = 0.9,
        lab=c(10, 6, 2), las = 1, mgp=c(3.5, 1, 0))

legend("right", c("Sains", "Infectés","resistants"),
       col = c("Green","Red","Blue"), lty = 1, lwd = 2)