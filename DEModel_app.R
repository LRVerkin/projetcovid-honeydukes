library(shiny)
library(rootSolve)
library(deSolve) 
library(phaseR)

ui <- fluidPage(
  titlePanel('Modèles épidémiologiques déterministes'),
  sidebarPanel(
    h2('Paramètres initiaux'),
    br(),
    p('La proportion initiale d\'infectés est fixée à 0.01%'),
    sliderInput(inputId = "pop", 
                label = "Total de la population", 
                value = 70*10^6, min = 10^6, max = 10^8),
    sliderInput(inputId = "vax", 
                label = "Pourcentage de vaccinés", 
                value = 0, min = 0, max = 99.99),
    sliderInput(inputId = "temps",
                label = "Temps de la simulation en jours",
                value = 60, min = 0, max = 2160),
    h2('Coefficients'),
    sliderInput(inputId = "alpha",
                label = "Taux d'infection alpha (.10^-8)",
                value = 0.6, min = 0, max = 3.9),
    sliderInput(inputId = "rho",
                label = "Taux de retrait",
                value = 0.17, min = 0, max = 1),
    sliderInput(inputId = "mu",
                label = "Taux de mortalité",
                value = 0.1, min = 0, max = 1)
  ),
  mainPanel(
    plotOutput("graph"),
    textOutput("R0"))
)

###################################################################################
server <- function(input, output) {

  SIRM <- function(t,y,parameters){ 
    ds <- -parameters[1]*y[1]*y[2]
    di <- parameters[1]*y[1]*y[2]-parameters[2]*y[2]-parameters[3]*y[2]
    dr <- parameters[2]*y[2]
    dm <- parameters[3]*y[2]
    list(c(ds,di,dr,dm))
  }
  
  output$graph <- renderPlot({
    # Paramètres :
    alpha <- input$alpha*10^(-8) # Taux d'infection
    rho <- input$rho # Taux de resistance
    mu <- input$mu # Taux de mort

    # Conditions initiales :
    Ri <- (input$vax/100)*input$pop # Taux de vax * pop
    Ii <- 0.01*input$pop # Taux d'infectés * pop
    Si <- input$pop - (Ri + Ii) # Reste de la population
    Mi <- 0 # Nombre de morts de la maladie
    
    # Temps de simulation :
    time <- seq(0,input$temps,by=0.01)
    
    initialisation <- c(Si, Ii, Ri, Mi)
    parameters <- c(alpha, rho, mu)
    
    result <- ode(initialisation, time, SIRM, parameters)
    par(mar=c(6,6,3,0.75))
    matplot(result[ ,1], result[ ,2:5],
            type = "l", lwd = 2, col = c("Green","Red","Blue","Black"), lty = 1,
            xlab = "Temps (jours)", ylab = "Nombre d'individus", main = "Modèle SIRM",
            cex.main = 1.5, cex.lab = 1.2, cex.axis = 0.9,
            lab=c(10, 6, 2), las = 1, mgp=c(3.5, 1, 0))

    legend("topright", c("Sains","Infectés","Résistants","Morts"),
           col = c("Green","Red","Blue","Black"), lty = 1, lwd = 2)
  })
  
  output$R0 <- renderText({
    paste(c("R0 = ", (input$pop - ((input$vax/100)*input$pop + 0.01*input$pop))*input$alpha*10^(-8)/(input$rho+input$mu)), collapse = " ")
  })

}

shinyApp(ui = ui, server = server)
