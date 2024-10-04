library(shiny)
source("helpers.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("IP Visualization by Tumor Size"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      helpText("Select a ligand-receptor pair."),
      width = 3,
      selectInput('lig_rec', label = NULL, choices = lig_rec_pairs, selected = "Ccl3=Ccr1"),
      helpText("It may take ~10 seconds to make the plots...")
    ),
    
    
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Small", plotOutput("sizeSmall")),
                  tabPanel("Medium", plotOutput("sizeMedium")),
                  tabPanel("Large", plotOutput("sizeLarge")),
                  tabPanel("Alluvial", plotOutput("AlluvialFlow"))
      ))
      

    
)
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  ipSmall <- reactive((sigIP$s%>%filter(lr_pair==input$lig_rec)%>%pull(name)))
  ipMedium <- reactive((sigIP$m%>%filter(lr_pair==input$lig_rec)%>%pull(name)))
  ipLarge <- reactive((sigIP$l%>%filter(lr_pair==input$lig_rec)%>%pull(name)))
  in_lrp <- reactive(input$lig_rec)
  output$sizeSmall <- renderPlot({
    ip <- ipSmall()
    IP_Size(data = smallD,ip,ip_details=sigIP$s)
  }, height = 1500)
  output$sizeMedium <- renderPlot({
    ip <- ipMedium()
    IP_Size(data = mediumD,ip,ip_details=sigIP$m)
  }, height = 1500)
  output$sizeLarge <- renderPlot({
    ip <- ipLarge()
    IP_Size(data = largeD,ip,ip_details=sigIP$l)
  }, height = 1500)
  output$AlluvialFlow <- renderPlot({
    lrp <- in_lrp()
    DynamicIP(lrp)
  }, height = 850,width=850)
}

shinyApp(ui = ui, server = server)

