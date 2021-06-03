library(shiny)
library(shinythemes)
library(readr)
library(markdown)
library(rsconnect)
library(googlesheets4)
library(DT)

pckgs <- c('dplyr','tibble',
           'cowplot','ggplot2',
           'survival','mvtnorm')
for (pp in pckgs) { library(pp,character.only = T)}

# Define UI for app that draws a histogram ----
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage("Re-operation after pyeloplasty",

                tabPanel("Survival Distribution",
                                           
                  sidebarLayout(
                  sidebarPanel(#HTML('Input the following patient information to output the individualized survival function after pyeloplasty.'),
                      HTML("<h4><b>Patient Information</b></h4>"),
                      selectInput("op_type", label = "Surgical Approach", 
                                  c("Open" = "open", "Minimally Invasive" = "MIS")),
                      selectInput("stent_type", label = "Type of Stent used", 
                                  c("JJ Stent" = "jj_s", "Salle Stent" = "salle", "Other" = "other")),
                      selectInput("use_cath", label = "Indwelling catheter left after surgery", 
                                  c("Yes" = "yes", "No" = "no")),
                      selectInput("use_drain", label = "JP drain Used", 
                                  c("Yes" = "yes", "No" = "no")),
                      selectInput("use_narcotic", label = "Narcotics used after surgery", 
                                  c("Yes" = "yes", "No" = "no")),
                      numericInput(inputId = "LOS", "Length of stay (days)", value = 2),
                      numericInput(inputId = "dur_IV", "Duration of IV fluids (hours)", value = 30),
                      HTML("<h4><b>Antero-Posterior Diameter (APD)</b></h4>"),
                      numericInput(inputId = "Pre_op_APD", "Pre-operative APD (mm)", value = 40),
                      numericInput(inputId = "Post_op_APD", "Post-operative APD (mm)", value = 20),
                      numericInput(inputId = "sec_APD", "Second Follow-up APD (mm)", value = 15),
                      actionButton("submitbutton", 
                                   "Determine Risk", 
                                   class = "btn btn-primary"),
                      HTML("<br> <br> Data may be collected for quality assurance purposes")
                    ),
                    mainPanel(plotOutput(outputId = "survPlot"),
                              htmlOutput("thirtymonth"),
                              htmlOutput("diagnosis"))
                  )),
                tabPanel("About", div(includeMarkdown("about.md"), align="justify"))
                
                )
)

server <- function(input, output, session) {
  
  source('funs_support.R')
  load('cox_mdl.RData')
  
  observeEvent(input$submitbutton,{
  output$survPlot <- renderPlot({
    
    if(input$use_cath == 'Yes'){has_catheter <- 1}else{has_catheter <- 0}
    if(input$use_drain == 'Yes'){has_drain <- 1}else{has_drain <- 0}
    if(input$use_narcotic == 'Yes'){has_narc_use <- 1}else{has_narc_use <- 0}
    duration_IV <- input$dur_IV
    post_op_APD <- input$Post_op_APD
    sec_APD <- input$sec_APD
    pct_improve_2nd <- (input$Pre_op_APD - input$sec_APD)/input$Pre_op_APD
    
    df = tibble(Who_indicated4=0, Surgeon6=0, Return_DietOther=0,
                Catheter1=has_catheter, Drain1=has_drain, 
                Narc_use=has_narc_use, Duration_IV=duration_IV,
                Post_op_APD=post_op_APD, sec_APD=sec_APD,
                percent_improve_2nd=pct_improve_2nd)
    stopifnot(all(colnames(df) %in% names(lst_cox$bhat)))
    
    X_row = rownames_to_column(data.frame(x=t(df)),'cn')
    dat_norm = tibble(cn=names(lst_cox$mu),bhat=lst_cox$bhat,
                      mu=lst_cox$mu,se=lst_cox$se)
    X_row = left_join(X_row,dat_norm,by='cn')
    X_row = X_row %>% mutate(x_s = (x-mu)/se)
    
    
    idx_align = match(X_row$cn,names(lst_cox$bhat))
    stopifnot(all(lst_cox$bhat[idx_align] == X_row$bhat))
    Sigma = lst_cox$Sigma[idx_align,idx_align]
    
    set.seed(1234)
    surv_dist = pm_surv(bhat=X_row$bhat, Sigma = Sigma, Eta=exp(lst_cox$eta),
                        Y = lst_cox$y, x=matrix(X_row$x_s,nrow=1),
                        nsim = 1000,alpha = 0.05)
    
      surv_24m <- surv_dist %>% filter(time==24)
      output$thirtymonth <- renderText({paste("<br>The probability of not requiring re-intervention within 24 months is <b>", round(100*surv_24m[[2]], digits=1), "% [95%CI:", round(100*surv_24m[[3]], digits=1), ",", round(100*surv_24m[[4]], digits=1), "]</b>.")})
    
      gg_km = ggplot(surv_dist,aes(x=time,y=mu)) +
      theme_bw(base_size=20) + geom_line() +
      labs(y='Survival probability',x='Months from pyeloplasty',subtitle = 'Shaded area is 95% confidence interval.') +
      geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5) +
      scale_y_continuous(limits=c(0,1))
    gg_km
  })
  
  output$diagnosis <- renderText({' <br> More information can be found under <b> About </b>.
    <br> <br> <b> Reference </b> <br> Application of Machine Learning Algorithms to Identify Patients at Risk for Recurrent UPJO after Dismembered Pyeloplasty. <br>
    <i> Drysdale E., Khondker A., Kim JK., Erdman L., Kwong JCC., Chua M., Keefe DT., Lolas M., Dos Santos J., Rickard M., Lorenzo AJ. (in preparation) </i> <br>
    <br> <b> Disclaimer </b> <br> This web application does not provide medical advice.
    Access to general information is provided for educational purposes only. Content is not recommended or endorsed 
    by any doctor or healthcare provider. The application provided is not a substitute for medical or professional 
    care, virtual care, consultation or the advice of your healthcare provider. You are liable or responsible for any advice, 
    course of treatment, diagnosis or any other information, services or product obtained through this web application.
    </h5>'})})
}

shinyApp(ui, server)
