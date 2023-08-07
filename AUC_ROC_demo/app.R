#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidymodels)
library(ggsci)

pal = pal_npg("nrc")(9)
named_pal = c("COVID-19" = pal[1],
              "Healthy" = pal[2])

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("ROC-AUC Demo"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("cutoff",
                        "Cutoff for COVID-19 diagnosis",
                        min = 0,
                        max = 1,
                        value = 0.5, 
                        step = 0.1),
            sliderInput("mean_shift",
                        "Difference in mean expression",
                        min = 0,
                        max = 500,
                        value = 100,
                        step = 50),
            sliderInput("n",
                        "Group n",
                        min = 0,
                        max = 100,
                        value = 10, 
                        step = 1)
            ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("model_plot"),
           fluidRow(
             column(width = 6,
                    plotOutput("roc_auc")
             ),
             column(width = 6,
                    plotOutput("confusion")
             )
             ),
           fluidRow(
             uiOutput("sensitivity")
           ),
           fluidRow(
             uiOutput("fpr")
           )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #midpoint between the two distributions
  division_point = 500
  
  #healthy is a normal distribution centered at division_point - (mean_shift / 2)
  #COVID is a normal distribution centered at division_point + (mean_shift / 2)
  set.seed(12345)
  get_df = reactive({
    set.seed(12345)
    healthy_dist = rnorm(n = input$n, mean = (division_point - (input$mean_shift / 2)), sd = 100)
    #set all negative values to 0 to better simulate cytokine data
    healthy_dist[healthy_dist < 0] = 0
    
    set.seed(12345)
    covid_dist = rnorm(n = input$n, mean = (division_point + (input$mean_shift / 2)), sd = 100)
    
    #join into a single df
    out = data.frame(diagnosis = factor(c(rep("Healthy", input$n), rep("COVID-19", input$n)),
                                       levels = c("COVID-19", "Healthy")),
                    CXCL10 = c(healthy_dist, covid_dist))
    return(out)
  })
  
    #generate ROC-AUC with highlighting of current threshold
    output$roc_auc <- renderPlot({
      df = get_df()
      
      #generate logistic model
      cxcl10_model = logistic_reg() %>%
        set_engine("glm") %>%
        set_mode("classification") %>%
        fit(diagnosis ~ CXCL10, data = df)
      
      #fine to predict on train for this purpose
      cxcl10_predictions = predict(cxcl10_model,
                                   new_data = df,
                                   type = "prob") %>% 
        cbind(df, .) %>% 
        dplyr::mutate(threshold_pred = factor(ifelse(`.pred_COVID-19` >= input$cutoff,
                                                     yes = "COVID-19",
                                                     no = "Healthy"),
                                              levels = c("COVID-19", "Healthy")))
      
      #generate data for ROC
      roc_data = roc_curve(cxcl10_predictions,
                           truth = diagnosis,
                           `.pred_COVID-19`)
      
      
      #get current values
      confusion_matrix = conf_mat(data = cxcl10_predictions,
                                  truth = diagnosis,
                                  estimate = threshold_pred)$table
      true_pos = confusion_matrix %>% 
        as_tibble() %>% 
        dplyr::filter(Prediction == "COVID-19" &
                        Truth == "COVID-19") %>% 
        .$n %>% 
        sum()
      false_neg = confusion_matrix %>% 
        as_tibble() %>% 
        dplyr::filter(Prediction == "Healthy" &
                        Truth == "COVID-19") %>% 
        .$n 
      sensitivity = true_pos / (true_pos + false_neg)
      true_neg = confusion_matrix %>% 
        as_tibble() %>% 
        dplyr::filter(Prediction == "Healthy" &
                        Truth == "Healthy") %>% 
        .$n %>% 
        sum()
      false_pos = confusion_matrix %>% 
        as_tibble() %>% 
        dplyr::filter(Prediction == "COVID-19" &
                        Truth == "Healthy") %>% 
        .$n 
      fpr = false_pos / (false_pos + true_neg)
      
      #finally, plot
      roc_data %>% 
        autoplot() +
        annotate(geom = "point", 
                 x = fpr, 
                 y = sensitivity,
                 size = 4,
                 color = "red") +
        labs(x = "1 - Specificity / False Positive Rate",
             y = "Sensitivity / True Positive Rate") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 18),
              axis.title.x = element_text(size = 18),
              legend.text = element_text(size = 12))
    })
    
      #Generate plot of logistic regression model with points; highlight current point
      output$model_plot <- renderPlot({
        #regenerate data and model
        df = get_df() %>% 
          dplyr::mutate(diagnosis = factor(diagnosis,
                                           levels = c("Healthy", "COVID-19")))
        cxcl10_model = logistic_reg() %>%
          set_engine("glm") %>%
          set_mode("classification") %>%
          fit(diagnosis ~ CXCL10, data = df)
        
        preds = data.frame(CXCL10 = seq(0, max(df$CXCL10) * 1.2, 10)) %>%
          dplyr::mutate(sum_beta = coef(cxcl10_model$fit)["(Intercept)"] + 
                          CXCL10 * coef(cxcl10_model$fit)["CXCL10"],
                        prediction = plogis(sum_beta))
        
        ggplot(NULL) +
          geom_point(data = NULL, aes(x = cxcl10_model$fit$data$CXCL10, y = cxcl10_model$fit$y,
                                      fill = cxcl10_model$fit$data$diagnosis),
                     color = "black",
                     pch=21,
                     size = 3) +
          geom_line(data = preds, aes(x = CXCL10, y = prediction), size = 1.1) +
          geom_hline(yintercept = input$cutoff, linetype = 2) +
          geom_vline(xintercept = max(preds$CXCL10[preds$prediction <= input$cutoff]),
                     linetype = 2) +
          theme_bw() +
          scale_fill_manual(values = named_pal,
                            name = "Diagnosis") +
          labs(x = "[CXCL10] (pg/mL)", y = "Probability of COVID-19") +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14))
        })
      
      #display confusion matrix   
      output$confusion <- renderPlot({
        df = get_df()
        
        #regenerate logistic model
        cxcl10_model = logistic_reg() %>%
          set_engine("glm") %>%
          set_mode("classification") %>%
          fit(diagnosis ~ CXCL10, data = df)
        #generate continuous predictions and break at new cutoff
        cxcl10_predictions = predict(cxcl10_model,
                                     new_data = df,
                                     type = "prob") %>% 
          cbind(df, .) %>% 
          dplyr::mutate(threshold_pred = factor(ifelse(`.pred_COVID-19` >= input$cutoff,
                                                       yes = "COVID-19",
                                                       no = "Healthy"),
                                                levels = c("COVID-19", "Healthy")))
        
        confusion_matrix = conf_mat(data = cxcl10_predictions,
                                    truth = diagnosis,
                                    estimate = threshold_pred)
        autoplot(confusion_matrix, type = "heatmap") +
          theme(axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                legend.text = element_text(size = 12))
      })
      
      #display sensitivity calculation
      output$sensitivity <- renderUI({
        df = get_df()
        
        #regenerate logistic model
        cxcl10_model = logistic_reg() %>%
          set_engine("glm") %>%
          set_mode("classification") %>%
          fit(diagnosis ~ CXCL10, data = df)
        #generate continuous predictions and break at new cutoff
        cxcl10_predictions = predict(cxcl10_model,
                                     new_data = df,
                                     type = "prob") %>% 
          cbind(df, .) %>% 
          dplyr::mutate(threshold_pred = factor(ifelse(`.pred_COVID-19` >= input$cutoff,
                                                       yes = "COVID-19",
                                                       no = "Healthy"),
                                                levels = c("COVID-19", "Healthy")))
        
        confusion_matrix = conf_mat(data = cxcl10_predictions,
                                    truth = diagnosis,
                                    estimate = threshold_pred)$table
        true_pos = confusion_matrix %>% 
          as_tibble() %>% 
          dplyr::filter(Prediction == "COVID-19" &
                          Truth == "COVID-19") %>% 
          .$n %>% 
          sum()
        false_neg = confusion_matrix %>% 
          as_tibble() %>% 
          dplyr::filter(Prediction == "Healthy" &
                          Truth == "COVID-19") %>% 
          .$n 
        sensitivity = round(true_pos / (true_pos + false_neg),
                            digits = 2)
        
        withMathJax(paste0("$$Sensitivity = \\frac{TP}{TP + FN} = \\frac{",
                           true_pos, 
                           "}{",
                           true_pos,
                           " + ",
                           false_neg,
                           "} = ",
                           sensitivity,
                           "$$"))
      })
      
      output$fpr <- renderUI({
        df = get_df()
        
        #regenerate logistic model
        cxcl10_model = logistic_reg() %>%
          set_engine("glm") %>%
          set_mode("classification") %>%
          fit(diagnosis ~ CXCL10, data = df)
        #generate continuous predictions and break at new cutoff
        cxcl10_predictions = predict(cxcl10_model,
                                     new_data = df,
                                     type = "prob") %>% 
          cbind(df, .) %>% 
          dplyr::mutate(threshold_pred = factor(ifelse(`.pred_COVID-19` >= input$cutoff,
                                                       yes = "COVID-19",
                                                       no = "Healthy"),
                                                levels = c("COVID-19", "Healthy")))
        
        confusion_matrix = conf_mat(data = cxcl10_predictions,
                                    truth = diagnosis,
                                    estimate = threshold_pred)$table
        true_neg = confusion_matrix %>% 
          as_tibble() %>% 
          dplyr::filter(Prediction == "Healthy" &
                          Truth == "Healthy") %>% 
          .$n %>% 
          sum()
        false_pos = confusion_matrix %>% 
          as_tibble() %>% 
          dplyr::filter(Prediction == "COVID-19" &
                          Truth == "Healthy") %>% 
          .$n 
        fpr = round(false_pos / (false_pos + true_neg),
                    digits = 2)
        
        withMathJax(paste0("$$FPR = \\frac{FP}{FP + TN} = \\frac{",
                           false_pos, 
                           "}{",
                           false_pos,
                           " + ",
                           true_neg,
                           "} = ",
                           fpr,
                           "$$"))
      })
}

# Run the application 
shinyApp(ui = ui, server = server)
