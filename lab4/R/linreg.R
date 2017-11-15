#' LAB 4
#'
#' @field formula formula.
#' @field data data.frame.
#' @field beta_hat matrix.
#' @field y_hat matrix.
#' @field e_hat matrix.
#' @field df numeric.
#' @field var_hat numeric.
#' @field p_value matrix.
#' @field t_value matrix.
#' @field var_beta_hat matrix.
#' @field data_name character.
#' @field res_val matrix.
#' 
#' @return print() should print out the coeﬃcients and coeﬃcient names,similar as done by the lm class,
#'         plot() should plot the following two plots using ggplot2,
#'         resid() should return the vector of residuals e_hat,
#'         pred() should return the predicted values y_hat,
#'         coef() should return the coeﬃcients as a named vector,
#'         summary() should return a similar printout as printed for lm objects.
#'         
#' @import methods
#' 
#' @references \href{https://en.wikipedia.org/wiki/Linear_regression}{Linear Regression Model}
#' @examples linreg_mod <- linreg$new(Petal.Length~Sepal.Width+Sepal.Length,data=iris)
#'  linreg_mod$print()
#'  linreg_mod$pred()
#'  linreg_mod$resid()
#'  linreg_mod$coef()
#'  linreg_mod$summary()
#'  linreg_mod$plot()
#'           
#' 

#Create the function:
linreg <- setRefClass(Class = "linreg",
                      fields = list(formula = "formula",data = "data.frame",beta_hat = "matrix",
                                    y_hat = "matrix",e_hat = "matrix",df = "numeric",
                                    var_hat = "numeric",p_value = "matrix",
                                    t_value = "matrix",var_beta_hat = "matrix",
                                    data_name = "character",res_val = "matrix"),
                      methods = list(initialize = function(formula, data){
                        formula <<- formula
                        data <<- data
                        data_name <<- deparse(substitute(data))
                        
                        #Create the matrix X
                        X <- model.matrix(formula, data)
                        
                        #pick out the dependent variable y
                        aux <- all.vars(expr = formula)[1]
                        y <- (data[,aux])
                        
                        #Regressions coefficients
                        beta_hat <<- solve((t(X) %*% X)) %*% t(X) %*% y
                        
                        
                        #The fitted values
                        y_hat <<- X %*% beta_hat
                        
                        #The residuals
                        e_hat <<- y - y_hat
                        
                        #The degrees of freedom
                        n <- nrow(X)
                        p <- ncol(X)
                        df <<- n - p
                        
                        #The Residual variance
                        var_hat <<- as.numeric((t(e_hat) %*% e_hat) / df)
                        
                        #The variance of the regression coefficients
                        var_beta_hat <<- var_hat * solve((t(X) %*% X))
                        
                        #The t-values for each coefficient
                        t_value <<- beta_hat / sqrt(diag(var_beta_hat))
                        
                        #The p-values for each regression coefficient
                        p_value <<- pt(abs(t_value), df = df,lower.tail=FALSE)

                        res_val <<- sqrt(abs((e_hat - mean(e_hat)) / sqrt(var_hat)))
                        },  

                        
                        # PRINT out the coefficients and coefficients names
                        print = function() {
                          
                          cat("\n","Call:","\n",
                              paste("linreg(", "formula = ", formula[2]," ", formula[1], " ", formula[3],", ", "data = ", data_name, ")",sep = "", collapse = "\n" ),
                              "\n","Coefficients:","\n",
                              paste(row.names(beta_hat),
                                    sep = "  ", collapse ="  " ),"\n",
                              format(round(beta_hat,2), justify = "centre",width = 10))
                        },
                        
                        # PLOT the following plots using ggplot2
                        plot=function() {
                          library(ggplot2) #necesary
                          library(ggThemeAssist) #necesary
                          
                          liu_blue <- "#54D8E0"
                          theme_liu <- theme(plot.margin = unit(c(1,1,1,1), "cm"),
                                             panel.background = element_rect(fill="white"),
                                             panel.grid.major.y = element_blank(),
                                             panel.grid.minor.y = element_blank(),
                                             panel.grid.major.x = element_blank(),
                                             panel.grid.minor.x = element_blank(),
                                             axis.line = element_line(color= "#58585b", size=0.1),
                                             axis.text.x = element_text(color="Black", size="10"),
                                             axis.text.y = element_text(color="Black", size="10"),
                                             axis.title.x = element_text(color="Black", size="10", face="bold"),
                                             axis.title.y = element_text(color="Black", size="10", face="bold"),
                                             axis.ticks.y = element_blank(),
                                             axis.ticks.x = element_line(color = "#58585b", size = 0.3),
                                             plot.title = element_text(color="Black", face="bold", size="14"),
                                             legend.position="bottom", legend.title = element_blank(),
                                             legend.key = element_blank(),
                                             legend.text = element_text(color="Black", size="10"))
                          
                          #Plot data frame
                          plot_df <- data.frame(df_resid = e_hat, df_fitted_values = y_hat)
                          
                          #FIRST PLOT: Residual vs Fitted
                          first_plot <-ggplot(data = plot_df, aes(x = df_fitted_values, y = df_resid)) +
                            geom_point(colour = liu_blue) +
                            geom_smooth(method = "loess",color = "red",se = FALSE) +
                            geom_abline(slope = 0,intercept = 0,linetype = "dotted") +
                            ggtitle("Residual vs Fitted") +
                            ylab("Residuals") +
                            xlab("Fitted Values") +
                            theme_liu
                          
                          #SECOND PLOT: Scale-Location
                          second_plot <- ggplot(data = plot_df,aes(x = df_fitted_values, y = sqrt(abs((df_resid - mean(df_resid)) / sqrt(var_hat))))) +
                            geom_point(colour = liu_blue) +
                            geom_smooth(method = "loess",color = "red",se = FALSE) +
                            ggtitle("Scale-Location") +
                            ylab(expression(sqrt(abs("Standardized Residuals")))) +
                            xlab("Fitted Values") +
                            theme_liu
                          
                          return(list(Residual_vs_Fitted = first_plot,Scale_Location = second_plot))
                        },
                        
                        #RESID should return the vector of residuals e_hat
                        resid = function() {
                          return((Residuals = round(e_hat, 2)))
                        },
                        
                        #PRED should return the predicted values of y_hat
                        pred = function() {
                          return((Fitted_values = round(y_hat, 2)))
                        },
                        
                        #COEF should return the coefficients as a named vector
                        coef = function() {
                          vector <- as.vector(beta_hat)
                          vector_names <- rownames(beta_hat)
                          names(vector) <- vector_names
                          return(vector)
                        },
                        
                        #SUMMARY should return a similar print out as printed of lm object
                        summary = function() {
                          coef_stand <- data.frame(var = rownames(beta_hat),stimate = round(beta_hat, 2),std.error = round(sqrt(diag(var_beta_hat)), 2),t_value = round(t_value, 2),p_value = round(p_value, 4))
                          coef_stand$var <- as.character(coef_stand$var)
                          # str(coef_stand)
                          row.names(coef_stand) <- NULL
                          
                          coef_stand <- rbind(c(" ", "Estimate", "Standard Error", "t-value", "P-Value"), coef_stand)
                          
                          #Accept or refuse the hypothesis
                          for(i in 2:nrow(coef_stand)){
                            if(coef_stand$p_value[i] == 0){
                              coef_stand$p_value[i] <- "***"
                            } else if(coef_stand$p_value[i] > 0 & coef_stand$p_value[i] <= 0.001){
                              coef_stand$p_value[i] <- "**"
                            } else if(coef_stand$p_value[i] > 0.001 & coef_stand$p_value[i] <= 0.01){
                              coef_stand$p_value[i] <- "*"
                            } else if(coef_stand$p_value[i] > 0.01 & coef_stand$p_value[i] <= 0.05){
                              coef_stand$p_value[i] <- "."
                            } else if(coef_stand$p_value[i] > 0.05 & coef_stand$p_value[i] <= 0.1){
                              coef_stand$p_value[i] <- " "
                            } else if(coef_stand$p_value[i] > 0.1){
                              coef_stand$p_value[i] <- " "
                            }
                          }
                          
                          for(c in 1:ncol(coef_stand)){
                            wdth <- max(nchar(as.character(coef_stand[, c])), na.rm = TRUE)
                            for(r in 1:nrow(coef_stand)){
                              coef_stand[r, c] <- format(coef_stand[r, c], width = wdth, justify = c("right"))
                            }
                          }
                          
                          #Print out the Console

                          cat("\n Call:\n\n",
                              paste("linreg(", "formula = ", formula[2]," ", formula[1], " ", formula[3],", ", "data = ", data_name, ")",sep = "", collapse = "\n" ),
                              "\n\n Coefficients:\n\n")
                          
                          for(i in 1:ncol(coef_stand)){
                            distance <- max(nchar(as.character(coef_stand[, i])), na.rm = TRUE)
                            for(j in 1:nrow(coef_stand)){
                              coef_stand[j, i] <- format(coef_stand[j, i], width=distance,justify = c("right"))
                            }
                          }
                          
                          for(i in 1:nrow(coef_stand)){
                            cat(paste(as.character(coef_stand[i, ]),collapse = " "),"\n")
                          }
                          cat("\n Residual standard error:",round(sqrt(var_hat), 2), "on", df, "degrees of freedom")
                          
                        }))
                          