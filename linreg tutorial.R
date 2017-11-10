# -------------------------------------------------------
# Constructing OLS-estimates using matrix algebra
# Storing results in a nicely looking regression table
# -------------------------------------------------------

# Function takes y-variable and any  number of x-variables as arguments
linreg <- function(y, ...) { 
  
  # Constructing Matrix for Dependent Variable
  y.matrix <- as.matrix(y)
  
  # Constructing Matrix for Independent Variable
  
    # The Object predictors is defined as a list of all predictor variables the user entered
    predictors <- list(...) # This command can accept any number of predictor variables

    # The list of predictors is bound together in a matrix
    x.matrix.initial <- do.call(cbind,predictors) 

    # A column of 1s for the intercept is added. 
    x.matrix <- cbind(x.matrix.initial, rep(1, nrow(x.matrix.initial)))
    # The number of 1s added depends on the number of cases (nrow(x.matrix.initial))
    # Any number will work. The adequate number will be added. 
  
  
  # Estimating Regression Model
  
    # Regression Coefficients
    # Inverse of (transpose of X multiplied with X)
    tX.X <-solve(t(x.matrix)%*%x.matrix) 
    
    # Multiplying with Transpose of X
    X <- tX.X%*%t(x.matrix) 
    
    # Multiplying with Y
    b <- X%*%y.matrix # ols coefficients
    
    # Predicted Values and Residuals
    yhat <- x.matrix%*%b
    res <- y.matrix-yhat
    
    # Standard Errors 
    # For the following steps, first, number of rows and number of columns has to be defined
    n = nrow(x.matrix)
    k = ncol(x.matrix)
    # The numbers assigned to n and k depend on how many predictor variables the user entered
    
    # Variance-Covariance Matrix is calculated
    VCV = 1/(n-k) * as.numeric(t(res)%*%res) * solve(t(x.matrix)%*%x.matrix)
    
    # Standard errors are defined as the square root of the diagonal of the Variance-Covariance Matrix
    StdErr <- sqrt(diag(VCV))
    
    # T-Statistic
    T.Statistic <- rbind(b/StdErr)
    
    # P-value (two-tailed)
    P.Value = rbind(2*pt(abs(b/StdErr), df=n-k,lower.tail= FALSE))
    
    # Confidence Interval
    # In the CI-formula, the standard error is multiplied with the adequate 
    # coefficient of the t-distribution. This coefficient depends on n. 
    # It cannot be multiplied with a fixed number. Thus, the coefficient depends on n here.
    coefficient <- qt(0.975,df=n-1) # Coefficient is determined
    CI_L <- rbind(b-coefficient*StdErr) # Lower Bound of 95% CI
    CI_H <- rbind(b+coefficient*StdErr) # Upper Bound of 95% CI
    
    # R Squared
      # Sum of Squares Explained by Model
      ssm <- sum((yhat-y.matrix)^2)
      
      # Total Sum of Squares
      sst <- sum((y.matrix - mean(y.matrix))^2)
      
      # 1- ratio of SSM and SST
      r2 <- 1-(ssm/sst)
    
    # R Squared Adjusted
    # Here, the R Squared is adjusted to the number of predictors entered by the user  
    r2adj <- 1-((1-r2)*(length(y.matrix)-1))/(length(y.matrix)-length(predictors)-1)
    
    # F-Statistic
    # Variance explained by model
    MSR <- sum((yhat - mean(y.matrix))^2)/(k-1)  
    
    # Error Variance
    MSE <- (sum((y.matrix - yhat)^2))/(n-k)
    
    # Ratio=F
    f <- MSR/MSE
  
  
  # Constructing the Regression Table
  
    # Binding elements together
    Table <- cbind(b, StdErr, T.Statistic, P.Value, CI_L, CI_H)
    
    # Renaming columns
    colnames(Table)[1] <- "Beta"
    colnames(Table)[2] <- "Std. Err"
    colnames(Table)[3] <- "T-Statistic"
    colnames(Table)[4] <- "P-Value"
    colnames(Table)[5] <- "95% CI Low"
    colnames(Table)[6] <- "95% CI High"
    
    # Renaming rows (= Names of Regression Parameters)
    object <- as.list(substitute(list(...)))[-1L] # List of predictor names that were entered by user is extracted
    my.names <- c(object, c="Intercept") # Name for the Intercept is added to that list
    rownames(Table) <- my.names # Rownames are assigned
    
    # Rounding every entry to three digits
    Table <- round(Table,3)
  
  
  # Constructing table for Fit Statistics
  
    # Binding elements together
    Fit <- cbind(r2, r2adj, f)
    
    # Renaming columns
    colnames(Fit)[1] <- "R Squared"
    colnames(Fit)[2] <- "R Squared Adj."
    colnames(Fit)[3] <- "F-Statistic"
    
    # Rounding all entries to 4 digits
    Fit <- round(Fit, digits=4)
  
  
  # Storing results in a list and displaying this list as output of the function
  output.list <- list("Regression Table"=Table, "Fit Statistics"=Fit)
  output.list
  
}


# ----------------------
# Executing function
# ----------------------

# Constructing dataframe of lowely correlated data
new.data <- matrix((rnorm(1000, mean = 0, sd = 1)), nrow=200, ncol=5)
data.names <- c("Var1", "Var2", "Var3", "Var4", "Var5")
colnames(new.data) <- data.names
data <- as.data.frame(new.data)

# Executing function
with(data, linreg(Var1, Var2, Var3))
