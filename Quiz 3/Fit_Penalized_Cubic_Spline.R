# ================================================
# Implementation - Penalized Cublic Splines
# 
# STAT 460/560 - Statistical Theory I
# 
# Anthony-A. Christidis
# Department of Statistics
# University of British Columbia
# ================================================

# Fitting penalized cubic spline function 
Fit_Penalized_Cubic_Spline <- function(x, y, alpha){
  
  # Sort the data
  data_order <- order(x)
  x <- x[data_order]
  y <- y[data_order]
  
  # Find new interpolation points
  delta <- x[-1] - x[-length(x)]
  R_Matrix <- Get_R_Matrix(x)
  Q_Matrix <- Get_Q_Matrix(x)
  s <- alpha*solve(alpha*diag(x = 1, nrow = length(y)) +
                     (1 - alpha)*t(Q_Matrix) %*% solve(R_Matrix) %*% Q_Matrix) %*% y
  gamma <- solve(R_Matrix) %*% Q_Matrix %*% s
  
  # Return interpolation points and gamma values
  return(list(x = x, s = s, delta = delta, gamma = c(0, gamma, 0)))
}




