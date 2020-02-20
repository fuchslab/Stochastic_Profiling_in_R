# Function to calculate the Difference of the two distirbutions

OVL_LN_LN <- function(mu_1,mu_2,sigma_1,sigma_2){
  f1 <- function(x){dlnorm(x, meanlog = mu_1, sdlog = sigma_1) }
  f2 <- function(x){dlnorm(x, meanlog = mu_2, sdlog = sigma_2) }  
  f3 <- function(x){pmin(f1(x), f2(x))}
  integrate( f3, lower = 0.000001, upper = Inf)$value 
  }



  
