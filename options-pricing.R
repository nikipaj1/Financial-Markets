library(ggplot2)
library(dplyr)
library(reshape2)


###############################################################
# =============================================================
# question 1 functions
# =============================================================
###############################################################


# function for matrix form binomial tree
bin_tree <- function(S0, u, d, N) {
  """ Function to calculate binomial trees based on the inputs
      
      The function calculates stock prices for each node in 
      the discrete model. The calculation involves going up i times
      and down j times for each entrance in the matrix. Then the value of the 
      stock is S0 * u^i * d^j.
     

      @parameter S0: initial price of the stock
      @parameter  u: the factor by which the price goes up
      @parameter  d: the factor by which the price goes down
      @parameter  N: the number of nodes

      Return Tree:  Binomial tree expressed in a matrix form.
  """
  tree <- matrix(NA, nrow = N+1, ncol = N+1)
  
  for (i in 1:(N+1)){
    for (j in 1:i) {
      # u and d come into action only at t1
      tree[i,j] = S0 * u^(j-1) * d^(i-j)
    }
  }
  return(tree)
}


q_calc <- function(r,u,d) { 
    """ Function to calculate the synthetic probability.

        @parameter r: risk-free rate
        @parameter  u: the factor by which the price goes up
        @parameter  d: the factor by which the price goes down

        Return q    synthetic probability
    """
  return((1.0+r-d)/(u-d))
}

# call/put return calculate
call_C <- function(St, K) {return(pmax(St - K,0))}
put_C <- function(St, K) { return(pmax(K-St,0))}


# function for the binomial tree, options price matrix and 
opt_price_EU <- function(S0, u, d, N, r, K, type = "European",ex_type = "Call"){
    """ Function to draw the binomial tree and calculate the stock prices.
        
        @parameter S0: initial price of the stock
        @parameter  u: the factor by which the price goes up
        @parameter  d: the factor by which the price goes down
        @parameter  N: the number of nodes
        @parameter  r: risk-free rate
        @parameter  K: strike price

        Return binomial_tree          -- binomial tree in the (N x N) matrix representation 
               opt_price_matrix       -- Options price matrix
               exer_mat               -- exercise matrix - matrix of ones and zeros showing whether the option is exercised at time t (American and European)
               q                      -- synthetic probability
    """
  
  # check for arbitrage
  if(!( (u > r + 1) & (r + 1 >= d) )){return("Your values of u, d and r permit arbitrage. Please check for u > 1+r > d.")}
  
  # draw binomial tree matrix
  binomial_tree <- bin_tree(S0 = S0, u = u, d = d, N = N)
  
  # initialise options price matrix
  opt_price_matrix <- matrix(0, nrow = nrow(binomial_tree), ncol = ncol(binomial_tree))

  if(ex_type=="Call"){
    # first we subtract the K from the last row and observe if its more than 0 (as usually)
    opt_price_matrix[nrow(opt_price_matrix),] = pmax(binomial_tree[nrow(binomial_tree),] - K, 0)
  }
  
  else if(ex_type=="Put"){
    opt_price_matrix[nrow(opt_price_matrix),] = pmax(K - binomial_tree[nrow(binomial_tree),],0)
  }
  
  else {return("Incorrect Option Call/Put type")}
  
  q = q_calc(r,u,d)
  
  # calculate the stock prices at each node
  for(i in (nrow(binomial_tree)-1):1){
    for(j in i:1){
      opt_price_matrix[i,j] = (q*opt_price_matrix[i+1,j+1] + (1-q)*opt_price_matrix[i+1,j])/(1+r)
    }
  }
      
      
  # initialise exercise matrix
  exer_mat <- matrix(0, nrow=nrow(opt_price_matrix), ncol = ncol(opt_price_matrix))

  if(type=="European") {
    # if EU return one only at the end
    
    exer_mat[N+1,] <- opt_price_matrix[N+1,]
    for(i in 0:N+1){
      if(exer_mat[N+1,i] > 0){ exer_mat[N+1,i] = 1}
    }
  }
  
  else if(type=="American") {
    # if American return ones at the nodes that are valuable
    
    exer_mat <- opt_price_matrix
    for(i in 1:N+1) {
      for (j in 1:i) {
        if(exer_mat[i,j] > 0) { exer_mat[i,j] = 1}
      }
    }
    exer_mat[1,1] <- 0
  }
  
  else {return("Please check the type of your option. You are only allowed to use European or American")}
  
  return(list(binomial_tree, opt_price_matrix, exer_mat, q))
}  


###############################################################
# =============================================================
# question 2 functions
# =============================================================
###############################################################

S_tph_EMM <- function(St, r, h, sigma, x_tph){
  return(St + r*St*h + sigma*St*x_tph*sqrt(h))
}

S_tph_MM <- function(St, r, h, sigma, x_tph){
  return(St + r*St*h + sigma*St*x_tph*sqrt(h) + 0.5*sigma^2*St*((x_tph*sqrt(h))^2-h))
}

C_AT <- function(S_T, S_avg){
  return(pmax(S_T - S_avg, 0))
}


asian_BSM <- function(M = 252, N = 10000, S0, r, sigma, k, Milstein = FALSE){
      """ Function to do Monte Carlo simulation of an Asian BSM. 
        
        @parameter  M: Number of trading days a year (inverse for step size)
        @parameter  N: Number of alterantive worlds for simulation
        @parameter S0: Initial price of the stock
        @parameter  r: risk-free rate
        @parameter  sigma: volatility
        @parameter  k: 
        @parameter  Milstein: If true: milstein method of simulation, otherwise Euler-Maruyama        

        Return return_vector    -- vector with returns in each world that are then averaged
    """
  
  # define step size h
  h = 1/M
  
  # initialise matrix
  matrix = matrix(0, N, M+1)
  matrix[,1] <- S0
  
  # initialise return vector
  return_vect <- rep(0, N)
  
  if(!Milstein){
    for (i in 1:N){
      for(j in 2:(M+1)){matrix[i,j] = rnorm(1)}
      
      
      # calculate stock price values for each node
      for (j in 2:(M+1)){
        matrix[i,j] = S_tph_EMM(matrix[i,j-1], r, h, sigma, matrix[i,j])
  
      }
      
      # average of stock prices
      A_T <- mean(matrix[i,])
  
      # calculate Asian-float returns
      return_vect[i] = C_AT(S_T = matrix[i,M+1], S_avg = A_T*k)
      
    }
  } else {
    
    for (i in 1:N){
      for(j in 2:(M+1)){matrix[i,j] = rnorm(1)}
      
      
      # calculate stock price values for each node
      for (j in 2:(M+1)){
        matrix[i,j] = S_tph_MM(matrix[i,j-1], r, h, sigma, matrix[i,j])
        
      }
      
      # average of stock prices
      A_T <- mean(matrix[i,])
      
      # calculate Asian-float returns
      return_vect[i] = C_AT(S_T = matrix[i,M+1], S_avg = A_T*k)
    }
  }
  
  x_vect <- seq(0,252,1)
  df <- data.frame("x" = x_vect)
  
  for(i in 1:100) {
    df[[paste("y", toString(i), sep = "")]] <- matrix[i,]
  } 
  
  df <- melt(df ,  id.vars = 'x', variable.name = 'series')
  
  ggplot(df, aes(x,value)) + geom_line(aes(colour = series)) +
  ggtitle("Asian Floating-Strike Call Option in Black-Sholes Model") +
  labs(y="Stock price / $", x = "Day") +
  theme(plot.title = element_text(hjust = 0.5))

  
  ggsave(paste('asian-floating-strike',"Melstein", toString(Milstein),".png",setp=""), device = "png", width = 40, height = 20, units = "cm")
  
  
  C_hat <- exp(-r)*mean(return_vect)
  

}

###############################################################
# =============================================================
# question 3 functions
# =============================================================
###############################################################


S_tph_hes <- function(St,Vt, r, h, sigma, x_tph){
  return(St + r*St*h + pmax(Vt,0)*St*x_tph*sqrt(h) + 0.5*pmax(Vt,0)^2*St*((x_tph*sqrt(h))^2-h))
}

V_tph <- function(Vt, theta, kappa, eta, h, epsilon) {
  return(Vt + (theta - kappa*pmax(Vt, 0))*h + eta*sqrt(pmax(Vt,0)*h)*epsilon)
}

# function to calculate Down-And-Out Barrier call
C_do <- function(vector, M, S_T, K, B) {
  for(i in 1:M){
    if(vector[i] < B)
      return(0)
  }
  
  return(pmax(S_T - K,0))
}


heston <- function(M,N,S0, V0, K, r, B, rho, kappa, theta, eta) {
    """ Function to do Monte Carlo simulation with Heston model
        
        @parameter  M: Number of trading days a year (inverse for step size)
        @parameter  N: Number of alterantive worlds for simulation
        @parameter S0: Initial price of the stock
        @parameter  r: risk-free rate
        @parameter  sigma: volatility
        @parameter  K: 
        @parameter rho: correlation between dWs, dWt
        @parameter kappa: param
        @parameter theta: param
        @parameter eta: param
        @parameter B  : boundary for down-and-out
        @parameter  Milstein: If true: milstein method of simulation, otherwise Euler-Maruyama        

        Return return_vector    -- vector with returns in each world that are then averaged
    """
  # define step size
  h = 1/M
  
  # initialise price matrix
  mat_S = matrix(0, N, (M+1))
  mat_S[,1] <- S0

  # initialise volatility matrix
  mat_V = matrix(0, N, (M+1))
  mat_V[,1] <- V0
  
  # initialise return vector
  return <- rep(0,N)
  
  for(i in 1:N) {
    for(j in 2:(M+1)) {
      mat_S[i,j] <- rnorm(1) 
      mat_V[i,j] <- rnorm(1) 
    }
    
    # calculate the returns and volatility
    for(j in 2:(M+1)) {
      epsilon <- rnorm(1)
      rnd <- rnorm(1)
      x_tph <- rho*epsilon + sqrt(1-rho^2)*rnd
      
      # calculate returns matrix with Misltein
      mat_S[i,j] = S_tph_hes(mat_S[i,j-1],mat_V[i, j-1], r, h, x_tph = x_tph)
      # calculate volatility matrix with Euler-Maruyama
      mat_V[i,j] = V_tph(mat_V[i, j-1], theta, kappa, eta, h, epsilon)
    }
    
    
    # compute payoff for each simulation with DO
    return[i] <- C_do(vector = mat_S[i,], M= M, S_T = mat_S[i,M], K = K, B=B)
    
  }
  answer <- mean(return)
  
  
  # Plot the simulation
  x_vect <- seq(0,252,1)
  df <- data.frame("x" = x_vect)
  
  for(i in 1:50) {
    df[[paste("y", toString(i), sep = "")]] <- mat_S[i,]
  } 
  
  df <- melt(df ,  id.vars = 'x', variable.name = 'series')
  
  ggplot(df, aes(x,value)) + geom_line(aes(colour = series)) +
    ggtitle("Heston") +
    labs(y="Stock price / $", x = "Day") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste('Heston',"rho", toString(rho), "theta", toString(theta),".png",sep=""), 
         device = "png", 
         width = 40, 
         height = 20, 
         units = "cm")
  
  return(list(answer))
}


###############################################################
# =============================================================
# Simulate
# =============================================================
###############################################################


# European call
europe_cal <- opt_price_EU(S0 = 100, u = 1.1, d = 0.9, N = 10, r = 0.05, K = 100, type = "European",ex_type = "Call")
europe_cal
print("===================================================================")
europe_put <- opt_price_EU(S0 = 100, u = 1.1, d = 0.9, N = 10, r = 0.05, K = 100, type = "European",ex_type = "Put")
europe_put
print("===================================================================")
america_call <- opt_price_EU(S0 = 100, u = 1.1, d = 0.9, N = 10, r = 0.05, K = 100, type = "American",ex_type = "Call")
america_call
print("===================================================================")
america_put <- opt_price_EU(S0 = 100, u = 1.1, d = 0.9, N = 10, r = 0.05, K = 100, type = "American",ex_type = "Put")
america_put

# ========================================================================
# question 2 simulation


# Euler Maruyama simulation
asian_BSM(S0 = 100, r = 0.04, sigma = 0.25, k=0.9, Milstein = FALSE)

# Milstein simulation
asian_BSM(S0 = 100, r = 0.04, sigma = 0.25, k=0.9, Milstein = TRUE)

# The Milstein is slightly more rigorous as it adds a few terms into the equation, tihs makes it longer to calculate. But in this case it does not bring us much value as the difference in prices is minute.

# ========================================================================
# question 3 simulation

heston(M = 252,N = 100, S0 = 90, V0 = (0.04/1.15), K = 75, B = 70, r = (0.02/0.58), rho = -0.64, kappa = 1.15, theta = 0.04, eta = 0.39)

heston(M = 252,N = 100, S0 = 90, V0 = (0.06/1.15), K = 75, B = 70, r = (0.02/0.58), rho = -0.64, kappa = 1.15, theta = 0.06, eta = 0.39)


# 18.68878 vs 17.71818 => the idea that as the volatility increases, the prices go down, holds. i.e. negative correlation between the volatility and stock prices.

heston(M = 252,N = 100, S0 = 90, V0 = (0.06/1.15), K = 75, B = 70, r = (0.02/0.58), rho = 0, kappa = 1.15, theta = 0.06, eta = 0.39)


# the result went even higher. this is because due to no correlation between volatility and the stock prices, the volatility is smaller (dW = sigma * h)
# if there is no correlation, we get a symmetric smiles
