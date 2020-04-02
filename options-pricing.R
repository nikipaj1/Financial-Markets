###############################################################
# =============================================================
# question 1 functions
# =============================================================
###############################################################


# function for matrix form binomial tree
bin_tree <- function(S0, u, d, N) {
  tree <- matrix(NA, nrow = N+1, ncol = N+1)
  
  for (i in 1:(N+1)){
    for (j in 1:i) {
      # u and d come into action only at t1
      tree[i,j] = S0 * u^(j-1) * d^(i-j)
    }
  }
    
  return(tree)
}


# function to calculate the risk-neutral synthetic probability
q_calc <- function(r,u,d) { 
  return((1.0+r-d)/(u-d))
}

# call/put return calculate
call_C <- function(St, K) {return(pmax(St - K,0))}
put_C <- function(St, K) { return(pmax(K-St,0))}


# function for the binomial tree, options price matrix and 
opt_price_EU <- function(S0, u, d, N, r, K, type = "European",ex_type = "Call"){
  
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
  # define step size h
  h = 1/M
  
  # initialise matrix
  matrix = matrix(0, N, M)
  matrix[,1] <- S0
  
  # initialise return vector
  return_vect <- rep(0, N)
  
  if(!Milstein){
    for (i in 1:N){
      for(j in 2:M){matrix[i,j] = rnorm(1)}
      
      
      # calculate stock price values for each node
      for (j in 2:M){
        matrix[i,j] = S_tph_EMM(matrix[i,j-1], r, h, sigma, matrix[i,j])
  
      }
      
      # average of stock prices
      A_T <- mean(matrix[i,])
  
      # calculate Asian-float returns
      return_vect[i] = C_AT(S_T = matrix[i,M], S_avg = A_T*k)
      
    }
  } else {
    
    for (i in 1:N){
      for(j in 2:M){matrix[i,j] = rnorm(1)}
      
      
      # calculate stock price values for each node
      for (j in 2:M){
        matrix[i,j] = S_tph_MM(matrix[i,j-1], r, h, sigma, matrix[i,j])
        
      }
      
      # average of stock prices
      A_T <- mean(matrix[i,])
      
      # calculate Asian-float returns
      return_vect[i] = C_AT(S_T = matrix[i,M], S_avg = A_T*k)
  }
  }
  
  
  C_hat <- exp(-r)*mean(return_vect)
  
  return(C_hat)
  

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
  # define step size
  h = 1/M
  
  # initialise price matrix
  mat_S = matrix(0, N, M)
  mat_S[,1] <- S0

  # initialise volatility matrix
  mat_V = matrix(0, N, M)
  mat_V[,1] <- V0
  
  # initialise return vector
  return <- rep(0,N)
  
  for(i in 1:N) {
    for(j in 2:M) {
      mat_S[i,j] <- rnorm(1) 
      mat_V[i,j] <- rnorm(1) 
    }
    
    # calculate the returns and volatility
    for(j in 2:M) {
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