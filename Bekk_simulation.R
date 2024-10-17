### Packages ###
library(readxl)
library(BEKKs)

### Data ###
data_path <- '../data/data_new.xlsx'
gold <- as.data.frame(read_excel(data_path, sheet = 'GOLD_long'))
msci <- as.data.frame(read_excel(data_path, sheet = 'MSCI_long'))
gold[2:length(gold[,2]),3] <- diff(log(gold[,2]))
msci[2:length(msci[,2]),3] <- diff(log(msci[,2]))
gold <- gold[-1,]
msci <- msci[-1,]
colnames(gold) <- c("Date", "Gold", "Ret")
colnames(msci) <- c("Date", "Msci", "Ret")
gold$Date <- as.Date(gold$Date)
msci$Date <- as.Date(msci$Date)
common_dates <- as.Date(intersect(gold$Date, msci$Date))
gold <- gold %>% filter(Date %in% common_dates)
msci <- msci %>% filter(Date %in% common_dates)
gold$Date <- as.POSIXct(gold$Date)
msci$Date <- as.POSIXct(msci$Date)

data <- data.frame(msci$Ret - mean(msci$Ret),gold$Ret - mean(gold$Ret))
data <- as.matrix(data)

### BEKK model ###

matroot2 <- function(A){
  return(t(chol(A)))
}

bekk_s1 <- bekk_spec(model = list(type="bekk",asymmetric = TRUE))
bekk_model <- bekk_fit(bekk_s1,data)

### BEKK return simulation ###
et <- bekk_model$e_t
C <- bekk_model$C0
A <- bekk_model$A
B <- bekk_model$B
G <- bekk_model$G

H_sim <- matrix(0, nrow(et), 4) 
H_sim[1,] <- bekk_model$H_t[1,] # Take H_t[1,] as first initialization for H_t
ret_sim <- matrix(0, nrow(et), 2)
ret_sim[1,] <- c(msci$Ret[1],gold$Ret[1]) # Take the first return for initialization 
# ret_sim[1,] <- matroot2(matrix(H_sim[1,],2,2)) %*% et[1,] # Use this for simulated et

for (i in 2:nrow(et)){
  current_returns <- ret_sim[i - 1,]
  eta <- c(0,0) # eta with indicator function
  if (current_returns[1] < 0){eta[1] <- current_returns[1]}
  if (current_returns[2] < 0){eta[2] <- current_returns[2]}
  H_sim[i,] <- as.vector(C %*% t(C) + t(A) %*% current_returns %*% t(current_returns) %*% A + t(B) %*% eta %*% t(eta) %*% B + t(G) %*% matrix(H_sim[i-1,],2,2) %*% G)
  ret_sim[i,] <- matroot2(matrix(H_sim[i,],2,2)) %*% et[i,]
}

plot(ret_sim[,1], col = "red", type = "l")
lines(data[,1], col = "black", type = "l")

#Problem: Warum sind die simulierten returns so unterschiedlich?
#Problem: So scheinen die simulierten results ja sehr unrealistisch, mit einem return von teilweise 30 %

all(H_sim[2,] == bekk_model$H_t[2,])
H_sim[2,] - bekk_model$H_t[2,]
#Problem: Selbst bei i = 2, also EINER Berechnung von mir weicht die H matrix schon minimal von der korrekten im BEKK model ab
#Problem: Wodran liegt das? Ist das auf numerische Ungenauigkeit zurückzuführen?
