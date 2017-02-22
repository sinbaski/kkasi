library(snowfall)     # Parallel computing
library(rlecuyer)     # RNG

################################################################################

sfInit(parallel = TRUE, cpus = 4, type = "SOCK")
# sfInit(parallel = FALSE)  # useful for trouble shooting
# load functions "CPT", "TS.SZ", "Teststatistik", "armasim_errorbreak", "archsim" below

EQ.estimate <- function(X, k.NUOS, prob){   # NUOS = number of upper order statistics
  n             <- length(X)
  Y             <- sort.int(X, decreasing = TRUE, method = "quick")
  Hill.estimate <- 1 / k.NUOS * sum(log(Y[1 : (k.NUOS + 1)] / Y[k.NUOS + 1]))
  return(log(Y[k.NUOS + 1] * ((n * prob) / k.NUOS)^( - Hill.estimate)))    # return extreme quantile estimate
}

sfExport("CPT", "TS.SZ", "TS.bootstrap", "Teststatistik", "armasim_errorbreak", "archsim", "EQ.estimate")
sfLibrary(matrixStats)
sfLibrary(TSA)
sfLibrary(boot)
sfClusterSetupRNG()   # type = "SPRNG", set.seed = 12345

wrapper <- function(p, d, n, comparison){
    CPT(AnzCPTs = 4000, AnzSim=n, k = c(0.08, 0.16, 0.24), d, p, comparison)
}

start <- Sys.time()
###############################

# for bootstrap only
p <- c(0.005, 0.01, 0.05, 0.1)

MA2_0907_2000 <- sfLapply(p, wrapper, d = 1, n=2000, comparison = "no")
MA2_0907_500  <- sfLapply(p, wrapper, d = 1, n=500,  comparison = "no")
ARCH1_08_2000 <- sfLapply(p, wrapper, d = 2, n=2000, comparison = "no")
ARCH1_08_500  <- sfLapply(p, wrapper, d = 2, n=500,  comparison = "no")

# Alternative:
MA2_break_05_2000 <- sfLapply(p, wrapper, d = 3, n=2000, comparison = "no")
MA2_break_05_500  <- sfLapply(p, wrapper, d = 3, n=500,  comparison = "no")
MA2_break_075_2000<- sfLapply(p, wrapper, d = 4, n=2000, comparison = "no")
MA2_break_075_500 <- sfLapply(p, wrapper, d = 4, n=500,  comparison = "no")

ARCH1a_break_05_2000 <- sfLapply(p, wrapper, d = 5, n=2000, comparison = "no")
ARCH1a_break_05_500  <- sfLapply(p, wrapper, d = 5, n=500,  comparison = "no")
ARCH1a_break_075_2000<- sfLapply(p, wrapper, d = 6, n=2000, comparison = "no")
ARCH1a_break_075_500 <- sfLapply(p, wrapper, d = 6, n=500,  comparison = "no")

ARCH1b_break_05_2000 <- sfLapply(p, wrapper, d = 7, n=2000, comparison = "no")
ARCH1b_break_05_500  <- sfLapply(p, wrapper, d = 7, n=500,  comparison = "no")
ARCH1b_break_075_2000<- sfLapply(p, wrapper, d = 8, n=2000, comparison = "no")
ARCH1b_break_075_500 <- sfLapply(p, wrapper, d = 8, n=500,  comparison = "no")


print(Sys.time() - start)

# end for bootstrap only

###############################

# Null:
p <- c(0.1, 0.05, 0.01, 0.005)    # (1-p)-quantiles of interest
MA2_0907_2000 <- sfLapply(p, wrapper, d = 1, n=2000, comparison = "no")
MA2_0907_500  <- sfLapply(p, wrapper, d = 1, n=500,  comparison = "yes")

ARCH1_08_2000 <- sfLapply(p, wrapper, d = 2, n=2000, comparison = "yes")
ARCH1_08_500  <- sfLapply(p, wrapper, d = 2, n=500,  comparison = "yes")

# Alternative:
MA2_break_05_2000 <- sfLapply(p, wrapper, d = 3, n=2000, comparison = "yes")
MA2_break_05_500  <- sfLapply(p, wrapper, d = 3, n=500,  comparison = "yes")
MA2_break_075_2000<- sfLapply(p, wrapper, d = 4, n=2000, comparison = "yes")
MA2_break_075_500 <- sfLapply(p, wrapper, d = 4, n=500,  comparison = "yes")

ARCH1a_break_05_2000 <- sfLapply(p, wrapper, d = 5, n=2000, comparison = "yes")
ARCH1a_break_05_500  <- sfLapply(p, wrapper, d = 5, n=500,  comparison = "yes")
ARCH1a_break_075_2000<- sfLapply(p, wrapper, d = 6, n=2000, comparison = "yes")
ARCH1a_break_075_500 <- sfLapply(p, wrapper, d = 6, n=500,  comparison = "yes")

ARCH1b_break_05_2000 <- sfLapply(p, wrapper, d = 7, n=2000, comparison = "yes")
ARCH1b_break_05_500  <- sfLapply(p, wrapper, d = 7, n=500,  comparison = "yes")
ARCH1b_break_075_2000<- sfLapply(p, wrapper, d = 8, n=2000, comparison = "yes")
ARCH1b_break_075_500 <- sfLapply(p, wrapper, d = 8, n=500,  comparison = "yes")

print(Sys.time() - start)

sfStop()

################################################################################

CPT <- function(AnzCPTs, AnzSim, k, d, p, comparison){       # Routine zur Ausführung der Change Point Tests mit AnzCPTs = Anzahl der Simulationen, Extremwertanteile fuer Tail Index Estimator
    AnzAntEx <- length(k)

    rej     <- array(0, dim = c(4, AnzAntEx))        # EQ & TI tests for 1%- and 5%-level (= 4)
    TeSta   <- array(0, dim = c(AnzAntEx, AnzCPTs))  # AnzAntEx different k's, AnzCPTs repetitions
    rej.SZ  <- matrix(0, nrow = 2, ncol = 2)          # no. of rejections (at 5% and 1%, columns) for SN and KS-test statistics (rows) from Shao & Zhang (JASA 2010)
    TeSta.SZ<- array(0, dim = c(AnzCPTs))             # test statistics for Shao & Zhang (JASA 2010) test

    for(i in 1:AnzCPTs){
        # Replication of tables
      
        # Results under the null

        # Table \ref{tab:SIM500a2000}
        if(d == 1){
            data <- abs(arima.sim(list(order = c(0, 0, 2), ma = c(0.9, 0.7)), n = AnzSim, innov = rt(n = AnzSim , df = 3)))
        }
        if(d == 2){
          data <- abs(garch.sim(alpha = c(0.01, 0.8), beta = 0.0, n = AnzSim, rnd = rnorm, ntrans = 100))
        }

        # Results under the alternative

        # Table \ref{tab:SIM500a2000_break}
        if(d == 3) {
          data <- abs(armasim_errorbreak(AR_coeff = 0, MA_coeff = c(0.9, 0.7), t_star = 0.5, n = AnzSim, ntrans = 100, df_pre = 3, varmult = 0.25))
        }
        if(d == 4) {
          data <- abs(armasim_errorbreak(AR_coeff = 0, MA_coeff = c(0.9, 0.7), t_star = 0.75, n = AnzSim, ntrans = 100, df_pre = 3, varmult = 0.25))
        }
        if(d == 5) {       
          data <- abs(archsim(alpha_pre = c(0.01, 0.8, 0.0), alpha_post = c(0.01, 0.3, 0.0), t_star = 0.5, n = AnzSim, ntrans = 100))
        }
        if(d == 6) {   
          data <- abs(archsim(alpha_pre = c(0.01, 0.8, 0.0), alpha_post = c(0.01, 0.3, 0.0), t_star = 0.75, n = AnzSim, ntrans = 100))
        }
        if(d == 7) {       
          data <- abs(archsim(alpha_pre = c(0.01, 0.8, 0.0), alpha_post = c(0.0211, 0.3, 0.0), t_star = 0.5, n = AnzSim, ntrans = 100))
        }
        if(d == 8) {   
          data <- abs(archsim(alpha_pre = c(0.01, 0.8, 0.0), alpha_post = c(0.0211, 0.3, 0.0), t_star = 0.75, n = AnzSim, ntrans = 100))
        }
        

        # Teststatistiken
#         Teststatistik <- Teststatistik(X = data, t0 = 0.2, k, p)
# 
#         TeSta[ , i]    <- Teststatistik$EQ
# 
#         for(w in 1:AnzAntEx){                         # loop over k's
#                 # rejections with TI test
#                 if(Teststatistik$TI[w] > 80.21){      # at 5%-level
#                     rej[1, w] <- rej[1, w] + 1
#                 }
#                 if(Teststatistik$TI[w] > 147.8){      # at 1%-level
#                     rej[3, w] <- rej[3, w] + 1
#                 }
#                 # rejections with EQ test
#                 if(Teststatistik$EQ[w] > 80.21){      # at 5%-level
#                     rej[2, w] <- rej[2, w] + 1
#                 }
#                 if(Teststatistik$EQ[w] > 147.8){      # at 1%-level
#                     rej[4, w] <- rej[4, w] + 1
#                 }
#         }
        #######################################
        # for bootstrap only
        
        Teststatistik <- TS.bootstrap(X = data, t0 = 0.2, k, p)
        TeSta[ , i]   <- Teststatistik$EQ
        
        for(w in 1:AnzAntEx){                         # loop over k's
          # rejections with TI test
          if(Teststatistik$TI[w] > 1.821){      # at 5%-level
            rej[1, w] <- rej[1, w] + 1
          }
          if(Teststatistik$TI[w] > 2.653){      # at 1%-level
            rej[3, w] <- rej[3, w] + 1
          }
          # rejections with EQ test
          if(Teststatistik$EQ[w] > 1.821){      # at 5%-level
            rej[2, w] <- rej[2, w] + 1
          }
          if(Teststatistik$EQ[w] > 2.653){      # at 1%-level
            rej[4, w] <- rej[4, w] + 1
          }
        }
        
        #######################################
        
        
        if(comparison == 'yes'){
            teststatistic.SZ <- TS.SZ(X = data, t0 = 0.2, p)
            TeSta.SZ[i]      <- teststatistic.SZ$SN
            if(teststatistic.SZ$SN > 80.21){        # at 5%-level
                rej.SZ[1, 1] <- rej.SZ[1, 1] + 1
            } 
            if(teststatistic.SZ$SN > 147.8){        # at 1%-level
                rej.SZ[1, 2] <- rej.SZ[1, 2] + 1
            }
            if(teststatistic.SZ$KS > 5.54){         # at 5%-level (critical values for KS from Hoga (2014, Table 1))
                rej.SZ[2, 1] <- rej.SZ[2, 1] + 1
            } 
            if(teststatistic.SZ$KS > 8.61){         # at 1%-level
                rej.SZ[2, 2] <- rej.SZ[2, 2] + 1
            }
             
        }
    }
    return(list(TI = matrix((1/AnzCPTs) * cbind(t(rej[1, ]), t(rej[3, ])), nrow=AnzAntEx, dimnames = list(k, c("5% H", "1% H"))), matrix(c(rowMeans(TeSta[ , ])), nrow = AnzAntEx, dimnames = list(k, c("H"))), EQ = matrix((1/AnzCPTs) * cbind(t(rej[2, ]), t(rej[4, ])), nrow=AnzAntEx, dimnames = list(k, c("5% H", "1% H"))), SZ = matrix((1/AnzCPTs) * rej.SZ, nrow=2, dimnames = list(c("SN", "KS"), c("5%", "1%"))), matrix(mean(TeSta.SZ), nrow = 1, dimnames = list(c("Test.statistic.SZ")))))
}               

# Routine nach Shao & Zhang (JASA 2010): Calculates test statistic
# with self-normalizer and subsampling estimator 

TS.SZ <- function(X, t0, p){
      n <- length(X)                   # t0\in(0,1/2) ist die strukturbruchlose Zeit, AntEx ein Vektor mit den extremen Stichprobenanteilen                                                    
      Tq_0t <- array(0, dim=c(n))      # array for forward quantile estimates
          
      for(i in (ceiling(n*t0) : n)){
          Tq_0t[i] <- quantile(X[1 : i], (1-p))  # you can insert any estimator of any population parameter here
      } 
      
      # Speichern der Werte für Tq(t,1)
      Tq_t1 <- array(0, dim=c(n))      # array for backward quantile estimates
      
      for(i in (floor(n * (1-t0)) : 1)){
          Tq_t1[i] <- quantile(X[i : n], (1-p))    
      }
      
      # Calculate self-normalized test statistic
      maxQ_Tq.SN <- 0                   # vector for suprema of self-normalized test statistic
      Q_Tq    <- array(0, dim = c(n))   # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden
      SN <- array(0, dim = c(n))        # Calculating the self-normalizer 

			for(i in ceiling(n*t0) : floor(n*(1-t0))){      # loop over tome
          t <- i/n				                            # discretize time
          Q_Tq[i] <- n * (t * (1-t) * (Tq_0t[i] - Tq_t1[i]))^2  # Q_1^{SN}
          SN[i] <- sum(((ceiling(n*t0) : i)/n * (Tq_0t[ceiling(n*t0) : i] - Tq_0t[i]))^2)    # calculate self-normalizer
          SN[i] <- SN[i] + sum((((n-(i : floor(n*(1-t0))) + 1) / n) * (Tq_t1[i : floor(n*(1-t0))] - Tq_t1[i]))^2)
      }
      maxQ_Tq.SN <- max(Q_Tq[ceiling(n*t0) : floor(n*(1-t0))] / SN[ceiling(n*t0) : floor(n*(1-t0))]) # elementwise division       
      
      # Calculate KS-test statistic
      l     <- 2 * n^(1/3)                    # suggested by Shao & Zhang (JASA 2010, p.1234)
      sn.l  <- floor(n / l)
      theta <- array(0, dim = c(sn.l))
      X.mat <- matrix(X, nrow = l)            # divide sample into subsamples of length l ...
      theta <- apply(X.mat, 2, quantile, 1-p) # ... and calculate (1-p)-quantile for all subsamples
      sigma2.med  <- l / sn.l * sum((theta - mean(theta))^2)
      maxQ_Tq.sub <- max(Q_Tq[ceiling(n*t0) : floor(n*(1-t0))]) / sigma2.med 
      
      return(list(SN = maxQ_Tq.SN, KS = maxQ_Tq.sub))
}

# Routine that calculates values of test statistic
   
Teststatistik <- function(X, t0, k, p){  # Routine, die zu gegebenem Datensatz X, die Teststatistiken ausrechnet
      n <- length(X)                                      # t0\in(0,1/2) ist die strukturbruchlose Zeit, AntEx ein Vektor mit den extremen Stichprobenanteilen
      AnzAntEx <- length(k)                               # a_n = n/fracN bezeichnet die Anzahl an Subsamples für die Schätzung der asympt. Var.
      k <- floor(k * n)                                   # (1-p)-Quantil soll auf Strukturbruch hin überprüft werden.
                                                    
      # Speichern der EQ Schätzwerte für T(0,t)
      T_0t  <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for tail index
      Tq_0t <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for extreme quantiles

      Y    <- sort.int(X[1 : floor(n * t0)], decreasing = TRUE, method = "quick")
          
      for(i in (floor(n * t0) + 1) : n){
          Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")      # sort X_1,...,X_[floor(n*t0)-1+i] in descending order

          for(j in (1 : AnzAntEx)){      
              # Hill-Schätzer
              T_0t[i, j]  <- 1/(k[j]*(i/n)) * sum(log(Y[1 : (floor(k[j]*(i/n))+1)] / Y[floor(k[j]*(i/n))+1])) # Hill estimator
              Tq_0t[i, j] <- Y[floor(k[j]*(i/n))+1] * ((n*p) / k[j])^(-T_0t[i, j])                            # Hill EQ estimator
          }
      } 
      
      T_t1  <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for tail index
      Tq_t1 <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for extreme quantiles
 
      Y    <- sort.int(X[(floor(n * (1 - t0)) + 1) : n], decreasing = TRUE, method = "quick")
          
      for(i in floor(n * (1 - t0)) : 1){      # loop over time
          Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")  # Die [ns]+1 bis [nt] ersten Werte der Daten X[ns]+1,...,X[nt] absteigend sortiert
          for(j in (1 : AnzAntEx)){           # loop over k's
              # Hill estimator
              T_t1[i, j] <- 1/(k[j]*(1 - (i-1)/n)) * sum(log(Y[1 : (floor(k[j]*(1 - (i-1)/n))+1)] / Y[floor(k[j]*(1 - (i-1)/n)) + 1]))
              Tq_t1[i, j] <- Y[floor(k[j]*(1 - (i-1)/n))+1] * ((n*p) / k[j])^(-T_t1[i, j])
          }
      }

      # Calculate SN test statistics
      maxQ_T  <- array(NA, dim = c(AnzAntEx))       # Save realised suprema for different no. of k's
      maxQ_Tq <- array(NA, dim = c(AnzAntEx)) 
      Q_T     <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer
      Q_Tq    <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer   
      SN      <- array(NA, dim = c(n, AnzAntEx))    # Calculating the self-normalizer "AnzAntEx" versch. k's für die "5" Schätzer
      SN_T    <- array(NA, dim = c(n, AnzAntEx))    # Calculating the self-normalizer "AnzAntEx" versch. k's für die "5" Schätzer
      l_sq    <- log(k / (n*p))^2

      
			for(i in (floor(n*t0)+1) : floor(n*(1-t0))){      # loop over time
    	  t <- i/n				                                # discretize time "t"
    		Q_Tq[i, ] <- k / l_sq * (t * (1-t) * log(Tq_0t[i, ] / Tq_t1[i, ]))^2 
    		Q_T[i, ]  <-       k  * (t * (1-t) * (T_0t[i, ] - T_t1[i, ]))^2
    				  
        for(m in 1:AnzAntEx){
            SN[i, m] <- 1/n * k[m] / l_sq[m] * sum((((floor(n*t0)+1) : i)/n * log(Tq_0t[(floor(n*t0)+1) : i, m] / Tq_0t[i, m]))^2)
            SN[i, m] <- SN[i, m] + 1/n * k[m] / l_sq[m] * sum((((n-(i : floor(n*(1-t0))) + 1) / n) * log(Tq_t1[i : floor(n*(1-t0)), m] / Tq_t1[i, m]))^2)
                  
            SN_T[i, m] <- 1/n * k[m] * sum((((floor(n*t0)+1) : i)/n * (T_0t[(floor(n*t0)+1) : i, m] - T_0t[i, m]))^2)
            SN_T[i, m] <- SN_T[i, m] + 1/n * k[m] * sum((((n-(i : floor(n*(1-t0))) + 1) / n) * (T_t1[i : floor(n*(1-t0)), m] - T_t1[i, m]))^2)
        }
			}
      
      maxQ_T  <- colMaxs( Q_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ] / SN_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ])   # elementwise division
      maxQ_Tq <- colMaxs(Q_Tq[(floor(n * t0) + 1) : floor(n * (1-t0)), ] /   SN[(floor(n * t0) + 1) : floor(n * (1-t0)), ])   # elementwise division

      
#############################################################
      
#       # leave this uncommented unless you want to reproduce Figure 1
#       # Figure 1 \ref{fig:TI_WTI}
#       # calculate rolling window estimates of tail index and std. dev.
#       T_t0  <- array(0, dim=c(n, AnzAntEx, 5))    # Hill estimates
#       T_01  <- array(0, dim=c(AnzAntEx, 5))       # Hill estimates based on whole sample
#       Sd <- rep.int(0, n)                         # standard deviation estimates
#       
#       Y    <- sort.int(X[1 : n], decreasing = TRUE, method = "quick")
#       for(j in (1 : AnzAntEx)){
#         # Hill-Schätzer
#         T_01[j, 1]    <- 1/k[j] * sum(log(Y[1 : (floor(k[j])+1)] / Y[floor(k[j])+1]))
#       }
#       
#       for(i in floor(n*t0) : n){
#         Sd[i] <- sd(X[(i - floor(n*t0) + 1) : i])
#       }
#       
#       for(i in floor(n*t0) : n){    # loop through observation period
#         Sd[i] <- sd(X[(i - floor(n*t0) + 1) : i])
#         Y <- sort.int(X[(i - floor(n*t0) + 1) : i], decreasing = TRUE, method = "quick")    # Die Werte der Daten X_1,...,X_[floor(n*t0)-1+i] werden absteigend sortiert
#         
#         for(j in (1 : AnzAntEx)){
#           # Hill-Schätzer
#           T_t0[i, j]    <- 1/(k[j]*t0) * sum(log(Y[1 : (floor(k[j]*t0)+1)] / Y[floor(k[j]*t0)+1]))
#         }
#       }
#     
#     
#       for(i in 1 : AnzAntEx){    # i bestimmt den Wert für k/n, hier k/n = 0.08 wegen i=2
#           split.screen(figs=c(3, 1))
#           par(oma=c(3,1,0,3))
#           screen(1)   # time series plot
#           par(mar = c(2,3,2,3))
#           plot(Datum, -WTI.logloss, type = "l", xlab = "", ylab = "Log-return", xaxs="i", yaxs="i", ylim=c(-0.02+min(-WTI.logloss), 0.02+max(-WTI.logloss)))    # Grafische Veranschaulichung für WTI
#           mtext("Log-return", side=2, line=2)
#           #plot(Datum, X, type = "l", xlab = "", ylab = "Log-return", xaxs="i", yaxs="i", ylim=c(-0.02+min(X), 0.02+max(X)))    # Grafische Veranschaulichung für WTI
#           #abline(h = c(-5, 0, 5, 10))
#           
#           screen(2)   # tests statistic plot
#           par(mar = c(2,3,2,3))
#           Q_H <- rep.int(0, n)
#           Q_MR <- rep.int(0, n)
#           Q_CV <- rep.int(0, n)
#           interval <- (floor(n * t0) + 1) : floor(n * (1-t0))
#           Q_H[interval]  <- Q_Tq[interval, i] / SN[interval, i]
#           
#           plot(Datum, Q_H, type = "l", xlab = "", ylab = "Test statistic", xaxs="i", yaxs="i", ylim = c(0, 1+maxQ_Tq[i, 1]))
#           mtext("Test statistic", side=2, line=2)
#           
#           screen(3)
#           par(mar = c(2,3,2,3))
#           plot(Datum, T_t0[ , i], axes=FALSE, xlab="", ylab="", type="n")
#           lines(Datum[floor(n*t0) : n], T_t0[floor(n*t0) : n, i])
#           abline(h = c(0, T_01[i, 1]), lwd = 0.5)
#           axis(2, ylim=c(0, max(T_t0[ , i])), las=1)  ## las=1 makes horizontal labels
#           mtext("Hill estimates", side=2, line=3)
#           
#           par(new=TRUE)       # Allow a second plot on the same graph
#           plot(Datum, Sd, pch=15,  xlab="", ylab="", ylim=c(0,max(Sd)), axes=FALSE, xaxs="i", type="n", col="red")
#           lines(Datum[floor(n*t0) : n], Sd[floor(n*t0) : n], lty="dotted")
#           
#           print(cbind(Sd, T_t0[, i]))
#           
#           ## a little farther out (line=3) to make room for labels
#           mtext("St.dev. estimates", side=4, line=3)
#           axis(4, ylim=c(0, max(Sd)), las=1)
#           
#           ## Draw the time axis
#           axis.Date(1, at = seq(min(Datum), max(Datum), "years"), pos=0)
#           mtext("Date", side=1, col="black", line=2.5)
#           abline(h = c(0, 0), lwd = 0.5)
#           
#           close.screen(all = TRUE)
#     }
      
#############################################################
       
      return(list(TI = maxQ_T, EQ = maxQ_Tq))
}


#############################################################

TS.bootstrap <- function(X, t0, k, p){  # Routine, die zu gegebenem Datensatz X, die Teststatistiken ausrechnet
  n <- length(X)                                      # t0\in(0,1/2) ist die strukturbruchlose Zeit, AntEx ein Vektor mit den extremen Stichprobenanteilen
  AnzAntEx <- length(k)                               # a_n = n/fracN bezeichnet die Anzahl an Subsamples für die Schätzung der asympt. Var.
  k <- floor(k * n)                                   # (1-p)-Quantil soll auf Strukturbruch hin überprüft werden.
  
  # Speichern der EQ Schätzwerte für T(0,t)
  T_0t  <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for tail index
  Tq_0t <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for extreme quantiles
  
  Y    <- sort.int(X[1 : floor(n * t0)], decreasing = TRUE, method = "quick")
  
  for(i in (floor(n * t0) + 1) : n){
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")      # sort X_1,...,X_[floor(n*t0)-1+i] in descending order
    
    for(j in (1 : AnzAntEx)){      
      # Hill-Schätzer
      T_0t[i, j]  <- 1/(k[j]*(i/n)) * sum(log(Y[1 : (floor(k[j]*(i/n))+1)] / Y[floor(k[j]*(i/n))+1])) # Hill estimator
      Tq_0t[i, j] <- Y[floor(k[j]*(i/n))+1] * ((n*p) / k[j])^(-T_0t[i, j])                            # Hill EQ estimator
    }
  } 
  
  T_t1  <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for tail index
  Tq_t1 <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for extreme quantiles
  
  Y    <- sort.int(X[(floor(n * (1 - t0)) + 1) : n], decreasing = TRUE, method = "quick")
  
  for(i in floor(n * (1 - t0)) : 1){      # loop over time
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")  # Die [ns]+1 bis [nt] ersten Werte der Daten X[ns]+1,...,X[nt] absteigend sortiert
    for(j in (1 : AnzAntEx)){           # loop over k's
      # Hill estimator
      T_t1[i, j] <- 1/(k[j]*(1 - (i-1)/n)) * sum(log(Y[1 : (floor(k[j]*(1 - (i-1)/n))+1)] / Y[floor(k[j]*(1 - (i-1)/n)) + 1]))
      Tq_t1[i, j] <- Y[floor(k[j]*(1 - (i-1)/n))+1] * ((n*p) / k[j])^(-T_t1[i, j])
    }
  }
  
  # Calculate SN test statistics
  maxQ_T  <- array(NA, dim = c(AnzAntEx))       # Save realised suprema for different no. of k's
  maxQ_Tq <- array(NA, dim = c(AnzAntEx)) 
  Q_T     <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer
  Q_Tq    <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer   
  l_sq    <- log(k / (n*p))^2
  
  
  for(i in (floor(n*t0)+1) : floor(n*(1-t0))){      # loop over time
    t <- i/n				                                # discretize time "t"
    Q_Tq[i, ] <- k / l_sq * (t * (1-t) * log(Tq_0t[i, ] / Tq_t1[i, ]))^2 
    Q_T[i, ]  <-       k  * (t * (1-t) * (T_0t[i, ] - T_t1[i, ]))^2
  }
  
  sigma.2 <- numeric(AnzAntEx)
  
  for(j in (1 : AnzAntEx)){
      sigma.2[j] <- k[j] / log(k[j] / (n*p))^2 * var(tsboot(X, statistic = EQ.estimate, k.NUOS = k[j], prob = p, R = 500, l = floor(n / 20), sim = "fixed")$t)
  }      
  
  maxQ_T  <- colMaxs( Q_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ] ) / sigma.2    # elementwise division
  maxQ_Tq <- colMaxs(Q_Tq[(floor(n * t0) + 1) : floor(n * (1-t0)), ] ) / sigma.2    # elementwise division
  
  return(list(TI = maxQ_T, EQ = maxQ_Tq))
}


################################################################################

# data generating processes
        
armasim_errorbreak <- function(AR_coeff, MA_coeff, t_star, n, ntrans, df_pre, varmult){      # ARMA(1,2) coefficients in "XX_coeff" (defined in ?arima)     
    m <- n + ntrans
    e <- numeric(m)                         # Innovationenvektor mit Strukuturbruch
    e[1 : floor(ntrans + n * t_star)]       <- rt(floor(ntrans + n * t_star), df = df_pre)
    e[(floor(ntrans + n * t_star) + 1) : m] <- sqrt(varmult) * rt(ceiling(n * (1 - t_star)), df = df_pre)
    y <- numeric(m)                         # intializes a vector of y's to be m long
    y[1:2] <- rt(n = 2, df = df_pre)        # start values
    for(i in 3:m){
      y[i] <- AR_coeff * y[i-1] + e[i] + MA_coeff[1] * e[i-1] + MA_coeff[2] * e[i-2]    # Simulate ARMA(1, 2) process  
    }
    y <- y[(ntrans + 1) : m]                # Drop the first ntrans and just call it y again
    return(y)
}

archsim <- function(alpha_pre, alpha_post, t_star, n, ntrans){   #ARCH(1) coefficients - alpha_pre = c(alpha0, alpha1)
    m <- n + ntrans
    e <- rnorm(m, mean = 0, sd = 1)  
    y <- numeric(m)                                             # intializes a vector of y's to be m long
    y[1] <- rnorm(n = 1, mean = 0, sd = 1)                      # start value
    for(i in 2:floor(ntrans + n * t_star)){
      y[i] <- e[i] * sqrt(alpha_pre[1] + alpha_pre[2] * y[i-1]^2 + alpha_pre[3] * e[i-1]^2)        # Simulate ARCH(1) process
    }
    for(i in (floor(ntrans + n * t_star)+1) : m){
      y[i] <- e[i] * sqrt(alpha_post[1] + alpha_post[2] * y[i-1]^2 + alpha_pre[3] * e[i-1]^2)  # Simulate ARCH(1) process
    }
    y <- y[(ntrans + 1) : m]                # Drop the first ntrans and just call it y again
    return(y)
}


################################################################################

# Calculate critical values for SZ test

quantilsup.KS <- function(AnzSim, AnzGitter, t0, Quantile){   # AnzSim = Anzahl Wiener-Pfade, AnzGitter = Diskretisierungspunkte der Wiener-Pfade, Quantile = Vektor der zu berechnenden Quantile
      RealZVseqc  <- rep.int(0, AnzSim)                       # Analog!
      GitterZVseqc<- rep.int(0, AnzGitter+1)
      for(j in 1:AnzSim){
            B <- rwiener(end = 1, frequency = AnzGitter)                # Simulation eines Standard-Wiener-Prozesses
            B <- c(0, B)                                                # Ergänzen um den Beginn in 0,...
            for(i in ceiling(AnzGitter*t0):floor(AnzGitter*(1-t0))){
                  t <- i/AnzGitter
                  GitterZVseqc[i+1] <- (min(1, (1-t)/t)*B[i+1] - min(1, t/(1-t))*(B[AnzGitter+1]-B[i+1]))^2  # analog! 
            }
            RealZVseqc[j]<- max(GitterZVseqc)             
      }
      return(quantile(RealZVseqc, probs = Quantile))          
}     # quantile-Funktion liefert die gewuenschten Quantile

start <- Sys.time()
Quantilvektor.KS <- quantilsup(AnzSim = 20000, AnzGitter = 10000, t0=0.2, Quantile = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99))
print(Sys.time() - start)
print(Quantilvektor.KS)

################################################################################

# Calculate critical values for SN-test with 1 possible breakpoint (others taken from Hoga, 2014)

quantilsup.SN1 <- function(AnzSim, AnzGitter, t0, Quantile){   # AnzSim = Anzahl Wiener-Pfade, AnzGitter = Diskretisierungspunkte der Wiener-Pfade, Quantile = Vektor der zu berechnenden Quantile
      RealZVseqc  <- rep.int(0, AnzSim)                    # Hier werden die "AnzSim"-vielen Realisationen der Grenz-Zufallsvariablen sup_{t\in[t0,1-t0]}(B(t)-tB(1))^{2} gespeichert
      GitterZVseqc<- rep.int(0, AnzGitter)                 # An Gitterpunkten (in t) wird (B(t)-tB(1))^{2} ausgewertet
      numerator   <- rep.int(0, AnzGitter)
      denom1      <- rep.int(0, AnzGitter)
      denom2      <- rep.int(0, AnzGitter)
      for(j in 1:AnzSim){
            B <- rwiener(end = 1, frequency = AnzGitter)                # Simulation eines Standard-Wiener-Prozesses
            for(i in ceiling(AnzGitter*t0) : floor(AnzGitter*(1-t0))){
                  t <- i/AnzGitter
                  numerator[i] <- (B[i] - t * B[AnzGitter])^2
                  denom1[i]    <- 1/AnzGitter * sum((B[ceiling(AnzGitter*t0) : i] - (ceiling(AnzGitter*t0) : i) / (AnzGitter*t) * B[i])^2)
                  denom2[i]    <- 1/AnzGitter * sum((B[AnzGitter] - B[i : floor(AnzGitter*(1-t0))] - (1 - (i : floor(AnzGitter*(1-t0))) / AnzGitter) / (1-t) * (B[AnzGitter] - B[i]))^2)
                  GitterZVseqc[i] <- numerator[i] / (denom1[i] + denom2[i])
            }
            RealZVseqc[j]<- max(GitterZVseqc)      # save realisation
      }
      return(quantile(RealZVseqc, probs = Quantile))
}     # quantile-function gives required quantiles

start <- Sys.time()
Quantilvektor.SN1 <- quantilsup.SN1(AnzSim = 30000, AnzGitter = 10000, t0=0.2, Quantile = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99))
print(Sys.time() - start)
print(Quantilvektor.SN)

################################################################################

# Calculate critical values for SN-test with 2 possible breakpoints

quantilsup.SN2 <- function(AnzSim, AnzGitter, t0, Quantile){   # AnzSim = Anzahl Wiener-Pfade, AnzGitter = Diskretisierungspunkte der Wiener-Pfade, Quantile = Vektor der zu berechnenden Quantile
  RealZV   <- rep.int(0, AnzSim)                                # Hier werden die "AnzSim"-vielen Realisationen der Grenz-Zufallsvariablen sup_{t\in[t0,1-t0]}(B(t)-tB(1))^{2} gespeichert
  GitterZV <- array(0, dim = c(AnzGitter, AnzGitter))              # An Gitterpunkten (in t) wird (B(t)-tB(1))^{2} ausgewertet
  t.0      <- floor(t0 * AnzGitter)
  one      <- AnzGitter
  denom1.1 <- rep.int(0, AnzGitter)
  denom2.1 <- rep.int(0, AnzGitter)
  denom1.2 <- rep.int(0, AnzGitter)
  denom2.2 <- rep.int(0, AnzGitter)
  numerator.1<- rep.int(0, AnzGitter)
  numerator.2<- rep.int(0, AnzGitter)
  
  for(rep in 1 : AnzSim){
    B <- rwiener(end = 1, frequency = AnzGitter)                # Simulation eines Standard-Wiener-Prozesses
    
    range.t1 <- (2*t.0) : (one - 4*t.0)
    range.t2 <- (4*t.0) : (one - 2*t.0)
    
    int.t1 <- t.0 : (one - 4*t.0)        # Integrationsbereich für t1
    int.t2 <- (4*t.0) : (one - t.0)      # Integrationsbereich für t2
    denom1.1[int.t1]   <- 1/AnzGitter * (cumsum(B[int.t1]^2) - 2 * B[int.t1] / (int.t1 / AnzGitter) * cumsum(B[int.t1] * (int.t1 / AnzGitter)) + (B[int.t1] / (int.t1 / AnzGitter))^2 * cumsum( (int.t1 / AnzGitter)^2 ))
    W.1mt2             <- ( B[one] - B[int.t2] ) / (one - int.t2) * AnzGitter
    denom2.2[int.t2] <- 1/AnzGitter * ( rev(cumsum( rev((B[one] - B[int.t2])^2) )) - 2* W.1mt2 * rev(cumsum( rev( (B[one] - B[int.t2]) * (1 - int.t2 / AnzGitter) ) )) + W.1mt2^2 * rev(cumsum( rev( (1 - int.t2 / AnzGitter)^2 ) )  ))
    
    for(t.1 in range.t1){
      subrange.t2 <- (t.1 + 2*t.0) : (one - 2*t.0) # range of t2 for specific value of t.1
      
      W.t2mt1 <- (B[(t.1 + t.0) : (one - 2*t.0)] - B[t.1]) / ((t.0 : (one - 2*t.0 - t.1)) / AnzGitter)   # = (W(t2) - W(t1)) / (t2 - t1) for all t2 e [2*t0, 1-t0]
      W.t2    <- B[(t.1 + t.0) : (one - 2*t.0)] # = W(t2) for all t2
      t2      <- ((t.1 + t.0) : (one - 2*t.0)) / AnzGitter # = t2 for all t2 e [t1 + t0, 1 - 2t0]
      int.1   <- 1/AnzGitter * (1 : (one - 3*t.0 - t.1 + 1)) 
      int.s   <- 1/AnzGitter * cumsum( (t.1 : (one - 3*t.0)) / AnzGitter )
      int.s.sq<- 1/AnzGitter * cumsum( ((t.1 : (one - 3*t.0)) / AnzGitter)^2 )
      int.W   <- 1/AnzGitter * cumsum( B[t.1 : (one - 3*t.0)] )
      int.W.sq<- 1/AnzGitter * cumsum( B[t.1 : (one - 3*t.0)]^2 )
      int.sW  <- 1/AnzGitter * cumsum( B[t.1 : (one - 3*t.0)] * (t.1 : (one - 3*t.0)) / AnzGitter )
      int.t2ms<- t2 * int.1 - int.s
      
      numerator.1[subrange.t2] <- (B[t.1] - t.1 / subrange.t2 * B[subrange.t2])^2
      denom2.1[(t.1 + t.0) : (one - 2*t.0)] <- W.t2^2 * int.1 - 2 * W.t2 * int.W - 2 * W.t2 * W.t2mt1 * int.t2ms + int.W.sq + 2 * W.t2mt1 * (t2 * int.W - int.sW) + W.t2mt1^2 * (t2^2 * int.1 - 2 * t2 * int.s + int.s.sq)
      
      # int. 1 ist wie oben
      int.smt1.sq <- 1/AnzGitter * cumsum( ((t.0 : (one - 2*t.0 - t.1)) / AnzGitter)^2 )
      int.smt1    <- 1/AnzGitter * cumsum( (t.0 : (one - 2*t.0 - t.1)) / AnzGitter )
      int.smt1W   <- 1/AnzGitter * cumsum( B[(t.1 + t.0) : (one - 2*t.0)] * (t.0 : (one - 2*t.0 - t.1)) / AnzGitter ) 
      int.W       <- 1/AnzGitter * cumsum( B[(t.1 + t.0) : (one - 2*t.0)] )
      int.W.sq    <- 1/AnzGitter * cumsum( B[(t.1 + t.0) : (one - 2*t.0)]^2 )
      
      numerator.2[(t.1 + t.0) : (one - 2*t.0)] <- (W.t2 - (1 - t2) / (1 - t.1 / AnzGitter) * B[t.1] - (t2 - t.1 / AnzGitter) / (1 - t.1 / AnzGitter) * B[one] )^2
      denom1.2[(t.1 + t.0) : (one - 2*t.0)] <- int.W.sq - 2 * B[t.1] * int.W - 2 * W.t2mt1 * int.smt1W + B[t.1]^2 * int.1 + 2 * B[t.1] * W.t2mt1 * int.smt1 + W.t2mt1^2 * int.smt1.sq
      GitterZV[t.1, subrange.t2] <- numerator.1[subrange.t2] / (denom1.1[t.1] + denom2.1[subrange.t2]) + numerator.2[subrange.t2] / (denom1.2[subrange.t2] + denom2.2[subrange.t2])
    }
    RealZV[rep]<- max(GitterZV)      # save realisation
  }
  return(quantile(RealZV, probs = Quantile))
}     # quantile-function gives required quantiles

start <- Sys.time()
Quantilvektor.SN2 <- quantilsup.SN2(AnzSim = 8000, AnzGitter = 10000, t0=0.1, Quantile = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99))
print(Sys.time() - start)
print(Quantilvektor.SN2)

################################################################################

# Calculate lower bound for Q_2^EQ:

LB.Q2 <- function(X, t0, t1, t2, k, p){
  n      <- length(X)
  k <- floor(k * n)  
  l_knp  <- log(k / (n*p))
  T_0t   <- array(NA, dim=c(n))     # 1 Schätzer für gamma: H
  Tq_0t  <- array(NA, dim=c(n))     # 1 Schätzer für EQ: H
  T_st2  <- array(NA, dim=c(n))     # 1 Schätzer für gamma: H
  Tq_st2 <- array(NA, dim=c(n))     # 1 Schätzer für EQ: H
  T_t1s  <- array(NA, dim=c(n))     # 1 Schätzer für gamma: H
  Tq_t1s <- array(NA, dim=c(n))     # 1 Schätzer für EQ: H
  T_s1   <- array(NA, dim=c(n))     # 1 Schätzer für gamma: H
  Tq_s1  <- array(NA, dim=c(n))     # 1 Schätzer für EQ: H
  
  # denom1
  Y    <- sort.int(X[1 : floor(n * t0)], decreasing = TRUE, method = "quick")
  for(i in (floor(n * t0) + 1) : floor(n * t1)){
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")      # Die Werte der Daten X_1,...,X_[floor(n*t0)-1+i] werden absteigend sortiert
    T_0t[i]  <- 1/(k*(i/n)) * sum(log(Y[1 : (floor(k*(i/n))+1)] / Y[floor(k*(i/n))+1]))     # Hill Schätzer
    Tq_0t[i] <- Y[floor(k*(i/n))+1] * ((n*p) / k)^(-T_0t[i]) 
  }
  # denom2
  Y <- sort.int(X[floor(n * (t2 - t0) + 1) : floor(n * t2)], decreasing = TRUE, method = "quick")
  for(i in (floor(n * (t2 - t0)) : floor(n * t1 + 1))){ 
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")    
    T_st2[i]  <- 1/(k*(length(Y)/n)) * sum(log(Y[1 : (floor(k*(length(Y)/n))+1)] / Y[floor(k*(length(Y)/n))+1]))     # Hill Schätzer
    Tq_st2[i] <- Y[floor(k*(length(Y)/n))+1] * ((n*p) / k)^(-T_st2[i]) 
  }
  # denom3
  Y <- sort.int(X[floor(n * t1 + 1) : floor(n * (t1 + t0))], decreasing = TRUE, method = "quick")
  for(i in (floor(n * (t1 + t0) + 1) : floor(n * t2))){ # 
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")    
    T_t1s[i]  <- 1/(k*(length(Y)/n)) * sum(log(Y[1 : (floor(k*(length(Y)/n))+1)] / Y[floor(k*(length(Y)/n))+1]))     # Hill Schätzer
    Tq_t1s[i] <- Y[floor(k*(length(Y)/n))+1] * ((n*p) / k)^(-T_t1s[i]) 
  }
  # denom4
  Y <- sort.int(X[floor(n * (1 - t0) + 1) : n], decreasing = TRUE, method = "quick")
  for(i in (floor(n * (1 - t0)) : floor(n * t2))){ # 
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")    
    T_s1[i]  <- 1/(k*(length(Y)/n)) * sum(log(Y[1 : (floor(k*(length(Y)/n))+1)] / Y[floor(k*(length(Y)/n))+1]))     # Hill Schätzer
    Tq_s1[i] <- Y[floor(k*(length(Y)/n))+1] * ((n*p) / k)^(-T_s1[i]) 
  }
  
  numerator1 <- ( t1 * (t2 - t1) / t2 * log( Tq_0t[floor(n * t1)] / Tq_t1s[floor(n * t2)] ) )^2
  numerator2 <- ( (t2 - t1) * (1 - t2) / (1 - t1) * log( Tq_t1s[floor(n * t2)] / Tq_s1[floor(n * t2)] ) )^2
  denom1     <- 1 / n * sum( ( ((floor(n * t0) + 1) : floor(n * t1)) / n * log( Tq_0t[(floor(n * t0) + 1) : floor(n * t1)] / Tq_0t[floor(n * t1)] ) )^2 )
  denom2     <- 1 / n * sum( ( (t2 - (floor(n * t1 + 1) : floor(n * (t2 - t0))) / n) * log( Tq_st2[floor(n * t1 + 1) : floor(n * (t2 - t0))] / Tq_st2[floor(n * t1 + 1)] ) )^2 )
  denom3     <- 1 / n * sum( ( ((floor(n * (t1 + t0) + 1) : floor(n * t2))/n - t1) * log( Tq_t1s[floor(n * (t1 + t0) + 1) : floor(n * t2)] / Tq_t1s[floor(n * t2)] )  )^2 )
  denom4     <- 1 / n * sum( ( 1 - (floor(n * t2) : floor(n * (1 - t0)))/n ) * (log( Tq_s1[floor(n * t2) : floor(n * (1 - t0))] / Tq_s1[floor(n * t2)] ) )^2 )
  
  return( numerator1 / (denom1 + denom2) + numerator2 / (denom3 + denom4) )
} 

################################################################################

# Figure \ref{fig:qsp}: Quantile stability plot

qsp <- function(X, t0, k, p = c(0.005, 0.01, 0.05, 0.1)){
  n    <- length(X) 
  # calculate rolling window estimates of tail index and extreme quantiles
  T_t0  <- array(0, dim=c(n))    # Hill estimates
  Tq_t0 <- array(0, dim=c(n, length(p)))
  
  for(i in floor(n*t0) : n){    # loop through observation period
    Y <- sort.int(X[(i - floor(n*t0) + 1) : i], decreasing = TRUE, method = "quick") - min(X)   # sort X_1,...,X_[floor(n*t0)-1+i] in descending order
    T_t0[i]    <- 1 / (k*t0) * sum(log(Y[1 : (floor(k*t0)+1)] / Y[floor(k*t0)+1])) # Hill estimate of \gamma
    for(j in 1 : length(p)){
      Tq_t0[i, j] <- Y[floor(k*t0)+1] * ((n*p[j]) / k)^(-T_t0[i]) + min(X)    # Weissman estimate (using Hill estimator)
    }
  }
  return(Tq_t0)
#   plot(Datum, Tq_t0[, 1], type = "n", xlab = "", ylab = "", xaxs="i", yaxs="i", ylim = c(0, 0.01+max(Tq_t0[, ])))
#   mtext("Quantiles", side=2, line=2)
#   lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 1], lty = "solid")     
#   lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 2], lty = "longdash")      
#   lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 3], lty = "dotdash") 
#   lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 4], lty = "twodash")
#   legend("bottomleft", c("0.005", "0.01", "0.05", "0.1"), lty=c("solid", "longdash", "dotdash", "twodash"), cex = 1)
}

################################################################################

# Section 4: Real world example WTI log-returns

WTI.data<- read.csv("WTI.csv", header=FALSE)
WTI     <- WTI.data$V2[251 : 2291]          # WTI prices from 1986-12-31 to 1994-12-30
Datum   <- WTI.data$V1[252 : 2291]          # read in dates and ...
Datum   <- as.Date(Datum)                   # ... convert dates.
WTI.logloss       <- -diff(log(WTI))        # Investigate losses
WTI.logloss.shift <- WTI.logloss - min(WTI.logloss)
n                 <- length(WTI.logloss)

# Replication of Figure 3 \ref{fig:qsp}
t0    <- 0.2
Tq_t0 <- qsp(X = WTI.logloss.shift, t0, k = 326) + min(WTI.logloss)  #  k = 326 corresponds to k/n = 0.16
par(mar = c(2,3,1,1))
plot(Datum, Tq_t0[, 1], type = "n", xlab = "", ylab = "", xaxs="i", yaxs="i", ylim = c(0, 0.01+max(Tq_t0[, ])))
mtext("Quantiles", side=2, line=2)
lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 1], lty = "solid")     
lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 2], lty = "longdash")      
lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 3], lty = "dotdash") 
lines(Datum[floor(n*t0) : n], Tq_t0[floor(n*t0) : n, 4], lty = "twodash")
legend("bottomleft", c("0.005", "0.01", "0.05", "0.1"), lty=c("solid", "longdash", "dotdash", "twodash"), cex = 1)

# Replication of Figure 1 \ref{fig:TI_WTI}
# uncomment last lines from function Teststatistik and issue command ...
Teststatistik(X = WTI.logloss.shift, t0 = 0.2, k = c(0.16, 0.16), p = 0.005)

# Result for Q_2^EQ
(LB <- LB.Q2(X = (WTI.logloss.shift), t0=0.1, t1=0.44, t2 = 0.57, k = 0.16, p = 0.005)) # = 149.72
(LB <- LB.Q2(X = (WTI.logloss.shift), t0=0.1, t1=0.45, t2 = 0.57, k = 0.16, p = 0.1))   # = 56.42

# Replication of Figure 2 \ref{fig:CPT_WTI}

# Test, if \gamma>0 with MoM plot

MoM <- function(X, k.max){
    n     <- length(X)
    Y     <- sort.int(X, decreasing = TRUE, method = "quick")
    Y     <- Y[Y[] > 0]          # Only consider positive values (bc. of log-transformation below)
    log.Y    <- log(Y)
    log.Y.sq <- log.Y^2
    Hill  <- cumsum(log.Y[1 : (k.max - 1)]) / (1 : (k.max - 1)) - log.Y[2:k.max]   # Hill estimates for k = 2, ..., k.max
    one   <- cumsum(log.Y.sq[1 : (k.max - 1)]) / (1 : (k.max - 1))
    two   <- 2 * log.Y[2 : k.max] * cumsum(log.Y[1 : (k.max - 1)]) / (1 : (k.max - 1))
    three <- log.Y.sq[2 : k.max]
    square<- one - two + three
    gamma_hat <- Hill[2 : (k.max-1)] + 1 - 0.5 / (1 - ((Hill[2 : (k.max-1)])^2) / square[2 : (k.max-1)])   # MoM estimates for k=3, ..., k.max-1
    return(cbind(gamma_hat, Hill[2 : (k.max-1)]))
}

k.min <- 10                 # minimal # of order statistics k for plot
k.max <- 700                # maximal # of order statistics k for plot
k <- k.min : (k.max-2)

X   <- WTI.logloss          # non-shifted time series WTI
n   <- length(X)
H.M <- MoM(X, k.max)        # returns Hill & MoM estimates for k=3, ..., k.max-1
H   <- H.M[, 2]             # Hill and ...
MM  <- H.M[, 1]             # ... MoM estimates for raw data
# For shifted data
X.s   <- X - min(X)           # geshiftete Zeitreihe
H.M.s <- MoM(X.s, k.max)
H.s   <- H.M.s[, 2]
MM.s  <- H.M.s[, 1]

# Calculate asymptotic variance for shifted data

H_asvar <- function(X, t0, k){  # calculate asymptotic variance according to Hoga(2016+, Theorem 2) 
  n <- length(X)
  AnzAntEx <- length(k)                               
  t_n.var <- 50 / n             # t_n.var = Stichprobenanteil für Varianzschätzung
  
  H_0t  <- array(0, dim=c(n, AnzAntEx))      # Hill estimator
  
  Y    <- sort.int(X[1 : floor(n * t_n.var)], decreasing = TRUE, method = "quick")
  
  for(i in (floor(n * t_n.var) + 1) : n){
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")      # sort X_1,...,X_[floor(n*t0)-1+i] in descending order
   
    for(j in (1 : AnzAntEx)){      
      # Hill estimator
      H_0t[i, j]  <- 1/(k[j]*(i/n)) * sum(log(Y[1 : (floor(k[j]*(i/n))+1)] / Y[floor(k[j]*(i/n))+1]))     # Hill estimator
    }
  }
  
  # calculate asymptotic variance according to Hoga(2016+, Theorem 2) 
  
  widehat_sigma_sq <- array(NA, dim=c(AnzAntEx))    
  Ti <- sweep(H_0t[(floor(n * t_n.var) + 1) : n, ], 2, H_0t[n, ], '-')^2
  widehat_sigma_sq <- colSums(Ti) 
  #widehat_sigma_sq <- 1/log(1/t_n) * k/n * widehat_sigma_sq[ ]   # without finite sample correction for log(1/t_n)
  widehat_sigma_sq <- 1/(sum(1 / ((floor(n * t_n.var) + 1) : n)) - (1 - t_n.var)) * k / n * widehat_sigma_sq  # with finite sample correction for log(1/t_n)   
  widehat_sigma_sq <- pmax(1 * H_0t[n, ]^2, widehat_sigma_sq)  # As.Var. = 1
  
  return(widehat_sigma_sq)
}

DepCoeff_H<- H_asvar(X.s, t0 = 0.2, k)
CI.min    <- H.s[k] - 1.96 * sqrt(DepCoeff_H / k)     # lower bound of c.i. for \gamma
CI.max    <- H.s[k] + 1.96 * sqrt(DepCoeff_H / k)     # upper bound for c.i. for \gamma

split.screen(figs=c(1, 2))
screen(1)
plot(k, MM[k], type = "l", xlab = "(a)", ylab = "", ylim=c(0.0, 0.6), asp=1000)
lines(k, H[k], lty = "dotted")
screen(2)
plot(k, MM.s[k], type = "l", xlab = "(b)", ylab = "", ylim=c(0.0, 0.6), asp=1000)
lines(k, H.s[k], lty = "dotted")
polygon(c(k, rev(k)), c(CI.max, rev(CI.min)), density = 20, angle = 45, lty = "dotted")
close.screen(all = TRUE)

################################################################################

# Plausibility check on results: Fit a GARCH(1,1)-model and inspect residuals!

library(fGarch)
set.seed(123)
GARCH.fit <- garchFit(data = WTI.logloss, cond.dist="std")   # fit a GARCH(1,1) with t_{n}-innovations 
summary(GARCH.fit)  # LB test: no correlation left in (raw and squared) standardized residuals
coef(GARCH.fit)
spec = garchSpec(model = list(omega = coef(GARCH.fit)[2], alpha = coef(GARCH.fit)[3], beta = coef(GARCH.fit)[4], shape = coef(GARCH.fit)[5]), cond.dist = "std")
plot(residuals(GARCH.fit) / volatility(GARCH.fit))  # standardized residuals look good!

sd.resid <- residuals(GARCH.fit) / volatility(GARCH.fit)
MoM(abs(sd.resid), k.max = 700)

GARCH.data <- garchSim(spec, n = length(WTI.logloss))
GARCH.fit <- garchFit(data = GARCH.data, cond.dist="std") 
Datum <- Sys.Date() - ((length(WTI.logloss) - 1) : 0)
Teststatistik(X = GARCH.data - min(GARCH.data), t0 = 0.2, k = c(0.04,0.04, 0.08, 0.16, 0.24), p = 0.005)

################################################################################

# Section 4: Asian crisis of 1997-98

# Load data from 3 stock exchanges: SET, KLCI, JCI

# SET: Stock Exchange of Thailand 
SET.data <- read.csv("SET_WSJ.csv", header=TRUE)   # make sure csv-file is in working directory!
SET      <- rev(SET.data$Close)                    # SET prices from 1995-01-02 to 1998-10-16
n.SET    <- length(SET) - 1                          # = 931
date.SET <- rev(SET.data$Date)[-1]                 # read in dates and ...
date.SET <- as.Date(date.SET, "%m/%d/%y")          # ... convert dates.
SET.logloss       <- -diff(log(SET))               # Investigate losses
SET.logloss.shift <- SET.logloss - min(SET.logloss)

# KLCI: Kuala Lumpur Composite Index
KLCI.data <- read.csv("KLCI_WSJ.csv", header=TRUE)  # make sure csv-file is in working directory!
KLCI      <- rev(KLCI.data$Close)                   # SET prices from 1995-01-02 to 1998-10-16
n.KLCI    <- length(KLCI) - 1                         # = 938
date.KLCI <- rev(KLCI.data$Date)[-1]                # read in dates and ...
date.KLCI <- as.Date(date.KLCI, "%m/%d/%y")         # ... convert dates.
KLCI.logloss       <- -diff(log(KLCI))              # Investigate losses
KLCI.logloss.shift <- KLCI.logloss - min(KLCI.logloss)

# JCI: Jakarta Composite Index
JCI.data <- read.csv("JCI_WSJ.csv", header=TRUE)  # make sure csv-file is in working directory!
JCI      <- rev(JCI.data$Close)                   # SET prices from 1995-01-02 to 1998-10-16
n.JCI    <- length(JCI) - 1                         # = 936
date.JCI <- rev(JCI.data$Date)[-1]                # read in dates and ...
date.JCI <- as.Date(date.JCI, "%m/%d/%y")         # ... convert dates.
JCI.logloss       <- -diff(log(JCI))              # Investigate losses
JCI.logloss.shift <- JCI.logloss - min(JCI.logloss)

# Replication of Figure 4 \ref{fig:qsp_Asian}
t0       <- 0.2
qsp.SET  <- qsp(X = SET.logloss,  t0, k = floor(0.16 * n.SET))   
qsp.KLCI <- qsp(X = KLCI.logloss, t0, k = floor(0.16 * n.KLCI))  
qsp.JCI  <- qsp(X = JCI.logloss,  t0, k = floor(0.16 * n.JCI))  

split.screen(figs=c(3, 1))
par(oma=c(1,1,0,1))

screen(1)   # qsp plot SET
par(mar = c(2,3,3,3))
plot(date.SET, qsp.SET[, 1], main = "SET",type = "n", xlab = "", ylab = "", xaxs="i", yaxs="i", ylim = c(0, 0.01+max(qsp.SET[, ])))
main = "SET"
lines(date.SET[floor(n.SET*t0) : n.SET], qsp.SET[floor(n.SET*t0) : n.SET, 1], lty = "solid")     
lines(date.SET[floor(n.SET*t0) : n.SET], qsp.SET[floor(n.SET*t0) : n.SET, 2], lty = "longdash")      
lines(date.SET[floor(n.SET*t0) : n.SET], qsp.SET[floor(n.SET*t0) : n.SET, 3], lty = "dotdash") 
lines(date.SET[floor(n.SET*t0) : n.SET], qsp.SET[floor(n.SET*t0) : n.SET, 4], lty = "twodash")
legend("bottomleft", c("0.005", "0.01", "0.05", "0.1"), lty=c("solid", "longdash", "dotdash", "twodash"), cex = 0.8)

screen(2)   # qsp plot KLCI
par(mar = c(2,3,3,3))
plot(date.KLCI, qsp.KLCI[, 1], main = "KLCI", type = "n", xlab = "", ylab = "", xaxs="i", yaxs="i", ylim = c(0, 0.01+max(qsp.KLCI[, ])))
lines(date.KLCI[floor(n.KLCI*t0) : n.KLCI], qsp.KLCI[floor(n.KLCI*t0) : n.KLCI, 1], lty = "solid")     
lines(date.KLCI[floor(n.KLCI*t0) : n.KLCI], qsp.KLCI[floor(n.KLCI*t0) : n.KLCI, 2], lty = "longdash")      
lines(date.KLCI[floor(n.KLCI*t0) : n.KLCI], qsp.KLCI[floor(n.KLCI*t0) : n.KLCI, 3], lty = "dotdash") 
lines(date.KLCI[floor(n.KLCI*t0) : n.KLCI], qsp.KLCI[floor(n.KLCI*t0) : n.KLCI, 4], lty = "twodash")
legend("bottomleft", c("0.005", "0.01", "0.05", "0.1"), lty=c("solid", "longdash", "dotdash", "twodash"), cex = 0.8)

screen(3)   # qsp plot JCI
par(mar = c(2,3,3,3))
plot(date.JCI, qsp.JCI[, 1], main = "JCI", type = "n", xlab = "", ylab = "", xaxs="i", yaxs="i", ylim = c(0, 0.01+max(qsp.JCI[, ])))
lines(date.JCI[floor(n.JCI*t0) : n.JCI], qsp.JCI[floor(n.JCI*t0) : n.JCI, 1], lty = "solid")     
lines(date.JCI[floor(n.JCI*t0) : n.JCI], qsp.JCI[floor(n.JCI*t0) : n.JCI, 2], lty = "longdash")      
lines(date.JCI[floor(n.JCI*t0) : n.JCI], qsp.JCI[floor(n.JCI*t0) : n.JCI, 3], lty = "dotdash") 
lines(date.JCI[floor(n.JCI*t0) : n.JCI], qsp.JCI[floor(n.JCI*t0) : n.JCI, 4], lty = "twodash")
legend("bottomleft", c("0.005", "0.01", "0.05", "0.1"), lty=c("solid", "longdash", "dotdash", "twodash"), cex = 0.8)

close.screen(all = TRUE)

# max to min

max(qsp.SET[floor(n.SET*t0) : n.SET, 1])/ min(qsp.SET[floor(n.SET*t0) : n.SET, 1]) # = 2.917543
max(qsp.SET[floor(n.SET*t0) : n.SET, 2])/ min(qsp.SET[floor(n.SET*t0) : n.SET, 2]) # = 2.87602
max(qsp.SET[floor(n.SET*t0) : n.SET, 3])/ min(qsp.SET[floor(n.SET*t0) : n.SET, 3]) # = 2.889368
max(qsp.SET[floor(n.SET*t0) : n.SET, 4])/ min(qsp.SET[floor(n.SET*t0) : n.SET, 4]) # = 3.078299

max(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 1])/ min(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 1]) # = 4.976796
max(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 2])/ min(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 2]) # = 5.004032
max(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 3])/ min(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 3]) # = 5.404087
max(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 4])/ min(qsp.KLCI[floor(n.KLCI*t0) : n.SET, 4]) # = 6.047208

max(qsp.JCI[floor(n.JCI*t0) : n.SET, 1])/ min(qsp.JCI[floor(n.JCI*t0) : n.SET, 1]) # = 5.384289
max(qsp.JCI[floor(n.JCI*t0) : n.SET, 2])/ min(qsp.JCI[floor(n.JCI*t0) : n.SET, 2]) # = 5.118554
max(qsp.JCI[floor(n.JCI*t0) : n.SET, 3])/ min(qsp.JCI[floor(n.JCI*t0) : n.SET, 3]) # = 4.867825
max(qsp.JCI[floor(n.JCI*t0) : n.SET, 4])/ min(qsp.JCI[floor(n.JCI*t0) : n.SET, 4]) # = 5.167637

# Replication of Figure 1 \ref{fig:TI_WTI}
# uncomment last lines from function Teststatistik and issue command ...
Teststatistik(X = SET.logloss,  t0 = 0.2, k = c(0.16, 0.16), p = 0.005)
Teststatistik(X = SET.logloss,  t0 = 0.2, k = c(0.16, 0.16), p = 0.01)
Teststatistik(X = SET.logloss,  t0 = 0.2, k = c(0.16, 0.16), p = 0.05)
Teststatistik(X = SET.logloss,  t0 = 0.2, k = c(0.16, 0.16), p = 0.1)

Teststatistik(X = KLCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.005)
Teststatistik(X = KLCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.01)
Teststatistik(X = KLCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.05)
Teststatistik(X = KLCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.1)

Teststatistik(X = JCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.005)
Teststatistik(X = JCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.01)
Teststatistik(X = JCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.05)
Teststatistik(X = JCI.logloss, t0 = 0.2, k = c(0.16, 0.16), p = 0.1)

Teststatistik(X = WTI.logloss.shift, t0 = 0.2, k = c(0.16, 0.16), p = 0.005)

library(matrixStats)
Teststatistik <- function(X, t0, k, p){  # Routine, die zu gegebenem Datensatz X, die Teststatistiken ausrechnet
  n <- length(X)                                      # t0\in(0,1/2) ist die strukturbruchlose Zeit, AntEx ein Vektor mit den extremen Stichprobenanteilen
  AnzAntEx <- length(k)                               # a_n = n/fracN bezeichnet die Anzahl an Subsamples für die Schätzung der asympt. Var.
  k <- floor(k * n)                                   # (1-p)-Quantil soll auf Strukturbruch hin überprüft werden.
  
  # Speichern der EQ Schätzwerte für T(0,t)
  T_0t  <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for tail index
  Tq_0t <- array(NA, dim=c(n, AnzAntEx))      # vector of subsample estimates for extreme quantiles
  
  Y    <- sort.int(X[1 : floor(n * t0)], decreasing = TRUE, method = "quick")
  
  for(i in (floor(n * t0) + 1) : n){
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")      # sort X_1,...,X_[floor(n*t0)-1+i] in descending order
    
    for(j in (1 : AnzAntEx)){      
      # Hill-Schätzer
      T_0t[i, j]  <- 1/(k[j]*(i/n)) * sum(log(Y[1 : (floor(k[j]*(i/n))+1)] / Y[floor(k[j]*(i/n))+1])) # Hill estimator
      Tq_0t[i, j] <- Y[floor(k[j]*(i/n))+1] * ((n*p) / k[j])^(-T_0t[i, j])                            # Hill EQ estimator
    }
  } 
  
  T_t1  <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for tail index
  Tq_t1 <- array(NA, dim=c(n, AnzAntEx))      # vector of (backward) subsample estimates for extreme quantiles
  
  Y    <- sort.int(X[(floor(n * (1 - t0)) + 1) : n], decreasing = TRUE, method = "quick")
  
  for(i in floor(n * (1 - t0)) : 1){      # loop over time
    Y   <- sort.int(c(Y, X[i]), decreasing = TRUE, method = "quick")  # Die [ns]+1 bis [nt] ersten Werte der Daten X[ns]+1,...,X[nt] absteigend sortiert
    for(j in (1 : AnzAntEx)){           # loop over k's
      # Hill estimator
      T_t1[i, j] <- 1/(k[j]*(1 - (i-1)/n)) * sum(log(Y[1 : (floor(k[j]*(1 - (i-1)/n))+1)] / Y[floor(k[j]*(1 - (i-1)/n)) + 1]))
      Tq_t1[i, j] <- Y[floor(k[j]*(1 - (i-1)/n))+1] * ((n*p) / k[j])^(-T_t1[i, j])
    }
  }
  
  # Calculate SN test statistics
  maxQ_T  <- array(NA, dim = c(AnzAntEx))       # Save realised suprema for different no. of k's
  maxQ_Tq <- array(NA, dim = c(AnzAntEx)) 
  Q_T     <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer
  Q_Tq    <- array(NA, dim = c(n, AnzAntEx))    # Initialisieren der Vektoren in den die Werte der Teststatistik gespeichert werden sollen: "5" Tests mit "AnzAntEx" versch. k's für die "5" Schätzer   
  SN      <- array(NA, dim = c(n, AnzAntEx))    # Calculating the self-normalizer "AnzAntEx" versch. k's für die "5" Schätzer
  SN_T    <- array(NA, dim = c(n, AnzAntEx))    # Calculating the self-normalizer "AnzAntEx" versch. k's für die "5" Schätzer
  l_sq    <- log(k / (n*p))^2
  
  
  for(i in (floor(n*t0)+1) : floor(n*(1-t0))){      # loop over time
    t <- i/n				                                # discretize time "t"
    Q_Tq[i, ] <- k / l_sq * (t * (1-t) * log(Tq_0t[i, ] / Tq_t1[i, ]))^2 
    Q_T[i, ]  <-       k  * (t * (1-t) * (T_0t[i, ] - T_t1[i, ]))^2
    
    for(m in 1:AnzAntEx){
      SN[i, m] <- 1/n * k[m] / l_sq[m] * sum((((floor(n*t0)+1) : i)/n * log(Tq_0t[(floor(n*t0)+1) : i, m] / Tq_0t[i, m]))^2)
      SN[i, m] <- SN[i, m] + 1/n * k[m] / l_sq[m] * sum((((n-(i : floor(n*(1-t0))) + 1) / n) * log(Tq_t1[i : floor(n*(1-t0)), m] / Tq_t1[i, m]))^2)
      
      SN_T[i, m] <- 1/n * k[m] * sum((((floor(n*t0)+1) : i)/n * (T_0t[(floor(n*t0)+1) : i, m] - T_0t[i, m]))^2)
      SN_T[i, m] <- SN_T[i, m] + 1/n * k[m] * sum((((n-(i : floor(n*(1-t0))) + 1) / n) * (T_t1[i : floor(n*(1-t0)), m] - T_t1[i, m]))^2)
    }
  }
  
  maxQ_T  <- colMaxs( Q_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ] / SN_T[(floor(n * t0) + 1) : floor(n * (1-t0)), ])   # elementwise division
  maxQ_Tq <- colMaxs(Q_Tq[(floor(n * t0) + 1) : floor(n * (1-t0)), ] /   SN[(floor(n * t0) + 1) : floor(n * (1-t0)), ])   # elementwise division
  
  
  #############################################################
  
        # leave this uncommented unless you want to reproduce Figure 1
        # Figure 1 \ref{fig:TI_WTI}
        # calculate rolling window estimates of tail index and std. dev.
        T_t0  <- array(0, dim=c(n, AnzAntEx))    # Hill estimates
        T_01  <- array(0, dim=c(AnzAntEx))       # Hill estimates based on whole sample
        Sd <- rep.int(0, n)                      # standard deviation estimates
        
        Y    <- sort.int(X[1 : n], decreasing = TRUE, method = "quick")
        for(j in (1 : AnzAntEx)){
          # Hill-Schätzer
          T_01[j]    <- 1/k[j] * sum(log(Y[1 : (floor(k[j])+1)] / Y[floor(k[j])+1]))
        }
        
        for(i in floor(n*t0) : n){
          Sd[i] <- sd(X[(i - floor(n*t0) + 1) : i])
        }
        
        for(i in floor(n*t0) : n){    # loop through observation period
          Sd[i] <- sd(X[(i - floor(n*t0) + 1) : i])
          Y <- sort.int(X[(i - floor(n*t0) + 1) : i], decreasing = TRUE, method = "quick")    # Die Werte der Daten X_1,...,X_[floor(n*t0)-1+i] werden absteigend sortiert
          
          for(j in (1 : AnzAntEx)){
            # rolling Hill estimates
            T_t0[i, j]    <- 1/(k[j]*t0) * sum(log(Y[1 : (floor(k[j]*t0)+1)] / Y[floor(k[j]*t0)+1]))
          }
        }
      
      
        for(i in 1 : AnzAntEx){    # i bestimmt den Wert für k/n, hier k/n = 0.08 wegen i=2
            split.screen(figs=c(3, 1))
            par(oma=c(3,1,0,3))
            screen(1)   # time series plot
            par(mar = c(2,3,2,3))
            plot(Datum, -WTI.logloss, type = "l", xlab = "", ylab = "Log-return", xaxs="i", yaxs="i", ylim=c(-0.02+min(-WTI.logloss), 0.02+max(-WTI.logloss)))    # Grafische Veranschaulichung für WTI
            mtext("Log-return", side=2, line=2)
            #plot(Datum, X, type = "l", xlab = "", ylab = "Log-return", xaxs="i", yaxs="i", ylim=c(-0.02+min(X), 0.02+max(X)))    # Grafische Veranschaulichung für WTI
            #abline(h = c(-5, 0, 5, 10))
            
            screen(2)   # tests statistic plot
            par(mar = c(2,3,2,3))
            Q_H <- rep.int(0, n)
            Q_MR <- rep.int(0, n)
            Q_CV <- rep.int(0, n)
            interval <- (floor(n * t0) + 1) : floor(n * (1-t0))
            Q_H[interval]  <- Q_Tq[interval, i] / SN[interval, i]
            
            plot(Datum, Q_H, type = "l", xlab = "", ylab = "Test statistic", xaxs="i", yaxs="i", ylim = c(0, 1+maxQ_Tq[i]))
            mtext("Test statistic", side=2, line=2)
            
            screen(3)
            par(mar = c(2,3,2,3))
            plot(Datum, T_t0[ , i], axes=FALSE, xlab="", ylab="", type="n")
            lines(Datum[floor(n*t0) : n], T_t0[floor(n*t0) : n, i])
            abline(h = c(0, T_01[i]), lwd = 0.5)
            axis(2, ylim=c(0, max(T_t0[ , i])), las=1)  ## las=1 makes horizontal labels
            mtext("Hill estimates", side=2, line=3)
            
            par(new=TRUE)       # Allow a second plot on the same graph
            plot(Datum, Sd, pch=15,  xlab="", ylab="", ylim=c(0,max(Sd)), axes=FALSE, xaxs="i", type="n", col="red")
            lines(Datum[floor(n*t0) : n], Sd[floor(n*t0) : n], lty="dotted")
            
            print(cbind(Sd, T_t0[, i]))
            
            ## a little farther out (line=3) to make room for labels
            mtext("St.dev. estimates", side=4, line=3)
            axis(4, ylim=c(0, max(Sd)), las=1)
            
            ## Draw the time axis
            axis.Date(1, at = seq(min(Datum), max(Datum), "years"), pos=0)
            mtext("Date", side=1, col="black", line=2.5)
            abline(h = c(0, 0), lwd = 0.5)
            
            close.screen(all = TRUE)
      }
  
  #############################################################
  
  return(list(TI = maxQ_T, EQ = maxQ_Tq))
}


################################################################################
