
# Load necessary libraries
library(plotly)
library(deSolve)
library(FBMS)
dt <- 0.0001
Time <- 50

repetitions <- 10

times <- seq(0, Time, by = dt)

finite_difference <- function(data, dt) {
  diff <- apply(data, 2, diff)
  diff / dt
}



# Define the Lorenz system
lorenz_system <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- sigma * (y - x)
    dy <- x * (rho - z) - y
    dz <- x * y - beta * z
    return(list(c(dx, dy, dz)))
  })
}

# Initial conditions and parameters for the Lorenz system
#state <- c(x = 1, y = 1, z = 1)  # Initial values for x, y, z
parameters <- c(sigma = 10, rho = 28, beta = 8/3)  # Lorenz parameters
x0_lorenz <- c(x = -0.5, y = -2, z = 3)
lorenz_data <- as.data.frame(ode(y = x0_lorenz, times = times, func = lorenz_system, parms = parameters))


# Extract the solution (x, y, z)
x <- lorenz_data[, "x"]
y <- lorenz_data[, "y"]
z <- lorenz_data[, "z"]

# Create the 3D plot using plotly
plot_ly(x = x, y = y, z = z, type = "scatter3d", mode = "lines", 
        line = list(color = "purple", width = 2)) %>%
  layout(title = "Lorenz 3D",
         scene = list(
           xaxis = list(title = "x"),
           yaxis = list(title = "y"),
           zaxis = list(title = "z")
         ))

y_lorenz <- finite_difference(as.matrix(lorenz_data[, c("x", "y", "z")]), dt)
X_lorenz <- cbind(1, lorenz_data$x, lorenz_data$y, lorenz_data$z, 
                  lorenz_data$x * lorenz_data$y, lorenz_data$x * lorenz_data$z, 
                  lorenz_data$y * lorenz_data$z)

#library(FBMS)

ground_truth_lorenz <- list()
ground_truth_lorenz[[1]] <- c("x.1","x.2")
ground_truth_lorenz[[2]] <- c("x.1","x.2","(x.1*x.3)","(x.3*x.1)")
ground_truth_lorenz[[3]] <- c("x.3","(x.1*x.2)","(x.2*x.1)")



noise_levels <- 0.1*2^c(0:7)  
# run FBMS and store results
results_lorenz <- data.frame()


for (noise_level in noise_levels) {
  for (i in 1:repetitions) {
    # Add noise to Lorenz data (assuming y_lorenz is a matrix with 3 columns)
    y_noisy <- y_lorenz + matrix(rnorm(n = length(y_lorenz), mean = 0, sd = noise_level), nrow = nrow(y_lorenz), ncol = ncol(y_lorenz))
    
    set.seed(i)
    
    # Initialize variables for TPR and FPR
    TPR_total <- 0
    FPR_total <- 0
    observed <- sample.int(449000,2000,replace = F) + 500
    predict.insample <- observed[1:1000]
    observed <- observed[1001:2000]
    predict.outsample <- c(1:500,449501:500000) 
    # Loop over each coordinate (x, y, z) in y_lorenz (which is 3D)
    
    r2.train <- 0
    r2.in.test  <- 0
    r2.out.test  <- 0
    r2.train.noisy  <- 0
    r2.in.test.noisy  <- 0
    r2.out.test.noisy <- 0
    
    
    for (coord in 1:3) {
      # Fit glmnet for each coordinate of y_lorenz
      data <- data.frame(yc=as.numeric(y_noisy[, coord]),x = X_lorenz[-1, 2:4, drop = FALSE])
      data.in.test <- data[predict.insample,]
      data.out.test <- data[predict.outsample,]
      data <- data[observed,]
      params <- gen.params.gmjmcmc(data)
      params$feat$pop.max <- 15
      params$feat$keep.org <- T
      params$feat$check.col <- T
      model <- fbms(formula = yc ~  x.1 + x.2 + x.3, family = "gaussian",data = data,method = "gmjmcmc.parallel",params = params,runs = 10, cores = 10,N.init = 500,
                    N.final = 500, transforms = c("sin_deg","cos_deg","p0","p2","p3","p05","pm1","pm2","pm05"),
                    P = 20)
      res = summary(model,tol = 0.5,labels = names(data)[-1])
      
      lentp <- ifelse(coord == 1, length(ground_truth_lorenz[[coord]]),length(ground_truth_lorenz[[coord]])-1)
      
      # Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
      TPR <- sum(res$feats.strings %in% ground_truth_lorenz[[coord]]) /  lentp
      FPR <- sum(! (res$feats.strings %in% ground_truth_lorenz[[coord]])) / length(res$feats.strings)
      
      # Accumulate TPR and FPR for each coordinate
      TPR_total <- TPR_total + TPR
      FPR_total <- FPR_total + FPR
      
      preds.in.test <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data.in.test)#predict(model,data.in.test,pop = "last")$aggr$mean
      preds.out.test <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data.out.test)#predict(model,data.out.test)$aggr$mean
      preds.train <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data)# predict(model,data)$aggr$mean
      
      r2.train <- r2.train + cor(preds.train,y_lorenz[observed,coord])^2
      r2.in.test <- r2.in.test + cor(preds.in.test,y_lorenz[predict.insample,coord])^2
      r2.out.test <-r2.out.test + cor(preds.out.test,y_lorenz[predict.outsample,coord])^2
      
      r2.train.noisy <- r2.train.noisy + cor(preds.train,y_noisy[observed,coord])^2
      r2.in.test.noisy <- r2.in.test.noisy + cor(preds.in.test,y_noisy[predict.insample,coord])^2
      r2.out.test.noisy <- r2.out.test.noisy + cor(preds.out.test,y_noisy[predict.outsample,coord])^2
    }
    
    # Average TPR and FPR across the 3 coordinates
    TPR_avg <- TPR_total / 3
    FPR_avg <- FPR_total / 3
    
    
    # Store results
    results_lorenz <- rbind(results_lorenz, data.frame(Noise = noise_level,  TPR = TPR_avg, FPR = FPR_avg,r2.train = r2.train/3,r2.in.test = r2.in.test/3,r2.out.test = r2.out.test/3,r2.train.noisy = r2.train.noisy/3,r2.in.test.noisy = r2.in.test.noisy/3,r2.out.test.noisy = r2.out.test.noisy/3))
  }
}

write.csv(results_lorenz,"lorenz.csv")
save(results_lorenz,file = "results_lorenz")

# Calculate the mean and 95% CI for TPR and FPR at each noise level
summary_results_lorenz <- results_lorenz %>%
  group_by(Noise) %>%
  summarise(
    TPR_mean = mean(TPR),
    TPR_lower = min(TPR),
    TPR_upper = max(TPR),
    FPR_mean = mean(FPR),
    FPR_lower = min(FPR),
    FPR_upper = max(FPR)
  )

# Plot the results with 95% confidence intervals
ggplot(summary_results_lorenz, aes(x = Noise)) +
  geom_line(aes(y = TPR_mean, color = "TPR"), size = 1) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper, fill = "TPR"), alpha = 0.1) +
  geom_line(aes(y = FPR_mean, color = "FPR"), size = 1) +
  geom_ribbon(aes(ymin = FPR_lower, ymax = FPR_upper, fill = "FPR"), alpha = 0.1) +
  labs(title = "TPR and FPR with 95% CI vs Noise", 
       x = "Noise", 
       y = "Value") +
  scale_color_manual(values = c("TPR" = "blue", "FPR" = "red")) +
  scale_fill_manual(values = c("TPR" = "blue", "FPR" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())



# Calculate the mean and 95% CI for R2s at each noise level
summary_results_lorenz_r2 <- results_lorenz %>%
  group_by(Noise) %>%
  summarise(
    r2.train_mean = mean(r2.train),
    r2.train_lower = min(r2.train),
    r2.train_upper = max(r2.train),
    r2.in.test_mean = mean(r2.in.test),
    r2.in.test_lower = min(r2.in.test),
    r2.in.test_upper = max(r2.in.test),
    r2.out.test_mean = mean(r2.out.test),
    r2.out.test_lower = min(r2.out.test),
    r2.out.test_upper = max(r2.out.test)
  )

# Plot the results with 95% confidence intervals
ggplot(summary_results_lorenz_r2, aes(x = Noise)) +
  geom_line(aes(y = r2.train_mean, color = "R2_train"), size = 1) +
  geom_ribbon(aes(ymin = r2.train_lower, ymax = r2.train_upper, fill = "R2_train"), alpha = 0.1) +
  geom_line(aes(y = r2.in.test_mean, color = "R2_in_test"), size = 1) +
  geom_ribbon(aes(ymin = r2.in.test_lower, ymax = r2.in.test_upper, fill = "R2_in_test"), alpha = 0.1) +
  geom_line(aes(y = r2.out.test_mean, color = "R2_out_test"), size = 1) +
  geom_ribbon(aes(ymin = r2.out.test_lower, ymax = r2.out.test_upper, fill = "R2_out_test"), alpha = 0.1) +
  labs(title = "R2s with 95% CI vs Noise", 
       x = "Noise", 
       y = "Value") +
  scale_color_manual(values = c("R2_train" = "blue","R2_in_test" = "green", "R2_out_test" = "red")) +
  scale_fill_manual(values = c("R2_train" = "blue", "R2_in_test" = "green","R2_out_test" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())

