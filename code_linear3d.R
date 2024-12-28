
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


# Define the Linear 3D system
linear_3d_system <- function(t, state, parameters) {
  with(as.list(state), {
    dx0 <- -1 * x0 + 20.0 * x1
    dx1 <- -20.0 * x0 - 1 * x1
    dx2 <- -3 * x2
    return(list(c(dx0, dx1, dx2)))
  })
}

# Initial conditions
state <- c(x0 = 2, x1 = 0, x2 = 1)
params <- c(a = -1, b = 20.0, c = -20.0, d = -1, e = -3)

# Solve the system with parameters
linear_3d_data <- as.data.frame(ode(
  y = state,
  times = times,
  func = linear_3d_system,
  parms = params
))
names(linear_3d_data) <- c("time", "x0", "x1", "x2")

# Plot the results
ggplot(linear_3d_data, aes(x = time)) +
  geom_line(aes(y = x0, color = "x0"), size = 1) +
  geom_line(aes(y = x1, color = "x1"), size = 1, alpha = 0.7) +
  geom_line(aes(y = x2, color = "x2"), size = 1, alpha = 0.7) +
  labs(title = "Linear 3D ODE - Time Series", x = "Time", y = "Values") +
  scale_color_manual(values = c("red", "blue", "green"), name = "Variables") +
  theme_minimal()



# 3D Trajectory Plot
plot_ly(
  data = linear_3d_data,
  x = ~x0,
  y = ~x1,
  z = ~x2,
  type = "scatter3d",
  mode = "lines",
  line = list(width = 4, color = "red")
) %>%
  layout(
    title = "Linear 3D ODE - Trajectory",
    scene = list(
      xaxis = list(title = "x0"),
      yaxis = list(title = "x1"),
      zaxis = list(title = "x2")
    )
  )

# Define ground truth for the Linear 3D ODE system
ground_truth_linear_3d <- list()
ground_truth_linear_3d[[1]] <- c("x.1", "x.2")
ground_truth_linear_3d[[2]] <- c("x.1", "x.2")
ground_truth_linear_3d[[3]] <- c("x.3")


noise_levels <- 0.1*2^c(0:7)  
# run FBMS and store results
results_linear3d <- data.frame()

y_linear_3d <- finite_difference(as.matrix(linear_3d_data[, c("x0", "x1", "x2")]), dt)
X_linear_3d <- cbind(1, linear_3d_data$x0, linear_3d_data$x1, linear_3d_data$x2)


for (noise_level in noise_levels) {
  for (i in 1:repetitions) {
    # Add noise to linear_3d data (assuming y_linear_3d is a matrix with 3 columns)
    y_noisy <- y_linear_3d + matrix(rnorm(n = length(y_linear_3d), mean = 0, sd = noise_level), nrow = nrow(y_linear_3d), ncol = ncol(y_linear_3d))
    
    set.seed(i)
    
    # Initialize variables for TPR and FPR
    TPR_total <- 0
    FPR_total <- 0
    observed <- sample.int(449000,20000,replace = F) + 500
    predict.insample <- observed[1:1000]
    observed <- observed[1001:2000]
    predict.outsample <- c(1:500,449501:500000) 
    # Loop over each coordinate (x, y, z) in y_linear_3d (which is 3D)
    
    r2.train <- 0
    r2.in.test  <- 0
    r2.out.test  <- 0
    r2.train.noisy  <- 0
    r2.in.test.noisy  <- 0
    r2.out.test.noisy <- 0
    
    
    for (coord in 1:3) {
      # Fit glmnet for each coordinate of y_linear_3d
      data <- data.frame(yc=as.numeric(y_noisy[, coord]),x =  X_linear_3d[-1, 2:4, drop = FALSE])
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
      
      lentp <- length(ground_truth_linear_3d[[coord]])
      
      # Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
      TPR <- sum(res$feats.strings %in% ground_truth_linear_3d[[coord]]) /  lentp
      FPR <- sum(! (res$feats.strings %in% ground_truth_linear_3d[[coord]])) / min(1,length(res$feats.strings))
      
      # Accumulate TPR and FPR for each coordinate
      TPR_total <- TPR_total + TPR
      FPR_total <- FPR_total + FPR
      
      if(length(res$feats.strings)==0)
        res[1,] = c("1",1)
      
      preds.in.test <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data.in.test)
      preds.out.test <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data.out.test)
      preds.train <- predict(lm(as.formula(paste0("yc ~",paste0(res$feats.strings,collapse = "+"),collapse = "")),data = data),newdata = data)
      
      
      r2.train <- r2.train + cor(preds.train,y_linear_3d[observed,coord])^2
      r2.in.test <- r2.in.test + cor(preds.in.test,y_linear_3d[predict.insample,coord])^2
      r2.out.test <-r2.out.test + cor(preds.out.test,y_linear_3d[predict.outsample,coord])^2
      
      r2.train.noisy <- r2.train.noisy + cor(preds.train,y_noisy[observed,coord])^2
      r2.in.test.noisy <- r2.in.test.noisy + cor(preds.in.test,y_noisy[predict.insample,coord])^2
      r2.out.test.noisy <- r2.out.test.noisy + cor(preds.out.test,y_noisy[predict.outsample,coord])^2
    }
    
    # Average TPR and FPR across the 3 coordinates
    TPR_avg <- TPR_total / 3
    FPR_avg <- FPR_total / 3
    
    # Store results
    results_linear3d  <- rbind(results_linear3d, data.frame(Noise = noise_level,  TPR = TPR_avg, FPR = FPR_avg,r2.train = r2.train/3,r2.in.test = r2.in.test/3,r2.out.test = r2.out.test/3,r2.train.noisy = r2.train.noisy/3,r2.in.test.noisy = r2.in.test.noisy/3,r2.out.test.noisy = r2.out.test.noisy/3))
  }
}

write.csv(results_linear3d,"linear_3d.csv")
save(results_linear3d,file = "results_linear_3d")

# Calculate the mean and 95% CI for TPR and FPR at each noise level
summary_results_linear_3d <- results_linear3d %>%
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
ggplot(summary_results_linear_3d, aes(x = Noise)) +
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
summary_results_linear_3d_r2 <- results_linear3d %>%
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
ggplot(summary_results_linear_3d_r2, aes(x = Noise)) +
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

