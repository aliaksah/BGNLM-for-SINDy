library(ggplot2)
library(dplyr)
library(tidyr)


load("results_lorenz")
load("results_linear_3d")
load("results_rossler")


# Combine data for space-efficient plotting
all_results <- bind_rows(
  results_lorenz %>% mutate(System = "Lorenz"),
  results_linear3d %>% mutate(System = "Linear 3D"),
  results_rossler %>% mutate(System = "Rössler-Lorenz Hybrid")
)

# Summarize results for TPR, FPR
summary_tpr_fpr <- all_results %>%
  group_by(System, Noise) %>%
  summarise(
    TPR_mean = mean(TPR,na.rm = T),
    TPR_lower = min(TPR,na.rm = T),
    TPR_upper = max(TPR,na.rm = T),
    FPR_mean = mean(FPR,na.rm = T),
    FPR_lower = min(FPR,na.rm = T),
    FPR_upper = max(FPR,na.rm = T),
  )

# TPR and FPR Plot
p1 <- ggplot(summary_tpr_fpr[summary_tpr_fpr$Noise<15,], aes(x = Noise)) +
  geom_line(aes(y = TPR_mean, color = "TPR"), size = 0.8) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper, fill = "TPR"), alpha = 0.1) +
  geom_line(aes(y = FPR_mean, color = "FPR"), size = 0.8) +
  geom_ribbon(aes(ymin = FPR_lower, ymax = FPR_upper, fill = "FPR"), alpha = 0.1) +
  facet_wrap(~ System, scales = "free_y") +
  labs(title = "TPR and FPR Across Systems", x = "Noise", y = "Value") +
  scale_color_manual(values = c("TPR" = "blue", "FPR" = "red")) +
  scale_fill_manual(values = c("TPR" = "blue", "FPR" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12, face = "bold")
  )

p1

# Summarize results for R2
summary_r2 <- all_results %>%
  group_by(System, Noise) %>%
  summarise(
    r2.train_mean = mean(r2.train,na.rm = T),
    r2.train_lower = min(r2.train,na.rm = T),
    r2.train_upper = max(r2.train,na.rm = T),
    r2.in.test_mean = mean(r2.in.test,na.rm = T),
    r2.in.test_lower = min(r2.in.test,na.rm = T),
    r2.in.test_upper = max(r2.in.test,na.rm = T),
    r2.out.test_mean = mean(r2.out.test,na.rm = T),
    r2.out.test_lower = min(r2.out.test,na.rm = T),
    r2.out.test_upper = max(r2.out.test,na.rm = T),
  )

# R2 Plot
p2 <- ggplot(summary_r2[summary_tpr_fpr$Noise<15,], aes(x = Noise)) +
  geom_line(aes(y = r2.train_mean, color = "Train"), size = 0.8) +
  geom_ribbon(aes(ymin = r2.train_lower, ymax = r2.train_upper, fill = "Train"), alpha = 0.1) +
  geom_line(aes(y = r2.in.test_mean, color = "In_test"), size = 0.8) +
  geom_ribbon(aes(ymin = r2.in.test_lower, ymax = r2.in.test_upper, fill = "In_test"), alpha = 0.1) +
  geom_line(aes(y = r2.out.test_mean, color = "Out_test"), size = 0.8) +
  geom_ribbon(aes(ymin = r2.out.test_lower, ymax = r2.out.test_upper, fill = "Out_test"), alpha = 0.1) +
  facet_wrap(~ System, scales = "free_y") +
  labs(title = "R2 Across Systems", x = "Noise", y = "Value") +
  scale_color_manual(values = c("Train" = "blue", "In_test" = "green", "Out_test" = "red")) +
  scale_fill_manual(values = c("Train" = "blue", "In_test" = "green", "Out_test" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12, face = "bold")
  )

p2


# Save plots
ggsave("result_tpr_fpr_plot.png", p1, width = 8, height = 4,bg = "white")
ggsave("result_r2_plot.png", p2, width = 8, height = 4,bg = "white")