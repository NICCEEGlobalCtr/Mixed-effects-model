library(nlme)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(writexl)
library(AICcmodavg)
##########Figure.3b analysis####
data <- read_excel("wheat_D_CK_f.xlsx")

treatments <- c("PK", "N", "NK", "NP", "NPK", "NPKm", "hNPKm", "NPKs")

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_NC - pred)^2))
  r2 <- 1 - sum((data$d_NC - pred)^2) / sum((data$d_NC - mean(data$d_NC))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_NC ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_NC <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_NC <- fit_mixed_model(data_subset)
  pred_mixed_NC <- predict(mixed_model_NC, newdata = data_subset)
  pred_mixed_NC_ci <- predictSE.lme(mixed_model_NC, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_NC_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_NC_ci$fit - 1.96 * pred_mixed_NC_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_NC_ci$fit + 1.96 * pred_mixed_NC_ci$se.fit

  goodness_mixed_NC <- calculate_goodness_of_fit(mixed_model_NC, data_subset, pred_mixed_NC)
  
  model_comparison_NC <- rbind(
    model_comparison_NC,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_NC)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_NC)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],
      AIC = goodness_mixed_NC$AIC,
      BIC = goodness_mixed_NC$BIC,
      RMSE = goodness_mixed_NC$RMSE,
      R2 = goodness_mixed_NC$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_NC
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_NC)["(Intercept)"],  
    Slope = summary(mixed_model_NC)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],  
    R2 = goodness_mixed_NC$R2  
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_NC, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-30, 90),
    breaks = seq(-30, 90, by = 30)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("wheat_mixed_effects_plot_NC.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

##################
##########Figure.S4 analysis####
calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_Yield - pred)^2))
  r2 <- 1 - sum((data$d_Yield - pred)^2) / sum((data$d_Yield - mean(data$d_Yield))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_Yield ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_Yield <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_Yield <- fit_mixed_model(data_subset)
  pred_mixed_Yield <- predict(mixed_model_Yield, newdata = data_subset)
  pred_mixed_Yield_ci <- predictSE.lme(mixed_model_Yield, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_Yield_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_Yield_ci$fit - 1.96 * pred_mixed_Yield_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_Yield_ci$fit + 1.96 * pred_mixed_Yield_ci$se.fit
 
  goodness_mixed_Yield <- calculate_goodness_of_fit(mixed_model_Yield, data_subset, pred_mixed_Yield)
  
  model_comparison_Yield <- rbind(
    model_comparison_Yield,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_Yield)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"],
      AIC = goodness_mixed_Yield$AIC,
      BIC = goodness_mixed_Yield$BIC,
      RMSE = goodness_mixed_Yield$RMSE,
      R2 = goodness_mixed_Yield$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_Yield
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_Yield)["(Intercept)"],  
    Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"], 
    R2 = goodness_mixed_Yield$R2  
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_Yield, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-100, 1300),
    breaks = seq(-100, 1300, by = 300)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("wheat_mixed_effects_plot_Yield.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

###########################
##########Figure.4b analysis####
data <- read_excel("maize_D_CK_f.xlsx")

treatments <- c("PK", "N", "NK", "NP", "NPK", "NPKm", "hNPKm", "NPKs")

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_NC - pred)^2))
  r2 <- 1 - sum((data$d_NC - pred)^2) / sum((data$d_NC - mean(data$d_NC))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_NC ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_NC <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_NC <- fit_mixed_model(data_subset)
  pred_mixed_NC <- predict(mixed_model_NC, newdata = data_subset)
  pred_mixed_NC_ci <- predictSE.lme(mixed_model_NC, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_NC_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_NC_ci$fit - 1.96 * pred_mixed_NC_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_NC_ci$fit + 1.96 * pred_mixed_NC_ci$se.fit
  
  goodness_mixed_NC <- calculate_goodness_of_fit(mixed_model_NC, data_subset, pred_mixed_NC)
  
  model_comparison_NC <- rbind(
    model_comparison_NC,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_NC)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_NC)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],
      AIC = goodness_mixed_NC$AIC,
      BIC = goodness_mixed_NC$BIC,
      RMSE = goodness_mixed_NC$RMSE,
      R2 = goodness_mixed_NC$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_NC
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_NC)["(Intercept)"], 
    Slope = summary(mixed_model_NC)$tTable["Time", "Value"], 
    PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],  
    R2 = goodness_mixed_NC$R2 
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_NC, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-30, 90),
    breaks = seq(-30, 90, by = 30)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("maize_mixed_effects_plot_NC.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

############################
##########Figure.S7 analysis####
calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_Yield - pred)^2))
  r2 <- 1 - sum((data$d_Yield - pred)^2) / sum((data$d_Yield - mean(data$d_Yield))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_Yield ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_Yield <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_Yield <- fit_mixed_model(data_subset)
  pred_mixed_Yield <- predict(mixed_model_Yield, newdata = data_subset)
  pred_mixed_Yield_ci <- predictSE.lme(mixed_model_Yield, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_Yield_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_Yield_ci$fit - 1.96 * pred_mixed_Yield_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_Yield_ci$fit + 1.96 * pred_mixed_Yield_ci$se.fit
  
  goodness_mixed_Yield <- calculate_goodness_of_fit(mixed_model_Yield, data_subset, pred_mixed_Yield)
  
  model_comparison_Yield <- rbind(
    model_comparison_Yield,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_Yield)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"],
      AIC = goodness_mixed_Yield$AIC,
      BIC = goodness_mixed_Yield$BIC,
      RMSE = goodness_mixed_Yield$RMSE,
      R2 = goodness_mixed_Yield$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_Yield
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_Yield)["(Intercept)"], 
    Slope = summary(mixed_model_Yield)$tTable["Time", "Value"], 
    PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"], 
    R2 = goodness_mixed_Yield$R2  
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_Yield, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-100, 4000),
    breaks = seq(-100, 4000, by = 1000)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("maize_mixed_effects_plot_Yield.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

###############################
##########Figure.5b analysis####
data <- read_excel("wheat_nue%.xlsx")

treatments <- c("N",  "NK", "NP", "NPK", "NPKm", "hNPKm", "NPKs") 

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$NUE_D_p - pred)^2))
  r2 <- 1 - sum((data$NUE_D_p - pred)^2) / sum((data$NUE_D_p - mean(data$NUE_D_p))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(NUE_D_p ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_nue <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_nue <- fit_mixed_model(data_subset)
  pred_mixed_nue <- predict(mixed_model_nue, newdata = data_subset)
  pred_mixed_nue_ci <- predictSE.lme(mixed_model_nue, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_nue_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_nue_ci$fit - 1.96 * pred_mixed_nue_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_nue_ci$fit + 1.96 * pred_mixed_nue_ci$se.fit

  goodness_mixed_nue <- calculate_goodness_of_fit(mixed_model_nue, data_subset, pred_mixed_nue)
  
  model_comparison_nue <- rbind(
    model_comparison_nue,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_nue)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_nue)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_nue)$tTable["Time", "p-value"],
      AIC = goodness_mixed_nue$AIC,
      BIC = goodness_mixed_nue$BIC,
      RMSE = goodness_mixed_nue$RMSE,
      R2 = goodness_mixed_nue$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_nue
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_nue)["(Intercept)"], 
    Slope = summary(mixed_model_nue)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_nue)$tTable["Time", "p-value"], 
    R2 = goodness_mixed_nue$R2  # R²
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = NUE_D_p, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-40, 100),
    breaks = seq(-40, 100, by = 30)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("wheat_mixed_effects_plot_nue.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

###############################
##########Figure.6b analysis####
data <- read_excel("maize_nue%.xlsx")

treatments <- c("N",  "NK", "NP", "NPK", "NPKm", "hNPKm", "NPKs")  

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$NUE_D_p - pred)^2))
  r2 <- 1 - sum((data$NUE_D_p - pred)^2) / sum((data$NUE_D_p - mean(data$NUE_D_p))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(NUE_D_p ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_nue <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)

  mixed_model_nue <- fit_mixed_model(data_subset)
  pred_mixed_nue <- predict(mixed_model_nue, newdata = data_subset)
  pred_mixed_nue_ci <- predictSE.lme(mixed_model_nue, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_nue_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_nue_ci$fit - 1.96 * pred_mixed_nue_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_nue_ci$fit + 1.96 * pred_mixed_nue_ci$se.fit
  
  goodness_mixed_nue <- calculate_goodness_of_fit(mixed_model_nue, data_subset, pred_mixed_nue)
  
  model_comparison_nue <- rbind(
    model_comparison_nue,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_nue)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_nue)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_nue)$tTable["Time", "p-value"],
      AIC = goodness_mixed_nue$AIC,
      BIC = goodness_mixed_nue$BIC,
      RMSE = goodness_mixed_nue$RMSE,
      R2 = goodness_mixed_nue$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_nue
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_nue)["(Intercept)"],  
    Slope = summary(mixed_model_nue)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_nue)$tTable["Time", "p-value"], 
    R2 = goodness_mixed_nue$R2 
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = NUE_D_p, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-50, 30),
    breaks = seq(-50, 30, by = 20)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("maize_mixed_effects_plot_nue.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 12, height = 8)

###############################
##########Figure.S5a analysis####
data <- read_excel("wheat_D_NPK_f.xlsx")

treatments <- c("NPKm", "hNPKm", "NPKs")

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_NC - pred)^2))
  r2 <- 1 - sum((data$d_NC - pred)^2) / sum((data$d_NC - mean(data$d_NC))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_NC ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_NC <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_NC <- fit_mixed_model(data_subset)
  pred_mixed_NC <- predict(mixed_model_NC, newdata = data_subset)
  pred_mixed_NC_ci <- predictSE.lme(mixed_model_NC, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_NC_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_NC_ci$fit - 1.96 * pred_mixed_NC_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_NC_ci$fit + 1.96 * pred_mixed_NC_ci$se.fit
 
  goodness_mixed_NC <- calculate_goodness_of_fit(mixed_model_NC, data_subset, pred_mixed_NC)
  
  model_comparison_NC <- rbind(
    model_comparison_NC,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_NC)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_NC)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],
      AIC = goodness_mixed_NC$AIC,
      BIC = goodness_mixed_NC$BIC,
      RMSE = goodness_mixed_NC$RMSE,
      R2 = goodness_mixed_NC$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_NC
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_NC)["(Intercept)"],  
    Slope = summary(mixed_model_NC)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],  
    R2 = goodness_mixed_NC$R2 
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_NC, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-30, 30),
    breaks = seq(-30, 30, by = 20)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("wheat_mixed_effects_plot_NC_manure.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 9, height = 4)

###############################
##########Figure.S6a analysis####
calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_Yield - pred)^2))
  r2 <- 1 - sum((data$d_Yield - pred)^2) / sum((data$d_Yield - mean(data$d_Yield))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_Yield ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_Yield <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_Yield <- fit_mixed_model(data_subset)
  pred_mixed_Yield <- predict(mixed_model_Yield, newdata = data_subset)
  pred_mixed_Yield_ci <- predictSE.lme(mixed_model_Yield, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_Yield_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_Yield_ci$fit - 1.96 * pred_mixed_Yield_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_Yield_ci$fit + 1.96 * pred_mixed_Yield_ci$se.fit
  
  goodness_mixed_Yield <- calculate_goodness_of_fit(mixed_model_Yield, data_subset, pred_mixed_Yield)
  
  model_comparison_Yield <- rbind(
    model_comparison_Yield,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_Yield)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"],
      AIC = goodness_mixed_Yield$AIC,
      BIC = goodness_mixed_Yield$BIC,
      RMSE = goodness_mixed_Yield$RMSE,
      R2 = goodness_mixed_Yield$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_Yield
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_Yield)["(Intercept)"],  
    Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"],  
    R2 = goodness_mixed_Yield$R2 
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_Yield, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = "dashed"),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-30, 250),
    breaks = seq(-30, 250, by = 50)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("wheat_mixed_effects_plot_Yield_manure.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 9, height = 4)

# ###############################
##########Figure.S5b analysis####
data <- read_excel("maize_D_NPK_f.xlsx")

treatments <- c("NPKm", "hNPKm", "NPKs")

calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_NC - pred)^2))
  r2 <- 1 - sum((data$d_NC - pred)^2) / sum((data$d_NC - mean(data$d_NC))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_NC ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_NC <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_NC <- fit_mixed_model(data_subset)
  pred_mixed_NC <- predict(mixed_model_NC, newdata = data_subset)
  pred_mixed_NC_ci <- predictSE.lme(mixed_model_NC, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_NC_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_NC_ci$fit - 1.96 * pred_mixed_NC_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_NC_ci$fit + 1.96 * pred_mixed_NC_ci$se.fit
  
  goodness_mixed_NC <- calculate_goodness_of_fit(mixed_model_NC, data_subset, pred_mixed_NC)
  
  model_comparison_NC <- rbind(
    model_comparison_NC,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_NC)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_NC)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],
      AIC = goodness_mixed_NC$AIC,
      BIC = goodness_mixed_NC$BIC,
      RMSE = goodness_mixed_NC$RMSE,
      R2 = goodness_mixed_NC$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_NC
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_NC)["(Intercept)"],  
    Slope = summary(mixed_model_NC)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_NC)$tTable["Time", "p-value"],  
    R2 = goodness_mixed_NC$R2 
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_NC, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-40, 40),
    breaks = seq(-40, 40, by = 20)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("maize_mixed_effects_plot_NC_manure.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 9, height = 4)

# ###############################
##########Figure.S6b analysis####
calculate_goodness_of_fit <- function(model, data, pred) {
  aic <- AIC(model)
  bic <- BIC(model)
  rmse <- sqrt(mean((data$d_Yield - pred)^2))
  r2 <- 1 - sum((data$d_Yield - pred)^2) / sum((data$d_Yield - mean(data$d_Yield))^2)
  return(data.frame(AIC = aic, BIC = bic, RMSE = rmse, R2 = r2))
}

fit_mixed_model <- function(data_subset) {
  lme(d_Yield ~ Time, random = ~ 1 | Sites, data = data_subset)
}

model_comparison_Yield <- data.frame()
predicted_data <- data.frame()
slope_p_values_mixed <- data.frame()

for (treatment in treatments) {
  data_subset <- subset(data, Treatment == treatment)
  data_subset <- na.omit(data_subset)
  
  mixed_model_Yield <- fit_mixed_model(data_subset)
  pred_mixed_Yield <- predict(mixed_model_Yield, newdata = data_subset)
  pred_mixed_Yield_ci <- predictSE.lme(mixed_model_Yield, newdata = data_subset, se.fit = TRUE, level = 0)
  data_subset$predicted_mixed <- pred_mixed_Yield_ci$fit
  data_subset$ci_lower_mixed <- pred_mixed_Yield_ci$fit - 1.96 * pred_mixed_Yield_ci$se.fit
  data_subset$ci_upper_mixed <- pred_mixed_Yield_ci$fit + 1.96 * pred_mixed_Yield_ci$se.fit
  
  goodness_mixed_Yield <- calculate_goodness_of_fit(mixed_model_Yield, data_subset, pred_mixed_Yield)
  
  model_comparison_Yield <- rbind(
    model_comparison_Yield,
    data.frame(
      Treatment = treatment,
      Model = "Mixed Effect",
      Time_Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],
      Time_SE = summary(mixed_model_Yield)$tTable["Time", "Std.Error"],
      Time_PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"],
      AIC = goodness_mixed_Yield$AIC,
      BIC = goodness_mixed_Yield$BIC,
      RMSE = goodness_mixed_Yield$RMSE,
      R2 = goodness_mixed_Yield$R2
    )
  )
  
  data_subset$predicted_mixed <- pred_mixed_Yield
  predicted_data <- rbind(predicted_data, data_subset)
  
  slope_p_values_mixed <- rbind(slope_p_values_mixed, data.frame(
    Treatment = treatment,
    Intercept = fixed.effects(mixed_model_Yield)["(Intercept)"],  
    Slope = summary(mixed_model_Yield)$tTable["Time", "Value"],  
    PValue = summary(mixed_model_Yield)$tTable["Time", "p-value"], 
    R2 = goodness_mixed_Yield$R2  
  ))
}

slope_p_values_mixed$LineType <- factor(
  ifelse(slope_p_values_mixed$PValue < 0.05, "solid", "dashed"),
  levels = c("solid", "dashed")
)

predicted_data_mixed <- predicted_data %>%
  left_join(slope_p_values_mixed[, c("Treatment", "LineType")], by = "Treatment")

predicted_data_mixed$Treatment <- factor(predicted_data_mixed$Treatment, levels = treatments)
slope_p_values_mixed$Treatment <- factor(slope_p_values_mixed$Treatment, levels = treatments)

plot_mixed <- ggplot(predicted_data_mixed, aes(x = Time, y = d_Yield, color = Treatment)) +
  geom_point(size = 4, alpha = 0.2,) +
  geom_abline(
    data = subset(slope_p_values_mixed),
    aes(slope = Slope, intercept = Intercept, color = Treatment, linetype = LineType),
    linewidth = 1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  facet_wrap(~ Treatment, nrow = 2, ncol = 4, scales = "fixed") +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 15,family = "Times New Roman"),
    strip.text = element_text(size = 15,family = "Times New Roman"),
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  scale_y_continuous(
    limits = c(-50, 300),
    breaks = seq(-50, 300, by = 50)
  ) +
  geom_text(data = slope_p_values_mixed,
            aes(x = Inf, y = Inf,
                label = sprintf("Slope = %.2f\nP = %.3f\nR² = %.2f", Slope, PValue, R2),
                color = Treatment),
            inherit.aes = FALSE,
            vjust = 1.1, hjust = 1.1, size = 5,
            family = "Times New Roman",
            lineheight = 0.8)

plot_mixed

ggsave("maize_mixed_effects_plot_Yield_manure.jpg", plot_mixed, device = "jpeg", dpi = 600, width = 9, height = 4)
