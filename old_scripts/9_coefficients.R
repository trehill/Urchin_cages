
# Load required libraries
library(glmmTMB)
library(broom)
library(broom.mixed)
library(ggplot2)
library(dplyr)

data <- read.csv("/Users/tessarehill/Dropbox/MSc/Exclosures/processed_data/cleaned_data.csv")

#kelp survival as proportion
data$kelp_prop <- data$kelp_survival/100

#transform response variable with small constant (beta regression can not handle zeros)
data$kelp_prop_adj <- (data$kelp_prop * (nrow(data) - 1) + 0.5) / nrow(data)
#standardizing predictors
data$mean_temperature_c <- scale(data$mean_temperature_c)
data$bottom_ph <- scale(data$bottom_ph)
data$bottom_salinity <- scale(data$bottom_salinity)
data$bottom_do <- scale(data$bottom_do)
data$mean_light_lux <- scale(data$mean_light_lux)
data$mean_depth_m <- scale(data$mean_depth_m)
data$urchin_density <- scale(data$urchin_density)
data$site_density <- scale(data$site_density)
data$simp_div <- scale(data$simp_div)
data$shannon_div <- scale(data$shannon_div)
data$species_richness <- scale(data$species_richness)
data$percent_occuring_algae <- scale(data$percent_occuring_algae)

#month needs to be turned back into a factor
data$site <- factor(data$site)
data$month <- factor(data$month)
data$treatment <- as.factor(data$treatment)
data$gravel_size <- as.factor(data$gravel_size)
data$plot_id <- as.factor(data$plot_id)


kelp_model <- glmmTMB(
  kelp_prop_adj ~ gravel_size + mean_temperature_c + simp_div + species_richness + percent_occuring_algae + (1 | site) + (1 | plot_id) + (1 | month), 
  data = data,
  family = beta_family(link = "logit")
)

urchin_plots  <- glmmTMB(urchin_density ~ treatment + mean_temperature_c + bottom_ph + (1 | plot_id) +(1 | site)  +(1 | month) , 
                                       family = tweedie(link = "log"), 
                                       data = data)

urchin_sites <- glm(site_density ~ mean_light_lux + 
                                   bottom_salinity + bottom_ph, 
                                 data = data, family = gaussian())

growth_model <- glmmTMB(kelp_growth ~ site_density + gravel_size + mean_temperature_c +  
                          bottom_ph + bottom_do + mean_depth_m + simp_div + (1 | site) + (1 | plot_id) + (1 | month), 
                      family = tweedie(), 
                      data = data)


# Extract coefficients and confidence intervals for each model
kelp_model_tidy <- tidy(kelp_model, conf.int = TRUE) %>%
  mutate(model = "Kelp Survival GLMM")

urchin_plots_tidy <- tidy(urchin_plots, conf.int = TRUE) %>%
  mutate(model = "Urchin Plot-Level Density GLMM")

growth_model_tidy <- tidy(growth_model, conf.int = TRUE) %>%
  mutate(model = "Kelp Growth GLMM")

# Correctly tidy the GLM (urchin_sites)
urchin_sites_tidy <- broom::tidy(urchin_sites, conf.int = TRUE) %>%
  mutate(model = "Urchin Site-Level Density GLM")

# Combine all models into one data frame
all_models <- bind_rows(kelp_model_tidy, urchin_plots_tidy, urchin_sites_tidy, growth_model_tidy)

# Filter out random effects (if you only want fixed effects)
all_models <- all_models %>%
  filter(effect == "fixed" | is.na(effect))  # GLM doesn't have an 'effect' column, so allow NA

# Define the color scheme for each model
model_colors <- c(
  "Kelp Survival GLMM" = "#eacf62",
  "Urchin Plot-Level Density GLMM" = "#62a4ea",
  "Urchin Site-Level Density GLM" = "#076fa6",
  "Kelp Growth GLMM" = "#FF5733"
)

# Create the coefficient plot with custom colors
plot <- ggplot(all_models, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = model)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  coord_flip() +
  scale_color_manual(values = model_colors) +
  labs(
    title = "Coefficient and Confidence Interval Plot",
    x = "Predictors",
    y = "Estimate",
    color = "Model"
  ) +
  theme_minimal()

plot

ggsave("/Users/tessarehill/Dropbox/MSc/Exclosures/plots/combined_dot_and_whisker_plot.png", plot = plot, width = 8, height = 6, dpi = 300)



