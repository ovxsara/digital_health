# Install if needed
install.packages(c("tidyverse", "haven", "naniar", "VIM", "mice", "survey", "MASS"))

# Load libraries
library(tidyverse)
library(haven)
library(naniar)
library(VIM)
library(mice)
library(survey)
library(MASS)

# Download NHANES files
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt",
              "DEMO_J.xpt", mode = "wb")
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PAQ_J.xpt",
              "PAQ_J.xpt", mode = "wb")
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.xpt",
              "BMX_J.xpt", mode = "wb")

# Read files
demo <- read_xpt("DEMO_J.xpt")
paq  <- read_xpt("PAQ_J.xpt")
bmx  <- read_xpt("BMX_J.xpt")

# Merge datasets
nhanes <- demo %>%
  inner_join(paq, by = "SEQN") %>%
  inner_join(bmx, by = "SEQN")

# Select variables
# PAQ665 = does respondent do moderate recreational activity? (1 yes, 2 no)
# PAQ670 = number of days moderate recreational activity
nhanes_subset <- nhanes %>%
  dplyr::select(SEQN, SDMVPSU, SDMVSTRA, WTMEC2YR,
         RIDAGEYR, RIAGENDR, BMXBMI,
         PAQ665, PAQ670) %>%
  rename(
    age = RIDAGEYR,
    gender = RIAGENDR,
    bmi = BMXBMI,
    modrec_yes = PAQ665,
    modrec_days = PAQ670
  )

# -------------------------------
# Explore missingness
# -------------------------------

# View missingness before cleaning
vis_miss(nhanes_subset)
gg_miss_var(nhanes_subset, show_pct = TRUE)

# Count missing values
nhanes_subset %>%
  summarise(across(everything(), ~sum(is.na(.))))

# Check skip-pattern missingness:
# If modrec_yes == 2 (No), then modrec_days should be 0
table(nhanes_subset$modrec_yes, useNA = "ifany")
table(nhanes_subset$modrec_days, useNA = "ifany")

# -------------------------------
# Clean outcome variable
# -------------------------------

nhanes_subset <- nhanes_subset %>%
  mutate(
    # Remove refusal/don't know if present
    modrec_yes = case_when(
      modrec_yes %in% c(7, 9) ~ NA_real_,
      TRUE ~ modrec_yes
    ),
    modrec_days = case_when(
      modrec_days %in% c(77, 99) ~ NA_real_,
      TRUE ~ modrec_days
    ),
    # Structural missingness: if respondent says "No" to activity, assign 0 days
    modrec_days = case_when(
      modrec_yes == 2 ~ 0,
      TRUE ~ modrec_days
    )
  )

# Re-check missingness after correcting structural missingness
vis_miss(nhanes_subset)
gg_miss_var(nhanes_subset, show_pct = TRUE)

# -------------------------------
# Impute predictor missingness
# -------------------------------

# I will use MICE for BMI, age, gender if needed.
# Not imputing the outcome if it is still missing; dropping remaining missing outcome rows.
impute_data <- nhanes_subset %>%
  dplyr::select(SEQN, SDMVPSU, SDMVSTRA, WTMEC2YR, age, gender, bmi, modrec_days)

md.pattern(impute_data)

set.seed(123)

imp <- mice(
  impute_data,
  m = 5,
  method = c("", "", "", "", "", "", "pmm", ""),
  seed = 123
)

# Inspect imputation
plot(imp)
densityplot(imp)

# Use first completed dataset
nhanes_imputed <- complete(imp, 1)

# Final cleaned dataset
nhanes_clean <- nhanes_imputed %>%
  drop_na(modrec_days, age, gender, bmi) %>%
  mutate(
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    modrec_days = as.numeric(modrec_days)
  )

# -------------------------------
# Poisson regression
# -------------------------------

poisson_model <- glm(
  modrec_days ~ age + gender + bmi,
  family = poisson(link = "log"),
  data = nhanes_clean
)

summary(poisson_model)
exp(coef(poisson_model))  # IRRs

# -------------------------------
# Check overdispersion
# -------------------------------

dispersion <- sum(residuals(poisson_model, type = "pearson")^2) /
  poisson_model$df.residual

dispersion

# If clearly > 1.5, fit Negative Binomial
if (dispersion > 1.5) {
  nb_model <- glm.nb(modrec_days ~ age + gender + bmi, data = nhanes_clean)
  summary(nb_model)
  exp(coef(nb_model))
}

# -------------------------------
# Predicted vs observed
# -------------------------------

nhanes_clean$predicted <- predict(poisson_model, type = "response")
nhanes_clean$residuals <- residuals(poisson_model, type = "response")

ggplot(nhanes_clean, aes(x = modrec_days, y = predicted)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Observed vs Predicted Moderate Recreational Activity Days",
    x = "Observed Days",
    y = "Predicted Days"
  )

ggplot(nhanes_clean, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Histogram of Residuals",
    x = "Residuals"
  )

ggplot(nhanes_clean, aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Predicted",
    x = "Predicted Days",
    y = "Residuals"
  )

# -------------------------------
# Survey-weighted Poisson model
# -------------------------------

nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  data = nhanes_clean,
  nest = TRUE
)

svy_poisson <- svyglm(
  modrec_days ~ age + gender + bmi,
  design = nhanes_design,
  family = poisson()
)

summary(svy_poisson)
exp(coef(svy_poisson))

# If dispersion is an issue in survey model, use quasi-Poisson
svy_quasi <- svyglm(
  modrec_days ~ age + gender + bmi,
  design = nhanes_design,
  family = quasipoisson()
)

summary(svy_quasi)
exp(coef(svy_quasi))
