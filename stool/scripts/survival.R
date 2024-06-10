library(tidyverse)
library(survminer)
require(survival)
library (readr)
survival <- read_csv("/Users/keithhazleton/Library/CloudStorage/Box-Box/Issa/aim1a_survival.csv", 
                     col_types = cols(Day = col_integer(), status = col_integer(), diet = col_character(), 
                     group = col_integer()))


# Make survival plot by diet
fit <- surv_fit(Surv(Day, status) ~ diet, data = survival)
ggsurvplot(fit, data = survival, ggtheme = theme_survminer (base_size = 14), 
           legend = "bottom", font.x = 14, font.y = 14, font.legend = 14,
           legend.title = 'diet:',
           pval = TRUE, pval.method =TRUE, conf.int = TRUE, xlab = 'Days Post-Infection')

# Cox PH for diet
surv_obj <- Surv(time = survival$Day, event = survival$status)
survival$diet_f = factor(survival$diet, levels = c("HF/LF", "LF+fiber","LF/LF", "chow", "HF+fiber"))
fit.cox <- coxph(surv_obj ~ diet_f, data = survival)
summary(fit.cox)

# Make survival plot by fiber content
fit <- surv_fit(Surv(Day, status) ~ fiber, data = survival)
ggsurvplot(fit, data = survival, ggtheme = theme_survminer (base_size = 14), 
           legend = "bottom", font.x = 14, font.y = 14, font.legend = 14,
           legend.title = 'diet:',
           pval = TRUE, pval.method =TRUE, conf.int = TRUE, xlab = 'Days Post-Infection')

# Cox PH for fiber content
surv_obj <- Surv(time = survival$Day, event = survival$status)
survival$fiber_f = factor(survival$fiber, levels = c("high_fiber", "low_fiber"))
fit.cox <- coxph(surv_obj ~ fiber_f, data = survival)
summary(fit.cox)

# Make survival plot by fat content
fit <- surv_fit(Surv(Day, status) ~ fat, data = survival)
ggsurvplot(fit, data = survival, ggtheme = theme_survminer (base_size = 14), 
           legend = "bottom", font.x = 14, font.y = 14, font.legend = 14,
           legend.title = 'diet:',
           pval = TRUE, pval.method =TRUE, conf.int = TRUE, xlab = 'Days Post-Infection')

# Cox PH for fat content
surv_obj <- Surv(time = survival$Day, event = survival$status)
survival$fat_f = factor(survival$fat, levels = c("low_fat", "high_fat"))
fit.cox <- coxph(surv_obj ~ fat_f, data = survival)
summary(fit.cox)

# repeat with chow removed
no_chow <- filter(survival, diet != "chow")

# Make survival plot by fiber content
fit <- surv_fit(Surv(Day, status) ~ fiber, data = no_chow)
ggsurvplot(fit, data = no_chow, ggtheme = theme_survminer (base_size = 14), 
           legend = "bottom", font.x = 14, font.y = 14, font.legend = 14,
           legend.title = 'diet:',
           pval = TRUE, pval.method =TRUE, conf.int = TRUE, xlab = 'Days Post-Infection')

# Cox PH for fiber content
surv_obj <- Surv(time = no_chow$Day, event = no_chow$status)
no_chow$fiber_f = factor(no_chow$fiber, levels = c("low_fiber", "high_fiber"))
fit.cox <- coxph(surv_obj ~ fiber_f, data = no_chow)
summary(fit.cox)

# Make survival plot by fat content
fit <- surv_fit(Surv(Day, status) ~ fat, data = no_chow)
ggsurvplot(fit, data = no_chow, ggtheme = theme_survminer (base_size = 14), 
           legend = "bottom", font.x = 14, font.y = 14, font.legend = 14,
           legend.title = 'diet:',
           pval = TRUE, pval.method =TRUE, conf.int = TRUE, xlab = 'Days Post-Infection')

# Cox PH for fat content
surv_obj <- Surv(time = no_chow$Day, event = no_chow$status)
no_chow$fat_f = factor(no_chow$fat, levels = c("low_fat", "high_fat"))
fit.cox <- coxph(surv_obj ~ fat_f, data = no_chow)
summary(fit.cox)
