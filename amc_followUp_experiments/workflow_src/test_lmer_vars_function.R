## running multiple linear mixed effects models and testing effects of fixed variables
## code courtesy of sofia burkhart
## 7-17-2025

## reading in needed libraries
library(dplyr)
library(broom)
library(magrittr)
library(lme4)


#### FUNCTION!! ####

# function for testing variables in models with random effects 
mixed_modeling <- function(predictors, 
                           response_vars, 
                           random_effect, 
                           data){ 
  for (response in response_vars ){ 
    # create full model 
    f <- paste0(predictors, collapse = "+")
    ff <- as.formula(paste0(response,"~", f, "+", random_effect)) 
    fm <- lmer(ff,
               data = data,
               REML = F)
    
    for(predictor in predictors){ 
      # create null models 
      pred <- predictors [! predictors %in% c(predictor)]
      pred <- paste0(pred, collapse = "+")
      nf <- as.formula((paste0(response, "~", pred, "+", random_effect)))
      nm <- lmer(nf, 
                 data = data, 
                 REML = F)
      
      # test for significance 
      res <- anova(nm, fm)
      
      ## need to pull out r/r-squared from this (or calculate it via sum of squares/total sum of squares = r-squared)
      # get beta 
      betas <- fm@beta
      cols <- c("Intercept", predictors)
      position <- which(cols == predictor)
      spec_beta <- as.character(fm@beta[position])
  
      ##print((paste0(response, " ~ ", predictor, "   p-val:", res$`Pr(>Chisq)`[2],"   beta:", spec_beta)))
      
      # if p-value is significant (alpha = 0.05)
      if (res$`Pr(>Chisq)`[2] < 0.05){ 
        # get beta 
        betas <- fm@beta
        cols <- c("Intercept", predictors)
        position <- which(cols == predictor)
        spec_beta <- as.character(fm@beta[position])
        
        print((paste0("significant results:", response, " ~ ", predictor, "   p-val:", res$`Pr(>Chisq)`[2],"   beta:", spec_beta)))
      }
      
    }
    
  }
}
