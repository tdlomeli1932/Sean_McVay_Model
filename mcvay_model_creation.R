# Load libraries
library(splines)
library(glmnet)
options(scipen = 999)

#Read and make final data modifications
mcvay <- read.csv("mcvay_analysis_data.csv")
mcvay$qtr <- factor(mcvay$qtr)
mcvay$down <- factor(mcvay$down)
mcvay$half_seconds_remaining <- NULL
mcvay$game_seconds_remaining <- NULL

#Double check for NAs
sapply(mcvay, function(x) {
  sum(is.na(x))
})


#Split into training and test sets. 10 Test games, 65 Training games
set.seed(9)
test_games <- sample(unique(mcvay$game_id), 10, replace = FALSE)
mcvay_test <- mcvay[mcvay$game_id %in% test_games, ]
mcvay_training <- mcvay[!(mcvay$game_id %in% test_games), ]

#Create knot selection function for natural splines
select_knot_count_for_each <- function(data, response, variables_to_test, fixed_vars = NULL, k_range = 3:6, verbose = TRUE) {
  results <- lapply(variables_to_test, function(var) {
    others <- setdiff(c(fixed_vars, variables_to_test), var)
    knot_results <- lapply(k_range, function(k) {
      df_spline <- k - 1
      spline_term <- paste0("ns(", var, ", df = ", df_spline, ")")
      rhs_terms <- c(spline_term, others)
      formula <- reformulate(rhs_terms, response = response)
      fit <- glm(formula, data = data, family = binomial)
      chi_sq <- fit$null.deviance - fit$deviance
      score <- chi_sq - 2 * k
      data.frame(
        variable = var,
        knots = k,
        df = df_spline,
        chi_sq = chi_sq,
        penalty = 2 * k,
        chi_sq_minus_2k = score,
        AIC = AIC(fit)
      )
    })
    do.call(rbind, knot_results)
  })
  all_results <- do.call(rbind, results)
  best_per_var <- do.call(rbind, lapply(split(all_results, all_results$variable), function(df) {
    df[which.max(df$chi_sq_minus_2k), ]
  }))
  if (verbose) {
    message("Top choice per variable (maximizing chi^2 - 2k):")
    print(best_per_var[order(best_per_var$variable), ])
  }
  return(invisible(list(all_results = all_results, best_k_per_var = best_per_var)))
}

#create a function to compare the spline function to the base one
compare_spline_to_linear <- function(data, formula, predictor, df_spline) {
  full_formula <- as.formula(formula)
  response <- all.vars(full_formula[[2]])
  rhs_terms <- attr(terms(full_formula, data = data), "term.labels")
  spline_terms <- ifelse(rhs_terms == predictor,
                         paste0("ns(", predictor, ", df = ", df_spline, ")"),
                         rhs_terms)
  spline_formula <- reformulate(spline_terms, response = response)
  model_linear <- glm(full_formula, data = data, family = binomial)
  model_spline <- glm(spline_formula, data = data, family = binomial)
  anova(model_linear, model_spline, test = "Chisq")
}

fixed_covariates <- c("isPost", "qtr", "shotgun", "no_huddle", "div_game", "isTurf", "isIndoor", "wind", "down", "goal_to_go")
test_vars <- c("yardline_100", "quarter_seconds_remaining", "fixed_drive", "ydstogo", "posteam_score",
               "score_differential", "temp", "def_adjusted_rush_ypp", "def_adjusted_pass_ypp")
knots <- select_knot_count_for_each(
  data = mcvay_training,
  response = "isPass",
  variables_to_test = test_vars,
  fixed_vars = fixed_covariates,
  k_range = 3:6
)

#Compare spline vs linear to decide which variables to spline and which to keep as linear
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "def_adjusted_pass_ypp", df_spline = 5)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "def_adjusted_rush_ypp", df_spline = 4)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "fixed_drive", df_spline = 2)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "posteam_score", df_spline = 3)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "quarter_seconds_remaining", df_spline = 5)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "score_differential", df_spline = 5)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "temp", df_spline = 3)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "yardline_100", df_spline = 3)
compare_spline_to_linear(data = mcvay_training[, -(1:5)], formula = isPass ~ ., predictor = "ydstogo", df_spline = 4)

#Construct splines
pass_ypp_spline <- ns(mcvay_training$def_adjusted_pass_ypp, df = 5)
colnames(pass_ypp_spline) <- paste0(rep("pass_ypp_spline", 5), 1:5)
posteam_score_spline <- ns(mcvay_training$posteam_score, df = 3)
colnames(posteam_score_spline) <- paste0(rep("posteam_score_spline", 3), 1:3)
score_differential_spline <- ns(mcvay_training$score_differential, df = 5)
colnames(score_differential_spline) <- paste0(rep("score_differential_spline", 5), 1:5)
yardline_100_spline <- ns(mcvay_training$yardline_100, df = 3)
colnames(yardline_100_spline) <- paste0(rep("yardline_100_spline", 3), 1:3)
ydstogo_spline <- ns(mcvay_training$ydstogo, df = 4)
colnames(ydstogo_spline) <- paste0(rep("ydstogo_spline", 4), 1:4)

#Build prediction dataset
mcvay_predict <- mcvay_training[,-c(1:5)]
mcvay_predict$ydstogo <- NULL
mcvay_predict$yardline_100 <- NULL
mcvay_predict$score_differential <- NULL
mcvay_predict$posteam_score <- NULL
mcvay_predict$def_adjusted_pass_ypp <- NULL
mcvay_predict <- cbind(mcvay_predict, pass_ypp_spline, posteam_score_spline, score_differential_spline,
                       yardline_100_spline, ydstogo_spline)

#Build model matrix to input to lasso model
training_matrix <- model.matrix(~ .*qtr*down*goal_to_go + .*isPost +.*shotgun +.*no_huddle +
                                  .*div_game +.*goal_to_go +.*isTurf + .* isIndoor +.* -
                                  qtr:fixed_drive - qtr : down : fixed_drive - qtr:goal_to_go:fixed_drive , data = mcvay_predict[,-1])[,-1]

col_counts <- colSums(training_matrix != 0)
training_matrix <- training_matrix[,col_counts > 20] # keep predictors used in at least ~1/3 of games

#Fit lasso model. This line will take a while to run
lasso_mcvay <- cv.glmnet(training_matrix, mcvay_predict$isPass, alpha = 1, family = "binomial", nfolds = 10)

beta <- coef(lasso_mcvay, s = "lambda.min")
mcvay_coefficients <- as.matrix(beta)
mcvay_coefficients <- mcvay_coefficients[mcvay_coefficients[, 1] != 0, , drop = FALSE]
mcvay_coefficients <- mcvay_coefficients[order(abs(mcvay_coefficients[, 1]), decreasing = TRUE), , drop = FALSE] #Sort by absolute value (descending) to examine more easily

training_predicted_probabilities <- predict(lasso_mcvay, newx = training_matrix, s = "lambda.min", type = "response")

training_predicted_probabilities <- as.numeric(training_predicted_probabilities)
training_predicted_probabilities[training_predicted_probabilities < 0.5] <- 0
training_predicted_probabilities[training_predicted_probabilities > 0.5] <- 1

training_success_rate <- 1 - sum(abs(training_predicted_probabilities - mcvay_predict$isPass)) / nrow(mcvay_predict)
training_baseline <- sum(mcvay_predict$isPass) / nrow(mcvay_predict)

apply_spline_from_train <- function(train_basis, test_x, varname) {
  test_basis <- predict(train_basis, newx = test_x)
  colnames(test_basis) <- paste0(varname, "_s", seq_len(ncol(test_basis)))
  return(test_basis)
}


#Apply training spline transformations to test data
test_spline_pass_ypp <- apply_spline_from_train(pass_ypp_spline, mcvay_test$def_adjusted_pass_ypp, "pass_ypp")
test_spline_posteam_score <- apply_spline_from_train(posteam_score_spline, mcvay_test$posteam_score, "posteam_score")
test_spline_score_differential <- apply_spline_from_train(score_differential_spline, mcvay_test$score_differential, "score_differential")
test_spline_yardline_100 <- apply_spline_from_train(yardline_100_spline, mcvay_test$yardline_100, "yardline_100")
test_spline_ydstogo <- apply_spline_from_train(ydstogo_spline, mcvay_test$ydstogo, "ydstogo")

#Prepare test data for modeling, removing non-splined variables and unused columns
test_predict <- mcvay_test[, -c(1:5)]
test_predict$ydstogo <- NULL
test_predict$yardline_100 <- NULL
test_predict$score_differential <- NULL
test_predict$posteam_score <- NULL
test_predict$def_adjusted_pass_ypp <- NULL

# Add spline basis variables
test_predict <- cbind(
  test_predict,
  test_spline_pass_ypp,
  test_spline_posteam_score,
  test_spline_score_differential,
  test_spline_yardline_100,
  test_spline_ydstogo
)

test_matrix <- model.matrix(
  ~ .*qtr*down*goal_to_go + .*isPost + .*shotgun + .*no_huddle +
    .*div_game + .*goal_to_go + .*isTurf + .*isIndoor + .* -
    qtr:fixed_drive - qtr:down:fixed_drive - qtr:goal_to_go:fixed_drive,
  data = test_predict[, -1]
)[, -1]


test_matrix <- test_matrix[, col_counts > 20] # Keep only features used in training

test_predicted_probabilities <- predict(
  lasso_mcvay,
  newx = test_matrix,
  s = "lambda.min",
  type = "response"
)

test_predicted_probabilities <- as.numeric(test_predicted_probabilities)
test_predicted_probabilities[test_predicted_probabilities < 0.5] <- 0
test_predicted_probabilities[test_predicted_probabilities > 0.5] <- 1

test_success_rate <- 1 - sum(abs(test_predicted_probabilities - test_predict$isPass)) / nrow(test_predict)
test_baseline <- sum(test_predict$isPass) / nrow(test_predict)

test_success_rate


