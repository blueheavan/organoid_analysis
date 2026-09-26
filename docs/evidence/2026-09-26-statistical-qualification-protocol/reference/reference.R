# Golden values for SG-5: lmerTest (Satterthwaite) and clubSandwich (CR2 + Satterthwaite).
# Rscript reference.R  -> reference_values.json
suppressMessages({library(lmerTest); library(clubSandwich); library(jsonlite)})
d_all <- read.csv("reference_datasets.csv", stringsAsFactors = FALSE)
out <- list(versions = list(R = R.version.string,
                            lme4 = as.character(packageVersion("lme4")),
                            lmerTest = as.character(packageVersion("lmerTest")),
                            clubSandwich = as.character(packageVersion("clubSandwich")),
                            sandwich = as.character(packageVersion("sandwich"))),
            datasets = list())
for (name in unique(d_all$dataset)) {
  d <- d_all[d_all$dataset == name, ]
  d$condition <- factor(d$condition)
  lev <- levels(d$condition); k <- length(lev)
  m <- lmer(y ~ 0 + condition + (1 | replicate), data = d, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5),
                                  check.conv.singular = "ignore"))
  vc <- as.data.frame(VarCorr(m))
  ols <- lm(y ~ 0 + condition, data = d)
  pairs <- list()
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    L <- matrix(0, 1, k); L[1, j] <- 1; L[1, i] <- -1
    ct <- contest(m, L, ddf = "Satterthwaite", joint = FALSE)
    cs <- linear_contrast(ols, vcov = "CR2", cluster = d$replicate, contrasts = L, test = "Satterthwaite")
    cs_vcov <- vcovCR(ols, cluster = d$replicate, type = "CR2")
    pairs[[length(pairs) + 1]] <- list(i = i - 1, j = j - 1,
      lmer = list(estimate = ct$Estimate, se = ct$`Std. Error`, df = ct$df, t = ct$`t value`, p = ct$`Pr(>|t|)`),
      cr2 = list(estimate = cs$Est, se = cs$SE, df = cs$df, ci_low = cs$CI_L, ci_high = cs$CI_U,
                 variance = as.numeric(L %*% as.matrix(cs_vcov) %*% t(L))))
  }
  out$datasets[[name]] <- list(
    singular = isSingular(m),
    sigma_e2 = vc$vcov[vc$grp == "Residual"], sigma_b2 = vc$vcov[vc$grp == "replicate"],
    reml_deviance = as.numeric(REMLcrit(m)), fixef = unname(fixef(m)), ols = unname(coef(ols)),
    pairs = pairs)
}
writeLines(toJSON(out, digits = NA, auto_unbox = TRUE, pretty = TRUE), "reference_values.json")
