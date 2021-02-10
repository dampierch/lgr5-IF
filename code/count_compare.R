##
## R Script for Comparing Counts of LGR5+ Cells in FAP vs HLT and LYN Crypts
##


library(readr)
library(dplyr)
library(geepack)
library(ggplot2)
library(cowplot)
library(scales) ## for muted colors


source("util.R")
source("themes.R")
source("interobs.R")
source("gee.R")
source("figs.R")


main <- function() {

    ## Main Analysis: Comparison of LGR5+ Cells in Healthy, Lynch, FAP Crypts

    ## Load Data
    l <- make_data()

    ## Check Inter-Observer Variability
    pl <- interobserver_plots(l)
    interobs_write(pl, "LGR5")
    interobs_write(pl, "Ectopic")

    ## Test Predictors of LGR5+ Cell Number
    ## Predictor: Side, Group: Healthy
    fits <- list(side=fit_gee_side(l$lgr5, "LGR5"))
    ## Predictor: Age, Group: Healthy
    fits$age <- fit_gee_age(l$lgr5, "LGR5")
    ## Predictor: Diagnosis, Group: All
    fits$diagttl <- fit_gee_diagnosis(l$lgr5, "LGR5")
    fits$diagect <- fit_gee_diagnosis(l$ectopic, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline
    fits$diagagettl <- fit_gee_diag_age(l$lgr5, "LGR5")
    fits$diagageect <- fit_gee_diag_age(l$ectopic, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Lynch As Baseline
    ss <- c("Lynch", "Healthy", "FAP")
    fits$diagagettllyn <- fit_gee_diag_age(l$lgr5, "LGR5", subset=ss)
    ## Check Fit and Residuals
    for (f in names(fits)) {
        print(summary(fits[[f]]))
        check_residuals(fits[[f]], f)
    }

    ## Make Figures
    ## Highlight Effect of Age
    # fit <- fit_gee_age(l$lgr5, "LGR5", subset=c("Healthy", "Lynch"))
    ggps <- list(age=make_ggp_age_predict(l$sumttl, fits$age))
    ## Show Main Effect per Subject
    ggps$mainsubj <- make_ggp_dot(l$sumttl)
    ## Show Main Effect per Crypt
    ggps$maincryp <- make_ggp_point(l$lgr5)
    ## Color Subjects
    ggps$colsubjhlt <- make_ggp_dot_colour(l$sumttl, "Healthy")
    ggps$colsubjlyn <- make_ggp_dot_colour(l$sumttl, "Lynch")
    ggps$colsubjfap <- make_ggp_dot_colour(l$sumttl, "FAP")
    ## Color Crypts
    ggps$colcryphlt <- make_ggp_point_colour(l$lgr5, "Healthy")
    ggps$colcryplyn <- make_ggp_point_colour(l$lgr5, "Lynch")
    ggps$colcrypfap <- make_ggp_point_colour(l$lgr5, "FAP")
    ## Highlight Crypts for Individual Subjects
    ggps$facsubj <- make_ggp_point_sub(l$lgr5)
    ## Main Figure 3
    ggps$main3 <- make_main_fig3(ggps)
}


main()
