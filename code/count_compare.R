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

    ## Main Analysis: Comparison of LGR5+ Cells in HLT, LYN, FAP Crypts

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
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline, Observer A
    fits$diagagettlA <- fit_gee_diag_age(l$lgr5, "LGR5", obs="A")
    fits$diagageectA <- fit_gee_diag_age(l$ectopic, "Ectopic", obs="A")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline, Observer B
    fits$diagagettlB <- fit_gee_diag_age(l$lgr5, "LGR5", obs="B")
    fits$diagageectB <- fit_gee_diag_age(l$ectopic, "Ectopic", obs="B")
    ## Predictor: Diagnosis and Age, Group: All, Lynch As Baseline
    ss <- c("Lynch", "Healthy", "FAP")
    fits$diagagettllyn <- fit_gee_diag_age(l$lgr5, "LGR5", subset=ss)
    fits$diagagettllyne <- fit_gee_diag_age(l$ectopic, "Ectopic", subset=ss)
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

}


final <- function() {

    ## Final Analysis For Publication

    ## Load Data
    l <- make_data()

    ## Check Inter-Observer Variability
    pl <- interobserver_plots(l)
    interobs_supp_fig(pl)

    ## Test Predictors of LGR5+ Cell Number
    ## Predictor: Side, Group: Healthy
    fits <- list(side=fit_gee_side(l$final, "LGR5"))
    fits$side_ect <- fit_gee_side(l$final, "Ectopic")
    ## Predictor: Age, Group: Healthy
    fits$age <- fit_gee_age(l$final, "LGR5")
    fits$age_ect <- fit_gee_age(l$final, "Ectopic")
    ## Predictor: Diagnosis, Group: All
    fits$diag <- fit_gee_diagnosis(l$final, "LGR5")
    fits$diag_ect <- fit_gee_diagnosis(l$final, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline
    fits$diagage <- fit_gee_diag_age(l$final, "LGR5")
    fits$diagage_ect <- fit_gee_diag_age(l$final, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline, Observer A
    fits$diagageA <- fit_gee_diag_age(l$final, "LGR5", obs="A")
    fits$diagageA_ect <- fit_gee_diag_age(l$final, "Ectopic", obs="A")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline, Observer B
    fits$diagageB <- fit_gee_diag_age(l$final, "LGR5", obs="B")
    fits$diagageB_ect <- fit_gee_diag_age(l$final, "Ectopic", obs="B")
    ## Predictor: Diagnosis and Age, Group: All, Lynch As Baseline
    ss <- c("Lynch", "Healthy", "FAP")
    fits$diagagelyn <- fit_gee_diag_age(l$final, "LGR5", subset=ss)
    fits$diagagelyn_ect <- fit_gee_diag_age(l$final, "Ectopic", subset=ss)
    ## Check Fit and Residuals
    for (f in names(fits)) {
        cat("\n", paste(
            paste(rep("#", 10), collapse=""), toupper(f), paste(rep("#", 10), collapse="")
        ), "\n")
        print(summary(fits[[f]]))
    }

    ## Make Figures
    ## Highlight Effect of Age
    ggps <- list(age=make_ggp_age_predict(l$sumfin, fits$age))
    ## Show Main Effect per Subject
    ggps$mainsubj <- make_ggp_dot(l$sumfin)
    ## Show Main Effect per Crypt
    ggps$maincryp <- make_ggp_point(l$final)
    ## Color Subjects
    ggps$colsubjhlt <- make_ggp_dot_colour(l$sumfin, "Healthy")
    ggps$colsubjlyn <- make_ggp_dot_colour(l$sumfin, "Lynch")
    ggps$colsubjfap <- make_ggp_dot_colour(l$sumfin, "FAP")
    ## Color Crypts
    ggps$colcryphlt <- make_ggp_point_colour(l$final, "Healthy")
    ggps$colcryplyn <- make_ggp_point_colour(l$final, "Lynch")
    ggps$colcrypfap <- make_ggp_point_colour(l$final, "FAP")
    ## Highlight Crypts for Individual Subjects
    ggps$facsubj <- make_ggp_point_sub(l$final)
    ## Main Figure 3
    ggps$main3 <- make_main_fig3(ggps)

}


args <- commandArgs(trailingOnly=TRUE)
NAME <- args[2]
cat(paste("NAME:", NAME, "\n"))


if (NAME == "main") {
    main()
} else {
    final()
}
