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
source("figs.R")
source("gee.R")


compare_counts <- function(dat, sum, version, pval=NULL) {
    fit <- gee_analysis(dat, sum, version)
    fit <- gee_analysis(dat, sum, version, pval=pval)
    ggp <- make_ggp_point(dat, version)
    ggp <- make_ggp_point2(dat, version)
    s <- unlist(
        unique(
            subset(dat, Phenotype == "Healthy")[ , "Subject_ID1"]
        ),
        use.names=FALSE
    )
    for (each in s) {
        ggp <- make_ggp_point_sub(dat, each, version)
    }
}


main <- function() {

    ## Main Analysis: Comparison of LGR5+ Cells in Healthy, Lynch, FAP Crypts

    ## Load Data
    l <- make_data()

    ## Check Inter-Observer Variability
    pl <- interobserver_plots(l)
    interobs_write(pl, "LGR5")
    interobs_write(pl, "Ectopic")

    ## Test Predictors of LGR5+ Cell Number
    fits <- list()
    ## Predictor: Side, Group: Healthy
    fits$side <- fit_gee_side(l$lgr5, "LGR5")
    ## Predictor: Age, Group: Healthy
    fits$age <- fit_gee_age(l$lgr5, "LGR5")
    ## Predictor: Diagnosis, Group: All
    fits$diagttl <- fit_gee_diagnosis(l$lgr5, "LGR5")
    fits$diagect <- fit_gee_diagnosis(l$ecto, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Healthy As Baseline
    fits$diagagettl <- fit_gee_diag_age(l$lgr5, "LGR5")
    fits$diagageect <- fit_gee_diag_age(l$ecto, "Ectopic")
    ## Predictor: Diagnosis and Age, Group: All, Lynch As Baseline
    ss <- c("Lynch", "Healthy", "FAP")
    fits$diagagettllyn <- fit_gee_diag_age(l$lgr5, "LGR5", subset=ss)
    ## Check Residuals
    for (f in names(fits)) {
        check_residuals(fits[[f]], f)
    }

    ## Make Figures
    ggps <- list()
    ## Highlight Effect of Age
    # fit <- fit_gee_age(l$lgr5, "LGR5", subset=c("Healthy", "Lynch"))
    ggps$age <- make_ggp_age_predict(l$sumttl, fits$age)
    ## Show Main Effect per Subject
    ggps$mainsubj <- make_ggp_dot(l$sumttl)
    ## Show Main Effect per Crypt
    ggps$maincryp <- make_ggp_point(l$lgr5)
    ## Color Subjects
    ## Color Crypts
    ## Highlight Healthy Subjects

}


#main()
