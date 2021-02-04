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
    ## Predictor: Diagnosis and Age, Group: All
    fits$diagagettl <- fit_gee_diag_age(l$lgr5, "LGR5")
    fits$diagageect <- fit_gee_diag_age(l$ecto, "Ectopic")

    ## Make Figures
    ggps <- list()
    ## Highlight Effect of Age
    # fit <- fit_gee_age(l$lgr5, "LGR5", subset=c("Healthy", "Lynch"))
    ggps$age <- make_ggp_age_predict(l$sumttl, fits$age)
    ## Show Main Effect per Sample
    ggps$X <- make_ggp_dot(l$sumttl)
    ## Show Main Effect per Crypt


    version <- "v4_1"
    dat <- l$data
    sum <- l$summary

    compare_counts(dat, sum, version)
    pval <- "==7.9e-05"
    compare_counts(dat, sum, version, pval=pval)

    version <- "v3_2"
    dat <- subset(dat, Subject_ID1 != "F4")
    sum <- subset(sum, Subject_ID1 != "F4")
    compare_counts(dat, sum, version)
    pval <- "==5.3e-09"
    compare_counts(dat, sum, version, pval=pval)


## compare age in healthy and lynch
fit_model_age(df)

}


sub_side <- function() {

    ## Sub-analysis: Comparison of Proximal vs Distal Healthy Crypts

    l <- make_data()
    dat <- l$data
    dat1 <- subset(dat, Phenotype=="Healthy")
    fit <- geepack::geeglm(
        LGR5_Count ~ factor(Region),
        family=gaussian(),
        data=dat1, id=factor(Subject_ID1),
        zcor=NULL, corstr="exchangeable", std.err="san.se"
    )
    summary(fit)

}


sub_age <- function() {

    ## Sub-analysis: Comparison of Age vs Count in Healthy Crypts

    l <- make_data()
    dat <- l$data
    dat1 <- subset(dat, Phenotype=="Healthy")
    fit <- geepack::geeglm(
        LGR5_Count ~ as.numeric(Age),
        family=gaussian(),
        data=dat1, id=factor(Subject_ID1),
        zcor=NULL, corstr="exchangeable", std.err="san.se"
    )
    summary(fit)

}
