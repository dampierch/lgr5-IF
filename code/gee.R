## source from count_compare.R


## generic comment: GEE model requires clustered obs to be in adjacent rows
## id: Data are assumed to be sorted so that observations on each cluster
## appear as contiguous rows in data. If data is not sorted this way, the
## function will not identify the clusters correctly.


fit_gee_side <- function(df, type, subset=c("Healthy")) {
    df <- subset(df, Diagnosis %in% subset)
    df <- df[order(df$Subject_ID), ]
    levels <- c("Right", "Left")
    if (type == "LGR5") {
        fit <- geepack::geeglm(
            LGR5_Count ~ factor(Side, levels=levels),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    } else {
        fit <- geepack::geeglm(
            Ectopic_Count ~ factor(Side, levels=levels),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    }
    return(fit)
}


fit_gee_age <- function(df, type, subset=c("Healthy")) {
    df <- subset(df, Diagnosis %in% subset)
    df <- df[order(df$Subject_ID), ]
    if (type == "LGR5") {
        fit <- geepack::geeglm(
            LGR5_Count ~ as.numeric(Age),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    } else {
        fit <- geepack::geeglm(
            Ectopic_Count ~ as.numeric(Age),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    }
    return(fit)
}


fit_gee_diagnosis <- function(df, type, subset=c("Healthy", "Lynch", "FAP")) {
    df <- subset(df, Diagnosis %in% subset)
    df <- df[order(df$Subject_ID), ]
    levels <- subset
    if (type == "LGR5") {
        fit <- geepack::geeglm(
            LGR5_Count ~ factor(Diagnosis, levels=levels),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    } else {
        fit <- geepack::geeglm(
            Ectopic_Count ~ factor(Diagnosis, levels=levels),
            family=gaussian(),
            data=df, id=factor(Subject_ID),
            zcor=NULL, corstr="exchangeable", std.err="san.se"
        )
    }
    return(fit)
}


fit_gee_diag_age <- function(df, type, obs=NULL, subset=c("Healthy", "Lynch", "FAP")) {
    df <- subset(df, Diagnosis %in% subset)
    df <- df[order(df$Subject_ID), ]
    levels <- subset
    if (type == "LGR5") {
        if (!is.null(obs) && obs == "A") {
            fit <- geepack::geeglm(
                LGR5_Count_A ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        } else if (!is.null(obs) && obs == "B") {
            fit <- geepack::geeglm(
                LGR5_Count_B ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        } else {
            fit <- geepack::geeglm(
                LGR5_Count ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        }
    } else {
        if (!is.null(obs) && obs == "A") {
            fit <- geepack::geeglm(
                Ectopic_Count_A ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        } else if (!is.null(obs) && obs == "B") {
            fit <- geepack::geeglm(
                Ectopic_Count_B ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        } else {
            fit <- geepack::geeglm(
                Ectopic_Count ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
                family=gaussian(),
                data=df, id=factor(Subject_ID),
                zcor=NULL, corstr="exchangeable", std.err="san.se"
            )
        }
    }
    return(fit)
}


check_residuals <- function(fit, version) {
    ## check model residuals
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "gee_fit_", version, ".pdf")
    pdf(file=target, width=5, height=5)
    plot(fit, which=1)
    dev.off()
    cat("plot saved to", target, "\n")
}
