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


fit_gee_diag_age <- function(df, type, subset=c("Healthy", "Lynch", "FAP")) {
    df <- subset(df, Diagnosis %in% subset)
    df <- df[order(df$Subject_ID), ]
    levels <- subset
    if (type == "LGR5") {
        fit <- geepack::geeglm(
            LGR5_Count ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
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
    return(fit)
}


check_residuals <- function(fit) {
    ## check model residuals
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "gee_fit_", version, ".pdf")
    pdf(file=target, width=5, height=5)
    plot(fit, which=1)
    dev.off()
    cat("plot saved to", target, "\n")
}


gee_analysis <- function(data, sum, version, pval=NULL) {
    fit <- fit_model(data)
    check_residuals(fit)

    ## annotate plots with p-value from GEE
    if (!is.null(pval)) {
        ggp <- make_ggp_dot(sum, version)
        ggp <- ggp + annotate(
            "text", x=1.5, y=1.1*max(sum$LGR5_Mean),
            label=paste0("italic(p)", pval), parse=TRUE, size=3
        )
        fp <- "~/projects/fap-lgr5/"
        fn <- paste0("fap_lgr5_count_dot_stat_", version, ".pdf")
        ggsave(paste0(fp, fn), ggp, device="pdf", height=3, width=2, unit="in")
        cat("plot saved to", paste0(fp, fn), "\n")
        ggp <- make_ggp_dot2(sum, "FAP", version)
        ggp <- make_ggp_dot2(sum, "Healthy", version)
    }
    return(fit)
}
