fit_model <- function(data) {
    ## statistical test with gee model
    levels <- c("Healthy", "Lynch", "FAP")
    fit1 <- geepack::geeglm(
        LGR5_Count ~ factor(Diagnosis, levels=levels),
        family=gaussian(),
        data=data, id=factor(Subject_ID),
        zcor=NULL, corstr="exchangeable", std.err="san.se"
    )
    fit2 <- geepack::geeglm(
        LGR5_Count ~ factor(Diagnosis, levels=levels) + as.numeric(Age),
        family=gaussian(),
        data=data, id=factor(Subject_ID),
        zcor=NULL, corstr="exchangeable", std.err="san.se"
    )
    cat("version:", version, "\n")
    print(summary(fit1))
    print(summary(fit2))
    return(fit2)
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
