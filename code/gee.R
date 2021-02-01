gee_analysis <- function(data, sum, version, pval=NULL) {
    ## statistical test with gee model
    fit <- geepack::geeglm(
        LGR5_Count ~ factor(Phenotype) + as.numeric(Age),
        family=gaussian(),
        data=data, id=factor(Subject_ID1),
        zcor=NULL, corstr="exchangeable", std.err="san.se"
    )
    cat("version:", version, "\n")
    print(summary(fit))
    ## check model residuals
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_gee_fit_", version, ".pdf")
    pdf(file=paste0(fp, fn))
    plot(fit, which=1)
    dev.off()
    cat("plot saved to", paste0(fp, fn), "\n")
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
