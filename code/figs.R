interobs_lm_eqn <- function(df) {
    m <- lm(y ~ x, df)
    eq <- substitute(
        italic(y) == a + b * italic(x) * "," ~~ italic(R)^2 == r2,
        list(
            a=format(unname(coef(m)[1]), digits=1),
            b=format(unname(coef(m)[2]), digits=2),
            r2=format(summary(m)$r.squared, digits=2)
        )
    )
    return(as.character(as.expression(eq)))
}


interobs_xbar_s <- function(df) {
    xbars <- substitute(
        italic(bar(x)) == a * "," ~~ italic(s) == b,
        list(
            a=format(mean(df$delta), digits=2),
            b=format(sd(df$delta), digits=2)
        )
    )
    return(as.character(as.expression(xbars)))
}


interobs_scatter <- function(df, observers) {
    ## figure works well as 3 in x 3 in
    ggp_xlab <- paste("Counts by Observer", observers[1])
    ggp_ylab <- paste("Counts by Observer", observers[2])
    ann_xval <- (max(df$x) - min(df$x)) / 2
    ann_yval <- max(df$y) * 0.95
    ggp <- ggplot(df, aes(x=x, y=y)) +
        geom_abline(intercept=0, slope=1, size=0.2, linetype="dashed") +
        geom_point(
            position=position_jitter(width=0.1, height=0.1),
            size=0.4, alpha=0.5
        ) +
        geom_smooth(method="lm", se=FALSE, color=muted("red"), size=0.8) +
        annotate(
            "text", x=ann_xval, y=ann_yval, label=interobs_lm_eqn(df),
            size=3, parse=TRUE
        ) +
        labs(x=ggp_xlab, y=ggp_ylab) +
        ggp_theme_default
    return(ggp)
}


interobs_blandaltman <- function(df, observers) {
    ## figure works well as 3 in x 3 in
    ggp_xlab <- paste(
        "Mean of Counts by Observers",
        paste(observers, collapse=", ")
    )
    ggp_ylab <- paste("Counts by Observer", paste(observers, collapse=" - "))
    ann_xval <- (max(df$mean) - min(df$mean)) / 2
    ann_yval <- 2.25
    ggp <- ggplot(df, aes(x=mean, y=delta)) +
        geom_point(
            position=position_jitter(width=0.1, height=0.05),
            size=0.4, alpha=0.5
        ) +
        geom_abline(
            intercept=mean(df$delta), slope=0, color=muted("red"), size=0.8
        ) +
        geom_abline(
            intercept=mean(df$delta) + 1.96 * sd(df$delta), slope=0, size=0.25,
            linetype="dashed"
        ) +
        geom_abline(
            intercept=mean(df$delta) - 1.96 * sd(df$delta), slope=0, size=0.25,
            linetype="dashed"
        ) +
        annotate(
            "text", x=ann_xval, y=ann_yval, label=interobs_xbar_s(df),
            size=3, parse=TRUE
        ) +
        labs(x=ggp_xlab, y=ggp_ylab) +
        scale_y_continuous(limits=c(-2.5, 2.5)) +
        ggp_theme_default
    return(ggp)
}


make_interobs_plots <- function(df, observers) {
    df$delta <- df$x - df$y
    df$mean <- (df$x + df$y) / 2
    pl <- list()
    pl$scatter <- interobs_scatter(df, observers)
    pl$blandalt <- interobs_blandaltman(df, observers)
    return(pl)
}


interobserver_plots <- function(l) {
    pl <- list()

    observers <- c("A", "B")
    type <- "LGR5"
    df <- data.frame(x=l$lgr5$LGR5_Count_A, y=l$lgr5$LGR5_Count_B)
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    observers <- c("A", "B")
    type <- "Ectopic"
    df <- data.frame(x=l$ecto$Ectopic_Count_A, y=l$ecto$Ectopic_Count_B)
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    return(pl)
}


interobs_write <- function(pl, type) {
    ## type is LGR5 or Ectopic
    cwp <- list()
    pl_scatter <- list(
        AB=pl[[paste(type, "AB", sep="_")]]$scatter
    )
    pl_blandalt <- list(
        AB=pl[[paste(type, "AB", sep="_")]]$blandalt
    )
    cwp$scatter <- cowplot::plot_grid(
        plotlist=pl_scatter, nrow=1, labels=c(NULL)
    )
    cwp$blandalt <- cowplot::plot_grid(
        plotlist=pl_blandalt, nrow=1, labels=c(NULL)
    )
    cwp$full <- cowplot::plot_grid(
        plotlist=cwp, nrow=2, labels=c("a", "b")
    )
    plot_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(plot_dir, "interobs_", tolower(type), ".pdf")
    pdf(target, height=6, width=3)
    print(cwp$full)
    dev.off()
    cat("plot written to", target, "\n")
}



make_ggp_dot <- function(data, version) {
    ggp <- ggplot(data, aes(x=Diagnosis, y=LGR5_Mean)) +
        geom_dotplot(
            binaxis="y", stackdir="center", stackratio=1, dotsize=0.25,
            binwidth=1, colour="grey", fill="grey"
        ) +
        stat_summary(
            fun.data=mean_sdl, fun.args=list(mult=1),
            geom="pointrange", colour="black", fill="black", shape=23,
            size=0.8
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count (mean)", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_count_dot_", version, ".pdf")
    ggsave(paste0(fp, fn), ggp, device="pdf", height=3, width=2, unit="in")
    cat("plot saved to", paste0(fp, fn), "\n")
    return(ggp)
}


make_ggp_dot2 <- function(data, group, version) {
    ## colour the Healthy or FAP dots
    ggp <- ggplot(data, aes(x=Diagnosis, y=LGR5_Mean)) +
        geom_dotplot(
            aes(colour=ID, fill=ID),
            binaxis="y", binwidth=1, binpositions="all",
            stackdir="center", stackratio=1, stackgroups=TRUE,
            dotsize=0.25
        ) +
        stat_summary(
            fun.data=mean_sdl, fun.args=list(mult=1),
            geom="pointrange", colour="black", fill="black", shape=23,
            size=0.8
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count (mean)", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default + theme(legend.position="right")
    if (group == "Healthy") {
        if (version == "v1") {
            ggp <- ggp +
                scale_colour_manual(
                    breaks=c("H1", "H2", "H3", "H4", "H5",
                        "H6", "H7", "H8", "H9", "H10"),
                    values=c(rep("grey", 4), scales::hue_pal()(10))
                ) +
                scale_fill_manual(
                    breaks=c("H1", "H2", "H3", "H4", "H5",
                        "H6", "H7", "H8", "H9", "H10"),
                    values=c(rep("grey", 4), scales::hue_pal()(10))
                ) +
                guides(colour=guide_legend(ncol=2))
        } else if (version == "v2") {
            ggp <- ggp +
                scale_colour_manual(
                    breaks=c("H1", "H2", "H3", "H4", "H5",
                        "H6", "H7", "H8", "H9", "H10"),
                    values=c(rep("grey", 3), scales::hue_pal()(10))
                ) +
                scale_fill_manual(
                    breaks=c("H1", "H2", "H3", "H4", "H5",
                        "H6", "H7", "H8", "H9", "H10"),
                    values=c(rep("grey", 3), scales::hue_pal()(10))
                ) +
                guides(colour=guide_legend(ncol=2))
        }
        w <- 3.5
    } else {
        if (version == "v1") {
            ggp <- ggp +
                scale_colour_manual(
                    breaks=c("F1", "F2", "F3", "F4"),
                    values=c(scales::hue_pal()(4), rep("grey", 10))
                ) +
                scale_fill_manual(
                    breaks=c("F1", "F2", "F3", "F4"),
                    values=c(scales::hue_pal()(4), rep("grey", 10))
                )
        } else if (version == "v2") {
            ggp <- ggp +
                scale_colour_manual(
                    breaks=c("F1", "F2", "F3"),
                    values=c(scales::hue_pal()(3), rep("grey", 10))
                ) +
                scale_fill_manual(
                    breaks=c("F1", "F2", "F3"),
                    values=c(scales::hue_pal()(3), rep("grey", 10))
                )
        }
        w <- 2.75
    }
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_count_dot2_", substr(group, 1, 1), version, ".pdf")
    ggsave(paste0(fp, fn), ggp, device="pdf", height=3, width=w, unit="in")
    cat("plot saved to", paste0(fp, fn), "\n")
    return(ggp)
}


make_ggp_point <- function(data, version) {
    ## color the FAP samples
    data[data$Phenotype=="Healthy", "Label"] <- "H"
    data[data$Phenotype=="FAP", "Label"] <- data[data$Phenotype=="FAP", "Subject_ID1"]
    ggp <- ggplot(data, aes(x=Phenotype, y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1) +
        geom_point(
            aes(colour=Label), size=1, position=position_jitter(width=0.2)
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default + theme(legend.position="right")
    if (version == "v1") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(scales::hue_pal()(4), "grey")
            )
    } else if (version == "v2") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3"),
                values=c(scales::hue_pal()(3), rep("grey", 10))
            )
    }
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_count_point_", version, ".pdf")
    ggsave(paste0(fp, fn), ggp, device="pdf", height=3.5, width=2.5, unit="in")
    cat("plot saved to", paste0(fp, fn), "\n")
    return(ggp)
}


make_ggp_point_sub <- function(data, subject, version) {
    ## color individual Healthy samples (the high outliers)
    ggp <- ggplot(data, aes(x=Phenotype, y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1) +
        geom_point(
            data=subset(dat, Subject_ID1 != subject),
            colour="grey", size=1, position=position_jitter(width=0.2)
        ) +
        geom_point(
            data=subset(dat, Subject_ID1 == subject),
            colour="black", size=2, position=position_jitter(width=0.2)
        ) +
        labs(title=paste("LGR5+ Cells:", subject), x=element_blank()) +
        scale_y_continuous(
            name="Count", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_count_point_", subject, "_", version, ".pdf")
    ggsave(paste0(fp, fn), ggp, device="pdf", height=3.5, width=1.75, unit="in")
    cat("plot saved to", paste0(fp, fn), "\n")
    return(ggp)
}


make_ggp_point2 <- function(data, version) {
    ## color the Healthy samples
    data[data$Phenotype=="Healthy", "Label"] <- data[data$Phenotype=="Healthy", "Subject_ID1"]
    data[data$Phenotype=="FAP", "Label"] <- "FAP"
    ggp <- ggplot(data, aes(x=Phenotype, y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1) +
        geom_point(
            aes(colour=Label), size=1, position=position_jitter(width=0.2)
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        scale_colour_manual(
            breaks=c("H1", "H2", "H3", "H4", "H5",
                "H6", "H7", "H8", "H9", "H10"),
            values=c("grey", scales::hue_pal()(10))
        ) +
        guides(colour=guide_legend(ncol=2)) +
        ggp_theme_default + theme(legend.position="right")
    fp <- "~/projects/fap-lgr5/"
    fn <- paste0("fap_lgr5_count_point2_", version, ".pdf")
    ggsave(paste0(fp, fn), ggp, device="pdf", height=3.5, width=3.5, unit="in")
    cat("plot saved to", paste0(fp, fn), "\n")
    return(ggp)
}
