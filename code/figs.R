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


make_ggp_age_predict <- function(df, fit) {
    ## figure works well as 3 in x 3 in
    ages <- 10:90
    dfpred <- data.frame(Ages=ages, LGR5pred=predict(fit, data.frame(Age=ages)))
    levels <- c("Healthy", "Lynch", "FAP")
    ggp <- ggplot(df, aes(
            x=Age, y=LGR5_Mean, colour=factor(Diagnosis, levels=levels)
        )) +
        geom_point(size=2) +
        geom_line(
            color=muted("red"), data=dfpred, aes(x=Ages, y=LGR5pred), size=0.8
        ) +
        labs(y="LGR5+ Cell Count") +
        scale_colour_grey(start=0.8, end=0.2) +
        ggp_theme_predict
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "age_predict.pdf")
    ggsave(target, ggp, device="pdf", height=3, width=3, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}


make_ggp_dot <- function(df) {
    levels <- c("Healthy", "Lynch", "FAP")
    ymax <- 1.05 * max(df$LGR5_Mean)
    dfseg1 <- data.frame(x1=1.1, x2=2.9, y=ymax)
    dfseg2 <- data.frame(x1=2.1, x2=2.9, y=ymax - 0.5)
    ggp <- ggplot(df, aes(x=factor(Diagnosis, levels=levels), y=LGR5_Mean)) +
        geom_dotplot(
            binaxis="y", stackdir="center", stackratio=1, dotsize=0.75,
            binwidth=0.5, colour="grey", fill="grey"
        ) +
        stat_summary(
            fun.data=mean_sdl, fun.args=list(mult=1),
            geom="pointrange", colour="black", fill="black", shape=23,
            size=0.8
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        geom_segment(aes(x=x1, xend=x2, y=y, yend=y), data=dfseg1) +
        geom_segment(aes(x=x1, xend=x2, y=y, yend=y), data=dfseg2) +
        annotate("text", x=2, y=ymax + 0.1, label="*", size=5) +
        annotate("text", x=2.5, y=ymax - 0.4, label="*", size=5) +
        scale_y_continuous(
            name="Mean Count Per Subject",
            breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "main_effect_dot.pdf")
    ggsave(target, ggp, device="pdf", height=3, width=3, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}


make_dot_legend <- function(df, group, levels) {
    ## colour the Healthy, Lynch, FAP dots
    ## Subject_ID levels must come from data.frame sorted for dotplot by x and y
    ggp <- ggplot(df, aes(x=Diagnosis, y=LGR5_Mean)) +
        geom_dotplot(
            aes(
                colour=factor(Subject_ID, levels=levels),
                fill=factor(Subject_ID, levels=levels)
            ),
            binaxis="y", binwidth=0.5, dotsize=0.75,
            stackdir="center", stackratio=1
        ) +
        guides(colour=guide_legend(ncol=2)) +
        ggp_theme_dotcolour
    if (group == "Healthy") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("H1", "H2", "H3", "H4", "H5",
                    "H6", "H7", "H8", "H9", "H10"),
                values=c(scales::hue_pal()(10), rep("grey", 7), rep("grey", 4))
            ) +
            scale_fill_manual(
                breaks=c("H1", "H2", "H3", "H4", "H5",
                    "H6", "H7", "H8", "H9", "H10"),
                values=c(scales::hue_pal()(10), rep("grey", 7), rep("grey", 4))
            )
    } else if (group == "Lynch") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("L1", "L2", "L3", "L4", "L5", "L6", "L7"),
                values=c(rep("grey", 10), scales::hue_pal()(7), rep("grey", 4))
            ) +
            scale_fill_manual(
                breaks=c("L1", "L2", "L3", "L4", "L5", "L6", "L7"),
                values=c(rep("grey", 10), scales::hue_pal()(7), rep("grey", 4))
            )
    } else {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(rep("grey", 10), rep("grey", 7), scales::hue_pal()(4))
            ) +
            scale_fill_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(rep("grey", 10), rep("grey", 7), scales::hue_pal()(4))
            )
    }
    leg <- cowplot::get_legend(ggp)
    return(leg)
}


make_ggp_dot_colour <- function(df, group) {
    ## colour the Healthy, Lynch, FAP dots
    ## fill inside aes() disrupts position of dotplot
    ## https://github.com/tidyverse/ggplot2/pull/1096#issuecomment-316443200
    levels <- c("Healthy", "Lynch", "FAP")
    df <- df[with(df, order(factor(Diagnosis, levels=levels), LGR5_Mean)), ]
    if (group == "Healthy") {
        df$colour <- c(scales::hue_pal()(10), rep("grey", 7), rep("grey", 4))
    } else if (group == "Lynch") {
        df$colour <- c(rep("grey", 10), scales::hue_pal()(7), rep("grey", 4))
    } else {
        df$colour <- c(rep("grey", 10), rep("grey", 7), scales::hue_pal()(4))
    }
    ggp <- ggplot(df, aes(x=factor(Diagnosis, levels=levels), y=LGR5_Mean)) +
        geom_dotplot(
            fill=df$colour, colour=df$colour,
            binaxis="y", binwidth=0.5, dotsize=0.75,
            stackdir="center", stackratio=1,
        ) +
        stat_summary(
            fun.data=mean_sdl, fun.args=list(mult=1),
            geom="pointrange", colour="black", fill="black", shape=23,
            size=0.8
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Mean Count Per Subject",
            breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_dotcolour
    leg <- make_dot_legend(df, group, df$Subject_ID)
    cwp <- cowplot::plot_grid(
        ggp, leg, ncol=2, labels=c(NULL, NULL), rel_widths=c(1, 0.4)
    )
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "subject_dot_", substr(group, 1, 1), ".pdf")
    pdf(target, height=3, width=3.25)
    print(cwp)
    dev.off()
    cat("plot saved to", target, "\n")
    return(cwp)
}


make_ggp_point <- function(df) {
    levels <- c("Healthy", "Lynch", "FAP")
    ymax <- 1.05 * max(df$LGR5_Count)
    dfseg1 <- data.frame(x1=1.1, x2=2.9, y=ymax)
    dfseg2 <- data.frame(x1=2.1, x2=2.9, y=ymax - 0.5)
    ggp <- ggplot(df, aes(x=factor(Diagnosis, levels=levels), y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1, width=0.5) +
        geom_point(
            size=1.25, shape=16, alpha=0.5, position=position_jitter(width=0.2)
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        geom_segment(aes(x=x1, xend=x2, y=y, yend=y), data=dfseg1) +
        geom_segment(aes(x=x1, xend=x2, y=y, yend=y), data=dfseg2) +
        annotate("text", x=2, y=ymax + 0.1, label="*", size=5) +
        annotate("text", x=2.5, y=ymax - 0.4, label="*", size=5) +
        scale_y_continuous(
            name="Count Per Crypt", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        ggp_theme_default +
        theme(
            axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1)
        )
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "main_effect_point.pdf")
    ggsave(target, ggp, device="pdf", height=3, width=2, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}


make_ggp_point_colour <- function(df, group) {
    ## color the Healthy, Lynch, FAP crypts
    levels <- c("Healthy", "Lynch", "FAP")
    i <- df$Diagnosis == group
    df[i, "Label"] <- df[i, "Subject_ID"]
    i <- df$Diagnosis != group
    df[i, "Label"] <- "O"
    ggp <- ggplot(df, aes(x=factor(Diagnosis, levels=levels), y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1, width=0.5) +
        geom_point(
            aes(colour=Label), size=1, position=position_jitter(width=0.2),
            alpha=0.75
        ) +
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count Per Crypt", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        guides(colour=guide_legend(ncol=2)) +
        ggp_theme_dotcolour
    if (group == "Healthy") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("H1", "H2", "H3", "H4", "H5",
                    "H6", "H7", "H8", "H9", "H10"),
                values=c(scales::hue_pal()(10), "grey")
            )
    } else if (group == "Lynch") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("L1", "L2", "L3", "L4", "L5", "L6", "L7"),
                values=c(scales::hue_pal()(7), "grey")
            )
    } else {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(scales::hue_pal()(4), "grey")
            )
    }
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "subject_point_", substr(group, 1, 1), ".pdf")
    ggsave(target, ggp, device="pdf", height=3, width=3.25, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}

## fix border, change subject order
make_ggp_point_sub <- function(df) {
    ## highlight crypts from each subject individually
    levels <- c("Healthy", "Lynch", "FAP")
    keep <- c("Diagnosis", "LGR5_Count", "Subject_ID")
    df <- df[ , keep]
    dfs <- list()
    for (s in unique(df$Subject_ID)) {
        tmp <- df
        tmp$sid <- s
        tmp$lab <- NA
        tmp[tmp$Subject_ID == s, "lab"] <- TRUE
        dfs[[s]] <- tmp
    }
    df <- do.call(rbind, dfs)
    ggp <- ggplot(df, aes(x=factor(Diagnosis, levels=levels), y=LGR5_Count)) +
        geom_boxplot(outlier.size=-1, width=0.5) +
        geom_point(
            aes(colour=lab), size=1, position=position_jitter(width=0.2),
            alpha=0.75
        ) +
        labs(title="LGR5+ Cells", x=element_blank()) +
        scale_y_continuous(
            name="Count", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        scale_colour_manual(
            values=c(muted("red"))
        ) +
        facet_wrap(vars(sid), nrow=6) +
        ggp_theme_dotcolour + theme(legend.position="none")
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "subject_point_facet.pdf")
    ggsave(target, ggp, device="pdf", height=10, width=6, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}
