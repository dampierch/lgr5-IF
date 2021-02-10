## source from count_compare.R


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
        labs(x=element_blank()) +
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
                values=c(scales::hue_pal()(10), rep("grey", 6), rep("grey", 4))
            ) +
            scale_fill_manual(
                breaks=c("H1", "H2", "H3", "H4", "H5",
                    "H6", "H7", "H8", "H9", "H10"),
                values=c(scales::hue_pal()(10), rep("grey", 6), rep("grey", 4))
            )
    } else if (group == "Lynch") {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("L1", "L2", "L3", "L5", "L6", "L7"),
                values=c(rep("grey", 10), scales::hue_pal()(6), rep("grey", 4))
            ) +
            scale_fill_manual(
                breaks=c("L1", "L2", "L3", "L5", "L6", "L7"),
                values=c(rep("grey", 10), scales::hue_pal()(6), rep("grey", 4))
            )
    } else {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(rep("grey", 10), rep("grey", 6), scales::hue_pal()(4))
            ) +
            scale_fill_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(rep("grey", 10), rep("grey", 6), scales::hue_pal()(4))
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
        df$colour <- c(scales::hue_pal()(10), rep("grey", 6), rep("grey", 4))
    } else if (group == "Lynch") {
        df$colour <- c(rep("grey", 10), scales::hue_pal()(6), rep("grey", 4))
    } else {
        df$colour <- c(rep("grey", 10), rep("grey", 6), scales::hue_pal()(4))
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
        labs(x=element_blank()) +
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
        labs(x=element_blank()) +
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
            alpha=0.9
        ) +
        labs(x=element_blank()) +
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
                breaks=c("L1", "L2", "L3", "L5", "L6", "L7"),
                values=c(scales::hue_pal()(6), "grey")
            )
    } else {
        ggp <- ggp +
            scale_colour_manual(
                breaks=c("F1", "F2", "F3", "F4"),
                values=c(rev(scales::hue_pal()(4)), "grey")
            )
    }
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "subject_point_", substr(group, 1, 1), ".pdf")
    ggsave(target, ggp, device="pdf", height=3, width=3.25, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}


make_ggp_point_sub <- function(df) {
    ## highlight crypts from each subject individually
    levels <- c("Healthy", "Lynch", "FAP")
    levs_h <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10")
    levs_l <- c("L1", "L2", "L3", "L5", "L6", "L7")
    levs_f <- c("F1", "F2", "F3", "F4")
    levels2 <- c(levs_h, levs_l, levs_f)
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
        labs(title="LGR5+ Cell Count", x=element_blank()) +
        scale_y_continuous(
            name="Count Per Crypt", breaks=seq(0, 9, 2), labels=seq(0, 9, 2)
        ) +
        scale_colour_manual(
            values=c(muted("red"))
        ) +
        facet_wrap(vars(factor(sid, levels=levels2)), nrow=5) +
        ggp_theme_pointsub
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "subject_point_facet.pdf")
    ggsave(target, ggp, device="pdf", height=10, width=6, unit="in")
    cat("plot saved to", target, "\n")
    return(ggp)
}


make_main_fig3 <- function(ggps) {
    dots <- list(hlt=ggps$colsubjhlt, lyn=ggps$colsubjlyn, fap=ggps$colsubjfap)
    pnts <- list(hlt=ggps$colcryphlt, lyn=ggps$colcryplyn, fap=ggps$colcrypfap)
    pl <- list(
        ad=cowplot::plot_grid(
            ggps$mainsubj, NULL, ggps$age, nrow=1, labels=c("a", "", "d"),
            rel_widths=c(1, 0.5, 1), scale=c(1, 1, 0.8)
        ),
        b=cowplot::plot_grid(plotlist=dots, nrow=1, labels=NULL),
        c=cowplot::plot_grid(plotlist=pnts, nrow=1, labels=NULL)
    )
    cwp <- cowplot::plot_grid(
        plotlist=pl, nrow=3, labels=c("", "b", "c"),
        rel_heights=c(1, 0.8, 0.8)
    )
    res_dir <- "~/projects/fap-lgr5/res/"
    target <- paste0(res_dir, "main_fig3.pdf")
    pdf(target, height=10, width=10)
    print(cwp)
    dev.off()
    cat("plot saved to", target, "\n")
    return(cwp)
}
