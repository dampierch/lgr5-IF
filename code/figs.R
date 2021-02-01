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
