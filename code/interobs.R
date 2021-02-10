## source from count_compare.R


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
            aes(colour=lab, size=lab, alpha=lab),
            position=position_jitter(width=0.1, height=0.1)
        ) +
        geom_smooth(method="lm", se=FALSE, color=muted("red"), size=0.8) +
        annotate(
            "text", x=ann_xval, y=ann_yval, label=interobs_lm_eqn(df),
            size=3, parse=TRUE
        ) +
        labs(x=ggp_xlab, y=ggp_ylab) +
        scale_size_manual(values=c(1, 0.5, 0.5)) +
        scale_colour_manual(values=c("red", "black", "black")) +
        scale_alpha_manual(values=c(0.75, 0.5, 0.5)) +
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
            aes(colour=lab, size=lab, alpha=lab),
            position=position_jitter(width=0.1, height=0.075)
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
        scale_size_manual(values=c(1, 0.5, 0.5)) +
        scale_colour_manual(values=c("red", "black", "black")) +
        scale_alpha_manual(values=c(0.75, 0.5, 0.5)) +
        scale_y_continuous(limits=c(-3.05, 3.05), breaks=-3:3) +
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
    df <- data.frame(
        x=l$lgr5$LGR5_Count_A, y=l$lgr5$LGR5_Count_B, lab=l$lgr5$Diagnosis
    )
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    observers <- c("A", "B")
    type <- "Ectopic"
    df <- data.frame(
        x=l$ectopic$Ectopic_Count_A, y=l$ectopic$Ectopic_Count_B,
        lab=l$ectopic$Diagnosis
    )
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    observers <- c("A", "C")
    type <- "LGR5"
    df <- data.frame(
        x=l$lgr5$LGR5_Count_A, y=l$lgr5$LGR5_Count_C, lab=l$lgr5$Diagnosis
    )
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    observers <- c("B", "C")
    type <- "LGR5"
    df <- data.frame(
        x=l$lgr5$LGR5_Count_B, y=l$lgr5$LGR5_Count_C, lab=l$lgr5$Diagnosis
    )
    lab <- paste(type, paste(observers, collapse=""), sep="_")
    pl[[lab]] <- make_interobs_plots(df, observers)

    return(pl)
}


interobs_write <- function(pl, type) {
    ## type is LGR5 or Ectopic
    cwp <- list()
    ht <- 6
    wd <- 3
    pl_scatter <- list(
        AB=pl[[paste(type, "AB", sep="_")]]$scatter
    )
    pl_blandalt <- list(
        AB=pl[[paste(type, "AB", sep="_")]]$blandalt
    )
    if (type == "LGR5") {
        pl_scatter$AC <- pl[[paste(type, "AC", sep="_")]]$scatter
        pl_blandalt$AC <- pl[[paste(type, "AC", sep="_")]]$blandalt
        pl_scatter$BC <- pl[[paste(type, "BC", sep="_")]]$scatter
        pl_blandalt$BC <- pl[[paste(type, "BC", sep="_")]]$blandalt
        wd <- 9
    }
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
    pdf(target, height=6, width=wd)
    print(cwp$full)
    dev.off()
    cat("plot written to", target, "\n")
}
