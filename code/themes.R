## source from count_compare.R


ggp_theme_default <- theme(
    panel.background=element_rect(fill="white"),
    panel.grid.major=element_line(color="white"),
    panel.grid.minor=element_line(color="white"),
    plot.background=element_rect(fill="white"),
    plot.margin=margin(t=0.75, r=0.75, b=0.75, l=0.75, unit="lines"),
    plot.title=element_text(size=12, face="plain", hjust=0.5),
    plot.subtitle=element_text(size=10, face="plain", hjust=0.5),
    axis.title.x=element_text(size=10, face="plain"),
    axis.title.y=element_text(size=10, face="plain"),
    axis.text.x=element_text(size=10, face="plain"),
    axis.text.y=element_text(size=10, face="plain"),
    axis.line.x.bottom=element_line(),
    axis.line.y.left=element_line(),
    legend.key=element_rect(fill="white"),
    legend.position="none",
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")
)


ggp_theme_predict <- ggp_theme_default + theme(
    plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
    legend.position=c(0.8, 0.8),
    legend.text=element_text(size=10),
    legend.key=element_blank(),
    legend.key.size=unit(0.5, "cm"),
    legend.key.width=unit(0.5, "cm"),
    legend.spacing.x=unit(0.001, "cm"),
    legend.spacing.y=unit(0.001, "cm"),
    legend.background=element_rect(fill="white", size=0.25, colour="grey")
)


ggp_theme_dotcolour <- ggp_theme_default + theme(
    plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
    axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1),
    legend.position="right",
    legend.spacing.x=unit(0.001, "cm"),
    legend.spacing.y=unit(0.001, "cm")
)


ggp_theme_pointsub <- ggp_theme_default +
    theme(
        panel.border=element_rect(colour="black", size=1.0, fill=NA),
        plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
        axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1),
        strip.text=element_text(colour="white"),
        axis.line.x.bottom=element_blank(),
        axis.line.y.left=element_blank(),
        legend.position="none"
    )
