## source from count_compare.R


read_inputs <- function() {
    data_dir <- "~/projects/fap-lgr5/data/"
    data <- list()
    target <- paste0(data_dir, "ann_subjects.tsv")
    data$subjects <- readr::read_tsv(target)
    target <- paste0(data_dir, "ann_crypts.tsv")
    data$crypts <- readr::read_tsv(target)
    target <- paste0(data_dir, "crypt_counts_blind_chd.tsv")
    data$counts_chd <- readr::read_tsv(target)
    target <- paste0(data_dir, "crypt_counts_blind_sjp.tsv")
    data$counts_sjp <- readr::read_tsv(target)
    return(data)
}


merge_inputs <- function(data) {
    fnames <- c("Crypt_ID", "Slide_Number")
    snames <- c("_A", "_B")
    df <- merge(data$crypts, data$subjects, by="Subject_ID")
    df <- merge(df, data$counts_chd, by=fnames)
    df <- merge(df, data$counts_sjp, by=fnames, suffixes=snames)
    return(df)
}


filter_df <- function(df, type, all_obs=FALSE) {
    singletons <- names(table(df$Subject_ID)[table(df$Subject_ID) < 2])
    fnames <- c(paste0(type, "_Count_A"), paste0(type, "_Count_B"))
    i1 <- !is.na(df[ , fnames[1]]) & !is.na(df[ , fnames[2]])
    if (all_obs) {
        i1 <- i1 & !is.na(df[ , paste0(type, "_Count_C")])
    }
    i2 <- !(df$Subject_ID %in% singletons)
    i <- i1 & i2
    df <- df[i, ]
    df$LGR5_Count <- (df$LGR5_Count_A + df$LGR5_Count_B) / 2
    df$Ectopic_Count <- (df$Ectopic_Count_A + df$Ectopic_Count_B) / 2
    if (all_obs) {
        df$LGR5_Count <- (df$LGR5_Count_A + df$LGR5_Count_B + df$LGR5_Count_C) / 3
        df$Ectopic_Count <- (df$Ectopic_Count_A + df$Ectopic_Count_B + df$Ectopic_Count_C) / 3
    }
    return(df)
}


sum_df <- function(df) {
    df1 <- df %>% dplyr::group_by(Subject_ID) %>%
        dplyr::summarise(
            Number=as.numeric(substring(dplyr::first(Subject_ID), 2)),
            ID=dplyr::first(Subject_ID),
            Sex=dplyr::first(Sex),
            Age=dplyr::first(Age),
            Race=dplyr::first(Race),
            Diagnosis=dplyr::first(Diagnosis),
            Crypts=n(),
            LGR5_Mean=mean(LGR5_Count),
            LGR5_SD=sd(LGR5_Count),
            Ectopic_Mean=mean(Ectopic_Count),
            Ectopic_SD=sd(Ectopic_Count)
        ) %>% dplyr::arrange(desc(Diagnosis), Number)
    return(df1)
}


make_data <- function() {
    data <- read_inputs()
    dfs <- list(full=merge_inputs(data))
    dfs$lgr5 <- filter_df(dfs$full, "LGR5")
    dfs$ectopic <- filter_df(dfs$full, "Ectopic")
    dfs$sumttl <- sum_df(dfs$lgr5)
    dfs$sumect <- sum_df(dfs$ectopic)
    dfs$final <- dfs$ectopic
    dfs$sumfin <- dfs$sumect
    return(dfs)
}
