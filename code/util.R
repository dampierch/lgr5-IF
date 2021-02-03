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
    dfs <- list()
    fnames <- c("Crypt_ID", "Slide_Number")
    snames <- c("_A", "_B")
    dfs$full <- merge(data$crypts, data$subjects, by="Subject_ID")
    dfs$full <- merge(dfs$full, data$counts_chd, by=fnames)
    dfs$full <- merge(dfs$full, data$counts_sjp, by=fnames, suffixes=snames)
    #dfs$full <- merge(dfs$full, data$counts_ltj, by=fnames)
    i <- !is.na(dfs$full$LGR5_Count_A) & !is.na(dfs$full$LGR5_Count_B)
    dfs$lgr5 <- dfs$full[i, ]
    i <- !is.na(dfs$full$Ectopic_Count_A) & !is.na(dfs$full$Ectopic_Count_B)
    dfs$ecto <- dfs$full[i, ]
    dfs$lgr5$LGR5_Count <- (dfs$lgr5$LGR5_Count_A + dfs$lgr5$LGR5_Count_B) / 2
    dfs$lgr5$Ectopic_Count <- (dfs$lgr5$Ectopic_Count_A + dfs$lgr5$Ectopic_Count_B) / 2
    dfs$ecto$LGR5_Count <- (dfs$ecto$LGR5_Count_A + dfs$ecto$LGR5_Count_B) / 2
    dfs$ecto$Ectopic_Count <- (dfs$ecto$Ectopic_Count_A + dfs$ecto$Ectopic_Count_B) / 2
    dfs$sumttl <- sum_df(dfs$lgr5)
    dfs$sumect <- sum_df(dfs$ecto)
    return(dfs)
}
