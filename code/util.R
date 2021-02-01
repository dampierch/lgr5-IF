read_inputs <- function() {
    data_dir <- "~/projects/fap-lgr5/data/"
    data <- list()
    target <- paste0(data_dir, "ann_subjects.tsv")
    data$subjects <- readr::read_tsv(target)
    target <- paste0(data_dir, "ann_crypts.tsv")
    data$crypts <- readr::read_tsv(target)
    target <- paste0(data_dir, "crypt_counts_blind_chd.tsv")
    data$counts_chd <- readr::read_tsv(target)
    return(data)
}


make_data <- function() {
    data <- read_inputs()
    df <- merge(data$crypts, data$subjects, by="Subject_ID")
    df <- merge(df, data$counts_chd, by=c("Crypt_ID", "Slide_Number"))
    dfsum <- df %>% dplyr::group_by(Subject_ID) %>%
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
    return(setNames(list(df, dfsum), c("data", "summary")))
}
