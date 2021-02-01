make_data <- function() {
    fp <- "~/projects/fap-lgr5/"
    fn <- "fap_lgr5_count_2.xlsx"
    data <- readxl::read_excel(paste0(fp, fn))
    datsum <- data %>% dplyr::group_by(Subject_ID1) %>%
        dplyr::summarise(
            Number=as.numeric(substring(dplyr::first(Subject_ID1), 2)),
            ID=dplyr::first(Subject_ID1),
            Sex=dplyr::first(Sex),
            Age=dplyr::first(Age),
            Diagnosis=dplyr::first(Phenotype),
            Crypts=n(),
            LGR5_Mean=mean(LGR5_Count),
            LGR5_SD = sd(LGR5_Count)
        ) %>% dplyr::arrange(desc(Diagnosis), Number)
    return(setNames(list(data, datsum), c("data", "summary")))
}
