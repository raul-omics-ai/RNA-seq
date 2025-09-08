########## 24/04/2025 ##########
# ==================================================== #
# ========== AUTOMATE CHECK MD5SUM FUNCTION ==========
# ==================================================== #

# This function will be perform the check of two md5files to check if the download was succesfull
# This function will not perform the data preparation of md5file, so the naming of columns and
# matching both dataframes are not included. So, I'm going to keep the data manipulation code 
# in this function in commented lines

# ======================================= #
# ========== DATA MANIPULATION ========== #
# ======================================= #

#md5sum_original <- read.table("/media/rfernandez/DISCO_2/06_RNAseq_Olink_TATTAA_102101/01_RawData/2025-04-12_716309240-102101.run_1.md5sum.txt")
#colnames(md5sum_original) <- c("md5sum", "path")
#md5sum_original$file_name <- basename(md5sum_original$path)

#md5sum_calculated <- read.table("/media/rfernandez/DISCO_2/06_RNAseq_Olink_TATTAA_102101/01_RawData/md5sum.txt",
#                                header = T)

# ============================== #
# ========== FUNCTION ========== #
# ============================== #
check_md5sum <- function(md5sum_original, md5sum_calculated, where_to_save = NULL) {
  # The preprocessing step is not included in this function, please check if that colnames are
  # the same and sorted in the same order.
  library(dplyr)
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  # Unir dataframes por nombre de archivo
  md5check_df <- md5sum_calculated %>%
    left_join(md5sum_original[, c("file_name", "md5sum")],
              by = "file_name", suffix = c("_calc", "_orig")) %>%
    mutate(md5check = ifelse(md5sum_calc == md5sum_orig, "OK", "NOT THE SAME"))
  
  write.csv(x = md5check_df, file = file.path(where_to_save, "md5sum_check.csv"))
  return(md5check_df)
}

#check_md5sum(md5sum_original = md5sum_original, md5sum_calculated = md5sum_calculated, where_to_save = "/media/rfernandez/DISCO_2/06_RNAseq_Olink_TATTAA_102101/01_RawData/01_RNAseq/")
