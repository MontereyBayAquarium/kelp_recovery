

rm(list=ls())


################################################################################
#Prep workspace

#required packages
require(librarian)
librarian::shelf(tidyverse, readxl, googledrive, readr, stringr, janitor, here,
                 googlesheets4)

#drive_auth()

#set directories 
gdir <- "1vsT-_TrHs0A3xWBG7Gd3sNPXxP3a7wxX"
outdir <- here::here("output")

#read scan metadata
scan_meta <- read_sheet("16oADjtgufnUgFXe1sEjKvXHFCJitIH5oLe7-ZOf8mtY", sheet=1)
scan_codes <- read_sheet("16oADjtgufnUgFXe1sEjKvXHFCJitIH5oLe7-ZOf8mtY", sheet=2)

#load scan data from Google drive
scan_raw <- drive_ls(as_id(gdir), type = "csv")

#download scan data and compile into a single file
scan_dat <- purrr::map_dfr(scan_raw$name, function(file_name) {
  # Get the file info row
  file_row <- scan_raw %>% filter(name == file_name)
  
  # Download temporarily
  temp_path <- tempfile(fileext = ".csv")
  drive_download(file = file_row$id, path = temp_path, overwrite = TRUE)
  
  # Extract date from filename
  date_str <- str_extract(file_name, "\\d{8}")
  date_parsed <- as.Date(date_str, format = "%m%d%Y")
  
  # Read and add date, force all columns to character to avoid type mismatch
  read_csv(temp_path, col_types = cols(.default = "c")) %>%
    mutate(date = date_parsed)
})

################################################################################
#Tidy up dataframe

# Step 1: Clean scan data
scan_build1 <- scan_dat %>%
  select(-1) %>%
  clean_names() %>%
  mutate(
    across(c(ind, large, small, temperatur), as.numeric),
    across(c(behav, canopy, kelptype, seafix, prey, prey2, prey3, prey4), as.character),
    temp = temperatur,
    year = year(date),
    quarter = quarter(date)
  ) %>%
  select(-temperatur)

# Step 2: Flatten scan_codes$code and standardize column names
scan_codes_clean <- scan_codes %>%
  mutate(
    code = sapply(code, as.character),     # flatten list column
    column_name = tolower(column_name)     # match scan_build1 column names
  )

# Step 3: Helper function to replace code with description
replace_code_with_desc <- function(df, codes, colname) {
  codes_sub <- codes %>%
    filter(column_name == colname) %>%
    mutate(code = as.character(code)) %>%
    select(code, value) %>%
    rename(code_key = code, value_desc = value)

  df %>%
    mutate(across(all_of(colname), as.character)) %>%
    rename(code_key = all_of(colname)) %>%
    left_join(codes_sub, by = "code_key") %>%
    select(-code_key) %>%
    rename(!!colname := value_desc)
}

# Step 4: Replace all relevant columns
cols_to_replace <- c("behav", "canopy", "kelptype", "prey", "prey2", "prey3", "prey4")

scan_build2 <- scan_build1
for (col in cols_to_replace) {
  cat("Replacing column with descriptions:", col, "\n")
  scan_build2 <- replace_code_with_desc(scan_build2, scan_codes_clean, col)
}

# Step 5: tidy and organize

scan_build3 <- scan_build2 %>%
                select(date, year, quarter, lat, long, ind, large, small, 
                       seafix, visibility, sky, wind, dir, seaop, swell,
                       temp, behav, canopy, kelptype, prey, prey2, prey3, prey4)


################################################################################
#Explore data structure

#glance at survey days per quarter
sampling_days_per_quarter <- scan_build3 %>%
  mutate(
    year = year(date),
    quarter = quarter(date)
  ) %>%
  distinct(date, year, quarter) %>%
  count(year, quarter, name = "n_sampling_days") %>%
  arrange(year, quarter)

ggplot(sampling_days_per_quarter, aes(x = factor(quarter), y = n_sampling_days, fill = factor(year))) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of Sampling Days per Quarter by Year",
    x = "Quarter",
    y = "Number of Unique Sampling Days",
    fill = "Year"
  ) +
  theme_minimal(base_size = 14)



################################################################################
#Export
write.csv(scan_build3, file = file.path(outdir, "scans", "scans_data.csv"), 
          row.names = F) #last write 31 July 2025


                       


















