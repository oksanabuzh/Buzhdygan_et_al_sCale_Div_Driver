# Purpose: Combine and clean raw data from 10 vegetation data sets
# Raw data are taken from iDiv cloud datasets_all_modified (tab1-tab11.csv)
# Raw data should be placed in a subfolder 1_prepare_data/data-raw
library(tidyverse)

# Read raw data -----------------------------------------------------------

# list all files and read them into a list
all_dataset_files <- list.files("1_prepare_data/data-raw/", full.names = TRUE)
all_dataset_names <- str_extract(
  list.files("1_prepare_data/data-raw/"),
  "[0-9]+"
)
all_datasets <- lapply(all_dataset_files, data.table::fread)
names(all_datasets) <- all_dataset_names

# transpose to long format and make everything a character
all_datasets <- lapply(all_datasets, function(x) {
  x %>%
    mutate_all(as.character) %>%
    pivot_longer(!c(layer, group, species)) %>%
    drop_na()
})

# combine all data sets
all_datasets <- data.table::rbindlist(all_datasets, idcol = "dataset")

# Clean data --------------------------------------------------------------

# remove special characters from species names
# \xd7 - hybrid
# \xa0 - cross sign
clean_data <- all_datasets %>% mutate(
  species =
    str_replace(species, "\xd7|\xa0", "")
)

clean_data <- mutate(clean_data, species = str_trim(species))

# Separate plot info: series&subplot&scale
clean_data <- clean_data %>%
  separate(name, into = c("series", "subplot", "scale"), sep = "&")

clean_data <- clean_data %>% mutate(
  response = case_when(
    scale == "10" | scale == "100" ~ "cover",
    TRUE ~ "p/a"
  )
)

# Convert all p/a to 1
clean_data <- clean_data %>%
  mutate(value = case_when(
    response == "p/a" & value != "1" ~ "1",
    TRUE ~ value
  ))

# write_csv(clean_data, "/data/complete_database.csv")

# Remove plots using different sampling design
to_remove <- c(
  "SB20138",
  "SB20140",
  "SB20143",
  "NFD21_19",
  "NFD21_22",
  "NFD21_25",
  "NFD21_29",
  "NFD21_32",
  "NFD21_41"
)
clean_data <- clean_data %>%
  filter(!(series %in% to_remove))

# Remove algae
clean_data <- clean_data %>%
  filter(group != "A")

# Remove plots with scale 1000
clean_data <- clean_data %>%
  filter(scale != "1000")

# Fix cover for 100m2 plots  ----------------------------------------------

# dataset 1, 2, and 3 only record presence/absence for 100 m2
# For common species, 100 m2 should be mean of the 2 10 m2 plots
# For rare species (that do not appear in 10 m2), we assume a cover of 0.1

# Problem: There are duplicates in the data -------------------------------
clean_data %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, value, response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# For now: Remove the duplicates from all datasets (including 7 and 8)
# This removes 17 duplicated rows
clean_data <- clean_data %>%
  distinct(dataset, layer, group, species, series, subplot, scale, value, response)

# Are there any duplicates with different cover values?
# There is one species duplicated for one series in 100m2 in dataset 8
clean_data %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

clean_data <- clean_data %>% filter(!(dataset == "8" &
  species == "Polytrichastrum sexangulare" & series == "NFD21_20" & scale == "100" &
    value == "0.5"))


# Calculate the correct cover values for 100m2 plots in dataset 1:3
correct_cover <- clean_data %>%
  filter(scale %in% c(10, 100) & dataset %in% 1:3) %>%
  unite(scale, c(subplot, scale)) %>%
  pivot_wider(names_from = scale, values_from = value) %>%
  mutate(across(c(x_100, NW_10, SE_10), .fns = as.numeric)) %>%
  mutate(NW_10 = case_when(
    is.na(NW_10) ~ 0,
    TRUE ~ NW_10
  )) %>% 
    mutate(SE_10 = case_when(
      is.na(SE_10) ~ 0,
      TRUE ~ SE_10
    )) %>%
  rowwise() %>%
  mutate(x_100 = case_when(
    NW_10==0 & SE_10==0 ~ 0.1,
    TRUE ~ mean(c(NW_10, SE_10), na.rm = TRUE)
  ))

correct_cover <- correct_cover %>%
  select(-NW_10, -SE_10) %>% 
  pivot_longer(x_100, names_to = "scale", values_to = "value") %>%
  separate(scale, into = c("subplot", "scale")) %>%
  mutate(across(everything(), as.character))

# check if now the data looks better
# view(correct_cover)

# Any duplicates?
correct_cover %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, value, response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# Combine the tables back together

# First check if the line numbers match
nrow(clean_data)
nrow(clean_data %>% filter(scale == 100 & dataset %in% 1:3)) +
  nrow(clean_data %>% filter(!(scale == 100 & dataset %in% 1:3)))

clean_data <- clean_data %>%
  filter(!(scale == 100 & dataset %in% 1:3)) %>%
  bind_rows(correct_cover)


clean_data %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, response, value) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# Unifiying layers -------------------------------------------------------

clean_data <- clean_data %>%
  mutate(layer = case_when(
    layer == "Juv" | layer == "Seedling" ~ "H",
    TRUE ~ layer
  ))
# Sum up all cover and p/a values for the different layers for all species
clean_data <- clean_data %>%
  mutate(value = as.numeric(value)) %>%
  group_by(dataset, layer, group, species, series, subplot, scale, response) %>%
  summarize(
    value = sum(value)
  ) %>%
  ungroup()

# change p/a to 1 (p/a will be bigger than 1 if species is present in multiple
# layers)
clean_data <- clean_data %>%
  mutate(value = case_when(
    response == "p/a" & value != 1 ~ 1,
    TRUE ~ value
  ))

# clean species names -----------------------------------------------------
# use species_match.csv to clean the species names
species_correct <- read_csv("data/species_match.csv")

species_correct <- species_correct %>%
  rename(species = species_data) %>%
  select(-group)

clean_data <- left_join(clean_data, species_correct, by = "species")

clean_data <- clean_data %>%
  select(-species) %>%
  rename(species = species_correct)

# Did the cleaning of the species names introduce duplicates?
clean_data %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, response, value) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# Remove duplicates that have the same value
clean_data <- clean_data %>%
  distinct(dataset, layer, group, species, series, subplot, scale, value, response)

# Some duplicates with different cover values in 10 m2 and 100 m2 remain
clean_data %>%
  dplyr::group_by(dataset, layer, group, species, series, subplot, scale, response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# For now, take the mean
clean_data <- clean_data %>%
  group_by(dataset, layer, group, species, series, subplot, scale, response) %>%
  summarize(value = mean(value)) %>%
  ungroup()

# write a unique species list
# write_csv(clean_data %>% select(species) %>% unique(), "data/unique_species_list.csv")

# Check clean data --------------------------------------------------------------

# Check all variables
clean_data %>%
  pull(dataset) %>%
  unique()
clean_data %>%
  pull(layer) %>%
  unique()
clean_data %>%
  pull(scale) %>%
  unique()
clean_data %>%
  pull(subplot) %>%
  unique()
clean_data %>%
  pull(series) %>%
  unique()

clean_data %>%
  count(series, scale) %>%
  arrange(-n)

clean_data %>%
  filter(response == "p/a") %>%
  pull(value) %>%
  unique()

# Check species list
# idea: check if all species in smaller plots are counted in the 100m2 plot
double_check <- list()
for (i in unique(clean_data$series)) {
  species_100 <- clean_data %>%
    filter(series == i & scale == 100) %>%
    pull(species)
  species_10 <- clean_data %>%
    filter(series == i & scale == 10) %>%
    pull(species) %>%
    unique()

  non_matching <- species_10[!(species_10 %in% species_100)]
  if (length(non_matching) > 0) {
    double_check[[i]] <- non_matching
  }
}

double_check <- lapply(double_check, function(x) {
  tibble(species = x)
})
# double_check should be empty, otherwise there are species in the smaller plots
# that are not counted in 100m2
double_check <- data.table::rbindlist(double_check, idcol = "series")

# Check values -> all numbers?
clean_data %>%
  mutate(value = as.numeric(value)) %>%
  summary()

# Match names of plots with header
headers <- read_csv("data/headers.csv")
headers <- headers %>%
  pull(series) %>%
  unique()
# This should all be TRUE, otherwise there are plot IDs in the data that are not
# present in the header
clean_data %>%
  pull(series) %>%
  unique() %in% headers

# Write clean data --------------------------------------------------------
write_csv(clean_data, "data/database_analysis.csv")
