# Equivocal results have been repeated
# The interpretation is an ordered factor so we group by rodent and slice the group to retain the non-equivocal (final) result

rodent_elisa <- read_rds(here("data", "input", "ELISA_output.rds")) %>%
  drop_na(rodent_uid) %>%
  group_by(rodent_uid) %>% 
  arrange(interpretation) %>%
  slice(1)

rodent_data <- read_csv(here("data", "input", "rodent_data.csv")) %>%
  left_join(rodent_elisa %>%
              select(rodent_uid, result, interpretation, valid)) %>%
  mutate(village = factor(village, levels = village_order, labels = str_to_sentence(village_order)),
         species = factor(clean_names, levels = str_replace_all(str_to_lower(all_species_order), " ", "_"), labels = all_species_order),
         # convert trap uid to first night as will be grouped
         trap_uid = paste0(str_split(trap_uid, pattern = "_", simplify = TRUE)[, 1], "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 2], "_",
                           1, "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 4], "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 5]),
         trap_uid = factor(trap_uid)) %>%
  select(rodent_uid, village, visit, trap_uid, species, genus, interpretation, sex, age_group, weight, head_body, hind_foot, ear, length_skull) %>%
  arrange(village, visit, rodent_uid) %>%
  mutate(rodent_uid = fct_inorder(rodent_uid))

trap_data <- read_csv(here("data", "input", "trap_data.csv")) %>%
  mutate(landuse = factor(case_when(str_detect(habitat_group, "agriculture") ~ "agriculture",
                                    str_detect(habitat_group, "village") ~ "village",
                                    !village %in% c("lambayama", "baiama") & str_detect(habitat_group, "forest") ~ "forest",
                                    village == "baiama" & grid_number == 2 ~ "agriculture",
                                    village == "baiama" & grid_number == 1 ~ "forest",
                                    village == "lambayama" & grid_number == 3 ~ "agriculture"), levels = c("forest", "agriculture", "village"))) %>%
  mutate(village = factor(village, levels = village_order, labels = str_to_sentence(village_order)),
         landuse = factor(landuse, levels = str_to_lower(names(landuse_palette)), labels = str_to_sentence(names(landuse_palette)))) %>%
  select(date_set, village, trap_uid, visit, grid_number, trap_number, landuse, geometry) %>%
  # Remove 'c(' and ')' and separate the longitude and latitude
  mutate(geometry = gsub("[c()]", "", geometry),  # Remove 'c(' and ')'
         geometry = strsplit(geometry, ", ")) %>%  # Split into list of coordinates
  # Convert the cleaned geometry into a simple feature column
  mutate(geometry = st_sfc(map(geometry, ~st_point(as.numeric(.))))) %>%
  # Convert to an sf object with the specified CRS
  st_as_sf(crs = default_CRS) %>%
  st_transform(crs = SL_UTM) %>%
  mutate(
    trap_easting = st_coordinates(.)[, 1],  # Extract easting (x-coordinate)
    trap_northing = st_coordinates(.)[, 2]  # Extract northing (y-coordinate)
  ) %>%
  as_tibble() %>%
  arrange(village, visit, grid_number, trap_number) %>%
  mutate(trap_uid = fct_inorder(trap_uid)) %>%
  select(-geometry)
