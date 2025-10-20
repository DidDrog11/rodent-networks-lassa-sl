source(here::here("R", "00_setup.R"))
source(here("R", "01_load_data.R"))

# Number of rodents trapped
nrow(rodent_data)

# Number of trapnights
43266

# Number of species
rodent_data %>%
  mutate(total = n()) %>%
  group_by(species) %>%
  mutate(n = n(),
         prop = n/total) %>%
  distinct(species, n, prop) %>%
  arrange(-prop)


# Describing communities --------------------------------------------------

left_join(rodent_data, trap_data) %>%
  group_by(landuse, visit) %>%
  summarise(richness = length(unique(species)),
            n = n()) %>%
  group_by(landuse) %>%
  summarise(median_richness = median(richness),
            IQR_richness = IQR(richness),
            median_n = median(n),
            IQR_n = IQR(n))

# Rarefaction based analysis
rodent_full <- rodent_data %>%
  left_join(trap_data, by = c("village", "visit", "trap_uid")) %>%
  mutate(site_id = paste(village, grid_number, visit, sep = "_"))

species_abundance <- rodent_full %>%
  count(site_id, species) %>%
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>%
  arrange(site_id) %>%
  column_to_rownames("site_id")

latlong <-  rodent_full %>%
  vect(geom = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
  project(y = default_CRS)

site_covariates <- rodent_full %>%
  mutate(x = crds(latlong)[ , 1],
         y = crds(latlong)[, 2]) %>%
  group_by(site_id) %>%
  summarise(village = first(village),
            grid_number = first(grid_number),
            visit = first(visit),
            landuse = first(landuse),
            x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE)) %>%
  arrange(site_id) %>%
  column_to_rownames("site_id")

mob_obj <- make_mob_in(species_abundance, site_covariates, coord_names = c("x", "y"), latlong = TRUE)

div_stats <- get_delta_stats(mob_in = mob_obj, env_var = "landuse", type = "discrete")

IBR_summary <- div_stats$S_df %>%
  filter(test == "SAD", sample == "indiv", effort == 10) %>%
  select(env, S = effect, low_CI = low_effect, median = med_effect, high_CI = high_effect)

# Number of antibody positive rodents
rodent_data %>%
  filter(interpretation == "Positive") %>%
  distinct(rodent_uid) %>%
  nrow()

# Positive rodents by species
table_1_df <- rodent_data %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0),
         total_positive = sum(n_positive)) %>%
  group_by(species) %>%
  summarise(n_individuals = n(),
            n_positive = sum(n_positive),
            total_positive = unique(total_positive)) %>%
  mutate(perc_positive = round((n_positive/n_individuals) * 100, 1),
         perc_all_positive = round((n_positive/total_positive) * 100, 1)) %>%
  arrange(-n_positive, -n_individuals) %>%
  mutate(Species = species,
         `Individuals (N)` = n_individuals,
         `LASV Antibody detected (%)` = paste0(n_positive, " (", perc_positive, "%)"),
         `Percentage of all positive individuals` = paste0(perc_all_positive, "%")) %>%
  select(Species, `Individuals (N)`, `LASV Antibody detected (%)`, `Percentage of all positive individuals`)


# Model seroprevalence ----------------------------------------------------

# Prevalence Odds Ratio (POR) by species 
por_df <- rodent_data %>%
  mutate(seropositive = case_when(interpretation == "Positive" ~ 1,
                                  TRUE ~ 0)) %>%
  select(species, seropositive) %>%
  group_by(species) %>%
  filter(n() >= 10) %>%
  ungroup()

# Set "Mastomys natalensis" as the reference category
por_df$species <- relevel(por_df$species, ref = "Mastomys natalensis")
por_df$species <- droplevels(por_df$species)

# Define the prior probability levels
prior_mas = 0.3       # Prior probability for Mastomys natalensis > 10%
prior_non_mas = 0.05  # Prior probability for non-Mastomys, non-Crocidura species (1-10%)
prio_croc = NULL      # Non-informative prior for Crocidura species

# Define the Bayesian logistic regression model
model <- brm(
  seropositive ~ species,
  data = por_df,
  family = bernoulli(link = "logit"),
  prior = set_prior("normal(0, 1)", class = "b"),
  sample_prior = TRUE,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

# Extract the posterior samples of the coefficients
post_samples <- as_draws_df(model, variable = "^b_species", regex = TRUE) %>%
  gather(Species, Estimate, starts_with("b_species"))

species_mapping <- c(
  "Crocidurabuettikoferi" = "Crocidura buettikoferi",
  "Crociduragrandiceps" = "Crocidura grandiceps",
  "Crociduraolivieri" = "Crocidura olivieri",
  "Lemniscomysstriatus" = "Lemniscomys striatus",
  "Lophuromyssikapusi" = "Lophuromys sikapusi",
  "Malacomysedwardsi" = "Malacomys edwardsi",
  "Musmusculus" = "Mus musculus",
  "Mussetulosus" = "Mus setulosus",
  "Praomysrostratus" = "Praomys rostratus",
  "Rattusrattus" = "Rattus rattus"
)

# Remove the "b_species" prefix and rename the species
summary_df <- post_samples %>%
  mutate(
    Species = gsub("^b_species", "", Species),
    Species = species_mapping[Species],
    Species = fct_rev(factor(Species, levels = species_order_n))
  )

# Calculate the mean and 95% credible interval
species_summary_stats <- summary_df %>%
  group_by(Species) %>%
  summarize(
    Mean_log_odds = mean(Estimate),
    CI_low_log_odds = quantile(Estimate, 0.025),
    CI_high_log_odds = quantile(Estimate, 0.975),
    .groups = "drop"
  ) %>%
  # Convert log-odds to odds ratios
  mutate(
    OR = exp(Mean_log_odds),
    CI_low_OR = exp(CI_low_log_odds),
    CI_high_OR = exp(CI_high_log_odds)
  )

# Create a forest plot using ggplot
odds_seropositivity <- ggplot() +
  geom_density_ridges(data = summary_df, 
                      aes(x = exp(Estimate), y = Species), 
                      scale = 1, fill = "lightblue", alpha = 0.6) +
  geom_point(data = species_summary_stats, 
             aes(x = OR, y = Species), 
             size = 3, color = "black") +
  geom_errorbarh(data = species_summary_stats, 
                 aes(xmin = CI_low_OR, xmax = CI_high_OR, y = Species), 
                 height = 0.2, color = "black", size = 1) +
  labs(x = expression('Odds Ratio of Seropositivity vs.'~italic('Mastomys natalensis')~'(Reference)'),
       y = "Species") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic")) +
  scale_x_continuous(breaks = c(0.5, 1, 2, 5)) +
  geom_vline(xintercept = 1, linetype = "dashed") + # Reference OR of 1
  coord_cartesian(xlim = c(0, 6))

save_plot(filename = here("output", "figures", "Figure_2.png"), plot = odds_seropositivity, base_width = 7, bg = "white")

# Add ORs to table
table_1_updated <- table_1_df %>%
  left_join(species_summary_stats %>%
              mutate(`OR (95% CrI)` = paste0(round(OR, 2), " (", round(CI_low_OR, 2), "-", round(CI_high_OR, 2), ")")) %>%
              select(Species, `OR (95% CrI)`),
            by = "Species")
# Add "ref" to table
table_1_updated[1, 5] <- "ref"

write_rds(table_1_updated, here("output", "tables", "Table_1.rds"))

table_1_ft <- flextable(table_1_updated) %>%
  italic(j = 1, part = "body")

write_rds(table_1_ft, here("output", "tables", "Table_1_flextable.rds"))

# Positive rodents by location and visit

rodent_detail <- rodent_data %>%
  select(rodent_uid, trap_uid, species, interpretation) %>%
  left_join(trap_data %>%
              select(date_set, village, visit, trap_uid, landuse),
            by = c("trap_uid"))

rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                         TRUE ~ 0)) %>%
  group_by(village) %>%
  summarise(tested_village = n(),
         n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_village = round(n_positive/tested_village * 100, 1))

rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0)) %>%
  group_by(landuse) %>%
  summarise(tested_landuse = n(),
            n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_landuse = round(n_positive/tested_landuse * 100, 1))

season_tab <- rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0)) %>%
  left_join(visit_season, by = c("visit")) %>%
  group_by(season) %>%
  summarise(tested_season = n(),
            n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_season = round(n_positive/tested_season * 100, 1))

season_tab_chi <- as.table(rbind(c(16, 444-16),
                                 c(23, 240-23)))
dimnames(season_tab_chi) <- list(season = c("Dry", "Rainy"),
                                 serostatus = c("Positive", "Negative"))

(seasonal <- chisq.test(season_tab_chi))


# Exploratory analysis ----------------------------------------------------

# Prevalence Odds Ratio (POR) by village
village_por_df <- rodent_detail %>%
  mutate(seropositive = case_when(interpretation == "Positive" ~ 1,
                                  TRUE ~ 0)) %>%
  select(village, seropositive) %>%
  group_by(village) %>%
  ungroup()

# Set one village (e.g., "Seilama") as the reference category
village_por_df$village <- relevel(as.factor(village_por_df$village), ref = "Seilama")

# Define the Bayesian logistic regression model
village_model <- brm(
  seropositive ~ village,
  data = village_por_df,
  family = bernoulli(link = "logit"),
  prior = set_prior("normal(0, 1)", class = "b"),
  sample_prior = TRUE,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

# Extract the posterior samples of the coefficients
village_post_samples <- as_draws_df(village_model, variable = "^b_village", regex = TRUE) %>%
  gather(Village, Estimate, starts_with("b_village"))

# Remove the "b_village" prefix for clarity
village_summary <- village_post_samples %>%
  mutate(
    Village = gsub("^b_village", "", Village),
    Village = fct_rev(factor(Village, levels = unique(village_por_df$village)))
  ) %>%
  group_by(Village) %>%
  summarize(
    Mean_log_odds = mean(Estimate),
    CI_low_log_odds = quantile(Estimate, 0.025),
    CI_high_log_odds = quantile(Estimate, 0.975),
    .groups = "drop"
  ) %>%
  # Convert log-odds to odds ratios
  mutate(
    OR = exp(Mean_log_odds),
    CI_low_OR = exp(CI_low_log_odds),
    CI_high_OR = exp(CI_high_log_odds)
  )

landuse_por_df <- rodent_detail %>%
  mutate(seropositive = case_when(interpretation == "Positive" ~ 1,
                                  TRUE ~ 0)) %>%
  select(landuse, seropositive) %>%
  group_by(landuse) %>%
  ungroup()

# Set one land use type (e.g., "Village") as the reference category
landuse_por_df$landuse <- relevel(as.factor(landuse_por_df$landuse), ref = "Village")

# Define the Bayesian logistic regression model
landuse_model <- brm(
  seropositive ~ landuse,
  data = landuse_por_df,
  family = bernoulli(link = "logit"),
  prior = set_prior("normal(0, 1)", class = "b"),
  sample_prior = TRUE,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

# Extract the posterior samples of the coefficients
landuse_post_samples <- as_draws_df(landuse_model, variable = "^b_landuse", regex = TRUE) %>%
  gather(LandUse, Estimate, starts_with("b_landuse"))

# Remove the "b_landuse" prefix for clarity
landuse_summary <- landuse_post_samples %>%
  mutate(
    LandUse = gsub("^b_landuse", "", LandUse),
    LandUse = fct_rev(factor(LandUse, levels = unique(landuse_por_df$landuse)))
  ) %>%
  group_by(LandUse) %>%
  summarize(
    Mean_log_odds = mean(Estimate),
    CI_low_log_odds = quantile(Estimate, 0.025),
    CI_high_log_odds = quantile(Estimate, 0.975),
    .groups = "drop"
  ) %>%
  # Convert log-odds to odds ratios
  mutate(
    OR = exp(Mean_log_odds),
    CI_low_OR = exp(CI_low_log_odds),
    CI_high_OR = exp(CI_high_log_odds)
  )

