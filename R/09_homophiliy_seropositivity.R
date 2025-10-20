source(here::here("R", "00_setup.R"))
source(here::here("R", "01_load_data.R"))

rodent_models <- read_rds(here("data", "temp", "rodent_models_2025-10-16.rds"))
included_networks <- c(11, 12, 13, 14, 15, 21, 23, 24, 25, 26, 29) # same as included models in 05_interpreting_network_models.R
rodent_homophily_models <- rodent_models$homophily[included_networks]
rm(rodent_models)
gc()

# Assessing homophily and seropositivty -----------------------------------

simulated_networks <- lapply(rodent_homophily_models, function(x) {
  simulate(x, nsim = 50, output = "network")
})

rm(rodent_homophily_models)
gc()

get_mnat_homophily <- function(net, species_name = "Mastomys natalensis") {
  
  # Extract vertex attributes
  network_tbl <- tibble(
    species = network::get.vertex.attribute(net, "Species"),
    node_name = network.vertex.names(net),
    observed_node = network::get.vertex.attribute(net, "Observed")
  )
  
  # Observed nodes
  observed_ids <- network_tbl %>%
    filter(observed_node) %>%
    pull(node_name)
  
  # Observed M. natalensis nodes
  mnat_ids <- network_tbl %>%
    filter(observed_node, species == species_name) %>%
    pull(node_name)
  
  # Initialise vectors
  homophily_scores <- numeric(length(mnat_ids))
  degree_scores <- integer(length(mnat_ids))
  
  for (j in seq_along(mnat_ids)) {
    i <- mnat_ids[j]  # focal node
    
    alter_ids <- which(net[i, ] == 1)  # All neighbours (by index)
    alter_names <- network_tbl$node_name[alter_ids]
    
    # Keep only alters that are observed
    observed_alter_ids <- alter_ids[alter_names %in% observed_ids]
    
    # If no observed alters, return NA
    if (length(observed_alter_ids) == 0) {
      homophily_scores[j] <- NA_real_
      degree_scores[j] <- 0L
    } else {
      # Count same-species observed alters
      same_species <- network_tbl$species[observed_alter_ids] == species_name
      homophily_scores[j] <- sum(same_species) / length(observed_alter_ids)
      degree_scores[j] <- length(observed_alter_ids)
    }
  }
  
  # Return data frame
  tibble(
    node = mnat_ids,
    homophily = homophily_scores,
    degree = degree_scores
  )
}

# Ensure network names are assigned
names(simulated_networks) <- as.character(included_networks)

# Extract homophily scores and tag with network ID and simulation ID
mnat_homophily_all_df <- purrr::imap_dfr(
  simulated_networks,
  function(sim_list, net_id) {
    purrr::map_dfr(sim_list, get_mnat_homophily, .id = "sim_id") %>%
      mutate(network = net_id)
  }
) %>%
  left_join(rodent_data %>%
              select(rodent_uid, species, interpretation), by = c("node" = "rodent_uid"))

homophily_positivity <- brm(
  formula = interpretation ~ homophily + (1 | network) + (1 | sim_id),
  family = bernoulli(),
  data = mnat_homophily_all_df,
  cores = 4, chains = 4,
  control = list(adapt_delta = 0.95)
)

homophily_degree_positivity <- brm(
  formula = interpretation ~ homophily + degree + (1 | network) + (1 | sim_id),
  family = bernoulli(),
  data = mnat_homophily_all_df,
  cores = 4, chains = 4,
  control = list(adapt_delta = 0.95)
)

homophily_degree_int_positivity <- brm(
  formula = interpretation ~ homophily * degree + (1 | network) + (1 | sim_id),
  family = bernoulli(),
  data = mnat_homophily_all_df,
  cores = 4, chains = 4,
  control = list(adapt_delta = 0.99)
)

loo1 <- loo(homophily_positivity, cores = 8)
loo2 <- loo(homophily_degree_positivity, cores = 8)
loo3 <- loo(homophily_degree_int_positivity, cores = 8)
loo_compare(loo1, loo2, loo3)

summary_df <- fixef(homophily_degree_int_positivity) %>%
  as.data.frame() %>%
  rownames_to_column("parameter")  %>%
  mutate(across(c(Estimate, Q2.5, Q97.5), exp)) %>%
  rename(OddsRatio = Estimate, OR_Lower = Q2.5, OR_Upper = Q97.5)
