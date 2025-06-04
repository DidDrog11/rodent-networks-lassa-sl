source(here::here("R", "00_setup.R"))

rodent_models_summary_15 <- read_rds(here("temp", "rodent_models_summary_s_15_2023-12-13.rds"))
rodent_network_15 <- read_rds(here("data", "rodent_network_site_visit_15.rds"))
rodent_data_15 <- read_rds(here("data", "expanded_assemblages_15_2023-12-13.rds"))

final_model_15 <- read_rds(here("temp", "final_model_s_15_2023-12-13.rds"))
final_model_summary_15 <- final_model_15$final_model
`%s%` <- network::`%s%`

homophily_term_15 <- final_model_summary_15

coefficients_15 <- list()

for(i in 1:length(homophily_term_15)) {
  
  if(is.list(homophily_term_15[[i]])) {
    
    coeff = homophily_term_15[[i]][["coefficients"]] %>%
      as.data.frame()
    
    coeff$coefficient = row.names(coeff) %>%
      str_replace(., "nodefactor.Species.", "Species - ") %>%
      str_replace(., "nodematch.Species.", "Match - ") %>%
      str_replace(., "edges", "Edges")
    
    coeff$Network = i
    coeff$Visit = unique(rodent_data_15$nodelist[[i]]$Visit)
    coeff$Landuse = unique(rodent_data_15$nodelist[[i]]$Landuse)
    coeff$`Observed M. natalensis` = rodent_data_15$nodelist[[i]] %>%
      filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
      nrow()
    coeff$`Observed all species` =  rodent_data_15$nodelist[[i]] %>%
      filter(Observed == TRUE) %>%
      nrow()
    coeff$`Unobserved M. natalensis` = rodent_data_15$nodelist[[i]] %>%
      filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
      nrow()
    coeff$N = nrow(rodent_data_15$nodelist[[i]])
    
    coefficients_15[[i]] <- coeff %>%
      tibble() %>%
      select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient = coefficient, Estimate, `Std. Error`)
    
  } else {
    
    coefficients_15[[i]] <- tibble(Network = i,
                                   Visit = unique(rodent_data_15$nodelist[[i]]$Visit),
                                   Landuse = unique(rodent_data_15$nodelist[[i]]$Landuse),
                                   `Observed M. natalensis` = rodent_data_15$nodelist[[i]] %>%
                                     filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
                                     nrow(),
                                   N = nrow(rodent_data_15$nodelist[[i]]),
                                   `Observed all species` =  rodent_data_15$nodelist[[i]] %>%
                                     filter(Observed == TRUE) %>%
                                     nrow(),
                                   `Unobserved M. natalensis` = rodent_data_15$nodelist[[i]] %>%
                                     filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
                                     nrow(),
                                   Coefficient = NA,
                                   Estimate = NA,
                                   `Std. Error` = NA) %>%
      tibble() %>%
      select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient, Estimate, `Std. Error`)
    
    
  }
  
}

coefficients_df_15 <- bind_rows(coefficients_15) %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village")),
         Estimate = case_when(Estimate == -Inf ~ as.numeric(NA),
                              TRUE ~ Estimate),
         Variance = `Std. Error`^2) %>%
  group_by(Landuse) %>%
  mutate(Network = paste0(Landuse, " ", Visit)) %>%
  filter(`Observed M. natalensis` > 1) %>% 
  filter(!Network %in% Network[is.na(Estimate)])

included_models_15 <- bind_rows(coefficients_15) %>%
  group_by(Network) %>% 
  filter(!Network %in% Network[is.na(Estimate) | is.infinite(Estimate)]) %>%
  distinct(Network)


# Speed up processing by removing networks
rm(rodent_data_15)
rm(final_model_15)
gc()
# Meta-analysis of edges --------------------------------------------------

edges_15 <- coefficients_df_15 %>%
  filter(Coefficient == "Edges") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .,
      slab = Network)

# for multi-level meta-analysis
#   rma.mv(yi = Estimate, V = Variance, data = .,
#   random = ~ 1 | Landuse/Visit,
#   slab = Network)

edges_rma_15 <- bind_rows(tibble(ES = edges_15$yi, SE = sqrt(edges_15$vi), Type = "Network", Network = factor(edges_15$data$Network), N = edges_15$data$`Observed M. natalensis`),
                          tibble(ES = edges_15$b[, 1], SE = edges_15$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

edges_combined <- ggplot(edges_rma_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed (all species)",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")


# Meta-analysis of edges stratified by landuse ----------------------------

edges_ag_15 <- coefficients_df_15 %>%
  filter(Coefficient == "Edges" & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_vil_15 <- coefficients_df_15 %>%
  filter(Coefficient == "Edges" & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_rma_df_15 <- bind_rows(tibble(ES = edges_ag_15$yi, SE = sqrt(edges_ag_15$vi), Type = "Network",
                                 Weight = weights(edges_ag_15),
                                 Network = factor(edges_ag_15$data$Network),
                                 N = edges_ag_15$data$`Observed M. natalensis`),
                          tibble(ES = edges_ag_15$b[, 1], SE = edges_ag_15$se, Type = "Summary",
                                 Weight = 100,
                                 Network = factor("Summary - Agriculture"),
                                 N = as.numeric(NA)),
                          tibble(ES = edges_vil_15$yi, SE = sqrt(edges_vil_15$vi), Type = "Network",
                                 Weight = weights(edges_vil_15),
                                 Network = factor(edges_vil_15$data$Network),
                                 N = edges_vil_15$data$`Observed M. natalensis`),
                          tibble(ES = edges_vil_15$b[, 1], SE = edges_vil_15$se, Type = "Summary",
                                 Weight = 100,
                                 Network = factor("Summary - Village"),
                                 N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

edges_stratified_15 <- ggplot(edges_rma_df_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::breaks_pretty()) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  guides(size = "none") +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "Weight",
       title = "Odds of a tie being observed (all species)")


# Meta-analysis of ties based on species type -----------------------------

species_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Species")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma_15 <- bind_rows(tibble(ES = species_15$yi, SE = sqrt(species_15$vi), Type = "Network", Network = factor(species_15$data$Network), N = species_15$data$`Observed M. natalensis`),
                         tibble(ES = species_15$b[, 1], SE = species_15$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

species_combined_15 <- ggplot(species_rma_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed for M. natalensis")

species_ag_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_vil_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma_df_15 <- bind_rows(tibble(ES = species_ag_15$yi, SE = sqrt(species_ag_15$vi), Type = "Network",
                                      Weight = weights(species_ag_15),
                                      Network = factor(species_ag_15$data$Network),
                                      N = species_ag_15$data$`Observed M. natalensis`),
                               tibble(ES = species_ag_15$b[, 1], SE = species_ag_15$se, Type = "Summary", Weight = 100, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                               tibble(ES = species_vil_15$yi, SE = sqrt(species_vil_15$vi), Type = "Network",
                                      Weight = weights(species_vil_15),
                                      Network = factor(species_vil_15$data$Network), N = species_vil_15$data$`Observed M. natalensis`),
                               tibble(ES = species_vil_15$b[, 1], SE = species_vil_15$se, Type = "Summary", Weight = 100, Network = factor("Summary - Village"), N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

species_stratified_15 <- ggplot(species_rma_df_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::pretty_breaks()) +
  theme_bw() +
  guides(size = "none") +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "Weight",
       title = "Odds of a tie being observed for M. natalensis")

# Meta-analysis of intra-specific ties ------------------------------------

match_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Match")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma_15 <- bind_rows(tibble(ES = match_15$yi, SE = sqrt(match_15$vi), Type = "Network", 
                                 Network = factor(match_15$data %>%
                                                    drop_na(Estimate) %>%
                                                    pull(Network)), 
                                 N = match_15$data %>%
                                   drop_na(Estimate) %>%
                                   pull(`Observed M. natalensis`)),
                          tibble(ES = match_15$b[, 1], SE = match_15$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

match_combined_15 <- ggplot(match_rma_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N M. natalensis",
       title = "Odds of a tie between two M. natalensis being observed",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis") +
  coord_cartesian(xlim = c(0, 40))

match_ag_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_vil_15 <- coefficients_df_15 %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma_df_15 <- bind_rows(tibble(ES = match_ag_15$yi, SE = sqrt(match_ag_15$vi), Type = "Network", 
                                    Weight = weights(match_ag_15),
                                    Network = factor(match_ag_15$data %>%
                                                       drop_na(Estimate) %>%
                                                       pull(Network)),
                                    N = match_ag_15$data %>%
                                      drop_na(Estimate) %>%
                                      pull(`Observed M. natalensis`)),
                             tibble(ES = match_ag_15$b[, 1], SE = match_ag_15$se, Type = "Summary", Weight = 100, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                             tibble(ES = match_vil_15$yi, SE = sqrt(match_vil_15$vi), Type = "Network",
                                    Weight = weights(match_vil_15),
                                    Network = factor(match_vil_15$data %>%
                                                       drop_na(Estimate) %>%
                                                       pull(Network)),
                                    N = match_vil_15$data %>%
                                      drop_na(Estimate) %>%
                                      pull(`Observed M. natalensis`)),
                             tibble(ES = match_vil_15$b[, 1], SE = match_vil_15$se, Type = "Summary", Weight = 100, Network = factor("Summary - Village"), N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

match_stratified_15 <- ggplot(match_rma_df_15) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::pretty_breaks()) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = element_blank(),
       title = "Odds of a tie between two M. natalensis being observed") +
  guides(size = "none") +
  coord_cartesian(xlim = c(0, 40))

p1_15 <- edges_stratified_15 + 
  theme(legend.position = "none") +
  labs(title = element_blank())
p2_15 <- species_stratified_15 +
  theme(legend.position = "none") +
  labs(title = element_blank())
p3_15 <- match_stratified_15 +
  theme(legend.position = "none") +
  labs(title = element_blank())
legend <- get_legend(match_stratified_15 +
                       theme(legend.position = "bottom"))

save_plot(plot = plot_grid(p1_15, p2_15, p3_15, legend, ncol = 1, labels = c("A", "B", "C", " "), rel_heights = c(1, 1, 1, 0.2)), filename = here("output", "Sensitivity_1_15m.png"), base_width = 7, base_height = 11)


combined_rma_15 <- list(edges_ag_15,
                        edges_vil_15,
                        species_ag_15,
                        species_vil_15,
                        match_ag_15,
                        match_vil_15)
write_rds(combined_rma_15, here("data", "Sensitivity_1_15m_results.rds"))

# Assessing GOF -----------------------------------------------------------
# If assessing GOF need to comment out the removal of final_model on line 138

models_in_rma <- c(included_models$Network)

gof_list <- final_model$final_model_gof[models_in_rma]

plot(gof_list[[1]])


# 50 meter sensitivity analysis -------------------------------------------

rodent_models_summary_50 <- read_rds(here("temp", "rodent_models_summary_s_50_2023-12-13.rds"))
rodent_network_50 <- read_rds(here("data", "rodent_network_site_visit_50.rds"))
rodent_data_50 <- read_rds(here("data", "expanded_assemblages_50_2023-12-13.rds"))

final_model_50 <- read_rds(here("temp", "final_model_s_50_2023-12-13.rds"))
final_model_summary_50 <- final_model_50$final_model
`%s%` <- network::`%s%`

homophily_term_50 <- final_model_summary_50

coefficients_50 <- list()

for(i in 1:length(homophily_term_50)) {
  
  if(is.list(homophily_term_50[[i]])) {
    
    coeff = homophily_term_50[[i]][["coefficients"]] %>%
      as.data.frame()
    
    coeff$coefficient = row.names(coeff) %>%
      str_replace(., "nodefactor.Species.", "Species - ") %>%
      str_replace(., "nodematch.Species.", "Match - ") %>%
      str_replace(., "edges", "Edges")
    
    coeff$Network = i
    coeff$Visit = unique(rodent_data_50$nodelist[[i]]$Visit)
    coeff$Landuse = unique(rodent_data_50$nodelist[[i]]$Landuse)
    coeff$`Observed M. natalensis` = rodent_data_50$nodelist[[i]] %>%
      filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
      nrow()
    coeff$`Observed all species` =  rodent_data_50$nodelist[[i]] %>%
      filter(Observed == TRUE) %>%
      nrow()
    coeff$`Unobserved M. natalensis` = rodent_data_50$nodelist[[i]] %>%
      filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
      nrow()
    coeff$N = nrow(rodent_data_50$nodelist[[i]])
    
    coefficients_50[[i]] <- coeff %>%
      tibble() %>%
      select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient = coefficient, Estimate, `Std. Error`)
    
  } else {
    
    coefficients_50[[i]] <- tibble(Network = i,
                                   Visit = unique(rodent_data_50$nodelist[[i]]$Visit),
                                   Landuse = unique(rodent_data_50$nodelist[[i]]$Landuse),
                                   `Observed M. natalensis` = rodent_data_50$nodelist[[i]] %>%
                                     filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
                                     nrow(),
                                   N = nrow(rodent_data_50$nodelist[[i]]),
                                   `Observed all species` =  rodent_data_50$nodelist[[i]] %>%
                                     filter(Observed == TRUE) %>%
                                     nrow(),
                                   `Unobserved M. natalensis` = rodent_data_50$nodelist[[i]] %>%
                                     filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
                                     nrow(),
                                   Coefficient = NA,
                                   Estimate = NA,
                                   `Std. Error` = NA) %>%
      tibble() %>%
      select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient, Estimate, `Std. Error`)
    
    
  }
  
}

coefficients_df_50 <- bind_rows(coefficients_50) %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village")),
         Estimate = case_when(Estimate == -Inf ~ as.numeric(NA),
                              TRUE ~ Estimate),
         Variance = `Std. Error`^2) %>%
  group_by(Landuse) %>%
  mutate(Network = paste0(Landuse, " ", Visit)) %>%
  filter(`Observed M. natalensis` > 1) %>% 
  filter(!Network %in% Network[is.na(Estimate)])

included_models_50 <- coefficients_df_50 %>%
  group_by(Network) %>% 
  filter(!Network %in% Network[is.na(Estimate) | is.infinite(Estimate)]) %>%
  distinct(Network)


# Speed up processing by removing networks
rm(rodent_data_50)
rm(final_model_50)
gc()
# Meta-analysis of edges --------------------------------------------------

edges_50 <- coefficients_df_50 %>%
  filter(Coefficient == "Edges") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .,
      slab = Network)

# for multi-level meta-analysis
#   rma.mv(yi = Estimate, V = Variance, data = .,
#   random = ~ 1 | Landuse/Visit,
#   slab = Network)

edges_rma_50 <- bind_rows(tibble(ES = edges_50$yi, SE = sqrt(edges_50$vi), Type = "Network", Network = factor(edges_50$data$Network), N = edges_50$data$`Observed M. natalensis`),
                          tibble(ES = edges_50$b[, 1], SE = edges_50$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

edges_combined_50 <- ggplot(edges_rma_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed (all species)",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")


# Meta-analysis of edges stratified by landuse ----------------------------

edges_ag_50 <- coefficients_df_50 %>%
  filter(Coefficient == "Edges" & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_vil_50 <- coefficients_df_50 %>%
  filter(Coefficient == "Edges" & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_rma_df_50 <- bind_rows(tibble(ES = edges_ag_50$yi, SE = sqrt(edges_ag_50$vi), Type = "Network",
                                    Weight = weights(edges_ag_50),
                                    Network = factor(edges_ag_50$data$Network),
                                    N = edges_ag_50$data$`Observed M. natalensis`),
                             tibble(ES = edges_ag_50$b[, 1], SE = edges_ag_50$se, Type = "Summary",
                                    Weight = 100,
                                    Network = factor("Summary - Agriculture"),
                                    N = as.numeric(NA)),
                             tibble(ES = edges_vil_50$yi, SE = sqrt(edges_vil_50$vi), Type = "Network",
                                    Weight = weights(edges_vil_50),
                                    Network = factor(edges_vil_50$data$Network),
                                    N = edges_vil_50$data$`Observed M. natalensis`),
                             tibble(ES = edges_vil_50$b[, 1], SE = edges_vil_50$se, Type = "Summary",
                                    Weight = 100,
                                    Network = factor("Summary - Village"),
                                    N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

edges_stratified_50 <- ggplot(edges_rma_df_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::breaks_pretty()) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  guides(size = "none") +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "Weight",
       title = "Odds of a tie being observed (all species)")


# Meta-analysis of ties based on species type -----------------------------

species_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Species")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma_50 <- bind_rows(tibble(ES = species_50$yi, SE = sqrt(species_50$vi), Type = "Network", Network = factor(species_50$data$Network), N = species_50$data$`Observed M. natalensis`),
                            tibble(ES = species_50$b[, 1], SE = species_50$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

species_combined_50 <- ggplot(species_rma_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed for M. natalensis")

species_ag_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_vil_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma_df_50 <- bind_rows(tibble(ES = species_ag_50$yi, SE = sqrt(species_ag_50$vi), Type = "Network",
                                      Weight = weights(species_ag_50),
                                      Network = factor(species_ag_50$data$Network),
                                      N = species_ag_50$data$`Observed M. natalensis`),
                               tibble(ES = species_ag_50$b[, 1], SE = species_ag_50$se, Type = "Summary", Weight = 100, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                               tibble(ES = species_vil_50$yi, SE = sqrt(species_vil_50$vi), Type = "Network",
                                      Weight = weights(species_vil_50),
                                      Network = factor(species_vil_50$data$Network), N = species_vil_50$data$`Observed M. natalensis`),
                               tibble(ES = species_vil_50$b[, 1], SE = species_vil_50$se, Type = "Summary", Weight = 100, Network = factor("Summary - Village"), N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

species_stratified_50 <- ggplot(species_rma_df_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::pretty_breaks()) +
  theme_bw() +
  guides(size = "none") +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "Weight",
       title = "Odds of a tie being observed for M. natalensis")

# Meta-analysis of intra-specific ties ------------------------------------

match_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Match")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma_50 <- bind_rows(tibble(ES = match_50$yi, SE = sqrt(match_50$vi), Type = "Network", 
                                 Network = factor(match_50$data %>%
                                                    drop_na(Estimate) %>%
                                                    pull(Network)), 
                                 N = match_50$data %>%
                                   drop_na(Estimate) %>%
                                   pull(`Observed M. natalensis`)),
                          tibble(ES = match_50$b[, 1], SE = match_50$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

match_combined_50 <- ggplot(match_rma_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N M. natalensis",
       title = "Odds of a tie between two M. natalensis being observed",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis") +
  coord_cartesian(xlim = c(0, 40))

match_ag_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_vil_50 <- coefficients_df_50 %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma_df_50 <- bind_rows(tibble(ES = match_ag_50$yi, SE = sqrt(match_ag_50$vi), Type = "Network", 
                                    Weight = weights(match_ag_50),
                                    Network = factor(match_ag_50$data %>%
                                                       drop_na(Estimate) %>%
                                                       pull(Network)),
                                    N = match_ag_50$data %>%
                                      drop_na(Estimate) %>%
                                      pull(`Observed M. natalensis`)),
                             tibble(ES = match_ag_50$b[, 1], SE = match_ag_50$se, Type = "Summary", Weight = 100, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                             tibble(ES = match_vil_50$yi, SE = sqrt(match_vil_50$vi), Type = "Network",
                                    Weight = weights(match_vil_50),
                                    Network = factor(match_vil_50$data %>%
                                                       drop_na(Estimate) %>%
                                                       pull(Network)),
                                    N = match_vil_50$data %>%
                                      drop_na(Estimate) %>%
                                      pull(`Observed M. natalensis`)),
                             tibble(ES = match_vil_50$b[, 1], SE = match_vil_50$se, Type = "Summary", Weight = 100, Network = factor("Summary - Village"), N = as.numeric(NA))) %>%
  mutate(Network = fct_inorder(Network))

match_stratified_50 <- ggplot(match_rma_df_50) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), colour = Weight, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.2, 0.5)) +
  scale_colour_viridis_c(direction = -1, breaks = scales::pretty_breaks()) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = element_blank(),
       title = "Odds of a tie between two M. natalensis being observed") +
  guides(size = "none") +
  coord_cartesian(xlim = c(0, 40))

p1_50 <- edges_stratified_50 + 
  theme(legend.position = "none") +
  labs(title = element_blank())
p2_50 <- species_stratified_50 +
  theme(legend.position = "none") +
  labs(title = element_blank())
p3_50 <- match_stratified_50 +
  theme(legend.position = "none") +
  labs(title = element_blank())
legend <- get_legend(match_stratified_50 +
                       theme(legend.position = "bottom"))

save_plot(plot = plot_grid(p1_50, p2_50, p3_50, legend, ncol = 1, labels = c("A", "B", "C", " "), rel_heights = c(1, 1, 1, 0.2)), filename = here("output", "Sensitivity_1_50m.png"), base_width = 7, base_height = 11)


combined_rma_50 <- list(edges_ag_50,
                        edges_vil_50,
                        species_ag_50,
                        species_vil_50,
                        match_ag_50,
                        match_vil_50)
write_rds(combined_rma_50, here("data", "Sensitivity_1_50m_results.rds"))
