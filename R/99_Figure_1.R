source(here::here("R", "00_setup.R"))

# A map of Sierra Leone in Africa -----------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

africa <- world %>%
  filter(str_detect(continent, "Africa")) %>%
  select(admin) %>%
  mutate(fill = case_when(str_detect(admin, "Sierra Leone") ~ "black",
                          TRUE ~ "grey")) %>%
  ms_filter_islands(min_area = 1e10)

sl_bbox <- africa %>%
  filter(str_detect(admin, "Sierra Leone")) %>%
  st_bbox() %>%
  st_as_sfc()

africa_map <- ggplot() +
  geom_sf(data = africa, aes(fill = fill), colour = "black", size = 0.2) +
  geom_sf(data = sl_bbox, fill = NA, colour = "black", size = 1.5) +
  scale_fill_manual(values = c("grey", "white")) +
  guides(fill = "none") +
  theme_void()

# Load Sierra Leone administrative boundaries
sle_sf <- geodata::gadm(country = "SLE", level = 1, path = here("data", "geodata")) %>%
  st_as_sf()

# Define Eastern Province
eastern_province <- sle_sf %>%
  filter(NAME_1 == "Eastern")

# Load administrative level 2 for finer-scale mapping
sle_sf2 <- geodata::gadm(country = "SLE", level = 2, path = here("data", "geodata")) %>%
  st_as_sf()

poi <- tibble(name = c("Kenema", "Baiama", "Lalehun", "Lambayama", "Seilama", "Freetown"),
              type = c("city", "site", "site", "site", "site", "city"),
              lat = c(7.8762, 7.8375, 8.1974, 7.8506, 8.1223, 8.4844),
              lon = c(-11.1908, -11.2684, -11.0803, -11.1969, -11.1936, -13.2299)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


# Bounding box for entire Sierra Leone
zoom_sl <- st_bbox(sle_sf)

# Main Sierra Leone map highlighting Eastern Province
study_map <- ggplot() + 
  # Full Sierra Leone map (lighter grey for non-Eastern Province)
  geom_sf(data = sle_sf, fill = "grey90", colour = "white") +
  # Highlight Eastern Province (darker grey)
  geom_sf(data = eastern_province, fill = "grey70", colour = "black", size = 1) +  
  # Plot trap sites as black circles and Kenema/Freetown with a different shape
  geom_sf(data = poi, aes(shape = type), size = 3, colour = "black") +
  # Label sites and cities
  ggrepel::geom_text_repel(data = poi, aes(label = name, geometry = geometry), 
                           stat = "sf_coordinates", min.segment.length = 0, size = 3) +
  # Set shapes: "+" for Kenema & Freetown, black circles for trap sites
  scale_shape_manual(values = c("city" = 3, "site" = 19)) +
  # Adjust bounding box to show all of Sierra Leone
  coord_sf(xlim = c(zoom_sl[1], zoom_sl[3]), ylim = c(zoom_sl[2], zoom_sl[4])) +
  # Map styling
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  # North arrow and scale bar
  ggspatial::annotation_north_arrow(style = north_arrow_minimal()) +
  ggspatial::annotation_scale(location = "br") +
  guides(shape = "none") +
  labs(x = NULL, y = NULL)

# Combine the main map and inset map
sl_inset_map <- ggdraw() +
  draw_plot(study_map) +
  draw_plot(africa_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)

# Example network in Figure 1 ---------------------------------------------
library(scico)

rodent_network <- read_rds(here("data", "processed_data", "rodent_network_site_visit.rds"))

# Extract network 24
network_24 <- rodent_network[[24]]

# Convert the network into an edgelist and create a tibble for nodes and edges
a <- as.edgelist(network_24)
vnames <- tibble(node = attr(a, "vnames"))
obs_edges <- tibble(from = attr(a, "vnames")[a[, 1]],
                    to = attr(a, "vnames")[a[, 2]])

# Create a tbl_graph object with node attributes
gg_graph_24 <- tbl_graph(nodes = vnames, edges = obs_edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(Landuse = network_24 %v% "Landuse",
         Observed = network_24 %v% "Observed",
         Species = network_24 %v% "Species",
         Village = network_24 %v% "Village",
         Visit = network_24 %v% "Visit") %>%
  mutate(Species = factor(Species, levels = c(species_order_n, "Other")))

# Generate the plot for network 24
network_plot_24 <- ggraph(gg_graph_24, layout = "kk") +
  geom_edge_fan() +
  geom_node_point(aes(colour = Species, alpha = Observed)) +
  scale_colour_scico_d(palette = "roma", drop = FALSE) +
  scale_alpha_discrete(range = c(0.1, 1)) +
  theme_graph() +
  guides(alpha = "none",
         size = "none")


save_plot(plot = plot_grid(plotlist = list(sl_inset_map,
                                           network_plot_24),
                           labels = c("A)", "B)"),
                           ncol = 1),
          filename = here("output", "figures", "Figure_1_combined.png"),
          base_height = 10,
          base_width = 8)
