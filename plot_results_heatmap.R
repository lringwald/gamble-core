# Heatmap visualizations of the model results.
#  1) composition delta: driver x subclass (D/O/F), faceted by species -> the ALLOCATION drivers
#  2) totals FE: driver x species -> the LEVEL effects (mostly ~0: drivers shape split, not amount)
# Output: output/plots/{composition_delta_heatmap,totals_level_heatmap}.png
suppressMessages({ library(data.table); library(ggplot2) })
dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)

# ---- 1. composition delta heatmap ----
p <- fread("output/composition/subclass_allocation_parameters.csv")
p[, type := substr(subclass, 4, 4)]                                   # D / O / F
p[, type := factor(type, levels = c("D","O","F"), labels = c("Dairy","Meat","Follower"))]
# order drivers by overall effect strength (max-min across cells), nicest at top
ord <- p[, .(strength = max(delta) - min(delta)), by = driver][order(strength), driver]
p[, driver := factor(driver, levels = ord)]
lim <- as.numeric(quantile(abs(p$delta), 0.95))                       # clip color so outliers don't wash it out
p[, dclip := pmax(pmin(delta, lim), -lim)]

g1 <- ggplot(p, aes(type, driver, fill = dclip)) +
  geom_tile(color = "grey92") +
  geom_point(data = p[visible == TRUE], aes(type, driver), shape = 20, size = 0.6, color = "black") +
  facet_wrap(~species, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = expression(delta)) +
  labs(title = "Subtype allocation drivers (composition delta)",
       subtitle = "+ = pushes that subtype up; dot = effect is 'visible' (90% CI excludes 0)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold"), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave("output/plots/composition_delta_heatmap.png", g1, width = 8, height = 8, dpi = 150)
cat("wrote output/plots/composition_delta_heatmap.png\n")

# ---- 2. totals level heatmap (FE per species) ----
rd <- function(sp) { f <- list.files("output", sprintf("parameter_summary_.*%s_.*csv", sp), full.names = TRUE)
  d <- fread(f[which.max(file.mtime(f))]); d[, species := sp][, .(species, driver = Covariate, beta = Post_Mean)] }
tot <- rbind(rd("BOV"), rd("SGT"))[driver != "intercept"]
tot[, driver := factor(driver, levels = tot[species=="BOV"][order(abs(beta)), driver])]
tl <- as.numeric(quantile(abs(tot$beta), 0.95)); tot[, bclip := pmax(pmin(beta, tl), -tl)]
g2 <- ggplot(tot, aes(species, driver, fill = bclip)) +
  geom_tile(color = "grey92") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-tl, tl), name = expression(beta)) +
  labs(title = "Totals model: driver effect on the LEVEL (how much livestock)",
       subtitle = "near-uniformly ~0 -> drivers shape composition, not the amount", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
ggsave("output/plots/totals_level_heatmap.png", g2, width = 5, height = 9, dpi = 150)
cat("wrote output/plots/totals_level_heatmap.png\n")

# ---- 3. gamma heatmap: FINAL merged coefficients (beta_total + delta), averaged over countries ----
gp <- fread("output/composition/subclass_country_parameters.csv")
gp <- gp[driver != "intercept", .(gamma = mean(gamma)), by = .(species, subclass, driver)]   # mean over countries
gp[, type := factor(substr(subclass, 4, 4), levels = c("D","O","F"), labels = c("Dairy","Meat","Follower"))]
gord <- gp[, .(strength = max(gamma) - min(gamma)), by = driver][order(strength), driver]
gp[, driver := factor(driver, levels = gord)]
gl <- as.numeric(quantile(abs(gp$gamma), 0.95)); gp[, gclip := pmax(pmin(gamma, gl), -gl)]
g3 <- ggplot(gp, aes(type, driver, fill = gclip)) +
  geom_tile(color = "grey92") +
  facet_wrap(~species, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-gl, gl), name = expression(gamma)) +
  labs(title = "Final merged allocation coefficients (gamma = beta_total + delta)",
       subtitle = "country-averaged; subclass weight per grid unit prop. exp(sum gamma*X + offset), renormalized to NUTS3",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold"), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave("output/plots/gamma_heatmap.png", g3, width = 8, height = 9, dpi = 150)
cat("wrote output/plots/gamma_heatmap.png\n")
