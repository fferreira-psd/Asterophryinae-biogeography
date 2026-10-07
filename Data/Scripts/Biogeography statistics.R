######### BioGeoBEARS visualization pipeline ############
# 
# Inputs:
#   DECJ.Rdata          BioGeoBEARS DEC+j results
#   BEAST_tree.trees    BEAST consensus tree
#   elevation.csv       Optional tip-elevation data
#
# Outputs:
#   Fig.2.pdf            Circular ancestral-state phylogeny
#   Tree with bars.pdf   Circular phylogeny + elevation
#   Fig3A.pdf            Connectivity through time
#   Fig.3B.pdf           Net dispersal balance
#   Fig. 3C.pdf          Within-area cladogenesis
#
# Criteria:
#   - ancestral-state with marginal probability >= 0.66
#   - single-area states only
#   - dispersal events restricted to branches <= 5 Myr
#   - time bins = 5 Myr
#
# R version used: 4.3.3


# ------------------------------------------------------------
# Packages and files
# ------------------------------------------------------------

library(BioGeoBEARS)
library(ape)
library(plotrix)
library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)

# SET THIS TO YOUR WORKING DIRECTORY #
data_dir <- "."


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

threshold <- 0.66
max_branch_length <- 5
time_bin <- 5

regs <- c(
  "Bird's Head", "Northern terranes", "Central Range",
  "Southern New Guinea", "EPCT", "New Britain",
  "Wallacea", "Philippines"
)

region_cols <- c(
  "olivedrab3", "goldenrod3", "red3", "slateblue3",
  "dodgerblue2", "orangered", "black", "seagreen"
)
names(region_cols) <- regs

ambig_col <- "grey70"

# Regions retained for downstream analyses
valid_regions <- c(
  "Bird's Head", "Northern terranes", "Central Range",
  "Southern New Guinea", "EPCT", "Wallacea"
)
region_cols_plot <- region_cols[valid_regions]

# Time bins
time_labels <- c("0–5", "5–10", "10–15", "15–20", "20–25", "25–30")
time_midpoints <- setNames(seq(2.5, 27.5, 5), time_labels)


# ------------------------------------------------------------
# Plotting theme
# ------------------------------------------------------------

theme_biogeo <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 13),
      axis.title.x = element_text(size = 16, margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, margin = margin(r = 15)),
      axis.text = element_text(size = 9)
    )
}


# ------------------------------------------------------------
# Load BioGeoBEARS results and tree
# ------------------------------------------------------------

load(file.path(data_dir, "DECJ.Rdata"))
resDECj <- res

tree <- read.tree(file.path(data_dir, "BEAST_tree.trees"))
tree$tip.label <- gsub("'", "", tree$tip.label)
tree <- ladderize(tree, right = TRUE)


# ------------------------------------------------------------
# Ancestral-state reconstruction
# ------------------------------------------------------------

anc <- resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
states <- resDECj$inputs$states_list

best_state <- max.col(anc)
best_prob <- apply(anc, 1, max)

node_age <- c(
  rep(0, Ntip(tree)),
  branching.times(tree)
)

tip_nodes <- seq_len(Ntip(tree))
internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))

# Retain only supported, single-area states
node_state <- rep(NA_integer_, nrow(anc))

for (i in seq_len(nrow(anc))) {
  best <- best_state[i]
  if (best_prob[i] < threshold) next
  
  st <- states[[best]]
  if (length(st) == 1)
    node_state[i] <- st
}

# Node and branch colours
node_cols <- rep(ambig_col, length(node_state))
node_cols[!is.na(node_state)] <-
  region_cols[node_state[!is.na(node_state)] + 1]

edge_cols <- node_cols[tree$edge[, 2]]

# Pie-chart colours
pie_cols <- sapply(states, function(x)
  if (length(x) == 1) region_cols[x + 1] else ambig_col
)


# ------------------------------------------------------------
# Infer dispersal events
# ------------------------------------------------------------

events_list <- vector("list", nrow(tree$edge))

for (i in seq_len(nrow(tree$edge))) {
  
  parent <- tree$edge[i, 1]
  daughter <- tree$edge[i, 2]
  
  from <- node_state[parent]
  to <- node_state[daughter]
  
  if (any(is.na(c(from, to))) || from == to) next
  
  branch_length <- node_age[parent] - node_age[daughter]
  if (branch_length > max_branch_length) next
  
  events_list[[i]] <- data.frame(
    from = regs[from + 1],
    to = regs[to + 1],
    midpoint_age = mean(c(node_age[parent], node_age[daughter]))
  )
}

events <- bind_rows(events_list) %>%
  filter(from %in% valid_regions, to %in% valid_regions)


# ------------------------------------------------------------
# Fig. 2  Circular ancestral-state phylogeny
# ------------------------------------------------------------

pdf("Fig.2.pdf", width = 20, height = 28)

plot(
  tree,
  type = "fan",
  edge.color = edge_cols,
  edge.width = 5,
  show.tip.label = FALSE,
  open.angle = 16,
  rotate.tree = 197
)

axisPhylo(side = 1, cex = 0.7, lwd = 5)

nodelabels(
  node = internal_nodes,
  pie = anc[internal_nodes, ],
  piecol = pie_cols,
  cex = 0.3
)

legend(
  "bottomright",
  legend = c(regs, "Ambiguous"),
  fill = c(region_cols, ambig_col),
  border = NA,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# Optional elevation plot
# ------------------------------------------------------------

elevation_file <- file.path(data_dir, "elevation.csv")

if (file.exists(elevation_file)) {
  
  elevation <- read.csv(elevation_file, sep = ";")
  elevation <- elevation[match(tree$tip.label, elevation$label), ]
  
  stopifnot(all(elevation$label == tree$tip.label))
  
  pdf("Tree with bars.pdf", width = 25, height = 33)
  
  plot(
    tree,
    type = "fan",
    edge.color = edge_cols,
    edge.width = 4,
    show.tip.label = FALSE,
    open.angle = 16,
    rotate.tree = 197
  )
  
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  tip_xy <- cbind(
    pp$xx[tip_nodes],
    pp$yy[tip_nodes]
  )
  
  gap <- 0.7
  scale_length <- 3
  max_alt <- 3000
  
  tip_r <- mean(sqrt(
    pp$xx[tip_nodes]^2 + pp$yy[tip_nodes]^2
  ))
  
  bar_length <- elevation$alt / max_alt * scale_length
  
  bar_cols <- ifelse(
    elevation$alt < 1000, "grey90",
    ifelse(elevation$alt <= 2000, "grey70", "grey20")
  )
  
  for (i in seq_along(tree$tip.label)) {
    
    if (is.na(elevation$alt[i])) next
    
    theta <- atan2(tip_xy[i, 2], tip_xy[i, 1])
    r0 <- tip_r + gap
    r1 <- r0 + bar_length[i]
    
    segments(
      r0 * cos(theta), r0 * sin(theta),
      r1 * cos(theta), r1 * sin(theta),
      col = bar_cols[i],
      lwd = 15,
      lend = "butt"
    )
  }
  
  # Reference circles
  theta <- seq(0, 2 * pi, length.out = 1000)
  
  for (alt in c(1000, 2000, 3000)) {
    
    r <- tip_r + gap + scale_length * alt / max_alt
    
    lines(
      r * cos(theta),
      r * sin(theta),
      col = "grey80",
      lwd = 1
    )
  }
  
  nodelabels(
    node = internal_nodes,
    pie = anc[internal_nodes, ],
    piecol = pie_cols,
    cex = 0.4
  )
  
  legend(
    "bottomright",
    legend = c(regs, "Ambiguous"),
    fill = c(region_cols, ambig_col),
    border = NA,
    bty = "n"
  )
  
  dev.off()
}


# ------------------------------------------------------------
# Add time bins to dispersal events for Fig. 3
# ------------------------------------------------------------

events_binned <- events %>%
  mutate(
    time_bin = cut(
      midpoint_age,
      breaks = seq(0, 30, by = time_bin),
      include.lowest = TRUE,
      right = FALSE,
      labels = time_labels
    )
  )


# ------------------------------------------------------------
# Fig. 3A  Connectivity through time
# ------------------------------------------------------------

region_order <- valid_regions

pdf("Fig3A.pdf", width = 18, height = 8)

par(mfrow = c(1, 5), mar = c(1, 1, 2, 1))

# Original figure order: oldest to youngest
windows <- c("20–25", "15–20", "10–15", "5–10", "0–5")

for (w in windows) {
  
  dat <- events_binned %>%
    filter(time_bin == w) %>%
    count(from, to, name = "n")
  
  if (!nrow(dat)) next
  
  mat <- matrix(
    0,
    nrow = length(region_order),
    ncol = length(region_order),
    dimnames = list(region_order, region_order)
  )
  
  link_col <- matrix(
    NA,
    nrow = length(region_order),
    ncol = length(region_order),
    dimnames = list(region_order, region_order)
  )
  
  for (i in seq_len(nrow(dat))) {
    
    mat[dat$from[i], dat$to[i]] <- dat$n[i]
    
    link_col[dat$from[i], dat$to[i]] <-
      add_transparency(
        region_cols_plot[dat$from[i]],
        transparency = 0.5
      )
  }
  
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.degree = 5
  )
  
  chordDiagram(
    mat,
    order = region_order,
    grid.col = region_cols_plot,
    col = link_col,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "triangle",
    link.arr.length = 0.25,
    link.arr.width = 0.15,
    annotationTrack = "grid"
  )
  
  title(paste0(w, " Ma"), cex.main = 3)
}

dev.off()
circos.clear()


# ------------------------------------------------------------
# Fig. 3B — Net source/sink dispersal balance
# ------------------------------------------------------------

balance <- events_binned %>%
  count(time_bin, region = from, name = "exports") %>%
  full_join(
    events_binned %>%
      count(time_bin, region = to, name = "imports"),
    by = c("time_bin", "region")
  ) %>%
  mutate(
    exports = replace_na(exports, 0),
    imports = replace_na(imports, 0),
    balance = exports - imports,
    total_events = exports + imports,
    time = time_midpoints[as.character(time_bin)],
    region = factor(region, levels = valid_regions)
  )

g_balance <- ggplot(
  balance,
  aes(x = time, y = balance, fill = region)
) +
  geom_col(width = 4.4, linewidth = 0.2) +
  geom_text(
    aes(
      y = ifelse(balance >= 0, balance + 0.7, balance - 0.7),
      label = total_events
    ),
    size = 3
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey40",
    linewidth = 0.4
  ) +
  facet_grid(. ~ region, switch = "x") +
  scale_fill_manual(values = region_cols_plot) +
  scale_x_reverse(
    limits = c(30, 0),
    breaks = seq(30, 0, -5),
    labels = c("30", "", "20", "", "10", "", "0"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = seq(-6, 6, 2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_blank(),
    axis.title.x = element_text(size = 10, margin = margin(t = 15)),
    axis.title.y = element_text(size = 10, margin = margin(r = 15)),
    axis.text = element_text(size = 9)
  ) +
  labs(
    x = "Time (Ma)",
    y = "Net dispersal balance\n(sink < 0 < source)"
  )

ggsave(
  "Fig.3B.pdf",
  g_balance,
  width = 11.5,
  height = 3
)


# ------------------------------------------------------------
# Fig. 3C — Within-area cladogenesis by genus
# ------------------------------------------------------------

genus_nodes <- c(
  Xenorhina = 442,
  Oreophryne = 253,
  Hylophorbus = 384
)

count_cladogenesis <- function(genus_name, genus_node) {
  
  clade_nodes <- intersect(
    c(genus_node, getDescendants(tree, genus_node)),
    internal_nodes
  )
  
  clado <- lapply(clade_nodes, function(node) {
    
    daughters <- tree$edge[tree$edge[, 1] == node, 2]
    
    # Only bifurcating nodes are considered
    if (length(daughters) != 2) return(NULL)
    
    states <- node_state[c(node, daughters)]
    
    # Ignore ambiguous ancestral or descendant states
    if (any(is.na(states))) return(NULL)
    
    # Within-area cladogenesis
    if (!all(states == states[1])) return(NULL)
    
    region <- regs[states[1] + 1]
    
    if (!region %in% valid_regions) return(NULL)
    
    data.frame(
      genus = genus_name,
      region = region,
      time = node_age[node]
    )
  })
  
  clado <- bind_rows(clado)
  
  if (!nrow(clado)) {
    return(data.frame(
      genus = character(),
      region = character(),
      time = numeric(),
      time_bin = factor(character(), levels = time_labels),
      n = integer()
    ))
  }
  
  clado %>%
    mutate(
      time_bin = cut(
        time,
        breaks = seq(0, 30, by = time_bin),
        include.lowest = TRUE,
        right = FALSE,
        labels = time_labels
      )
    ) %>%
    count(genus, region, time_bin, name = "n")
}

clado_all <- bind_rows(
  lapply(
    names(genus_nodes),
    \(g) count_cladogenesis(g, genus_nodes[g])
  )
)

clado_all <- clado_all %>%
  mutate(
    time = time_midpoints[as.character(time_bin)],
    region = factor(region, levels = valid_regions),
    genus = factor(
      genus,
      levels = c("Xenorhina", "Oreophryne", "Hylophorbus")
    )
  )

### Plot ###

g_clado_genus <- ggplot(
  clado_all,
  aes(x = time, y = n, fill = region)
) +
  geom_col(width = 4.4, linewidth = 0.2) +
  facet_grid(genus ~ region) +
  scale_fill_manual(values = region_cols_plot) +
  geom_hline(
    yintercept = 0,
    colour = "grey40",
    linewidth = 0.4
  ) +
  scale_x_reverse(
    limits = c(30, 0),
    breaks = seq(30, 0, -5),
    labels = c("30", "", "20", "", "10", "", "0")
  ) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, 3)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_blank(),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    axis.title.x = element_text(
      size = 10,
      margin = margin(t = 15)
    ),
    axis.title.y = element_text(
      size = 10,
      margin = margin(r = 15)
    ),
    axis.text = element_text(size = 9)
  ) +
  labs(
    x = "Time (Ma)",
    y = "Within-area cladogenesis"
  )

ggsave(
  "Fig.3C.pdf",
  g_clado_genus,
  width = 11.5,
  height = 6
)
