# Studentized permutation test (Hahn, 2012) to compare the 2D interaction between neurons 
# with and without aNP-driven LI
# Supp. Fig.14

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gridExtra)
library(tidyverse)
library(spatstat)
library(latex2exp)

# GGPLOT COLOR PALETTE
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

dir.create("neurons")

# STATISTICAL ANALYSIS
setEPS()
n.sim1 = 499
n.sim2 = 500
n.sim = n.sim1 + n.sim2

matlabObject = R.matlab::readMat(paste0("data/Neurons_positions_20200426_withLI_f059.mat"))
matlabList = matlabObject[[1]]

neuronsWithLI1 = map_dfr(1L:18L, function(index) {
    matlabList[[index]] %>% 
        as_tibble() %>% 
        rename(x = V1, y = V2, z = V3) %>% 
        mutate(simulation = as.character(index)) %>% 
        select(simulation, everything())
})

matlabObject = R.matlab::readMat(paste0("data/Neurons_positions_20200426_withLI_f059_more_simulations.mat"))
matlabList = matlabObject[[1]]

neuronsWithLI2 = map_dfr(1L:16L, function(index) {
    matlabList[[index]] %>% 
        as_tibble() %>% 
        rename(x = V1, y = V2, z = V3) %>% 
        mutate(simulation = as.character(index)) %>% 
        select(simulation, everything())
}) %>% mutate(simulation = as.character(as.numeric(simulation) + 18L))

neuronsWithLI = bind_rows(neuronsWithLI1, neuronsWithLI2)

patternsWithLI = solapply(1L:34L, function(index) {
    cells = neuronsWithLI %>% filter(simulation == index)
    cellWindow = ripras(x = cells$x, y = cells$y, shape = "rectangle")
    p = ppp(x = cells$x, y = cells$y, window = cellWindow)
    return(p)
})

matlabObject = R.matlab::readMat(paste0("data/Neurons_positions_20200426_withoutLI_f1.mat"))
matlabList = matlabObject[[1]]

neuronsWithoutLI1 = map_dfr(1L:17L, function(index) {
    matlabList[[index]] %>% 
        as_tibble() %>% 
        rename(x = V1, y = V2, z = V3) %>% 
        mutate(simulation = as.character(index)) %>% 
        select(simulation, everything())
})

matlabObject = R.matlab::readMat(paste0("data/Neurons_positions_20200426_withoutLI_f1_more_simulations.mat"))
matlabList = matlabObject[[1]]

neuronsWithoutLI2 = map_dfr(1L:16L, function(index) {
    matlabList[[index]] %>% 
        as_tibble() %>% 
        rename(x = V1, y = V2, z = V3) %>% 
        mutate(simulation = as.character(index)) %>% 
        select(simulation, everything())
}) %>% mutate(simulation = as.character(as.numeric(simulation) + 17L))

neuronsWithoutLI = bind_rows(neuronsWithoutLI1, neuronsWithoutLI2)

patternsWithoutLI = solapply(1L:33L, function(index) {
    cells = neuronsWithoutLI %>% filter(simulation == index)
    cellWindow = ripras(x = cells$x, y = cells$y, shape = "rectangle")
    p = ppp(x = cells$x, y = cells$y, window = cellWindow)
    return(p)
})

neurons = hyperframe(
    patterns = c(patternsWithLI, patternsWithoutLI), 
    group = c(
        rep("withLI", length(patternsWithLI)), 
        rep("withoutLI", length(patternsWithoutLI))
        )
    )
neurons$L = with(neurons, Lest(patterns, rmax = 20, ratio = TRUE, correction = "best"))
Lpool = anylapply(split(neurons$L, neurons$group), pool)
plot(Lpool, cbind(pooliso, pooltheo, loiso, hiiso) - r ~ r, shade=c("loiso", "hiiso"), legend = FALSE, equal.scales = TRUE)

LpoolTib = map_dfr(names(Lpool), function(mgroup) {
    tibble(group = mgroup, L = list(Lpool[[mgroup]]))
})

LpoolUnnest = LpoolTib %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

ggplot(LpoolUnnest, 
       aes(x=r, y=pooliso-r, ymin=loiso-r, ymax=hiiso-r, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$L(r) - r$")) +
    theme_bw()
ggsave("./neurons/Neurons_positions_20200426_2d_both.eps", device=cairo_ps, width=14, height=8, units="in")

Leach = anylapply(split(neurons$L, neurons$group), collapse.fv, same="theo", different="iso")
plot(Leach, . - r ~ r, legend=FALSE, equal.scales = TRUE)

test <- studpermu.test(neurons, patterns ~ group, summaryfunction = Lest, 
                       rinterval = c(0, 20), correction = "best", use.Tbar = TRUE)
test
plot(test, . - r ~ r)

postscript(file="./neurons/Neurons_positions_20200329_2d_pool.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(Lpool, cbind(pooliso, pooltheo, loiso, hiiso) - r ~ r, shade=c("loiso", "hiiso"), legend = FALSE, equal.scales = TRUE)
dev.off()

postscript(file="./neurons/Neurons_positions_20200329_2d.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(test, . - r ~ r)
dev.off()

p1 = ggplot(neuronsWithLI %>% filter(as.numeric(simulation) <= 3), aes(x=x, y=y)) +
    geom_point(colour = "dodgerblue1") +
    facet_wrap(~ simulation, ncol = 1) +
    coord_equal() +
    theme_grey()

p2 = ggplot(neuronsWithoutLI %>% filter(as.numeric(simulation) <= 3), aes(x=x, y=y)) +
    geom_point(colour = "firebrick2") +
    facet_wrap(~ simulation, ncol = 1) +
    coord_equal() +
    theme_grey()

cowplot::plot_grid(p1, p2)
