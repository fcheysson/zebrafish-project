# Studentized permutation test (Hahn, 2012) to compare the 3D interaction between neurons 
# with and without aNP-driven LI
# Fig.4G

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

patternsWithLI = anylapply(1L:34L, function(index) {
    cells = neuronsWithLI %>% filter(simulation == index)
    box = box3(range(cells$x), range(cells$y), range(cells$z), unitname=c("micron", "microns"))
    p = pp3(x = cells$x, y = cells$y, z = cells$z, box)
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

patternsWithoutLI = anylapply(1L:33L, function(index) {
    cells = neuronsWithoutLI %>% filter(simulation == index)
    box = box3(range(cells$x), range(cells$y), range(cells$z), unitname=c("micron", "microns"))
    p = pp3(x = cells$x, y = cells$y, z = cells$z, box)
    return(p)
})

neurons = hyperframe(
    patterns = c(patternsWithLI, patternsWithoutLI), 
    group = c(
        rep("withLI", length(patternsWithLI)), 
        rep("withoutLI", length(patternsWithoutLI))
    )
)
neurons$K = with(neurons, K3est(patterns, rmax = 60, ratio = TRUE, correction = "best"))
Kpool = anylapply(split(neurons$K, neurons$group), pool)
plot(Kpool, cbind(pooliso, pooltheo, loiso, hiiso) - 4/3*pi*r^3 ~ r, shade=c("loiso", "hiiso"), legend = FALSE, equal.scales = TRUE)

Keach = anylapply(split(neurons$K, neurons$group), collapse.fv, different="iso")
plot(Keach, . - 4/3*pi*r^3 ~ r, legend=FALSE, equal.scales = TRUE)

sptest = function(hf) {
    Keach = anylapply(split(hf$K, hf$group), collapse.fv, different = "iso")
    Ks = anylapply(Keach, function(fv) fv %>% as.matrix() %>% {.[,-1]} )
    Kmean = anylapply(Ks, function(mat) apply(mat, 1, mean))
    Kvar = anylapply(Ks, function(mat) apply(mat, 1, var))
    T = (Keach[[1]]$r[2] - Keach[[1]]$r[1]) * 
        sum( (Kmean[[1]] - Kmean[[2]])^2 / (Kvar[[1]] / ncol(Ks[[1]]) + Kvar[[2]] / ncol(Ks[[2]])) , na.rm = TRUE)
    return(T)
}

Ktib = map_dfr(1:nrow(neurons), function(index) {
    neurons$K[[index]] %>% as_tibble() %>% mutate(sim = index, group = neurons$group[index])
})

ggplot(Ktib %>% filter(r <= 30), aes(x=r, y=iso-4/3*pi*r^3, colour=group, group=sim)) +
    geom_line() +
    ylab(TeX("$K_3(r) - 4/3\\pi r^3$")) +
    theme_bw()

KpoolTib = map_dfr(names(Kpool), function(mgroup) {
    tibble(group = mgroup, K = list(Kpool[[mgroup]]))
})

KpoolUnnest = KpoolTib %>% mutate(K = lapply(K, as.data.frame)) %>% unnest()

ggplot(KpoolUnnest %>% filter(r <= 30), 
       aes(x=r, y=pooliso-4/3*pi*r^3, ymin=loiso-4/3*pi*r^3, ymax=hiiso-4/3*pi*r^3, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$K_3(r) - 4/3\\pi r^3$")) +
    theme_bw()

nperm = 999
Tvalues = c()
for (k in 1:nperm) {
    sampled = neurons
    sampled$K = neurons$K[sample(nrow(neurons), nrow(neurons))]
    Tvalues[k] = sptest(sampled)
}

1 - sum(Tvalues < sptest(neurons)) / 1000

postscript(file="./neurons/Neurons_positions_20200426_3d_pool.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(Kpool, cbind(pooliso, pooltheo, loiso, hiiso) - 4/3*pi*r^3 ~ r, shade=c("loiso", "hiiso"), legend = FALSE, equal.scales = TRUE)
dev.off()

ggplot(Ktib %>% filter(r <= 30), aes(x=r, y=iso-4/3*pi*r^3, colour=group, group=sim)) +
    geom_line() +
    ylab(TeX("$K_3(r) - 4/3\\pi r^3$")) +
    theme_bw()
ggsave("./neurons/Neurons_positions_20200426_3d.eps", device="ps", width=14, height=8, units="in")

ggplot(KpoolUnnest %>% filter(r <= 30), 
       aes(x=r, y=pooliso-4/3*pi*r^3, ymin=loiso-4/3*pi*r^3, ymax=hiiso-4/3*pi*r^3, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$K_3(r) - 4/3\\pi r^3$")) +
    theme_bw()
ggsave("./neurons/Neurons_positions_20200426_3d_both.eps", device=cairo_ps, width=14, height=8, units="in")

ggplot(KpoolUnnest %>% filter(r <= 30), 
       aes(x=r, y=sqrt(pooliso)-sqrt(4/3*pi*r^3), 
           ymin=sqrt(loiso)-sqrt(4/3*pi*r^3), 
           ymax=sqrt(hiiso)-sqrt(4/3*pi*r^3), 
           fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$\\sqrt{K_3(r)} - \\sqrt{4/3\\pi r^3}$")) +
    theme_bw()
ggsave("./neurons/Neurons_positions_20200426_3d_both_sqrt.eps", device=cairo_ps, width=14, height=8, units="in")
