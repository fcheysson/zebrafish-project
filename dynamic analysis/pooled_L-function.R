# Pooled L-function for the experimental fish

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

names = c("bibi", "mimi", "titi")
deltas = 1L:5L

hyperMC = tibble(fish = factor(character(0L), levels = names),
                 delta = integer(0L),
                 Leach = list(),
                 L = list())

for (fish in names) {
    
    display = paste0("LOADING FISH '", fish, "'")
    cat("\n", rep("=", nchar(display)), "\n", rep("=", nchar(display)), "\n", display, "\n\n", sep = "")
    
    source(paste0(fish, "/import_inhom.R"))
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    
    # Rescale to cell diameter
    # cells.diameter = plist %>% sapply(function(p) mean(nndist(p))) %>% mean()
    cells.diameter = plist %>% sapply(function(p) 2 * sqrt(area(Window(p)) / (pi * npoints(p)))) %>% mean()
    plist = plist %>% solapply(function(p) rescale(X=p, s=cells.diameter, unitname="cell diameter"))
    
    # Delta 0
    cellPatterns = solapply(plist, function(p) split(p)$MC)
    
    Leach = anylapply(cellPatterns, Lest, r = seq(0, 6, by=0.1), ratio = TRUE, correction = "best")
    Lpool = pool(Leach)
    
    hyperMC = rbind(hyperMC, 
                    tibble(fish = fish, 
                           delta = 0L,
                           Leach = list(Leach),
                           L = list(Lpool)))
    
    # Other deltas
    for (delta in deltas) {
        
        cells.delta = map_dfr(1:(8-delta), function(image) {
            cells %>% filter(cliche == image) %>% bind_rows(
                cells %>% filter(cliche == image + delta, type == "MC") %>% mutate(type = "MCfuture")) %>% 
                mutate(cliche = paste0(image, image + delta), type = as.factor(type))
        })
        
        cliches = unique(cells.delta$cliche)
        plist.delta = cliches %>% solapply(function(image) {
            cells.delta %>% filter(cliche == image) %>% {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = .$type)} %>% rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        Leach = anylapply(plist.delta, Lcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
        Lpool = pool(Leach)
        
        hyperMC = rbind(hyperMC, 
                        tibble(fish = fish,
                               delta = delta,
                               Leach = list(Leach),
                               L = list(Lpool)))
        
    }
    
}

hyperMC1 = hyperMC

## Run the simulation comparisons now

LpoolTib = map_dfr(0L:5L, function(delay) { 
    map_dfr(groups, function(mgroup) {
        hf = hyperMC %>% filter(delta == delay, group == mgroup) %>% select(-L) %>% unnest(Leach)
        Lpool = pool(as.anylist(hf$Leach))
        tibble(group = mgroup, delta = delay, L = list(Lpool))
    })
})

LpoolUnnest = LpoolTib %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

LpoolTib1 = map_dfr(0L:5L, function(delay) { 
    hf = hyperMC1 %>% filter(delta == delay) %>% select(-L) %>% unnest(Leach)
    Lpool = pool(as.anylist(hf$Leach))
    tibble(group = "fish", delta = delay, L = list(Lpool))
})

LpoolUnnest1 = LpoolTib1 %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

ggplot(LpoolUnnest %>% bind_rows(LpoolUnnest1), 
       aes(x=r, y=pooliso-r, ymin=loiso-r, ymax=hiiso-r, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$L(r) - r$")) +
    facet_wrap(~ delta) +
    theme_bw()

ggsave("comparison-with-without-LI-ribbon.eps", device=cairo_ps, width=14, height=8, units="in")

ggplot(hyperMC1 %>% select(-Leach) %>% unnest(L), 
       aes(x=r, y=pooliso-r, colour=fish, group=fish)) +
    geom_line() +
    ylab(TeX("$L(r) - r$")) +
    scale_color_manual(name = "Modèles", 
                       values = gg_color_hue(5)) +
    facet_wrap(~ delta) +
    theme_bw()

ggsave("comparison-with-without-LI-ribbon-individual.eps", device=cairo_ps, width=14, height=8, units="in")
