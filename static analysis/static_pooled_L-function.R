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

zebra.dir = list.dirs(recursive=FALSE)[!grepl("./_images", list.dirs(recursive=FALSE))] # removes '_images' dir
imports.dir = list.dirs(zebra.dir)[!(list.dirs(zebra.dir) %in% zebra.dir)]
imports.path = paste0(imports.dir, "/import.R")
imports.dir = imports.dir %>% substring(2)

hyperMC = tibble(fish = factor(character(0L), levels = imports.dir),
                 type = character(0L),
                 L = list())

for (import in seq_along(imports.dir)) {
    
    import.display=paste0("Import #", import, ": ", imports.dir[import])
    cat("\n", rep("=", nchar(import.display)), "\n", import.display, "\n\n", sep = "")
    
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    source(imports.path[import])
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    
    # Rescale to cell diameter
    # cells.diameter = mean(nndist(p))
    cells.diameter = 2 * sqrt(area(Window(p)) / (pi * npoints(p)))
    p = rescale(X=p, s=cells.diameter, unitname="cell diameter")
    cells.mindist = min(nndist(p))
    
    ## If no_db
    # p = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNP")) %>% 
    #     mutate(type = factor(ifelse(type == "aNP", "aNP", ifelse(type == "qNSC", "qNSC", "aNSC")))) %>% 
    #     {ppp(.$x, .$y, window=cells.owin, marks=.$type)} %>% 
    #     rescale(s=cells.diameter, unitname="cell diameter")
    
    ## Permutation cellules type 2
    L = Lest(split(p)$aNSC, r = seq(0, 6, by=0.1), ratio = TRUE, correction = "best")
    
    hyperMC = rbind(hyperMC, 
                    tibble(fish = imports.dir[import],
                           type = "22",
                           L = list(L)))
    
    ## Permutation cellules type 3
    L = Lest(split(p)$aNP, r = seq(0, 6, by=0.1), ratio = TRUE, correction = "best")
    
    hyperMC = rbind(hyperMC, 
                    tibble(fish = imports.dir[import],
                           type = "33",
                           L = list(L)))
    
    ## Permutation cellules type 32
    L = Lcross(p, i="aNP", j="aNSC", r = seq(0, 6, by=0.1), ratio = TRUE, correction = "best")
    
    hyperMC = rbind(hyperMC, 
                    tibble(fish = imports.dir[import],
                           type = "32",
                           L = list(L)))
    
}

hyperMC1 = hyperMC %>% separate(fish, c(NA, "group", "fish"), sep="/")

## Run the simulation comparisons now

LpoolTib = map_dfr(c("22", "33", "32"), function(test) { 
    map_dfr(groups, function(mgroup) {
        hf = hyperMC %>% filter(type == test, group == mgroup) %>% select(-L) %>% unnest(Leach)
        Lpool = pool(as.anylist(hf$Leach))
        tibble(group = mgroup, type = test, L = list(Lpool))
    })
})

LpoolUnnest = LpoolTib %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

LpoolTib1 = map_dfr(c("22", "33", "32"), function(test) { 
    map_dfr(c("DMSO24h", "LY24h", "WT40x"), function(mgroup) {
        hf = hyperMC1 %>% filter(type == test, group == mgroup)
        Lpool = pool(as.anylist(hf$L))
        tibble(group = mgroup, type = test, L = list(Lpool))
    })
})

LpoolUnnest1 = LpoolTib1 %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

col2A = gg_color_hue(3)[1]
col2B = gg_color_hue(3)[3]
colG = gg_color_hue(3)[2]

ggplot(LpoolUnnest %>% bind_rows(LpoolUnnest1) %>% filter(group %in% c("2A", "WT40x", "2B")), 
       aes(x=r, y=pooliso-r, ymin=loiso-r, ymax=hiiso-r, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$L(r) - r$")) +
    scale_colour_manual(values = c(col2A, col2B, "grey50")) +#gg_color_hue(3)[c(1, 3, 2)]) +
    scale_fill_manual(values = c(col2A, col2B, "grey50")) +#gg_color_hue(3)[c(1, 3, 2)]) +
    facet_wrap(~ type) +
    theme_bw()

ggsave("comparison-with-without-LI-ribbon-2A-WT40x-2B.eps", device=cairo_ps, width=14, height=8, units="in")
