## Spatial analysis of aNSCs and aNPs
## Fig.1E,F, Fig.2D, Supp. Fig.3B-E, Supp. Fig.4

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)
library(clue)

zebra.dir = list.dirs(recursive=FALSE)[!grepl("./_results", list.dirs(recursive=FALSE))] # removes '_images' dir
imports.dir = list.dirs(zebra.dir)[!(list.dirs(zebra.dir) %in% zebra.dir)]
imports.path = paste0(imports.dir, "/import.R")
imports.dir = imports.dir %>% substring(2)

pval = vector("list", length(imports.dir))
diam = vector("list", length(imports.dir))

setEPS()

for (import in seq_along(imports.dir)) {
    
    dir.create(paste0("_images_alpha01/_eps", imports.dir[import]), recursive = TRUE)
    dir.create(paste0("_images_alpha05/_eps", imports.dir[import]), recursive = TRUE)

    import.display=paste0("Import #", import, ": ", imports.dir[import])
    cat("\n", rep("=", nchar(import.display)), "\n", import.display, "\n\n", sep = "")
    
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    source(imports.path[import])
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

    # Nice plot
    g = ggplot(cells, aes(x=x, y=y)) +
        geom_point(aes(colour=type), size=3) +
        geom_path(data=bind_rows(bdry, bdry[1,])) +
        scale_color_manual(name = "Cells", 
                           labels = c("qNSC", "aNSC_nodb", "aNSC_db", "aNSC_db2sg", "aNP"),
                           values = c("#C3D69B", "#7030A0", "red", "blue","#FF8000")) +
        coord_equal() +
        theme_grey()
    
    ggsave(filename=paste0("./_images_alpha01/_eps", imports.dir[import], "/pattern.eps"),
           plot = g, device="ps", width=12, height=8, units="in")
    ggsave(filename=paste0("./_images_alpha05/_eps", imports.dir[import], "/pattern.eps"),
           plot = g, device="ps", width=12, height=8, units="in")

    # Rescale to cell diameter
    # cells.diameter = mean(nndist(p))
    cells.diameter = 2 * sqrt(area(Window(p)) / (pi * npoints(p)))
    p = rescale(X=p, s=cells.diameter, unitname="cell diameter")
    cells.mindist = min(nndist(p))
    
    # "Dernier plus proche" voisin des types 3
    # Le voisin direct le plus éloigné
    # Nb de voisins directs des aNP = 4 en moyenne
    win.eroded = erosion(Window(p), 3)
    snn.quantiles = quantile(nndist(p, k=4)[marks(p) == "aNP" & inside.owin(p$x, p$y, win.eroded)], probs = .95)
    
    ## Permutation cellules "aNSC (w/ db2sg)" vs "aNSC (w/ db2sg)"
    source("aNSC_all vs aNSC_all.R")
    
    ## Permutation cellules "aNSC_nodb" vs "aNSC_nodb"
    source("aNSC_nodb vs aNSC_nodb.R")
    
    ## Permutation cellules "aNP" vs "aNP"
    source("aNP vs aNP.R")
    
    ## Permutation cellules "aNP" vs "aNSC w/ db2sg"
    source("aNP vs aNSC_all.R")
    
    ## Permutation cellules "aNP" vs "aNSC_nodb"
    source("aNP vs aNSC_nodb.R")
    
    ## Permutation cellules "aNSC_db2sg" vs "aNSC_nodb"
    source("aNSC_db2sg vs aNSC_nodb.R")
    
    ## Permutation cellules type2 vs type2 AND type 3 vs type 2
    source("M and Kd functions.R")
    
    pval[[import]] = tibble(
        fish = imports.dir[import],
        test = c("aNSC-aNSC", "aNSC_nodb-aNSC_nodb", "aNP-aNP", "aNP-aNSC", "aNP-aNSC_nodb", "aNSC_db2sg-aNSC_nodb"),
        pval = c(
            mad.test(EL2.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL2_nodb.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL3.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL32.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL32_nodb.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL2ab.simultaneous, rinterval = c(cells.mindist, 2))$p.value
        )
    )
    
    diam[[import]] = tibble(
        fish = imports.dir[import],
        diam = cells.diameter
    )
    
}

save(pval, file = "pval.RData")
save(diam, file = "cells_diameter.RData")

pdf("pval.pdf", height=6, width=5)
gridExtra::grid.table(
    pval %>% bind_rows() %>% separate(fish, c(NA, "fishType", "fish"), sep = "/") %>% 
        group_by(fishType, test) %>% add_tally() %>% 
        summarise(pval.pooled = 1 - pchisq(-2 * sum(log(pval)), df = 2 * n()))
)
dev.off()

pdf("pval_all.pdf", height=4, width=15)
gridExtra::grid.table(
    pval %>% bind_rows() %>% separate(fish, c(NA, "fishType", "fish"), sep = "/") %>% 
        spread(test, pval)
)
dev.off()  
