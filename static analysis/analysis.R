## Spatial analysis of aNSCs and aNPs
## Fig.1E,F, Fig.2D, Supp. Fig.3B-E, Supp. Fig.4

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)
library(clue)

zebra.dir = list.dirs(recursive=FALSE)[!grepl("./_images", list.dirs(recursive=FALSE))] # removes '_images' dir
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
    n.sim = 999
    cells12 = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg"))
    cells2 = cells %>% filter(type %in% c("aNSC_nodb", "aNSC_db2sg"))
    ppsim = list()
    
    for (k in 1:n.sim) {
        permut = cells12[ sample(x = n1 + n2, size = n2, replace = FALSE), ]
                                                     # rbinom(1, size = n1 + n2, prob = n2 / (n1 + n2))
        ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    }
    
    p2 = ppp(cells2$x, cells2$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    
    EL2.simultaneous = envelope(p2, Lest, simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL2.simultaneous)
    dev.off()
    
    EL2.simultaneous2 = envelope(EL2.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL2.simultaneous2)
    dev.off()
    
    Eg2.pointwise = envelope(p2, pcf, simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg2pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg2.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
    dev.off()
    
    Eg2.pointwise2 = envelope(Eg2.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg2pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg2.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
    dev.off()
    
    ## Permutation cellules "aNP" vs "aNP"
    n.sim = 999
    ppsim = list()
    
    for (k in 1:n.sim) {
        permut = cells[ sample(x = n1 + n2 + n3, size = n3, replace = FALSE), ]
                                                        # rbinom(1, size = n1 + n2 + n3, prob = n3 / (n1 + n2 + n3))
        ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    }
    
    p3 = ppp(cells3$x, cells3$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    
    EL3.simultaneous = envelope(p3, Lest, simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL3.simultaneous)
    dev.off()
    
    EL3.simultaneous2 = envelope(EL3.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL3.simultaneous2)
    dev.off()
    
    Eg3.pointwise = envelope(p3, pcf, simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg3pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg3.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
    dev.off()
    
    Eg3.pointwise2 = envelope(Eg3.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg3pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg3.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
    dev.off()
    
    ## Permutation cellules "aNP" vs "aNSC w/ db2sg"
    ## K_i,j is 1/int_j times the expected number of points of type j within a distance r of a typical point in i
    ## (chapter 31.3.1+, p.190+ of Analysing spatial point patterns in R, 2010)
    n.sim = 999
    cells12 = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg"))
    cells23 = cells %>% filter(type %in% c("aNSC_nodb", "aNSC_db2sg", "aNP")) %>% 
        mutate(type = factor(ifelse(type == "aNP", "aNP", "aNSC")))
    ppsim = list()
    
    for (k in 1:n.sim) {
        permut = cells12[ sample(x = n1 + n2, size = n2, replace = FALSE), ] %>% 
                                                     # rbinom(1, size = n1 + n2, prob = n2 / (n1 + n2))
            mutate(type = "aNSC") %>% 
            bind_rows(cells3) %>% 
            mutate(type = as.factor(type))
        ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin, marks=permut$type) %>% rescale(s=cells.diameter, unitname="cell diameter")
    }
    
    p23 = ppp(cells23$x, cells23$y, window=cells.owin, marks=cells23$type) %>% rescale(s=cells.diameter, unitname="cell diameter")
    
    EL32.simultaneous = envelope(p23, Lcross, i="aNP", j="aNSC", simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL32.simultaneous)
    dev.off()
    
    EL32.simultaneous2 = envelope(EL32.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL32.simultaneous2)
    dev.off()
    
    Eg32.pointwise = envelope(p23, pcfcross, i="aNP", j="aNSC", simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg32pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg32.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    Eg32.pointwise2 = envelope(Eg32.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg32pointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg32.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    ## Permutation cellules "aNSC_db2sg" vs "aNSC_nodb"
    ## K_i,j is 1/int_j times the expected number of points of type j within a distance r of a typical point in i
    ## (chapter 31.3.1+, p.190+ of Analysing spatial point patterns in R, 2010)
    n.sim = 999
    cells12a = cells %>% filter(type %in% c("qNSC", "aNSC_nodb"))
    cells2ab = cells %>% filter(type %in% c("aNSC_nodb", "aNSC_db2sg"))
    ppsim = list()
    
    for (k in 1:n.sim) {
        permut = cells12a[ sample(x = n1 + n2a, size = n2a, replace = FALSE), ] %>% 
                                                       # rbinom(1, size = n1 + n2a, prob = n2a / (n1 + n2a))
            mutate(type = "aNSC_nodb") %>% 
            bind_rows(cellsB) %>% 
            mutate(type = as.factor(type))
        ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin, marks=permut$type) %>% rescale(s=cells.diameter, unitname="cell diameter")
    }
    
    p2ab = ppp(cells2ab$x, cells2ab$y, window=cells.owin, marks=cells2ab$type) %>% rescale(s=cells.diameter, unitname="cell diameter")
    
    EL2ab.simultaneous = envelope(p2ab, Lcross, i="aNSC_db2sg", j="aNSC_nodb", simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL2ab.simultaneous)
    dev.off()
    
    EL2ab.simultaneous2 = envelope(EL2ab.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(EL2ab.simultaneous2)
    dev.off()
    
    Eg2ab.pointwise = envelope(p2ab, pcfcross, i="aNSC_db2sg", j="aNSC_nodb", simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg2abpointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg2ab.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    Eg2ab.pointwise2 = envelope(Eg2ab.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg2abpointwise.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(Eg2ab.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    ## Permutation cellules type2 vs type2 AND type 3 vs type 2
    ## M(r) = proportion de points d'intérêt dans un voisinage à celle observée sur l'ensemble de la fenêtre
    mcells = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg", "aNP")) %>% 
        mutate(type = factor(ifelse(type == "aNP", "aNP", ifelse(type == "qNSC", "qNSC", "aNSC")))) %>% 
        rename(PointType=type) %>% select(-z) %>% as.data.frame()
    mp = wmppp(mcells, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    attr(mp, "class") = c("wmppp", "ppp")
    
    cells12 = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg")) %>% 
        mutate(type = factor(ifelse(type == "qNSC", "qNSC", "aNSC")))
    mcells12 = cells12 %>% rename(PointType=type) %>% select(-z) %>% as.data.frame()
    mp12 = wmppp(mcells12, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    attr(mp12, "class") = c("wmppp", "ppp")
    
    cells12ab = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg"))
    mcells12ab = cells12ab %>% rename(PointType=type) %>% select(-z) %>% as.data.frame()
    mp12ab = wmppp(mcells12ab, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
    attr(mp12ab, "class") = c("wmppp", "ppp")
    
    ME2.simultaneous = MEnvelope(mp12, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType="aNSC", Global=TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/ME2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME2.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE2.simultaneous = KdEnvelope(mp12, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType = "aNSC", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/KdE2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE2.simultaneous)
    dev.off()
    
    ME3.simultaneous = MEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType="aNP", Global=TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/ME3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME3.simultaneous, ylim=c(0, max(ME3.simultaneous$obs, na.rm=TRUE)))
    dev.off()
    
    KdE3.simultaneous = KdEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType = "aNP", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/KdE3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE3.simultaneous)
    dev.off()
    
    ME32.simultaneous = MEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType="aNP", NeighborType="aNSC", Global=TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/ME32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME32.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE32.simultaneous = KdEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType = "aNP", NeighborType = "aNSC", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/KdE32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE32.simultaneous)
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    ME2ab.simultaneous = MEnvelope(mp12ab, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType="aNSC_db2sg", NeighborType="aNSC_nodb", Global=TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/ME2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME2ab.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE2ab.simultaneous = KdEnvelope(mp12ab, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.01, ReferenceType="aNSC_db2sg", NeighborType="aNSC_nodb", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/KdE2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE2ab.simultaneous)
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    ME2.simultaneous = MEnvelope(mp12, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType="aNSC", Global=TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/ME2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME2.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE2.simultaneous = KdEnvelope(mp12, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType = "aNSC", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/KdE2.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE2.simultaneous)
    dev.off()
    
    ME3.simultaneous = MEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType="aNP", Global=TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/ME3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME3.simultaneous, ylim=c(0, max(ME3.simultaneous$obs, na.rm=TRUE)))
    dev.off()
    
    KdE3.simultaneous = KdEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType = "aNP", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="RandomLocation")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/KdE3.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE3.simultaneous)
    dev.off()
    
    ME32.simultaneous = MEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType="aNP", NeighborType="aNSC", Global=TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/ME32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME32.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE32.simultaneous = KdEnvelope(mp, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType = "aNP", NeighborType = "aNSC", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/KdE32.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE32.simultaneous)
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    ME2ab.simultaneous = MEnvelope(mp12ab, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType="aNSC_db2sg", NeighborType="aNSC_nodb", Global=TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/ME2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(ME2ab.simultaneous, ylim=c(0, 2))
    dev.off()
    
    KdE2ab.simultaneous = KdEnvelope(mp12ab, r=seq(0, 5, by=.1), NumberOfSimulations = 500, Alpha=0.05, ReferenceType="aNSC_db2sg", NeighborType="aNSC_nodb", Original = FALSE, StartFromMinR = TRUE, Global = TRUE, SimulationType="PopulationIndependence")
    postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/KdE2ab.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
    plot(KdE2ab.simultaneous)
    abline(v=snn.quantiles, col="seagreen")
    dev.off()
    
    pval[[import]] = tibble(
        fish = imports.dir[import],
        test = c("aNSC-aNSC", "aNP-aNP", "aNP-aNSC", "aNSC_db2sg-aNSC_nodb"),
        pval = c(
            mad.test(EL2.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL3.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
            mad.test(EL32.simultaneous, rinterval = c(cells.mindist, 2))$p.value,
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

pdf("pval.pdf", height=4, width=5)
gridExtra::grid.table(
    pval %>% bind_rows() %>% separate(fish, c(NA, "fishType", "fish"), sep = "/") %>% 
        group_by(fishType, test) %>% add_tally() %>% 
        summarise(pval.pooled = 1 - pchisq(-2 * sum(log(pval)), df = 2 * n()))
)
dev.off()
    