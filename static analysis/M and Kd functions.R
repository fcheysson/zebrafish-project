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