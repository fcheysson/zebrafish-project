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