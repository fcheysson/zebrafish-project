## Permutation cellules "aNSC_nodb" vs "aNSC_nodb"
n.sim = 999
cells12 = cells %>% filter(type %in% c("qNSC", "aNSC_nodb"))
cells2 = cells %>% filter(type %in% c("aNSC_nodb"))
ppsim = list()

for (k in 1:n.sim) {
    permut = cells12[ sample(x = n1 + n2a, size = n2a, replace = FALSE), ]
    # rbinom(1, size = n1 + n2, prob = n2 / (n1 + n2))
    ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")
}

p2 = ppp(cells2$x, cells2$y, window=cells.owin) %>% rescale(s=cells.diameter, unitname="cell diameter")

EL2_nodb.simultaneous = envelope(p2, Lest, simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL2_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(EL2_nodb.simultaneous)
dev.off()

EL2_nodb.simultaneous2 = envelope(EL2_nodb.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL2_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(EL2_nodb.simultaneous2)
dev.off()

Eg2_nodb.pointwise = envelope(p2, pcf, simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg2pointwise_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(Eg2_nodb.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
dev.off()

Eg2_nodb.pointwise2 = envelope(Eg2_nodb.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg2pointwise_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(Eg2_nodb.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
dev.off()