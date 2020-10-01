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