## Permutation cellules "aNP" vs "aNSC_nodb"
## K_i,j is 1/int_j times the expected number of points of type j within a distance r of a typical point in i
## (chapter 31.3.1+, p.190+ of Analysing spatial point patterns in R, 2010)
n.sim = 999
cells12 = cells %>% filter(type %in% c("qNSC", "aNSC_nodb"))
cells23 = cells %>% filter(type %in% c("aNSC_nodb", "aNP"))
ppsim = list()

for (k in 1:n.sim) {
    permut = cells12[ sample(x = n1 + n2a, size = n2a, replace = FALSE), ] %>% 
        # rbinom(1, size = n1 + n2, prob = n2 / (n1 + n2))
        mutate(type = "aNSC_nodb") %>% 
        bind_rows(cells3) %>% 
        mutate(type = as.factor(type))
    ppsim[[k]] = ppp(permut$x, permut$y, window=cells.owin, marks=permut$type) %>% rescale(s=cells.diameter, unitname="cell diameter")
}

p23 = ppp(cells23$x, cells23$y, window=cells.owin, marks=cells23$type) %>% rescale(s=cells.diameter, unitname="cell diameter")

EL32_nodb.simultaneous = envelope(p23, Lcross, i="aNP", j="aNSC_nodb", simulate=ppsim, nsim=499, nsim2=500, nrank=5, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/EL32_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(EL32_nodb.simultaneous)
dev.off()

EL32_nodb.simultaneous2 = envelope(EL32_nodb.simultaneous, nsim=499, nsim2=500, nrank=25, global=TRUE, ginterval=c(cells.mindist, 5), funargs=list(rmax = 5, correction = "best"))
postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/EL32_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(EL32_nodb.simultaneous2)
dev.off()

Eg32_nodb.pointwise = envelope(p23, pcfcross, i="aNP", j="aNSC_nodb", simulate=ppsim, nsim=999, nrank=5, global=FALSE, funargs=list(rmax = 5, correction = "best"), savefuns=TRUE)
postscript(file=paste0("./_images_alpha01/_eps", imports.dir[import], "/Eg32pointwise_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(Eg32_nodb.pointwise, xlim=c(cells.mindist, 5), ylim=c(0,2))
abline(v=snn.quantiles, col="seagreen")
dev.off()

Eg32_nodb.pointwise2 = envelope(Eg32_nodb.pointwise, nsim=999, nrank=25, global=FALSE, funargs=list(rmax = 5, correction = "best"))
postscript(file=paste0("./_images_alpha05/_eps", imports.dir[import], "/Eg32pointwise_nodb.eps"), horizontal = FALSE, onefile = FALSE, paper = "special", width=14, height=8)
plot(Eg32_nodb.pointwise2, xlim=c(cells.mindist, 5), ylim=c(0,2))
abline(v=snn.quantiles, col="seagreen")
dev.off()