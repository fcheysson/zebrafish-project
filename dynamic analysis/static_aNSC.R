# Spatial (non temporal) analysis of aNSCs in movies
# Fig.4A,B

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gridExtra)
library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)

# GGPLOT COLOR PALETTE
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

names = c("bibi", "mimi", "titi")

n.sim1 = 499
n.sim2 = 500
n.sim = n.sim1 + n.sim2

n.test = length(names)
pval = tibble(fish = character(n.test), 
              L.simult2 = numeric(n.test),
              L.simult6 = numeric(n.test))
pval.it = 1

for (fish in names) {
    
    display = paste0("LOADING FISH '", fish, "'")
    cat("\n", rep("=", nchar(display)), "\n", rep("=", nchar(display)), "\n", display, "\n\n", sep = "")
    
    source(paste0(fish, "/import_inhom_beforediv.R"))
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    
    # Rescale to cell diameter
    # cells.diameter = plist %>% sapply(function(p) mean(nndist(p))) %>% mean()
    cells.diameter = plist %>% sapply(function(p) 2 * sqrt(area(Window(p)) / (pi * npoints(p)))) %>% mean()
    plist = plist %>% solapply(function(p) rescale(X=p, s=cells.diameter, unitname="cell diameter"))
    
    path.to.image = paste0("_images_static/aNSC/", fish, "/")
    dir.create(path.to.image, recursive = TRUE)
    pval$fish[pval.it] = fish
    
    ggplot(cells, aes(x=dim1, y=dim2)) +
        geom_point(aes(colour=type), size=3) +
        geom_path(aes(x=x, y=y), data=bind_rows(bdry, bdry[1,])) +
        facet_wrap(~ cliche, ncol=2) +
        # scale_colour_manual(values = c(gg_color_hue(3), "black")) +
        scale_colour_manual(values = c("#7030A0", "#00A0C6", "#C3D69B")) +
        coord_equal() +
        theme_grey()
    
    ggsave(filename=paste0(path.to.image, "patterns.eps"), device="ps", width=12, height=8, units="in")

    cliches = unique(cells$cliche)
    cells.tosample = cells
    cells.aNSC = cells %>% filter(type == "aNSC_beforediv")
    ppsim = anylapply(1:n.sim, function(sim) {
        permut = map_dfr(cliches, function(image) {
            cells.permut = cells.tosample %>% filter(cliche==image)
            n.aNSC = cells.aNSC %>% filter(cliche == image) %>% nrow()
            n.tot = cells.permut %>% nrow()
            cells.permut[ sample(x = n.tot, size = n.aNSC, replace=FALSE),] %>% 
                mutate(type="aNSC_beforediv")
        })
        
        solapply(cliches, function(image) {
            permut %>% filter(cliche == image) %>% 
            {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = as.factor(.$type))} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
    })
    
    plist.aNSC = solapply(cliches, function(image) {
        cells.aNSC %>% filter(cliche == image) %>% 
        {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = as.factor(.$type))} %>% 
            rescale(s=cells.diameter, unitname="cell diameter")
    })
    
    ##############
    ### L FUNCTION
    cat("L. ")
    
    Lsim = anylapply(1:n.sim, function(sim) {
        Leach = anylapply(ppsim[[sim]], Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
        pool(Leach)
    })
    
    Leach = lapply(plist.aNSC, Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
    L = pool(as.anylist(Leach))
    plot(L)
    
    # No function for envelopes of replicated point patterns
    # -> separate Lsim into the simulations for mean (n.sim2) and the simulations for deviation (n.sim1)
    # -> compute mean
    # -> compute maximum deviation
    Lmean = map_dfc((n.sim1+1):n.sim, function(sim) {
        Lsim[[sim]]$pooliso
    }) %>% apply(1, mean)
    
    # Go for bivariate testing
    Ldeviation = sapply(c(2.0, 6.0), function(rmax) {
        imax = 1 + 10 * rmax
        sapply(1:n.sim1, function(sim) {
            max( abs(Lsim[[sim]]$pooliso[1L:imax] - Lmean[1L:imax]) )
        })
    })
    
    Ldeviation.sorted = apply(Ldeviation, 2, sort, decreasing = TRUE)
    Ldeviation.observed = sapply(c(2.0, 6.0), function(rmax) {
        imax = 1 + 10 * rmax
        max( abs(L$pooliso[1L:imax] - Lmean[1L:imax]) )
    })
    pval$L.simult2[pval.it] = 1 - sum(Ldeviation.sorted[,1] < Ldeviation.observed[1]) / (n.sim1 + 1)
    pval$L.simult6[pval.it] = 1 - sum(Ldeviation.sorted[,2] < Ldeviation.observed[2]) / (n.sim1 + 1)
    
    EL.simultaneous = list(
        r = L$r,
        obs = L$pooliso,
        mmean = Lmean,
        lo = Lmean - Ldeviation.sorted[5, 2],
        hi = Lmean + Ldeviation.sorted[5, 2]
    ) %>% as_tibble()
    
    ggplot(EL.simultaneous, aes(x=r)) +
        geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
        geom_line(aes(y=mmean), col="red", lty="dashed") +
        geom_line(aes(y=obs)) +
        xlab(expression(paste(italic("r"), " (cell diameter)"))) +
        ylab(expression(italic(L["aNSC, aNSC"](r)))) +
        ggtitle("EL.simultaneous") +
        theme_bw() +
        theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
    
    ggsave(paste0(path.to.image, "ELsimultaneous01.eps"), device="ps", width=12, height=8, units="in")
    
    EL.simultaneous = list(
        r = L$r,
        obs = L$pooliso,
        mmean = Lmean,
        lo = Lmean - Ldeviation.sorted[25, 2],
        hi = Lmean + Ldeviation.sorted[25, 2]
    ) %>% as_tibble()
    
    ggplot(EL.simultaneous, aes(x=r)) +
        geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
        geom_line(aes(y=mmean), col="red", lty="dashed") +
        geom_line(aes(y=obs)) +
        xlab(expression(paste(italic("r"), " (cell diameter)"))) +
        ylab(expression(italic(L["aNSC, aNSC"](r)))) +
        ggtitle("EL.simultaneous") +
        theme_bw() +
        theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
    
    ggsave(paste0(path.to.image, "ELsimultaneous05.eps"), device="ps", width=12, height=8, units="in")
    
    rm(list=c("ppsim", "Lsim"))
    pval.it = pval.it + 1
    
}

pdf(paste0("pval_aNSC-aNSC.pdf"), height=2, width=3)
gridExtra::grid.table(
    bind_rows(pval, pval %>% 
        gather("Test", "pval", L.simult2:L.simult6) %>% 
        group_by(Test) %>% 
        summarise(sum.pval = sum(log(pval))) %>% 
        mutate(pooled = 1 - pchisq(-2*sum.pval, df = 2*3)) %>% 
        select(-sum.pval) %>% 
        spread(Test, pooled) %>% 
        mutate(fish = "pooled")
    ) %>% 
        mutate_at(vars(starts_with("L")), ~ format(round(., 3), nsmall = 3))
)
dev.off()
