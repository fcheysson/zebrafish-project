# Spatio-temporal analysis of MCs in movies
# Fig.3E-E", Supp. Fig.8, Supp. Fig.9

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
deltas = 1L:5L

n.test = length(names) * length(deltas)
pval = tibble(fish = character(n.test), 
              delta = integer(n.test), 
              L.simult = numeric(n.test),
              L.pointwise1 = numeric(n.test),
              L.pointwise2 = numeric(n.test),
              L.pointwise3 = numeric(n.test),
              L.pointwise4 = numeric(n.test),
              L.pointwise5 = numeric(n.test),
              L.pointwise6 = numeric(n.test),
              L.dclf1 = numeric(n.test),
              L.dclf2 = numeric(n.test),
              L.dclf3 = numeric(n.test),
              L.dclf4 = numeric(n.test),
              L.dclf5 = numeric(n.test),
              L.dclf6 = numeric(n.test),
              M.simult = numeric(n.test),
              M.pointwise1 = numeric(n.test),
              M.pointwise2 = numeric(n.test),
              M.pointwise3 = numeric(n.test),
              M.pointwise4 = numeric(n.test),
              M.pointwise5 = numeric(n.test),
              M.pointwise6 = numeric(n.test),
              M.dclf1 = numeric(n.test),
              M.dclf2 = numeric(n.test),
              M.dclf3 = numeric(n.test),
              M.dclf4 = numeric(n.test),
              M.dclf5 = numeric(n.test),
              M.dclf6 = numeric(n.test),
              g.simult = numeric(n.test),
              g.pointwise1 = numeric(n.test),
              g.pointwise2 = numeric(n.test),
              g.pointwise3 = numeric(n.test),
              g.pointwise4 = numeric(n.test),
              g.pointwise5 = numeric(n.test),
              g.pointwise6 = numeric(n.test),
              g.dclf1 = numeric(n.test),
              g.dclf2 = numeric(n.test),
              g.dclf3 = numeric(n.test),
              g.dclf4 = numeric(n.test),
              g.dclf5 = numeric(n.test),
              g.dclf6 = numeric(n.test))
pval.it = 1

for (fish in names) {
    
    display = paste0("LOADING FISH '", fish, "'")
    cat("\n", rep("=", nchar(display)), "\n", rep("=", nchar(display)), "\n", display, "\n\n", sep = "")
    
    source(paste0(fish, "/import_inhom.R"))
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    
    # Rescale to cell diameter
    cells.diameter = plist %>% sapply(function(p) mean(nndist(p))) %>% mean()
    plist = plist %>% solapply(function(p) rescale(X=p, s=cells.diameter, unitname="cell diameter"))
    
    # # Compute density of whole pattern
    # p = rescale(X=p, s=cells.diameter, unitname="cell diameter")
    # p.dens.MC = density(p) * (1.0 / length(plist)) * (sum(cells$type == "MC") / nrow(cells))
    
    for (delta in deltas) {
        
        display = paste0("Starting computations for delta '", delta, "'")
        cat("\n", rep("=", nchar(display)), "\n", display, "\n\n", sep = "")
        
        path.to.image = paste0("_images/", fish, "/delta", delta, "/")
        dir.create(path.to.image, recursive = TRUE) 
        pval$fish[pval.it] = fish
        pval$delta[pval.it] = delta
     
        cells.delta = map_dfr(1:(8-delta), function(image) {
            cells %>% filter(cliche == image) %>% bind_rows(
                cells %>% filter(cliche == image + delta, type == "MC") %>% mutate(type = "MCfuture")) %>% 
                mutate(cliche = paste0(image, image + delta), type = as.factor(type))
        })
        
        cliches = unique(cells.delta$cliche)
        plist.delta = cliches %>% solapply(function(image) {
            cells.delta %>% filter(cliche == image) %>% {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = .$type)} %>% rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        ggplot(cells.delta, aes(x=dim1, y=dim2)) +
            geom_point(aes(colour=type), size=3) +
            geom_path(aes(x=x, y=y), data=bind_rows(bdry, bdry[1,])) +
            facet_wrap(~ cliche, ncol=2) +
            # scale_colour_manual(values = c(gg_color_hue(3), "black")) +
            scale_colour_manual(values = c("#7030A0", "#00A0C6", "#FF0000", "#C3D69B")) +
            coord_equal() +
            theme_grey()
        
        ggsave(filename=paste0(path.to.image, "patterns.eps"), device="ps", width=12, height=8, units="in")
        
        n.sim1 = 499
        n.sim2 = 500
        n.sim = n.sim1 + n.sim2
        cells.delta.tosample = cells.delta %>% filter(type == "qNSC" | type == "aNSC" | type == "MCfuture")
        cells.delta.MC = cells.delta %>% filter(type == "MC")
        ppsim = anylapply(1:n.sim, function(sim) {
            permut = map_dfr(cliches, function(image) {
                cells.delta.permut = cells.delta.tosample %>% filter(cliche==image)
                n.MCfuture = cells.delta.permut %>% filter(type == "MCfuture") %>% nrow()
                n.tot = cells.delta.permut %>% nrow()
                cells.delta.permut[ sample(x = n.tot, size = n.MCfuture, replace=FALSE),] %>% 
                    mutate(type="MCfuture") %>% 
                    bind_rows(cells.delta.MC %>% filter(cliche==image)) %>% 
                    mutate(type=as.factor(type))
            })
            
            solapply(cliches, function(image) {
                permut %>% filter(cliche == image) %>% 
                {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = as.factor(.$type))} %>% 
                    rescale(s=cells.diameter, unitname="cell diameter")
            })
        })
        
        ##############
        ### L FUNCTION
        cat("L. ")
        
        Lsim = anylapply(1:n.sim, function(sim) {
            Leach = anylapply(ppsim[[sim]], Lcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
            pool(Leach)
        })
        
        Leach = lapply(plist.delta, Lcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
        L = pool(as.anylist(Leach))
        plot(L)
        
        # No function for envelopes of replicated point patterns
        # -> separate Lsim into the simulations for mean (n.sim2) and the simulations for deviation (n.sim1)
        # -> compute mean
        # -> compute maximum deviation
        Lmean = map_dfc((n.sim1+1):n.sim, function(sim) {
            Lsim[[sim]]$pooliso
        }) %>% apply(1, mean)
        
        # Go for univariate testing
        Ldeviation = sapply(1:n.sim1, function(sim) {
            min( pmin(Lsim[[sim]]$pooliso - Lmean, 0) )
        })
        
        Ldeviation.sorted = sort(Ldeviation, decreasing = FALSE)
        Ldeviation.observed = min( pmin(L$pooliso - Lmean, 0) )
        pval$L.simult[pval.it] = base::Position(function(bool) bool == FALSE, Ldeviation.sorted < Ldeviation.observed) / (n.sim1 + 1)
        
        Lpointwise = map_dfc(1:n.sim, function(sim) {
            Lsim[[sim]]$pooliso 
        })
        Lpointwise.sorted = apply(Lpointwise, 1, sort, decreasing = FALSE)
        Lpointwise.pval = sapply(1:ncol(Lpointwise.sorted), function(index) {
            base::Position(function(bool) bool == TRUE,
                           Lpointwise.sorted[,index] > L$pooliso[index]
            ) / (n.sim + 1)
        })
        pval$L.pointwise1[pval.it] = Lpointwise.pval[11]  # r starts at 0
        pval$L.pointwise2[pval.it] = Lpointwise.pval[21]  # r starts at 0
        pval$L.pointwise3[pval.it] = Lpointwise.pval[31]  # r starts at 0
        pval$L.pointwise4[pval.it] = Lpointwise.pval[41]  # r starts at 0
        pval$L.pointwise5[pval.it] = Lpointwise.pval[51]
        pval$L.pointwise6[pval.it] = Lpointwise.pval[61]
        
        EL.simultaneous = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lmean,
            lo = Lmean + Ldeviation.sorted[5],
            hi = Inf
        ) %>% as_tibble()
        
        EL.pointwise = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lpointwise %>% apply(1, mean),
            lo = Lpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[10]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(EL.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(L["3, 2"](r)))) +
            ggtitle("EL.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "ELsimultaneous01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(EL.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(L["3, 2"](r)))) +
            ggtitle("EL.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "ELpointwise01.eps"), device="ps", width=12, height=8, units="in")
        
        EL.simultaneous = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lmean,
            lo = Lmean + Ldeviation.sorted[25],
            hi = Inf
        ) %>% as_tibble()
        
        EL.pointwise = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lpointwise %>% apply(1, mean),
            lo = Lpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[50]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(EL.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(L["3, 2"](r)))) +
            ggtitle("EL.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "ELsimultaneous05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(EL.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(L["3, 2"](r)))) +
            ggtitle("EL.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "ELpointwise05.eps"), device="ps", width=12, height=8, units="in")
        
        Ldclf.deviation = sapply(1:n.sim1, function(sim) {
            pmin(Lsim[[sim]]$pooliso - Lmean, 0)^2 %>% cumsum()
        })
        
        Ldclf.deviation.sorted = apply(Ldclf.deviation, 1, sort, decreasing = TRUE)
        Ldclf.observed = pmin(L$pooliso - Lmean, 0)^2 %>% cumsum()
        Ldclf.pval = sapply(1:ncol(Ldclf.deviation.sorted), function(index) {
            base::Position(function(bool) bool == FALSE,
                           Ldclf.deviation.sorted[,index] >= Ldclf.observed[index]
            ) / (n.sim1 + 1)
        })
        pval$L.dclf1[pval.it] = Ldclf.pval[11]  # r starts at 0
        pval$L.dclf2[pval.it] = Ldclf.pval[21]  # r starts at 0
        pval$L.dclf3[pval.it] = Ldclf.pval[31]  # r starts at 0
        pval$L.dclf4[pval.it] = Ldclf.pval[41]  # r starts at 0
        pval$L.dclf5[pval.it] = Ldclf.pval[51]
        pval$L.dclf6[pval.it] = Ldclf.pval[61]
        
        Ldclf.progress = list(
            r = L$r,
            obs = Ldclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(Ldclf.deviation.sorted, 2, function(x) x[5])
        ) %>% as_tibble()
        
        Ldclf.sigtrace = list(
            r = L$r,
            obs = Ldclf.pval,
            risk = 0.01
        ) %>% as_tibble()
        
        ggplot(Ldclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(L["3, 2"](s)*ds, 0, r)))) +
            ggtitle("Ldclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Ldclfprogress01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Ldclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("Ldclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Ldclfsigtrace01.eps"), device="ps", width=12, height=8, units="in")
        
        Ldclf.progress = list(
            r = L$r,
            obs = Ldclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(Ldclf.deviation.sorted, 2, function(x) x[25])
        ) %>% as_tibble()
        
        Ldclf.sigtrace = list(
            r = L$r,
            obs = Ldclf.pval,
            risk = 0.05
        ) %>% as_tibble()
        
        ggplot(Ldclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(L["3, 2"](s)*ds, 0, r)))) +
            ggtitle("Ldclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Ldclfprogress05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Ldclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("Ldclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Ldclfsigtrace05.eps"), device="ps", width=12, height=8, units="in")
        
        ##############
        ### M FUNCTION
        cat("M. ")
        
        mplist.delta = cliches %>% solapply(function(image) {
            mp = cells.delta %>% filter(cliche == image) %>% 
                select(dim1, dim2, type) %>% 
                rename(x = dim1, y = dim2, PointType = type) %>% 
                as.data.frame() %>% 
                wmppp(window=cells.owin) %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
            attr(mp, "class") = c("wmppp", "ppp")
            mp
        })
        
        mppsim = anylapply(1:n.sim, function(sim) {
            permut = map_dfr(cliches, function(image) {
                cells.delta.permut = cells.delta.tosample %>% filter(cliche==image)
                n.MCfuture = (cells.delta.permut$type == "MCfuture") %>% sum()
                n.tot = cells.delta.permut %>% nrow()
                MC.sampled = sample(x = n.tot, size = n.MCfuture, replace=FALSE)
                cells.delta.permut[ MC.sampled,] %>% 
                    mutate(type="MCfuture") %>% 
                    bind_rows(
                        cells.delta.permut[-MC.sampled,] %>% mutate(type=plyr::mapvalues(type, "MCfuture", "qNSC"))
                    ) %>% 
                    bind_rows(
                        cells.delta.MC %>% filter(cliche==image)
                    ) %>% 
                    mutate(type=as.factor(type))
            })
            
            solapply(cliches, function(image) {
                mp = permut %>% filter(cliche == image) %>% 
                    select(dim1, dim2, type) %>% 
                    rename(x = dim1, y = dim2, PointType = type) %>% 
                    as.data.frame() %>% 
                    wmppp(window=cells.owin) %>% 
                    rescale(s=cells.diameter, unitname="cell diameter")
                attr(mp, "class") = c("wmppp", "ppp")
                mp
            })
        })
        
        Msim = anylapply(1:n.sim, function(sim) {
            Meach = anylapply(mppsim[[sim]], Mhat, r=seq(0, 6, by=0.1), ReferenceType="MC", NeighborType="MCfuture", CheckArguments = FALSE)
            pool(Meach)
        })
        
        Meach = anylapply(mplist.delta, Mhat, r=seq(0, 6, by=0.1), ReferenceType="MC", NeighborType="MCfuture", CheckArguments = FALSE)
        M = pool(as.anylist(Meach))
        plot(M$r, M$poolM, type="l")
        lines(M$r, M$loM, lty="dashed")
        lines(M$r, M$hiM, lty="dashed")
        abline(h=1, col="grey")
        
        # No function for envelopes of replicated point patterns
        # -> separate Msim into the simulations for mean (n.sim2) and the simulations for deviation (n.sim1)
        # -> compute mean
        # -> compute maximum deviation
        Mmean = map_dfc((n.sim1+1):n.sim, function(sim) {
            Msim[[sim]]$poolM
        }) %>% apply(1, mean)
        
        # Go for univariate testing
        Mdeviation = sapply(1:n.sim1, function(sim) {
            min( pmin(Msim[[sim]]$poolM - Mmean, 0, na.rm = TRUE) )
        })
        
        Mdeviation.sorted = sort(Mdeviation, decreasing = FALSE)
        Mdeviation.observed = min( pmin(M$poolM - Mmean, 0, na.rm = TRUE) )
        pval$M.simult[pval.it] = base::Position(function(bool) bool == FALSE, Mdeviation.sorted <= Mdeviation.observed) / (n.sim1 + 1)
        
        Mpointwise = map_dfc(1:n.sim, function(sim) {
            Msim[[sim]]$poolM 
        })
        Mpointwise.sorted = apply(Mpointwise, 1, sort, decreasing = FALSE, na.last = FALSE)
        Mpointwise.pval = sapply(1:ncol(Mpointwise.sorted), function(index) {
            if (any(is.na(Mpointwise.sorted[,index]))) return(NA)
            base::Position(function(bool) bool == TRUE,
                           Mpointwise.sorted[,index] > M$poolM[index]
            ) / (n.sim + 1)
        })
        pval$M.pointwise1[pval.it] = Mpointwise.pval[11]  # r starts at 0
        pval$M.pointwise2[pval.it] = Mpointwise.pval[21]  # r starts at 0
        pval$M.pointwise3[pval.it] = Mpointwise.pval[31]  # r starts at 0
        pval$M.pointwise4[pval.it] = Mpointwise.pval[41]  # r starts at 0
        pval$M.pointwise5[pval.it] = Mpointwise.pval[51]
        pval$M.pointwise6[pval.it] = Mpointwise.pval[61]
        
        EM.simultaneous = list(
            r = M$r,
            obs = M$poolM,
            mmean = Mmean,
            lo = Mmean + Mdeviation.sorted[5],
            hi = Inf
        ) %>% as_tibble()
        
        EM.pointwise = list(
            r = M$r,
            obs = M$poolM,
            mmean = Mpointwise %>% apply(1, mean),
            lo = Mpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[10]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(EM.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(M["3, 2"](r)))) +
            ggtitle("EM.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "EMsimultaneous01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(EM.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(M["3, 2"](r)))) +
            ggtitle("EM.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "EMpointwise01.eps"), device="ps", width=12, height=8, units="in")
        
        EM.simultaneous = list(
            r = M$r,
            obs = M$poolM,
            mmean = Mmean,
            lo = Mmean + Mdeviation.sorted[25],
            hi = Inf
        ) %>% as_tibble()
        
        EM.pointwise = list(
            r = M$r,
            obs = M$poolM,
            mmean = Mpointwise %>% apply(1, mean),
            lo = Mpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[50]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(EM.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(M["3, 2"](r)))) +
            ggtitle("EM.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "EMsimultaneous05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(EM.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(M["3, 2"](r)))) +
            ggtitle("EM.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "EMpointwise05.eps"), device="ps", width=12, height=8, units="in")
        
        Mdclf.deviation = sapply(1:n.sim1, function(sim) {
            pmin(Msim[[sim]]$poolM - Mmean, 0, na.rm = TRUE)^2 %>% cumsum()
        })
        
        Mdclf.deviation.sorted = apply(Mdclf.deviation, 1, sort, decreasing = TRUE)
        Mdclf.observed = pmin(M$poolM - Mmean, 0, na.rm = TRUE)^2 %>% cumsum()
        Mdclf.pval = sapply(1:ncol(Mdclf.deviation.sorted), function(index) {
            base::Position(function(bool) bool == FALSE,
                           Mdclf.deviation.sorted[,index] >= Mdclf.observed[index]
            ) / (n.sim1 + 1)
        })
        pval$M.dclf1[pval.it] = Mdclf.pval[11]  # r starts at 0
        pval$M.dclf2[pval.it] = Mdclf.pval[21]  # r starts at 0
        pval$M.dclf3[pval.it] = Mdclf.pval[31]  # r starts at 0
        pval$M.dclf4[pval.it] = Mdclf.pval[41]  # r starts at 0
        pval$M.dclf5[pval.it] = Mdclf.pval[51]
        pval$M.dclf6[pval.it] = Mdclf.pval[61]
        
        Mdclf.progress = list(
            r = M$r,
            obs = Mdclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(Mdclf.deviation.sorted, 2, function(x) x[5])
        ) %>% as_tibble()
        
        Mdclf.sigtrace = list(
            r = M$r,
            obs = Mdclf.pval,
            risk = 0.01
        ) %>% as_tibble()
        
        ggplot(Mdclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(M["3, 2"](s)*ds, 0, r)))) +
            ggtitle("Mdclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Mdclfprogress01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Mdclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("Mdclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Mdclfsigtrace01.eps"), device="ps", width=12, height=8, units="in")
        
        Mdclf.progress = list(
            r = M$r,
            obs = Mdclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(Mdclf.deviation.sorted, 2, function(x) x[25])
        ) %>% as_tibble()
        
        Mdclf.sigtrace = list(
            r = M$r,
            obs = Mdclf.pval,
            risk = 0.05
        ) %>% as_tibble()
        
        ggplot(Mdclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(M["3, 2"](s)*ds, 0, r)))) +
            ggtitle("Mdclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Mdclfprogress05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Mdclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("Mdclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Mdclfsigtrace05.eps"), device="ps", width=12, height=8, units="in")
        
        ##############
        ### g FUNCTION
        cat("g. ")
        
        gsim = anylapply(1:n.sim, function(sim) {
            geach = anylapply(ppsim[[sim]], pcfcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
            pool(geach)
        })
        
        geach = lapply(plist.delta, pcfcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
        g = pool(as.anylist(geach))
        plot(g)
        
        # No function for envelopes of replicated point patterns
        # -> separate gsim into the simulations for mean (n.sim2) and the simulations for deviation (n.sim1)
        # -> compute mean
        # -> compute maximum deviation
        gmean = map_dfc((n.sim1+1):n.sim, function(sim) {
            gsim[[sim]]$pooliso
        }) %>% apply(1, mean)
        
        # Go for univariate testing
        gdeviation = sapply(1:n.sim1, function(sim) {
            min( pmin(gsim[[sim]]$pooliso - gmean, 0, na.rm = TRUE) )
        })
        
        gdeviation.sorted = sort(gdeviation, decreasing = FALSE)
        gdeviation.observed = min( pmin(g$pooliso - gmean, 0, na.rm = TRUE) )
        pval$g.simult[pval.it] = base::Position(function(bool) bool == FALSE, gdeviation.sorted <= gdeviation.observed) / (n.sim1 + 1)
        
        gpointwise = map_dfc(1:n.sim, function(sim) {
            gsim[[sim]]$pooliso 
        })
        gpointwise.sorted = apply(gpointwise, 1, sort, decreasing = FALSE, na.last = FALSE)
        gpointwise.pval = sapply(1:ncol(gpointwise.sorted), function(index) {
            if (any(is.na(Mpointwise.sorted[,index]))) return(NA)
            base::Position(function(bool) bool == TRUE,
                           gpointwise.sorted[,index] > g$pooliso[index]
            ) / (n.sim + 1)
        })
        pval$g.pointwise1[pval.it] = gpointwise.pval[11]  # r starts at 0
        pval$g.pointwise2[pval.it] = gpointwise.pval[21]  # r starts at 0
        pval$g.pointwise3[pval.it] = gpointwise.pval[31]  # r starts at 0
        pval$g.pointwise4[pval.it] = gpointwise.pval[41]  # r starts at 0
        pval$g.pointwise5[pval.it] = gpointwise.pval[51]
        pval$g.pointwise6[pval.it] = gpointwise.pval[61]
        
        Eg.simultaneous = list(
            r = g$r,
            obs = g$pooliso,
            mmean = gmean,
            lo = gmean + gdeviation.sorted[5],
            hi = Inf
        ) %>% as_tibble()
        
        Eg.pointwise = list(
            r = g$r,
            obs = g$pooliso,
            mmean = gpointwise %>% apply(1, mean),
            lo = gpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[10]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(Eg.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(g["3, 2"](r)))) +
            ggtitle("Eg.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Egsimultaneous01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Eg.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(g["3, 2"](r)))) +
            ggtitle("Eg.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Egpointwise01.eps"), device="ps", width=12, height=8, units="in")
        
        Eg.simultaneous = list(
            r = g$r,
            obs = g$pooliso,
            mmean = gmean,
            lo = gmean + gdeviation.sorted[25],
            hi = Inf
        ) %>% as_tibble()
        
        Eg.pointwise = list(
            r = g$r,
            obs = g$pooliso,
            mmean = gpointwise %>% apply(1, mean),
            lo = gpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[50]}),
            hi = Inf
        ) %>% as_tibble()
        
        ggplot(Eg.simultaneous, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(g["3, 2"](r)))) +
            ggtitle("Eg.simultaneous") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Egsimultaneous05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(Eg.pointwise, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(g["3, 2"](r)))) +
            ggtitle("Eg.pointwise") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "Egpointwise05.eps"), device="ps", width=12, height=8, units="in")
        
        gdclf.deviation = sapply(1:n.sim1, function(sim) {
            pmin(gsim[[sim]]$pooliso - gmean, 0, na.rm = TRUE)^2 %>% cumsum()
        })
        
        gdclf.deviation.sorted = apply(gdclf.deviation, 1, sort, decreasing = TRUE)
        gdclf.observed = pmin(g$pooliso - gmean, 0, na.rm = TRUE)^2 %>% cumsum()
        gdclf.pval = sapply(1:ncol(gdclf.deviation.sorted), function(index) {
            base::Position(function(bool) bool == FALSE,
                           gdclf.deviation.sorted[,index] >= gdclf.observed[index]
            ) / (n.sim1 + 1)
        })
        pval$g.dclf1[pval.it] = gdclf.pval[11]  # r starts at 0
        pval$g.dclf2[pval.it] = gdclf.pval[21]  # r starts at 0
        pval$g.dclf3[pval.it] = gdclf.pval[31]  # r starts at 0
        pval$g.dclf4[pval.it] = gdclf.pval[41]  # r starts at 0
        pval$g.dclf5[pval.it] = gdclf.pval[51]
        pval$g.dclf6[pval.it] = gdclf.pval[61]
        
        gdclf.progress = list(
            r = g$r,
            obs = gdclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(gdclf.deviation.sorted, 2, function(x) x[5])
        ) %>% as_tibble()
        
        gdclf.sigtrace = list(
            r = g$r,
            obs = gdclf.pval,
            risk = 0.01
        ) %>% as_tibble()
        
        ggplot(gdclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(g["3, 2"](s)*ds, 0, r)))) +
            ggtitle("gdclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "gdclfprogress01.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(gdclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("gdclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "gdclfsigtrace01.eps"), device="ps", width=12, height=8, units="in")
        
        gdclf.progress = list(
            r = g$r,
            obs = gdclf.observed,
            mmean = 0,
            lo = 0,
            hi = apply(gdclf.deviation.sorted, 2, function(x) x[25])
        ) %>% as_tibble()
        
        gdclf.sigtrace = list(
            r = g$r,
            obs = gdclf.pval,
            risk = 0.05
        ) %>% as_tibble()
        
        ggplot(gdclf.progress, aes(x=r)) +
            geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
            geom_line(aes(y=mmean), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic(integral(g["3, 2"](s)*ds, 0, r)))) +
            ggtitle("gdclf.progress") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "gdclfprogress05.eps"), device="ps", width=12, height=8, units="in")
        
        ggplot(gdclf.sigtrace, aes(x=r)) +
            geom_hline(aes(yintercept=risk), col="red", lty="dashed") +
            geom_line(aes(y=obs)) +
            xlab(expression(paste(italic("r"), " (cell diameter)"))) +
            ylab(expression(italic("p-value"))) +
            ggtitle("gdclf.sigtrace") +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
        
        ggsave(paste0(path.to.image, "gdclfsigtrace05.eps"), device="ps", width=12, height=8, units="in")
        
        rm(list=c("ppsim", "Lsim", "Msim", "gsim"))
        pval.it = pval.it + 1
        
    }
    
}

logna0 = function(x) ifelse(!is.na(x), log(x), 0)
pval.byfish = gather(pval, "stat", "pvalue", L.simult:g.dclf6) %>% spread(fish, pvalue) %>% 
    mutate(pooled = 1 - pchisq(-2*(logna0(bibi) + logna0(mimi) + logna0(titi)), 
                               df = 6 - 2 * is.na(bibi) - 2 * is.na(mimi) - 2 * is.na(titi)))
save(pval.byfish, file = "pval.RData")

pdf("pval.pdf", height=25, width=6)
grid.table(
    pval.byfish %>%
        filter(
            !str_detect(stat, "[gM]\\.dclf"),
            !str_detect(stat, "[gM]\\.simult"),
            !str_detect(stat, "\\.pointwise[5,6]"),
            !str_detect(stat, ".dclf[1,2]")
        )
)
dev.off()
