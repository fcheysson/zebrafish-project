# Spatio-temporal analysis of MCs in the simulations from the NSC lattice model

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gridExtra)
library(tidyverse)
library(spatstat)
library(dbmss)

# GGPLOT COLOR PALETTE
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

n.sim1 = 99
n.sim2 = 100
n.sim = n.sim1 + n.sim2

simulations = 1L:34L
deltas = 1L:5L
n.test = length(deltas)

# System path and folder manipulation
file = "Positions_withLI_f059_simulation"
dir.create(file)

for (simulation in simulations) {
    
    dir.create(paste0(file, "/", file, simulation))
    dir.create(paste0(file, "/", file, simulation, "/images"))
    
    matlabObject = R.matlab::readMat(paste0("data/", file, simulation, ".mat"))
    matlabList = matlabObject[[1]]
    cellTypes = attr(matlabList, "dimnames")[[1]]
    
    cells = map_dfr(1:dim(matlabList)[3], function(timeStep) {
        map_dfr(1:length(cellTypes), function(index) {
            
            if (is.null(matlabList[,,timeStep][[index]])) return()
            if (!nrow(matlabList[,,timeStep][[index]])) return()
            
            matlabList[,,timeStep][[index]] %>% 
                as_tibble() %>% 
                setNames(c("x", "y")) %>% 
                mutate(
                    type = factor(cellTypes[index], cellTypes),
                    timeStep = timeStep
                )
            
        })
    })
    
    timeStepCount = cells %>% pull(timeStep) %>% range() %>% {.[2] - .[1] + 1}
    meanPointsPerTimeStep = nrow(cells) / timeStepCount
    
    cellWindow = ripras(x = cells$x, y = cells$y, shape = "rectangle", f = (meanPointsPerTimeStep + 1) / (meanPointsPerTimeStep - 1))
    
    # Pool each 3 consecutive time steps
    timeStepMin = cells %>% pull(timeStep) %>% min()
    timeStepMax = cells %>% pull(timeStep) %>% max()
    cells = cells %>% 
        mutate(nts = (timeStep - timeStepMin) %/% 3) %>% 
        filter(type == "MC" | (type %in% c("qNSC", "aNSC") & (timeStep - timeStepMin) %% 3 == 0)) %>% 
        select(-timeStep) %>% 
        rename(timeStep = nts)
    
    # CREATE SOLIST
    availableTimeSteps = cells %>% pull(timeStep) %>% unique()
    timeStepCount = availableTimeSteps %>% range() %>% {.[2] - .[1] + 1}
    cellPatterns = solapply(availableTimeSteps, function(image) {
        cells %>% filter(timeStep == image) %>% {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)}
    })
    
    ggplot(cells %>% filter(timeStep %in% 0:8), aes(x = x, y = y, colour = type)) + geom_point() + 
        geom_hline(yintercept = cellWindow$yrange) +
        geom_vline(xintercept = cellWindow$xrange) +
        scale_colour_manual(values = c("grey", gg_color_hue(length(cellTypes)-1))) +
        facet_wrap(~ timeStep)
    
    # STATISTICAL ANALYSIS
    plist = cellPatterns # Easier to work with plist since below is imported from other file
    
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
                  L.dclf6 = numeric(n.test))
    
    pval.it = 1
    
    # Rescale to cell diameter
    # cells.diameter = plist %>% sapply(function(p) mean(nndist(p))) %>% mean()
    cells.diameter = cellPatterns %>% sapply(function(p) 2 * sqrt(area(Window(p)) / (pi * npoints(p)))) %>% mean()
    plist = plist %>% solapply(function(p) rescale(X=p, s=cells.diameter, unitname="cell diameter"))
    
    for (delta in deltas) {
        
        display = paste0("Starting computations for delta '", delta, "'")
        cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")
        
        path.to.image = paste0(file, "/", file, simulation, "/images/delta", delta, "/")
        dir.create(path.to.image)
        pval$fish[pval.it] = "fish"
        pval$delta[pval.it] = delta
        
        cells.delta = map_dfr(availableTimeSteps[1:(timeStepCount-delta)], function(image) {
            cells %>% filter(timeStep == image) %>% bind_rows(
                cells %>% filter(timeStep == image + delta, type == "MC") %>% mutate(type = "MCfuture")) %>% 
                mutate(timeStep = paste0(image, image + delta))
        }) %>% mutate(type = as.factor(type))
        
        # Remove timeSteps without either MC or MCfuture
        hasMC = cells.delta %>% 
            group_by(timeStep, type) %>% 
            count() %>% 
            spread(type, n) %>% 
            filter(is.na(MC) | is.na(MCfuture)) %>% 
            pull(timeStep)
        cells.delta = cells.delta %>% filter(!timeStep %in% hasMC)
        
        timeSteps = unique(cells.delta$timeStep)
        plist.delta = timeSteps %>% solapply(function(image) {
            cells.delta %>% filter(timeStep == image) %>% {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        ggplot(cells.delta %>% filter(timeStep %in% timeSteps[1:9]), aes(x=x, y=y)) +
            geom_point(aes(colour=type), size=3) +
            geom_hline(yintercept = cellWindow$yrange) +
            geom_vline(xintercept = cellWindow$xrange) +
            facet_wrap(~ timeStep) +
            scale_colour_manual(values = c(gg_color_hue(length(cellTypes)), "grey")) +
            coord_equal() +
            theme_grey()
        
        # ggsave(filename=paste0(path.to.image, "patterns.png"), width=12, height=8, units="in")
        ggsave(filename=paste0(path.to.image, "patterns.eps"), device="ps", width=12, height=8, units="in")
        
        cells.delta.tosample = cells.delta %>% filter(type == "qNSC" | type == "aNSC" | type == "MCfuture") # | type == "aNP" | type == "div.aNP")
        cells.delta.MC = cells.delta %>% filter(type == "MC")
        ppsim = anylapply(1:n.sim, function(sim) {
            permut = map_dfr(timeSteps, function(image) {
                cells.delta.permut = cells.delta.tosample %>% filter(timeStep==image)
                n.MCfuture = cells.delta.permut %>% filter(type == "MCfuture") %>% nrow()
                n.tot = cells.delta.permut %>% nrow()
                cells.delta.permut[ sample(x = n.tot, size = n.MCfuture, replace=FALSE),] %>% 
                    mutate(type="MCfuture") %>% 
                    bind_rows(cells.delta.MC %>% filter(timeStep==image)) %>% 
                    mutate(type=as.factor(type))
            })
            
            solapply(timeSteps, function(image) {
                permut %>% filter(timeStep == image) %>% 
                {ppp(x = .$x, y = .$y, window = cellWindow, marks = as.factor(.$type))} %>% 
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
        pval$L.simult[pval.it] = 1 - sum(Ldeviation.sorted > Ldeviation.observed) / (n.sim1 + 1)
        
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
            lo = Lmean + Ldeviation.sorted[ceiling((n.sim1+1)/100)],
            hi = Inf
        ) %>% as_tibble()
        
        EL.pointwise = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lpointwise %>% apply(1, mean),
            lo = Lpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[ceiling((n.sim+1)/100)]}),
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
            lo = Lmean + Ldeviation.sorted[ceiling((n.sim1+1)/20)],
            hi = Inf
        ) %>% as_tibble()
        
        EL.pointwise = list(
            r = L$r,
            obs = L$pooliso,
            mmean = Lpointwise %>% apply(1, mean),
            lo = Lpointwise %>% apply(1, function(v) {sort(v, decreasing=FALSE)[ceiling((n.sim+1)/20)]}),
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
            hi = apply(Ldclf.deviation.sorted, 2, function(x) x[ceiling((n.sim1+1)/100)])
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
            hi = apply(Ldclf.deviation.sorted, 2, function(x) x[ceiling((n.sim1+1)/20)])
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
        
        rm(list=c("ppsim", "Lsim"))
        pval.it = pval.it + 1
        
    }
    
    save(pval, file = paste0(file, "/", file, simulation, "/pval.RData"))
    
    pdf(paste0(file, "/", file, simulation, "/pval.pdf"), height=2, width=8)
    grid.table(
        pval %>%
            select(delta, L.simult, matches("^L\\.pointwise[1-3]$"), matches("^L\\.dclf[3-6]$"))
    )
    dev.off()
    
}