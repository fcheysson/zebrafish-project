# Spatial (non temporal) analysis of the simulations from the NSC lattice model
# Fig.4E, Supp. Fig.12

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

# STATISTICAL ANALYSIS
setEPS()
n.sim1 = 499
n.sim2 = 500
n.sim = n.sim1 + n.sim2

model = "Positions_withLI_f059_simulation"
simulations = 1L:18L
n.test = length(simulations)
pval = tibble(simulation = numeric(n.test), 
              L22.simult = numeric(n.test),
              L33.simult = numeric(n.test),
              L32.simult = numeric(n.test))

dir.create("static_images")

pval.it = 1

for (simulation in simulations) {
    
    import.display = paste0("Simulation #", simulation)
    cat("\n", rep("=", nchar(import.display)), "\n", import.display, "\n\n", sep = "")
    
    pval$simulation[pval.it] = simulation

    matlabObject = R.matlab::readMat(paste0("data/", model, simulation, ".mat"))
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
    
    # CREATE SOLIST
    availableTimeSteps = cells %>% pull(timeStep) %>% unique() %>% {.[.>=500]}
    timeStepCount = availableTimeSteps %>% range() %>% {.[2] - .[1] + 1}
    cellPatterns = solapply(availableTimeSteps, function(image) {
        cells %>% 
            filter(timeStep == image) %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)}
    })
    
    # Rescale to cell diameter
    # cells.diameter = cellPatterns %>% sapply(function(p) mean(nndist(p))) %>% mean()
    cells.diameter = cellPatterns %>% sapply(function(p) 2 * sqrt(area(Window(p)) / (pi * npoints(p)))) %>% mean()
    
    # Which timestep to use
    timeSteps = c(400, 450, 500, 550, 600)
    cells = cells %>% filter(timeStep %in% timeSteps)
    
    ggplot(cells, aes(x=x, y=y)) +
        geom_point(aes(colour=type), size=3) +
        geom_hline(yintercept = cellWindow$yrange) +
        geom_vline(xintercept = cellWindow$xrange) +
        facet_wrap(~ timeStep) +
        scale_colour_manual(values = c("grey", gg_color_hue(length(cellTypes)-1))) +
        coord_equal() +
        theme_grey()
    
    ggsave(paste0("static_images/", model, simulation, "_pattern.eps"), device="ps", width=12, height=8, units="in")

    ##############################
    ## Permutation cellules type 2
    cat("L22.")
    
    plist = timeSteps %>% solapply(function(image) {
        cells  %>% 
            filter(timeStep == image, type == "aNSC") %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
            rescale(s=cells.diameter, unitname="cell diameter")
    })
    
    cells.tosample = cells %>% filter(type == "qNSC" | type == "aNSC")
    ppsim = anylapply(1:n.sim, function(sim) {
        permut = map_dfr(timeSteps, function(image) {
            cells.permut = cells.tosample %>% filter(timeStep == image)
            n.aNSC = cells.permut %>% filter(type == "aNSC") %>% nrow()
            n.tot = cells.permut %>% nrow()
            cells.permut[ sample(x = n.tot, size = n.aNSC, replace=FALSE),] %>% 
                                                   #rbinom(1, size = n.tot, prob = n.aNSC / n.tot)
                mutate(type="aNSC")
        })
        
        solapply(timeSteps, function(image) {
            permut %>% filter(timeStep == image) %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = as.factor(.$type))} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
    })
    
    Lsim = anylapply(1:n.sim, function(sim) {
        Leach = anylapply(ppsim[[sim]], Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
        pool(Leach)
    })
    
    Leach = lapply(plist, Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
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
    Ldeviation = sapply(1:n.sim1, function(sim) {
        max( abs(Lsim[[sim]]$pooliso - Lmean) )
    })
    
    Ldeviation.sorted = sort(Ldeviation, decreasing = FALSE)
    Ldeviation.observed = max( abs(L$pooliso - Lmean) )
    pval$L22.simult[pval.it] = 1 - base::Position(function(bool) bool == TRUE, Ldeviation.sorted < Ldeviation.observed, right = TRUE) / (n.sim1 + 1)
    
    EL.simultaneous = list(
        r = L$r,
        obs = L$pooliso,
        mmean = Lmean,
        lo = Lmean - Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)], # 100 for alpha = 0.01
        hi = Lmean + Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)]  # 100 for alpha = 0.01
    ) %>% as_tibble()
    
    ggplot(EL.simultaneous, aes(x=r)) +
        geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
        geom_line(aes(y=mmean), col="red", lty="dashed") +
        geom_line(aes(y=obs)) +
        xlab(expression(paste(italic("r"), " (cell diameter)"))) +
        ylab(expression(italic(L["2, 2"](r)))) +
        ggtitle("EL.simultaneous") +
        theme_bw() +
        theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
    
    ggsave(paste0("static_images/", model, simulation, "_EL22simultaneous05.eps"), device="ps", width=12, height=8, units="in")
    rm(list=c("ppsim", "Lsim"))
    
    ##############################
    ## Permutation cellules type 3
    cat("L33.")
    
    plist = timeSteps %>% solapply(function(image) {
        cells  %>% 
            filter(timeStep == image, type %in% c("aNP", "div.aNP")) %>% 
            mutate(type = "aNP") %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
            rescale(s=cells.diameter, unitname="cell diameter")
    })
    
    cells.tosample = cells %>% filter(type == "qNSC" | type == "aNSC" | type == "aNP" | type == "div.aNP")
    ppsim = anylapply(1:n.sim, function(sim) {
        permut = map_dfr(timeSteps, function(image) {
            cells.permut = cells.tosample %>% filter(timeStep == image)
            n.aNP = cells.permut %>% filter(type %in% c("aNP", "div.aNP")) %>% nrow()
            n.tot = cells.permut %>% nrow()
            cells.permut[ sample(x = n.tot, size = n.aNP, replace=FALSE),] %>% 
                #rbinom(1, size = n.tot, prob = n.aNP / n.tot)
                mutate(type="aNP")
        })
        
        solapply(timeSteps, function(image) {
            permut %>% filter(timeStep == image) %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = as.factor(.$type))} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
    })
    
    Lsim = anylapply(1:n.sim, function(sim) {
        Leach = anylapply(ppsim[[sim]], Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
        pool(Leach)
    })
    
    Leach = lapply(plist, Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
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
    Ldeviation = sapply(1:n.sim1, function(sim) {
        max( abs(Lsim[[sim]]$pooliso - Lmean) )
    })
    
    Ldeviation.sorted = sort(Ldeviation, decreasing = FALSE)
    Ldeviation.observed = max( abs(L$pooliso - Lmean) )
    pval$L33.simult[pval.it] = 1 - base::Position(function(bool) bool == TRUE, Ldeviation.sorted < Ldeviation.observed, right = TRUE) / (n.sim1 + 1)
    
    EL.simultaneous = list(
        r = L$r,
        obs = L$pooliso,
        mmean = Lmean,
        lo = Lmean - Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)], # 100 for alpha = 0.01
        hi = Lmean + Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)]  # 100 for alpha = 0.01
    ) %>% as_tibble()
    
    ggplot(EL.simultaneous, aes(x=r)) +
        geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey80") +
        geom_line(aes(y=mmean), col="red", lty="dashed") +
        geom_line(aes(y=obs)) +
        xlab(expression(paste(italic("r"), " (cell diameter)"))) +
        ylab(expression(italic(L["3, 3"](r)))) +
        ggtitle("EL.simultaneous") +
        theme_bw() +
        theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5))
    
    ggsave(paste0("static_images/", model, simulation, "_EL33simultaneous05.eps"), device="ps", width=12, height=8, units="in")
    rm(list=c("ppsim", "Lsim"))
    
    ###############################
    ## Permutation cellules type 32
    cat("L32.")
    
    plist = timeSteps %>% solapply(function(image) {
        cells  %>% 
            filter(timeStep == image) %>% 
            mutate(type = recode(type, div.aNP = "aNP")) %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
            rescale(s=cells.diameter, unitname="cell diameter")
    })
    
    cells.tosample = cells %>% filter(type == "qNSC" | type == "aNSC")
    cells.aNP = cells %>% 
        filter(type %in% c("aNP", "div.aNP")) %>% 
        mutate(type = recode(type, div.aNP = "aNP")) 
    ppsim = anylapply(1:n.sim, function(sim) {
        permut = map_dfr(timeSteps, function(image) {
            cells.permut = cells.tosample %>% filter(timeStep == image)
            n.aNSC = cells.permut %>% filter(type == "aNSC") %>% nrow()
            n.tot = cells.permut %>% nrow()
            cells.permut[ sample(x = n.tot, size = n.aNSC, replace=FALSE),] %>% 
                #rbinom(1, size = n.tot, prob = n.aNSC / n.tot)
                mutate(type = recode(type, qNSC = "aNSC", .default = "aNSC"))
        }) %>% bind_rows(cells.aNP)
        
        solapply(timeSteps, function(image) {
            permut %>% filter(timeStep == image) %>% 
            {ppp(x = .$x, y = .$y, window = cellWindow, marks = as.factor(.$type))} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
    })
    
    Lsim = anylapply(1:n.sim, function(sim) {
        Leach = anylapply(ppsim[[sim]], Lcross, r=seq(0, 6, by=0.1), i="aNP", j="aNSC", ratio=TRUE, correction="best")
        pool(Leach)
    })
    
    Leach = lapply(plist, Lcross, r=seq(0, 6, by=0.1), i="aNP", j="aNSC", ratio=TRUE, correction="best")
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
    Ldeviation = sapply(1:n.sim1, function(sim) {
        max( abs(Lsim[[sim]]$pooliso - Lmean) )
    })
    
    Ldeviation.sorted = sort(Ldeviation, decreasing = FALSE)
    Ldeviation.observed = max( abs(L$pooliso - Lmean) )
    pval$L32.simult[pval.it] = 1 - base::Position(function(bool) bool == TRUE, Ldeviation.sorted < Ldeviation.observed, right = TRUE) / (n.sim1 + 1)
    
    EL.simultaneous = list(
        r = L$r,
        obs = L$pooliso,
        mmean = Lmean,
        lo = Lmean - Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)], # 100 for alpha = 0.01
        hi = Lmean + Ldeviation.sorted[ceiling((n.sim1+1) - (n.sim1+1)/20)]  # 100 for alpha = 0.01
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
    
    ggsave(paste0("static_images/", model, simulation, "_EL32simultaneous05.eps"), device="ps", width=12, height=8, units="in")
    rm(list=c("ppsim", "Lsim"))
    
    pval.it = pval.it + 1
    save(pval, file = "static_pval.RData")
    
}

pval %>% summarise_at(vars(starts_with("L")), ~ 1 - pchisq(- 2 * sum(log(.)), df = 2 * length(.)))

pdf("static_pval.pdf", height=5.5, width=4.5)
grid.table(
    pval %>% mutate_at(vars(starts_with("L")), ~ format(round(., 3), nsmall = 3))
)
dev.off()