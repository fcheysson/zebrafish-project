# Studentized permutation test (Hahn, 2012) to compare the interaction between MCs 
# with and without aNP-driven LI
# Fig.7B

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gridExtra)
library(tidyverse)
library(spatstat)
library(latex2exp)

# GGPLOT COLOR PALETTE
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

groups = c("2A", "2B")
files = c("Positions_withLI_f059_simulation", 
          "Positions_withoutLI_f1_simulation")
simulations = list(1L:18L, 1L:18L)

hyperMC = tibble(group = factor(character(0L), levels = groups),
                     simulation = integer(0L), 
                     delta = integer(0L),
                     Leach = anylist(),
                     L = list())

deltas = 1L:5L

for (gr in 1L:length(groups)) {

    display = paste0("Starting computations for '", files[gr], "'")
    cat("\n", rep("=", nchar(display)), "\n", display, "\n", sep = "")
    
    for (simulation in simulations[[gr]]) {
        
        cat(simulation, ". ", sep = "")
        
        matlabObject = R.matlab::readMat(paste0("data/", files[gr], simulation, ".mat"))
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
        availableTimeSteps = cells %>% pull(timeStep) %>% unique()
        
        # Pattern observation window
        cellWindow = ripras(x = cells$x, y = cells$y, shape = "rectangle", f = (meanPointsPerTimeStep + 1) / (meanPointsPerTimeStep - 1))
        
        # Pattern mean cell diameter
        cellPatterns = solapply(availableTimeSteps, function(image) {
            cells %>% filter(timeStep == image) %>% {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)}
        })
        # cells.diameter = cellPatterns %>% sapply(function(p) mean(nndist(p))) %>% mean()
        cells.diameter = cellPatterns %>% sapply(function(p) 2 * sqrt(area(Window(p)) / (pi * npoints(p)))) %>% mean()
        
        # Pool each 3 consecutive time steps and keep only MC
        timeStepMin = cells %>% pull(timeStep) %>% min()
        timeStepMax = cells %>% pull(timeStep) %>% max()
        cells = cells %>% 
            mutate(nts = (timeStep - timeStepMin) %/% 3) %>% 
            filter(type == "MC") %>% 
            select(-timeStep) %>% 
            rename(timeStep = nts) %>% 
            add_count(timeStep)
        
        availableTimeSteps = cells %>% pull(timeStep) %>% unique()
        timeStepCount = availableTimeSteps %>% range() %>% {.[2] - .[1] + 1}
        
        # Delta 0
        availableTS = cells %>% filter(n > 1) %>% pull(timeStep) %>% unique()
        cellPatterns = solapply(availableTS, function(image) {
            cells %>% filter(timeStep == image) %>% 
                {ppp(x = .$x, y = .$y, window = cellWindow)} %>% 
                rescale(s=cells.diameter, unitname=c("cell diameter", "cell diameters"))
        })
        
        Leach = anylapply(cellPatterns, Lest, r = seq(0, 6, by=0.1), ratio = TRUE, correction = "best")
        Lpool = pool(Leach)
        
        hyperMC = rbind(hyperMC, 
                        tibble(group = groups[gr], 
                                   simulation = simulation, 
                                   delta = 0L,
                                   Leach = Leach,
                                   L = list(Lpool)))
        
        # Other deltas
        for (delta in deltas) {
            
            cells.delta = map_dfr(availableTimeSteps[1:(timeStepCount-delta)], function(image) {
                cells %>% filter(timeStep == image) %>% bind_rows(
                    cells %>% filter(timeStep == image + delta) %>% mutate(type = "MCfuture")) %>% 
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
                cells.delta %>% filter(timeStep == image) %>% 
                    {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
                    rescale(s=cells.diameter, unitname=c("cell diameter", "cell diameters"))
            })
            
            Leach = lapply(plist.delta, Lcross, r=seq(0, 6, by=0.1), i="MC", j="MCfuture", ratio=TRUE, correction="best")
            Lpool = pool(as.anylist(Leach))
            
            hyperMC = rbind(hyperMC, 
                            tibble(group = groups[gr], 
                                   simulation = simulation, 
                                   delta = delta,
                                   Leach = Leach,
                                   L = list(Lpool)))
            
        }
            
    }

}

save(hyperMC, file = "pooledFV.RData")

# Studentized permutation test
spstat = function(hf, delay, group1, group2) {
    Leach = hf %>% filter(delta == delay, group %in% c(group1, group2)) %>% 
        {anylapply(split(.$L, .$group), collapse.fv, same="pooltheo", different="pooliso")}
    Ls = anylapply(Leach, function(fv) fv %>% as.matrix() %>% {.[,-1:-2]} )
    Lmean = anylapply(Ls, function(mat) apply(mat, 1, mean))
    Lvar = anylapply(Ls, function(mat) apply(mat, 1, var))
    T = (Leach[[1]]$r[2] - Leach[[1]]$r[1]) * 
        sum( (Lmean[[1]] - Lmean[[2]])^2 / (Lvar[[1]] / ncol(Ls[[1]]) + Lvar[[2]] / ncol(Ls[[2]])) , na.rm = TRUE)
    return(T)
}

sptest = function(hf, delays = 0L:5L, group1, group2, nperm = 999) {
    spvalues = tibble(delta = integer(0), pval = numeric(0))
    for (delay in delays) {
        cat(delay, ". ", sep = "")
        Tvalues = c()
        for (k in 1:nperm) {
            sampled = hf %>% filter(delta == delay)
            sampled$L = sampled$L[sample(nrow(sampled), nrow(sampled))]
            Tvalues[k] = spstat(sampled, delay, group1, group2)
        }
        
        spvalues = bind_rows(spvalues, 
                             tibble(delta = delay, 
                                    pval = 1 - sum(Tvalues < spstat(hf, delay, group1, group2)) / 1000))
    }
    return(spvalues)
}

spval = map_dfr("2B", function(group2) {
    sptest(hyperMC, group1 = "2A", group2 = group2) %>% 
        mutate(comparison = paste("2A", group2, sep = "-"))
})

Ltib = function(delay) {
    hf = hyperMC %>% filter(delta == delay)
    map_dfr(1:nrow(hf), function(index) {
        hf$L[[index]] %>% as_tibble() %>% mutate(sim = index, group = hf$group[index], delta = delay)
    })
}

Ltibble = map_dfr(0L:5L, Ltib)

ggplot(Ltibble, aes(x=r, y=pooliso-r, colour=group, group=sim)) +
    geom_line() +
    ylab(TeX("$L(r) - r$")) +
    scale_color_manual(name = "Modèles", 
                       values = gg_color_hue(5)) +
    facet_wrap(~ delta) +
    theme_bw()

pdf("comparison-with-without-LI-pval.pdf", height=2, width=1.5)
grid.table(spval %>% mutate(pval = round(pval, digits = 3)) %>% spread(comparison, pval))
dev.off()

ggsave("comparison-with-without-LI.eps", device="ps", width=12, height=8, units="in")

# Plots of the L-function assessing the interaction between MCs 
# with and without aNP-driven LI
# Fig.7B
LpoolTib = map_dfr(0L:5L, function(delay) { 
    map_dfr(groups, function(mgroup) {
        hf = hyperMC %>% filter(delta == delay, group == mgroup)
        Lpool = pool(hf$Leach)
        tibble(group = mgroup, delta = delay, L = list(Lpool))
    })
})

LpoolUnnest = LpoolTib %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

ggplot(LpoolUnnest, 
       aes(x=r, y=pooliso-r, ymin=loiso-r, ymax=hiiso-r, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$L(r) - r$")) +
    facet_wrap(~ delta) +
    theme_bw()

ggsave("comparison-with-without-LI-ribbon.eps", device=cairo_ps, width=14, height=8, units="in")
