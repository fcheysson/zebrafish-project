# Plots of the L-function assessing the interaction between MCs 
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
simulations = list(1L:34L, 1L:33L)

hyperMC = tibble(group = factor(character(0L), levels = groups),
                 simulation = integer(0L), 
                 delta = integer(0L),
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

        hyperMC = rbind(hyperMC, 
                        tibble(group = groups[gr], 
                               simulation = simulation, 
                               delta = 0L,
                               L = Leach))
        
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

            hyperMC = rbind(hyperMC, 
                            tibble(group = groups[gr], 
                                   simulation = simulation, 
                                   delta = delta,
                                   L = Leach))
            
        }
        
    }
    
}

LpoolTib = map_dfr(0L:5L, function(delay) { 
    map_dfr(groups, function(mgroup) {
        hf = hyperMC %>% filter(delta == delay, group == mgroup)
        Lpool = pool(hf$L)
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
