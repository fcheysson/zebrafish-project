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
                 type = character(0L),
                 Leach = anylist(),
                 L = list())

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
        
        cellWindow = ripras(x = cells$x, y = cells$y, shape = "rectangle", f = (meanPointsPerTimeStep + 1) / (meanPointsPerTimeStep - 1))
        
        # CREATE SOLIST
        availableTimeSteps = cells %>% pull(timeStep) %>% unique() %>% {.[.>=400]}
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
        timeSteps = c(400, 433, 467, 500, 533, 567, 600) %>% 
            sapply(function(ts) {
                if (nrow(filter(cells, timeStep == ts, type == "aNSC")) < 2)
                    return(c())
                return(ts)
            }) %>% unlist()  # Remove timesteps with less than 2 aNSCs
        cells = cells %>% filter(timeStep %in% timeSteps)
        
        # ggplot(cells, aes(x=x, y=y)) +
        #     geom_point(aes(colour=type), size=3) +
        #     geom_hline(yintercept = cellWindow$yrange) +
        #     geom_vline(xintercept = cellWindow$xrange) +
        #     facet_wrap(~ timeStep) +
        #     scale_colour_manual(values = c("grey", gg_color_hue(length(cellTypes)-1))) +
        #     coord_equal() +
        #     theme_grey()
        
        ## Permutation cellules type 2
        plist = timeSteps %>% solapply(function(image) {
            cells  %>% 
                filter(timeStep == image, type == "aNSC") %>% 
                {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        Leach = anylapply(plist, Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
        Lpool = pool(Leach)
        
        hyperMC = rbind(hyperMC, 
                        tibble(group = groups[gr], 
                               simulation = simulation, 
                               type = "22",
                               Leach = Leach,
                               L = list(Lpool)))
        
        ## Permutation cellules type 3
        plist = timeSteps %>% solapply(function(image) {
            cells  %>% 
                filter(timeStep == image, type %in% c("aNP", "div.aNP")) %>% 
                mutate(type = "aNP") %>% 
                {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        Leach = anylapply(plist, Lest, r=seq(0, 6, by=0.1), ratio=TRUE, correction="best")
        Lpool = pool(Leach)
        
        hyperMC = rbind(hyperMC, 
                        tibble(group = groups[gr], 
                               simulation = simulation, 
                               type = "33",
                               Leach = Leach,
                               L = list(Lpool)))
        
        ## Permutation cellules type 32
        plist = timeSteps %>% solapply(function(image) {
            cells  %>% 
                filter(timeStep == image) %>% 
                mutate(type = recode(type, div.aNP = "aNP")) %>% 
                {ppp(x = .$x, y = .$y, window = cellWindow, marks = .$type)} %>% 
                rescale(s=cells.diameter, unitname="cell diameter")
        })
        
        Leach = anylapply(plist, Lcross, r=seq(0, 6, by=0.1), i="aNP", j="aNSC", ratio=TRUE, correction="best")
        Lpool = pool(Leach)
        
        hyperMC = rbind(hyperMC, 
                        tibble(group = groups[gr], 
                               simulation = simulation, 
                               type = "32",
                               Leach = Leach,
                               L = list(Lpool)))
        
    }
    
}

save(hyperMC, file = "pooledFV.RData")

# Studentized permutation test
spstat = function(hf, test, group1, group2) {
    Leach = hf %>% filter(type == test, group %in% c(group1, group2)) %>% 
        {anylapply(split(.$L, .$group), collapse.fv, same="pooltheo", different="pooliso")}
    Ls = anylapply(Leach, function(fv) fv %>% as.matrix() %>% {.[,-1:-2]} )
    Lmean = anylapply(Ls, function(mat) apply(mat, 1, mean))
    Lvar = anylapply(Ls, function(mat) apply(mat, 1, var))
    T = (Leach[[1]]$r[2] - Leach[[1]]$r[1]) * 
        sum( (Lmean[[1]] - Lmean[[2]])^2 / (Lvar[[1]] / ncol(Ls[[1]]) + Lvar[[2]] / ncol(Ls[[2]])) , na.rm = TRUE)
    return(T)
}

sptest = function(hf, tests = c("22", "33", "32"), group1, group2, nperm = 999) {
    spvalues = tibble(type = character(0), pval = numeric(0))
    for (test in tests) {
        cat(test, ". ", sep = "")
        Tvalues = c()
        for (k in 1:nperm) {
            sampled = hf %>% filter(type == test)
            sampled$L = sampled$L[sample(nrow(sampled), nrow(sampled))]
            Tvalues[k] = spstat(sampled, test, group1, group2)
        }
        
        spvalues = bind_rows(spvalues, 
                             tibble(type = test, 
                                    pval = 1 - sum(Tvalues < spstat(hf, test, group1, group2)) / 1000))
    }
    return(spvalues)
}

spval = map_dfr("2B", function(group2) {
    sptest(hyperMC, group1 = "2A", group2 = group2) %>% 
        mutate(comparison = paste("2A", group2, sep = "-"))
})

Ltib = function(test) {
    hf = hyperMC %>% filter(type == test)
    map_dfr(1:nrow(hf), function(index) {
        hf$L[[index]] %>% as_tibble() %>% mutate(sim = index, group = hf$group[index], type = test)
    })
}

Ltibble = map_dfr(c("22", "33", "32"), Ltib)

ggplot(Ltibble, aes(x=r, y=pooliso-r, colour=group, group=sim)) +
    geom_line() +
    ylab(TeX("$L(r) - r$")) +
    scale_color_manual(name = "Modèles", 
                       values = gg_color_hue(5)) +
    facet_wrap(~ type) +
    theme_bw()

pdf("comparison-with-without-LI-pval.pdf", height=2, width=1.5)
grid.table(spval %>% mutate(pval = round(pval, digits = 3)) %>% spread(comparison, pval))
dev.off()

ggsave("comparison-with-without-LI.eps", device="ps", width=12, height=8, units="in")

# Plots of the L-function assessing the interaction between cells 
# with and without aNP-driven LI (pooled L ribbon graphs)
# Fig.7B
LpoolTib = map_dfr(c("22", "33", "32"), function(test) { 
    map_dfr(groups, function(mgroup) {
        hf = hyperMC %>% filter(type == test, group == mgroup)
        Lpool = pool(hf$Leach)
        tibble(group = mgroup, type = test, L = list(Lpool))
    })
})

LpoolUnnest = LpoolTib %>% mutate(L = lapply(L, as.data.frame)) %>% unnest()

ggplot(LpoolUnnest, 
       aes(x=r, y=pooliso-r, ymin=loiso-r, ymax=hiiso-r, fill=group)) +
    geom_line(aes(colour=group)) +
    geom_ribbon(alpha=.4) +
    ylab(TeX("$L(r) - r$")) +
    facet_wrap(~ type) +
    theme_bw()

ggsave("comparison-with-without-LI-ribbon.eps", device=cairo_ps, width=14, height=8, units="in")
