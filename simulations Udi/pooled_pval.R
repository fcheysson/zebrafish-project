# Pooled p-values for the analysis of MCs in the simulations from the NSC lattice model
# Run 'dynamic_consecutive3timesteps.R' first
# Supp. Fig.13

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

ptab = tibble(simulation = integer(0), delta = integer(0), statistic = character(0), pval = numeric(0))

file = "Positions_withLI_f059_simulation"
simulations = 1L:34L

for (simulation in simulations) {
    load(paste0(file, "/", file, simulation, "/pval.RData"))
    ptab = ptab %>% bind_rows(
        pval %>% select(delta, matches("L\\.\\w*")) %>% 
            gather("statistic", "pval", matches("L\\.\\w*")) %>% 
            mutate(simulation = simulation)
    )
}

Ltab = ptab %>% filter(statistic == "L.simult") %>%
    spread(simulation, pval) %>% 
    mutate(pooled = 1 - pchisq(
        q = -2 * rowSums(log(select(., matches('\\d+')))), 
        df = 2 * length(simulations))) %>% 
    gather("simulation", "pval", `1`:pooled) %>% 
    spread(delta, pval) %>% 
    arrange(match(simulation, c(1:19, "pooled")))

ptab = ptab %>% mutate(log.pval = log(pval)) %>% 
    group_by(delta, statistic) %>% 
    summarise(sum.pval = sum(log.pval)) %>% 
    mutate(pooled = 1 - pchisq(-2*sum.pval, df = 2*length(unique(ptab$simulation)))) %>% 
    select(-sum.pval) %>% 
    spread(statistic, pooled)

pdf(paste0("pval_", file, ".pdf"), height=2, width=8)
gridExtra::grid.table(
    ptab %>%
        mutate_at(vars(starts_with("L")), ~ format(round(., 3), nsmall = 3)) %>% 
        select(delta, L.simult, matches("^L\\.pointwise[1-3]$"), matches("^L\\.dclf[3-6]$"))
)
dev.off()

pdf(paste0("pval_Lsimult_", file, ".pdf"), height=6, width=4.5)
gridExtra::grid.table(
    Ltab %>%
        mutate_at(vars(matches("\\d+")), ~ format(round(., 3), nsmall = 3)) %>% 
        select(-statistic)
)
dev.off()