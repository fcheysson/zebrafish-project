# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)
library(clue)

insertRow <- function(existingDF, newrow, r) {
    if (r < nrow(existingDF))
        existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
}

cells1 = read.table("qNSC.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type="qNSC")
cells2 = read.table("aNSC_nodb.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type="aNSC_nodb")
cells3 = read.table("aNP.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type="aNP")

# Replace doublets by mean singlets
cellsA = read.table("out.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type = "aNSC_db")
distMat = pairdist(cellsA[, c('x', 'y')])
diag(distMat) = rep(1e6, nrow(cellsA))
doublets = as.numeric(solve_LSAP(distMat))
cellsB = cellsA %>%
    mutate(x = .5 * (x + x[doublets]), 
           y = .5 * (y + y[doublets]), 
           z = .5 * (z + z[doublets])) %>% 
    distinct() %>% 
    mutate(type = "aNSC_db2sg")

cells = cells1 %>% bind_rows(cells2) %>% bind_rows(cells3) %>% 
    bind_rows(cellsA) %>% bind_rows(cellsB) %>% 
    mutate(type=factor(type, levels = c("qNSC", "aNSC_nodb", "aNSC_db", "aNSC_db2sg", "aNP")))
geomCells = cells %>% filter(type %in% c("qNSC", "aNSC_nodb", "aNSC_db2sg", "aNP")) %>% 
    mutate(type = factor(ifelse(type == "aNP", "aNP", ifelse(type == "qNSC", "qNSC", "aNSC")))) 

# plot_ly(x=cells$x, y=cells$y, z=cells$z, type="scatter3d", mode="markers", color=cells$type) %>% 
#     layout(scene=list(aspectmode="data"))
# plot_ly(x=cells$x, y=cells$y, type="scatter", mode="markers", color=cells$type) %>% 
#     layout(scene=list(aspectmode="data"))

n1 = nrow(cells1)
n2 = nrow(cells2) + nrow(cellsB)
n3 = nrow(cells3)
n2a = nrow(cells2)
n2b = nrow(cellsB)

xrange = range(cells$x)
yrange = range(cells$y)
zrange = range(cells$z)

p = ppp(geomCells$x, geomCells$y, xrange, yrange, marks=geomCells$type)
Window(p) = ripras(p)
unitname(p) = "µm"

bdry = Window(p)$bdry[[1]] %>% as.data.frame()
bdry = insertRow(bdry, c('x'=164, 'y'=327), 11)
bdry = insertRow(bdry, c('x'=131, 'y'=320), 12)
bdry = insertRow(bdry, c('x'=82, 'y'=182), 16)
bdry = insertRow(bdry, c('x'=100, 'y'=108), 17)
bdry = insertRow(bdry, c('x'=116, 'y'=58), 18)
bdry = insertRow(bdry, c('x'=110, 'y'=-14), 19)
bdry = insertRow(bdry, c('x'=152, 'y'=-67), 22)
bdry = insertRow(bdry, c('x'=194, 'y'=-106), 23)
bdry = insertRow(bdry, c('x'=278, 'y'=-162), 31)
bdry = insertRow(bdry, c('x'=277, 'y'=-89), 32)
bdry = insertRow(bdry, c('x'=289, 'y'=-36), 33)
bdry = insertRow(bdry, c('x'=313, 'y'=92), 34)
bdry = insertRow(bdry, c('x'=318, 'y'=181), 35)
# cells[which( abs(cells$x - 350) < 5 & abs(cells$y - 100) < 5 ), c("x", "y")]

# which(!inside.owin(bdry$x, bdry$y, cells.owin))
# bdry[11,] = bdry[11,] + c(-1e-6, +1e-6) # Prevents one point from getting outed by the window

cells.owin = bdry %>% as.list() %>% owin(poly=.)
Window(p) = cells.owin

# cells.density = density(p)
# plot(cells.density)
# plot(p, add=TRUE, cex=.5)

# ggplot(cells, aes(x=x, y=y)) +
#     geom_point(aes(colour=type)) +
#     # geom_path(data=bind_rows(bdry, bdry[1,])) +
#     geom_point(data=bdry, shape=2) +
#     geom_path(data=bdry) +
#     scale_x_continuous(minor_breaks = seq(-600, 600, by=10), breaks = seq(-600, 600, by=50)) +
#     scale_y_continuous(minor_breaks = seq(-600, 600, by=10), breaks = seq(-600, 600, by=50)) +
#     theme(panel.grid.major = element_line(size=1.5)) +
#     theme_grey()
