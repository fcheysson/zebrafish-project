# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)
library(FactoMineR)
library(factoextra)

insertRow <- function(existingDF, newrow, r) {
    if (r < nrow(existingDF))
        existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
}

# GGPLOT COLOR PALETTE
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

# IMPORT DATA
cells = map_dfr(1:8, function(image) {
    
    cells1 = read.table(paste0("mimi_", image, "_qNSC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
        as_tibble() %>% mutate(type="qNSC")
    cells2 = read.table(paste0("mimi_", image, "_aNSC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
        as_tibble() %>% mutate(type="aNSC")
    cells3 = read.table(paste0("mimi_", image, "_MC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
        as_tibble() %>% mutate(type="MC")
    
    cells1 %>% bind_rows(cells2) %>% bind_rows(cells3) %>% mutate(type = as.factor(type), cliche = image)
    
})

# FIND PLAN OF MOST VARIATION
cells.pca = PCA(cells[,1:3], scale.unit = FALSE)
cells = cells %>% mutate(dim1 = cells.pca$ind$coord[,1], 
                         dim2 = cells.pca$ind$coord[,2], 
                         dim3 = cells.pca$ind$coord[,3])

# Remove duplicated points (aNSC -> MC)
cells = cells %>% select(cliche, dim1, dim2) %>% duplicated(fromLast = TRUE) %>% {cells[!.,]}

# plot_ly(x=cells$x, y=cells$y, z=cells$z, type="scatter3d", mode="markers", color=cells$type) %>%
#     layout(scene=list(aspectmode="data"))
# plot_ly(x=cells$dim1, y=cells$dim2, type="scatter", mode="markers", color=cells$type) %>%
#     layout(scene=list(aspectmode="data"))

# CREATE APPROPRIATE BOUNDARY
xrange = range(cells$dim1)
yrange = range(cells$dim2)
zrange = range(cells$dim3)

p = ppp(cells$dim1, cells$dim2, xrange, yrange, marks=cells$type)
Window(p) = ripras(p, f = 1 / sqrt(1 - 240 / nrow(cells)))
unitname(p) = "µm"

bdry = Window(p)$bdry[[1]] %>% as.data.frame()
bdry = insertRow(bdry, c('x'=219, 'y'=60), 7)
bdry = insertRow(bdry, c('x'=201, 'y'=53), 8)
bdry = insertRow(bdry, c('x'=162, 'y'=58), 9)
bdry = insertRow(bdry, c('x'=113, 'y'=68), 11)
bdry = insertRow(bdry, c('x'=94, 'y'=47), 12)
bdry = insertRow(bdry, c('x'=68, 'y'=60), 13)
bdry = insertRow(bdry, c('x'=52, 'y'=61), 14)
bdry = insertRow(bdry, c('x'=-28, 'y'=56), 15)
bdry = insertRow(bdry, c('x'=-85, 'y'=47), 16)
bdry = insertRow(bdry, c('x'=-106, 'y'=55), 17)
bdry = insertRow(bdry, c('x'=-209, 'y'=-1), 23)
bdry[26,] = bdry[26,] + c(4, -1) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=-169, 'y'=-40), 27)
bdry = insertRow(bdry, c('x'=-158, 'y'=-55), 28)
bdry[30,] = bdry[30,] + c(6, 0) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=-18, 'y'=-50), 31)
bdry = insertRow(bdry, c('x'=33, 'y'=-58), 32)
bdry = insertRow(bdry, c('x'=61, 'y'=-61), 33)
bdry = insertRow(bdry, c('x'=92, 'y'=-54), 34)
bdry = insertRow(bdry, c('x'=160, 'y'=-60), 35)
bdry[36,] = bdry[36,] + c(-6, -1) # Prevents one point from getting outed by the window
# which(!inside.owin(cells$x, cells$y, cells.owin))

cells.owin = bdry %>% as.list() %>% owin(poly=.)
Window(p) = cells.owin

# ggplot(cells) +
#     geom_point(aes(x=dim1, y=dim2, colour=type)) +
#     # geom_path(data=bind_rows(bdry, bdry[1,])) +
#     geom_point(aes(x=x, y=y), data=bdry, shape=2) +
#     geom_path(aes(x=x, y=y), data=bdry) +
#     scale_x_continuous(minor_breaks = seq(-300, 300, by=10), breaks = seq(-300, 300, by=50)) +
#     scale_y_continuous(minor_breaks = seq(-300, 300, by=10), breaks = seq(-300, 300, by=50)) +
#     theme(panel.grid.major = element_line(size=1.5)) +
#     theme_grey()

# Rescale to account for inhomogeneity
rho = rhohat(p, "x")
transf = rho %>% mutate(int = cumsum(rho)) %>% {approxfun(.$x, .$int)}

cells = cells %>% mutate(dim1 = transf(dim1) * (xrange[2] - xrange[1]) / transf(xrange[2]))

# CREATE APPROPRIATE BOUNDARY
xrange = range(cells$dim1)
yrange = range(cells$dim2)
zrange = range(cells$dim3)

p = ppp(cells$dim1, cells$dim2, xrange, yrange, marks=cells$type)
Window(p) = ripras(p, f = 1 / sqrt(1 - 240 / nrow(cells)))
unitname(p) = "µm"

bdry = Window(p)$bdry[[1]] %>% as.data.frame()
bdry = insertRow(bdry, c('x'=463, 'y'=60), 8)
bdry = insertRow(bdry, c('x'=420, 'y'=60), 9)
bdry = insertRow(bdry, c('x'=395, 'y'=75), 11)
bdry = insertRow(bdry, c('x'=379, 'y'=68), 12)
bdry = insertRow(bdry, c('x'=360, 'y'=49), 13)
bdry = insertRow(bdry, c('x'=338, 'y'=62), 14)
bdry = insertRow(bdry, c('x'=257, 'y'=56), 15)
bdry = insertRow(bdry, c('x'=185, 'y'=47), 16)
bdry = insertRow(bdry, c('x'=158, 'y'=54), 17)
bdry = insertRow(bdry, c('x'=22, 'y'=18), 23)
bdry = insertRow(bdry, c('x'=23, 'y'=-7), 24)
bdry = insertRow(bdry, c('x'=27, 'y'=-44), 28)
bdry = insertRow(bdry, c('x'=73, 'y'=-41), 29)
bdry = insertRow(bdry, c('x'=96, 'y'=-57), 30)
bdry = insertRow(bdry, c('x'=146, 'y'=-63), 33)
bdry = insertRow(bdry, c('x'=242, 'y'=-52), 34)
bdry = insertRow(bdry, c('x'=340, 'y'=-61), 35)
bdry = insertRow(bdry, c('x'=349, 'y'=-54), 36)
bdry = insertRow(bdry, c('x'=425, 'y'=-63), 37)
bdry = insertRow(bdry, c('x'=437, 'y'=-76), 38)
# which(!inside.owin(cells$x, cells$y, cells.owin))

cells.owin = bdry %>% as.list() %>% owin(poly=.)
Window(p) = cells.owin

# ggplot(cells) +
#     geom_point(aes(x=dim1, y=dim2, colour=type)) +
#     # geom_path(data=bind_rows(bdry, bdry[1,])) +
#     geom_point(aes(x=x, y=y), data=bdry, shape=2) +
#     geom_path(aes(x=x, y=y), data=bdry) +
#     scale_x_continuous(minor_breaks = seq(0, 600, by=10), breaks = seq(0, 600, by=50)) +
#     scale_y_continuous(minor_breaks = seq(-300, 300, by=10), breaks = seq(-300, 300, by=50)) +
#     theme(panel.grid.major = element_line(size=1.5)) +
#     theme_grey()

# CREATE SOLIST
plist = solapply(1:8, function(image) {
    cells %>% filter(cliche == image) %>% {ppp(x = .$dim1, y = .$dim2, window = cells.owin, marks = .$type)}
})

# ggplot(cells, aes(x=dim1, y=dim2)) +
#     geom_point(aes(colour=type)) +
#     geom_path(aes(x=x, y=y), data=bind_rows(bdry, bdry[1,])) +
#     facet_wrap(~ cliche, nrow=2) +
#     scale_colour_manual(values = c(gg_color_hue(2), "black")) +
#     coord_equal() +
#     theme_grey()