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
    
    cells1 = read.table(paste0("titi_", image, "_qNSC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
        as_tibble() %>% mutate(type="qNSC")
    cells2 = read.table(paste0("titi_", image, "_aNSC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
        as_tibble() %>% mutate(type="aNSC")
    cells3 = read.table(paste0("titi_", image, "_MC.csv"), header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
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
bdry = insertRow(bdry, c('x'=92, 'y'=68), 6)
bdry = insertRow(bdry, c('x'=7, 'y'=63), 7)
bdry = insertRow(bdry, c('x'=-22, 'y'=53), 8)
bdry = insertRow(bdry, c('x'=-86, 'y'=61), 9)
bdry[10,] = bdry[10,] + c(6, 0) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=-126, 'y'=63), 11)
bdry = insertRow(bdry, c('x'=-140, 'y'=50), 12)
bdry[20,] = bdry[20,] + c(8, -1) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=-70, 'y'=-53), 21)
bdry = insertRow(bdry, c('x'=2, 'y'=-68), 22)
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
transf = rho %>% as_tibble() %>% mutate(int = cumsum(rho)) %>% {approxfun(.$x, .$int)}

cells = cells %>% mutate(dim1 = transf(dim1) * (xrange[2] - xrange[1]) / transf(xrange[2]))

# CREATE APPROPRIATE BOUNDARY
xrange = range(cells$dim1)
yrange = range(cells$dim2)
zrange = range(cells$dim3)

p = ppp(cells$dim1, cells$dim2, xrange, yrange, marks=cells$type)
Window(p) = ripras(p, f = 1 / sqrt(1 - 240 / nrow(cells)))
unitname(p) = "µm"

bdry = Window(p)$bdry[[1]] %>% as.data.frame()
bdry = insertRow(bdry, c('x'=198, 'y'=58), 6)
bdry[7,] = bdry[7,] + c(6, 2) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=113, 'y'=65), 8)
bdry = insertRow(bdry, c('x'=106, 'y'=51), 9)
bdry = insertRow(bdry, c('x'=14, 'y'=47), 10)
bdry[18,] = bdry[18,] + c(10, -3) # Prevents one point from getting outed by the window
bdry = insertRow(bdry, c('x'=195, 'y'=-48), 19)
bdry = insertRow(bdry, c('x'=257, 'y'=-69), 20)
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