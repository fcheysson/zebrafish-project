# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(spatstat)
library(plotly)
library(dbmss)

insertRow <- function(existingDF, newrow, r) {
    if (r < nrow(existingDF))
        existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
}

cells1 = read.table("qNSC.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type=1)
cells2 = read.table("aNSC_nodb.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type=2)
cells3 = read.table("aNP.csv", header=FALSE, sep=",", col.names=c("x", "y", "z")) %>% 
    as_tibble() %>% mutate(type=3)

cells = cells1 %>% bind_rows(cells2) %>% bind_rows(cells3) %>% mutate(type=as.factor(type))

# plot_ly(x=cells$x, y=cells$y, z=cells$z, type="scatter3d", mode="markers", color=cells$type) %>% 
#     layout(scene=list(aspectmode="data"))
# plot_ly(x=cells$x, y=cells$y, type="scatter", mode="markers", color=cells$type) %>% 
#     layout(scene=list(aspectmode="data"))

n1 = nrow(cells1)
n2 = nrow(cells2)
n3 = nrow(cells3)

xrange = range(cells$x)
yrange = range(cells$y)
zrange = range(cells$z)

p = ppp(cells$x, cells$y, xrange, yrange, marks=cells$type)
Window(p) = ripras(p)
unitname(p) = "µm"

bdry = Window(p)$bdry[[1]] %>% as.data.frame()
bdry = insertRow(bdry, c('x'=406, 'y'=528), 15)
bdry = insertRow(bdry, c('x'=393, 'y'=479), 16)
bdry = insertRow(bdry, c('x'=355, 'y'=443), 17)
bdry = insertRow(bdry, c('x'=355, 'y'=443), 18)
bdry = insertRow(bdry, c('x'=337, 'y'=421), 19)
bdry = insertRow(bdry, c('x'=246, 'y'=308), 20)
bdry = insertRow(bdry, c('x'=191, 'y'=257), 21)
bdry = insertRow(bdry, c('x'=155, 'y'=249), 22)
bdry = insertRow(bdry, c('x'=118, 'y'=230), 23)
bdry = insertRow(bdry, c('x'=333, 'y'=62), 32)
bdry = insertRow(bdry, c('x'=365, 'y'=131), 33)

# cells[which( abs(cells$x - 350) < 5 & abs(cells$y - 100) < 5 ), c("x", "y")]

# which(!inside.owin(cells$x, cells$y, cells.owin))
# bdry[35,] = bdry[35,] + c(+1e-6, -1e-6) # Prevents one point from getting outed by the window

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
#     scale_x_continuous(minor_breaks = seq(50, 600, by=10), breaks = seq(0, 600, by=50)) +
#     scale_y_continuous(minor_breaks = seq(0, 600, by=10), breaks = seq(0, 600, by=50)) +
#     theme(panel.grid.major = element_line(size=1.5)) +
#     theme_grey()
