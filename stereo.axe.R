library(here)
library(StereoMorph)
library(geomorph)
library(tidyverse)
library(wesanderson)

# read data and define number of sLMs
shapes <- readShapes("shapes")
coords<- readland.shapes(shapes, nCurvePts = c(5,10,10,10,3,5,5))

# read qualitative data
qdata <- read.csv("qdata.csv", header = TRUE, row.names = 1)

# gpa ----
Y.gpa <- gpagen(coords, print.progress = FALSE)
plot(Y.gpa)

# geomorph data frame ----
gdf <- geomorph.data.frame(shape = Y.gpa$coords, 
                           size = Y.gpa$Csize,
                           mark = qdata$mark)

# add centroid size to qdata ----
qdata$csz <- Y.gpa$Csize

# attributes for boxplots ----
csz <- qdata$csz # centroid size
cask <- qdata$cask # cask
mark <- qdata$mark  # mark

# boxplot of axe centroid size by mark ----
csz.mark <- ggplot(qdata, aes(x = mark, y = csz, color = mark)) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3) +
  scale_colour_manual(values = wes_palette("Moonrise2")) +
  theme(legend.position = "none") +
  labs(x = 'Cask', y = 'Centroid Size')
# render plot
csz.mark

# principal components analysis ----
pca<-gm.prcomp(Y.gpa$coords)
summary(pca)

# models ----
## general allometry ----
fit.size <- procD.lm(shape ~ size, 
                     data = gdf, 
                     print.progress = FALSE, 
                     iter = 9999)

# general allometry - does french trade axe shape change with size? 
anova(fit.size)

# hypothesis.test ----

# set plot parameters to plot by mark
pch.gps <- c(15,17)[as.factor(mark)]
col.gps <- wes_palette("Moonrise2")[as.factor(mark)]
col.hull <- c("#798E87","#C27D38")

# plot pca by mark
pc.plot1 <- plot(pca, 
                 asp = 1,
                 pch = pch.gps,
                 col = col.gps)
shapeHulls(pc.plot1, 
           groups = mark,
           group.cols = col.hull)

## size as a function of mark ----
fit.size.mark <- procD.lm(size ~ mark, 
                          data = gdf, 
                          print.progress = FALSE, 
                          iter = 9999)

# differences in size by mark?
anova(fit.size.mark)

## shape as a function of mark ----
fit.shape.mark <- procD.lm(shape ~ mark, 
                           data = gdf, 
                           print.progress = FALSE, 
                           iter = 9999)

## differences in shape by mark? ----
anova(fit.shape.mark)


# mean shapes ----
new.coords<-coords.subset(A = Y.gpa$coords, 
                          group = qdata$mark)
names(new.coords)

# group shape means
mean <- lapply(new.coords, mshape)

# plot mean shapes
plot(mean$asterisk)
plot(mean$DG)

# comparison plots
plotRefToTarget(mean$asterisk, 
                mean$DG, 
                method = "vector",
                mag = 2)

# end of code ----