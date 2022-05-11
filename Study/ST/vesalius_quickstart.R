suppressMessages(library(vesalius))
data(vesalius,package = "vesalius")
# First create a copy of raw data for later use
# For sake of simplicity we seperate counts from Vesalius images.
counts <- vesalius

# Preprocessing Spatial transcriptomic data with Seurat
vesalius <- NormalizeData(vesalius,verbose =F)
vesalius <- FindVariableFeatures(vesalius, nfeatures = 2000,verbose =F)
vesalius <- ScaleData(vesalius,verbose =F)
# PCA
vesaliusPCA <- rgbPCA(vesalius,slices = 1,verbose =F)

# UMAP
vesaliusUMAP <- rgbUMAP(vesalius,pcs = 30,verbose =F)

vesalius <- buildImageArray(vesaliusPCA, sliceID = 1, invert =T,verbose =F)


# Let's have a look shall we?
Vesalius_image <- imagePlot(vesalius,as.cimg =FALSE) + theme_void()


print(Vesalius_image)

# Histogram EQ
vesalius <- equalizeHistogram(vesalius,sleft = 2.5, sright=2.5,invert =T,verbose =F)

# Variance regularisation
vesalius <- regulariseImage(vesalius, lambda = 10,
                            niter = 200, normalise=T,verbose =F)

# Smoothing
# !!!Optional!!! -
vesalius <- smoothArray(vesalius,method = "box", box = 1,verbose =F)

# Segmentation
# This function provides smoothing internally.
vesalius <- iterativeSegmentation.array(vesalius,
                                        colDepth = 6,
                                        smoothIter = 20,
                                        method = c("iso","median"),
                                        sigma=1.5,
                                        box = 10,
                                        useCenter = F,
                                        invert =T,
                                        verbose =F)

# Isolating territories from colour segments
vesalius <- isolateTerritories.array(vesalius,
                                     captureRadius = 0.1,
                                     minBar = 10,
                                     verbose =F)
# Let's have a look!
imgTerritory <- territoryPlot(vesalius, randomise = TRUE,cex =15 , cex.pt=3.5)
print(imgTerritory)

# extract terriotry from seurat object
ter <- extractTerritories(vesalius,counts,seedID= 1,verbose =F)

# get all markers related to all territories
all <- extractAllMarkers(vesalius,counts,verbose =F)


# get markers associated with on or more territories
one <- extractMarkers(vesalius, counts, seed = c(1),verbose =F)

# comparing territories
two <- extractMarkers(vesalius, counts, seed = c(1), query = c(2,3),verbose =F)
