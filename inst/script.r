

### Read data and plot

set.cty.seas()
# We assume that the working directory is in folder "../inst" relative to "R/Data/"
xAbs = x.read(fn.path = "../R/Data/", 2017:2019)
xAvl = x.normalize(xAbs, 2017:2019, "../R/SolarWindcapacity.csv")

# plot time series of wind (DE: off-, onshore), solar, and optionally load
x.plot(xAvl,"CH","2017") 
x.plot(xAvl,"AT","2017")
x.plot(xAvl,"DE","2017")
x.plot(xAvl,"FR","2017")
x.plot(xAvl,"IT","2017")

x = sel.wind(xAvl, "DE", windType = 3) # 3 is total wind (on- and offshore)
colnames(x)

cor(x)
x.cor(x) # Correlations for each year separately

plot24(x, "ATsolar")
plot24(x, "DEsolar")
plot24(x, "CHsolar")
plot24(x, "ITsolar")
plot24(x, "FRsolar")
plot24(x, "DEwind", xpd=TRUE) # wind needs more plot space
plot24(x, "ATwind", xpd=TRUE)
plot24(x, "CHwind", xpd=TRUE)
plot24(x, "ITwind", xpd=TRUE)
plot24(x, "FRwind", xpd=TRUE)


l.wide24 = x.24wide(x, yr='2017/2019')
# l.wide24 is a list over seasons. List elements: data-frame of realizations in a wide format: Hours of days are in columns

# List over seasons of stacked data of all countries for external use (e.g. other clustering):

l.All = l.stacked.All(x)
for (s in seasNames) write.csv(l.All[[s]], paste0("SolarWindAvlCountry_",s,".csv"))

contour.wide24(l.wide24, "DE", "SU", "2017-2019")
contour.wide24(l.wide24
               ,cty = "DE|IT"
               ,seas= c("SU","FA")
               ,yrRange = "2017-2019"
               ,upper.triangle = FALSE
               ,subrange = FALSE
               ,old.version = FALSE)




write.csv(agggregate.avl.BEM(l.wide24), "avlBEM.csv")




## Hierarchical Clustering


scnBEM = scn.clust(l.wide24, nclus=20, cty = "DE", s = "SU", interactive = TRUE)
# scnBEM is a data.frame of scenarios in stacked format (suitable for import in electricity market model BEM)

## Cross-country scenarios using hierarchical clustering
scnAll = scn.clust.cross(l.wide24, nclust=20)
write.csv(scnAll, "clustAllCty.csv") 



### Generae senarios with PCA factor model



## PCA for country and season
fac.pcaloads = PCAfac(l.wide24, "DE","SU")
fac.pcaloads = PCAfac(l.wide24, "IT","SU")
fac.pcaloads = PCAfac(l.wide24, "DE|IT","SU")
fac.pcaloads = PCAfac(l.wide24, "DE","WI|SP|SU|FA")

fac = fac.pcaloads[["fac"]] # series of factor model of PCA
pcaloads = fac.pcaloads[["pcaloads"]] # Loadings of PCA
m.sw =  fac.pcaloads[["m.sw"]]  # shift from de-meaning

# Fit the distributions of the factors:
fitParam = facAnalysis(fac)

## Generate scenarios for the country and season
scnData = scn.PCAfac(fitParam, pcaloads, m.sw)
plotScnData(scnData, "DE", "SU")

## Generate scenarios over all countries and all seasons
J.cases = list(c(5,2,1))
scnAll = scn.All(l.wide24, J.cases)
## Sensitivity analysis for single country
J.cases = list(c(5,0), c(5,2,0), c(5,2,1))
scnAllComp = scn.All(l.wide24, J.cases, compareCty = "DE", comparePCA = TRUE)


writeScnAll(scnAll) #write scenarios into files

plotScn(scnAll, "DE", "SU")



### Quality check of scenarios:



## Differences in correlations:


cov.H = plotCor.Emp(l.wide24, "DE", "SU") # plot empirical correlations
# generate cluster scenarios and plot correlation
cov.C = plotCor.genClust(l.wide24, "DE", "SU", n=20) 
# generate PCA scenarios and plot correlations
cov.P = plotCor.genPCA(l.wide24,"DE","WI",J = c(2,2,1,0), textArg = "Scenarios")

## Error in correlations to empirical distribution
## under different matrix norms 
err = c(
   PCA = norm((cov.H-cov.P)/4,"F"),
   Clu = norm((cov.H-cov.C)/4,"F"),
   PCA.x = norm(cov.H[25:48,1:24]-cov.P[25:48,1:24],"F"),
   Clu.x = norm(cov.H[25:48,1:24]-cov.C[25:48,1:24],"F"),
   PCA.O = norm((cov.H-cov.P)/4,"O"),
   Clu.O = norm((cov.H-cov.C)/4,"O"),
   PCA.xO = norm(cov.H[25:48,1:24]-cov.P[25:48,1:24],"O"),
   Clu.xO = norm(cov.H[25:48,1:24]-cov.C[25:48,1:24],"O")
  )
print(signif(err,2))

diffcorr.P = cov2cor(cov.H) - cov2cor(cov.P) # weighted.corr$cor
diffcorr.C = cov2cor(cov.H) - cov2cor(cov.C) # weighted.corr$cor
diffcorr.P[is.na(diffcorr.P)]=0
diffcorr.C[is.na(diffcorr.C)]=0
windows()
lbl = paste0("DE", ", ", "SU", ", Difference PCA")
contour.plot(diffcorr.P[1:24,25:48], lbl, "DE", "SU", triang=F, oldVer=F)
#dev.print(pdf, paste0("corDiffPCAContour_",cty,"_",s,"_rev1.pdf"))
windows()
lbl = paste0("DE", ", ", "SU", ", Difference Clustering")
contour.plot(diffcorr.C[1:24,25:48], lbl, cty, s, triang=F, oldVer=F)
#dev.print(pdf, paste0("corDiffClusteringContour_",cty,"_",s,"_rev1.pdf"))



## Differences in ECDF:


plotECDF(scnAll, "DE", "SU") # for summer only
plotECDF(scnAll, "DE") # for all seasons
# Plot sensitivity (works for single country and season only):
plotECDF(scnAllComp, "DE", "SU", comparePCA = TRUE,
         J.cases = list(c(5,0), c(5,2,0), c(5,2,1)))



## Differences in availabilities (analysis is for region DE only):


scnBEM = scnAll[["DE"]]
avlPCA = avl.PCA(scnBEM) # Seasonal availability factors implied by PCA
# yearly availability factor of scenarios:
mean(avlPCA[,"wind","mean"])
mean(avlPCA[,"solar","mean"])
avlEmp = avl.emp(x)
diff.avl(avlPCA, avlEmp) #difference for DE  (hard-coded)
# Correlations of availabilities:
corScn = cov.wt(scnBEM[,c("windAvl", "solarAvl")],
                wt = scnBEM[,"prob"] , cor = TRUE)
corScn$cor
cor(x['2017/2019'])[paste0("DE","wind"),paste0("DE","solar")]





### Cross-regional scenarios using Copulas




c.d = daily.cor.cross(x,"SU")  # summer season
c.d.all = daily.cor.cross(x)  # all seasons
c.d.small = daily.cor.cross.simple(x,"FA")

# Lambda
x.d = x.daily(x,"FA") # get daily time series from hourly
lambda.nonpar(x.d)
lambda.t(x.d) # counterfactual


#Consider only an averaged solar column:
x.d.small = avg.Cols(x.d)

# t-copula
rtC = fit.tC.and.sample(x.d.small, sample.size = 20)
plot(rtC)

flo1 = fitLambda(rtC, method="Schmid.Schmidt")
flo2 = fitLambda(rtC, method="t", verbose=TRUE)
corrplot.mixed(f.display(flo1,names(x.d.small)))
corrplot.mixed(f.display(flo2$Lambda,names(x.d.small))) 


# Gaussian (normal) copula
rnC = fit.nC.and.sample(c.d.small, sample.size = 20)
plot(rnC)


## Cross-regional scenarios for a season:
# input: scnAll (from single region analysis): List of single region scenarios of all countries all seasons
scnCross = scn.cross(scnAll, "SU", rtC, sample.size = 20)
scn.cross.plot(scnCross, "SU", sample.size = 20)

## Cross-regional scenarios over all seasons:
# using random samples from t-copula:
scnCrossAll = scn.cross.year(x, scnAll, sample.size = 20, tCop = TRUE) ### HERE IS STILL AN ERROR
scn.cross.write(scnCrossAll, sample.size = 20, "_tcopula")
##  using random samples from gaussian copula:
scnCrossAll = scn.cross.year(x, scnAll, sample.size = 20, tCop = FALSE)
scn.cross.write(scnCrossAll, sample.size = 20, "_gCopula")




####### Read also the load

# We assume that the working directory is in folder "../inst" relative to "R/Data/"
xAbs = x.read(fn.path = "../R/Data/", 2017:2019, load = TRUE)
xNorm = x.normalize(xAbs, 2017:2019, "../R/SolarWindcapacity.csv", load = TRUE)

# plot time series of wind (DE: off-, onshore), solar, and load
x.plot(xNorm,"CH","2017") 


cor(x)
x.cor(x) # Correlations for each year separately


# Only available if load was imported before:
plot24(x, "ATload")
plot24(x, "DEload")
plot24(x, "CHload")
plot24(x, "FRload")
plot24(x, "ITload")

x = sel.wind(xAvl, "DE", windType = 3) # 3 is total wind (on- and offshore)
colnames(x)

l.wide24 = x.24wide(x, yr='2017/2019')

contour.wide24(l.wide24, "DE", "SU", "2017-2019")

fac.pcaloads = PCAfac(l.wide24, "DE","SU")



#------------------------------------------------
# Comparison: Simple t-copula, dim=3
# -----------------------------------------------

tC1.3 = tCopula(c(-.6,0.75,0), dim=3, dispstr="un")
set.seed(5640)
r.tC1.3 = rCopula(500, tC1.3)
## A t copula with variable df (df.fixed=FALSE)
## Maximum likelihood (start = (rho[1:3], df))
tc.ml = fitCopula(tCopula(dim=3, dispstr="un"),
                  r.tC1.3, method="ml",
                  start = c(0,0,0,10))
summary(tc.ml)

fitLambda(r.tC1.3, method="t", lower.tail = FALSE)
fitLambda(r.tC1.3, method="t")
