# Model code associated to the paper entitled: "Fine-scale interplay between decline and growth determines the spatial recovery of coral communities within a reef" paper 
# Created by Julie Vercelloni and Murray Logan 
# Contact: j.vercelloni@aims.gov.au 

source("R/packages.R")

# Read predictive locations sampled from Heron Reef 
pred_loc <- read.csv("data/pred_locations.csv") 

# Generate synthetic data - beta responses longitudinal data with varying observation and interval times
def <- defData(varname = "xbase", dist = "normal", formula = 20, variance = 3)

#nCount defines the number of measurements for an individual
def <- defData(def, varname = "nCount", dist = "noZeroPoisson", formula = 6)

#mInterval specifies the average time between intervals for a subject
def <- defData(def, varname = "mInterval", dist = "gamma", formula = 2, variance = 0.01)

# Create the dataset
dt <- genData(50, def)
dtPeriod <- addPeriods(dt)

# Add response variable 

def2 <- defDataAdd(varname = "g.mean", formula = "0.1 + 0.3*time",
                   variance = 2.6, dist="beta", link = "logit")

df <- addColumns(def2, dtPeriod)%>%
  mutate(idyear = time +1) %>%
  mutate(idsub = id) %>%
  dplyr::select(idsub,idyear,g.mean)

# Add random locations and associated Geomorphic zones from the predictive locations - for the example only and merge for final dataset

pred_sample <- pred_loc %>% sample_n(length(unique(df$idsub))) %>%
  dplyr::select(Long,Lat,Geomorph) %>%
  mutate(idsub = unique(df$idsub)) 

df <- inner_join(df,pred_sample) %>%
  data.frame()%>%
  mutate(g.mean = ifelse(g.mean == 0, 0.001, g.mean))%>%
  mutate(g.mean = ifelse(g.mean == 1, 0.999, g.mean))

# Create the INLA mesh

coords<-unique(df[c("idsub","Long","Lat")])

max.edge = diff(range(df$Long))/15
bound.outer = diff(range(df$Long))/15

mesh <- inla.mesh.2d(loc = cbind(df$Long, df$Lat),
                     max.edge = c(1,5)*max.edge,
                     cutoff = max.edge/20,
                     offset = NULL) 

plot(mesh, main="", lwd=0.5); points(df$Long, df$Lat, 
  col="red")                     

# Create the spde, index set and projection matrix 
spde = inla.spde2.pcmatern(mesh, alpha=2, prior.range = c(.02, .01), prior.sigma = c(3, 0.01))

timesn <- length(unique(df$idyear))
indexs <- inla.spde.make.index ("s", n.spde = spde$n.spde, n.group = timesn)

A <- inla.spde.make.A(mesh = mesh, loc=cbind(df$Long, df$Lat),group= as.numeric(as.factor(df$idyear)),n.group=timesn) 

# Create the stack

stk.e <- inla.stack(data=list(y=df$g.mean), A=list(A,1,1,1),
                    effects=list(c(indexs, list(b0=1)),
                                 list(Geomorph=df[ ,"Geomorph"]), list(iidx=as.numeric(as.factor(df$idsub))),
                                 list(year=as.numeric(as.factor(df$idyear)))), tag="est")

# Run the INLA model (~8min on local machine)

rprior<- list(theta = list(prior = "pccor1", param = c(0,0.9)))

formula <- y ~ -1 + b0 + Geomorph + f(s, model=spde, group=s.group,
                                      control.group = list(model="ar1",hyper=rprior)) + 
  f(iidx, model = "iid") 

start.time <- Sys.time()

mod <-  inla(formula,data = inla.stack.data(stk.e),
                  family= 'beta',
                  control.predictor = list(compute = TRUE,
                  link = 1,A = inla.stack.A(stk.e)),
                  control.compute = list(return.marginals.predictor=TRUE,
                  config = TRUE, dic= TRUE), 
                  verbose = FALSE)

end.time <- Sys.time()

# Time model took to run 
end.time - start.time

# Model checks 

plot_fixed_marginals(mod, priors=F)
plot_hyper_marginals(mod)

# Model fit (data)

index.est <- inla.stack.index(stack = stk.e,tag="est")$data

df$pred <- mod$summary.fitted.values$mean[index.est]
df$lower <- mod$summary.fitted.values[index.est,"0.025quant"]
df$upper <- mod$summary.fitted.values[index.est,"0.975quant"]

# Retrieve predictions 
# An example of 10 draws from predictive posterior distributions
# Use 2000 draws in the real study (~30 min on a local machine)

draws <- inla.posterior.sample(10, result=mod, add.names=FALSE)

# Retrieve the spatial.fields posteriors

proj.grid <- inla.mesh.projector(mesh, loc=as.matrix(cbind(pred_loc$Long, pred_loc$Lat)))
cellmeans = sapply(draws, function(x) x[['latent']])

i.mod <- lapply(c('APredictor','^Predictor','s','Geomorph','iidx', 'b'),
                function(x) grep(x, draws[[1]]$latent %>% rownames))

cellmeans.full <- cellmeans[i.mod[[3]],] %>%          
  as.data.frame %>%                                 
  mutate(idyear = rep(as.numeric(unique(df$idyear)),
                      each = which(indexs$s.group == 1) %>% length)) %>%
  group_by(idyear) %>%
  nest() %>%
  #   ## project onto spatial field
  mutate(Spatial = map(.x = data,
                       .f = function(x)
                         as.matrix(inla.mesh.project(proj.grid, x))))%>%
  mutate(geometry = list(pred_loc %>%
                           dplyr::select(Long, Lat)))

# Retrieve the fixed effects 
Xmat <- cbind(1, model.matrix(~ -1 + Geomorph, data=pred_loc)) 

wch <- c(6,4)
ii = unlist(i.mod[wch])
cellmeans.full.1 <- t(cellmeans[ii,]) %*% t(Xmat)

cellmeans.fixed <- pred_loc %>%
  dplyr::select(Long, Lat, Geomorph) %>%
  cbind(V = t(cellmeans.full.1)) %>%
  as.data.frame %>%                           
  dplyr::select(starts_with("V"))%>%
  slice(rep(row_number(), length(unique(df$idyear))))%>%
  mutate(idyear = rep(unique(df$idyear),each=nrow(pred_loc)))%>%
  group_by(idyear) %>%
  nest() 

# Add the posteriors together

cellmeans.full.c <-
  cellmeans.full %>%
  full_join(cellmeans.fixed %>%
              rename(data1 = data)) %>%
  mutate(value = map2(.x = Spatial, .y = data1,
                      .f = function(.x, .y) as.data.frame(.x + .y))) %>%
  dplyr::select(idyear, geometry, value) %>%
  unnest(cols = c(geometry, value)) %>%
  pivot_longer(c = starts_with("V"), names_to = "Rep") %>%
  mutate(Rep = gsub('\\.','',Rep)) %>%
  ungroup()

cellmeans.full.cc <- cellmeans.full.c %>%
  mutate(value = plogis(value))

# Summarized posterior distributions for each year and predictive locations 

p <- c(.25,.5,.75)
p_names <- paste0(p*100, "%")
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

pred.sum <- cellmeans.full.cc %>%
  group_by(idyear, Long, Lat) %>% 
  summarize_at(vars(value), funs(!!!p_funs))