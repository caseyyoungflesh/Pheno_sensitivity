######################
# arr ~ lat + green-up SVC + sensitivity
# Intercept as function of lat, spatially varying betas (effect of gr)
# 
######################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/pheno_trends/'

#Xanadu
#dir <- '/labs/Tingley/'


# model dir ------------------------------------------------------------

#IAR data
arr_date <- '2020-07-21'

#run date
run_date <- '2020-08-25'

#env date
env_date <- '2020-08-06'


# create dir and copy script ----------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date)),
       dir.create(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date)),
       FALSE)

system(paste0('cp ', dir, 'Scripts/3-arr-gr-SVC-sens.R ', dir, 'Results/arr-gr-SVC-sens-', run_date, '/3-arr-gr-SVC-sens-', run_date, '.R'))


# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)


# Filter data ------------------------------------------------------------------

#master arrival data (from IAR output)
setwd(paste0(dir, 'Data/arrival_master_', arr_date))
arr_master <- readRDS(paste0('arrival_master_', arr_date, '.rds'))
#invalid GAM cells as NA
na.idx <- which(arr_master$VALID_GAM == FALSE)
arr_master$arr_GAM_mean[na.idx] <- NA

#greenup
setwd(paste0(dir, 'Data/environment/processed/', env_date))
gr <- readRDS(paste0('MidGreenup-', env_date, '-forest.rds'))

#merge with greenup data
df_gr <- which(colnames(gr) %in% c('cell_lat', 'cell_lng'))
mrg_f <- dplyr::left_join(arr_master, gr[,-df_gr], by = c('cell', 'year'))

#only rows that have greenup data and cells that had > 10000 cells with valid greenup data (from forest pixels) - 2500 km^2 greenup data
mrg_f2 <- dplyr::filter(mrg_f, !is.na(gr_mn), gr_ncell > 10000, per_ovr >= 0.05)
sp_k <- unique(mrg_f2$species)

#what proportion of species range are there data available for in more recent year
fd <- dplyr::filter(mrg_f2, !is.na(arr_GAM_mean))
#number total cells
ntc <- aggregate(cell ~ species, mrg_f2, function(x) length(unique(x)))
#total lat range
ntl <- aggregate(cell_lat ~ species, mrg_f2, function(x) max(x) - min(x))
#number obs each cell/species
cnt_csp <- plyr::count(fd, c('species', 'cell', 'cell_lat'))
ge3_csp <- dplyr::filter(cnt_csp, freq >= 3)
#number cells greater than or equal to 3
ncd <- aggregate(cell ~ species, ge3_csp, function(x) length(unique(x)))
#range of cells greater than or equal to 3
ncl <- aggregate(cell_lat ~ species, ge3_csp, function(x) max(x) - min(x))
#merge
jj <- dplyr::left_join(ntc, ncd, by = 'species')
jj2 <- dplyr::left_join(jj, ntl, by = 'species')
jj3 <- dplyr::left_join(jj2, ncl, by = 'species')
#proportion cells with >= 3 data points
jj3$pcd <- jj3[,3] /  jj3[,2]
#proportion lat range where cells have >= 3 data points
jj3$lrng <- jj3[,5] /  jj3[,4]
colnames(jj3) <- c('species', 'n_cells', 'n_cells_ge3', 'rng_lat', 
                   'rng_lat_ge3', 'prop_cell_ge3', 'prop_lat_ge3')
tsid <- which(jj3$prop_cell_ge3 > 0.4)
sp_k2 <- sp_k[tsid]

#filter species
mrg_f3 <- dplyr::filter(mrg_f2, species %in% sp_k2)
mrg_f4 <- mrg_f3

mrg_f4$sp_idx <- as.numeric(factor(mrg_f4$species))
usp <- unique(mrg_f4$sp_idx)

cells <- sort(unique(mrg_f4$cell))
ncell <- length(cells)


#add sp/cell id and scaled greenup
mrg_f4$cn_id <- NA
gr_df <- data.frame(sp_idx = rep(NA, NROW(mrg_f4)), 
                    year = NA, 
                    cell = NA, 
                    sc_gr = NA)
counter <- 0
counter_gr <- 1
for (i in 1:length(usp))
{
  #i <- 1
  t_idx <- which(mrg_f4$sp_id == usp[i])
  temp <- mrg_f4[t_idx,]
  na_idx <- which(is.na(temp$arr_GAM_mean))
  temp$gr_mn[na_idx] <- NA
  
  tj <- as.numeric(factor(temp$cell))
  mrg_f4$cn_id[t_idx] <- counter + tj
  
  counter <- counter + max(tj)
  
  u_cell <- unique(temp$cell)
  #anomaly for each cell across years
  for (j in 1:length(u_cell))
  {
    #j <- 2
    temp2 <- dplyr::filter(temp, cell == u_cell[j])
    sc_t <- scale(temp2$gr_mn, scale = FALSE)[,1]
    
    nyc <- counter_gr + length(sc_t) - 1
    gr_df$sp_idx[counter_gr:nyc] <- usp[i]
    gr_df$cell[counter_gr:nyc] <- u_cell[j]
    gr_df$year[counter_gr:nyc] <- temp2$year
    gr_df$sc_gr[counter_gr:nyc] <- sc_t
    
    counter_gr <- nyc + 1
  }
}

#add scaled greenup to dataframe
mrg_f5 <- dplyr::left_join(mrg_f4, gr_df, by = c('sp_idx', 'year', 'cell'))

# #data for Shiny app
# setwd("~/Google_Drive/R/pheno_trends/Data/For_Bruna/shiny_data/input")
# saveRDS(mrg_f5, paste0('mrg_f5-', run_date, '.rds'))

#Csp = for each species what the cnids are (must be in same order as lambdas)
#species in rows, cells in columns - pad with -999
Csp <- matrix(NA, nrow = length(usp), ncol = ncell)
#number of cells for each species
clen <- rep(NA, length(usp))
#hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)
#neighbors from adj matrices (# cell * 6 is > max number of neighbors)
#pad with -999
node1 <- matrix(NA, nrow = length(usp), ncol = ncell*6)
node2 <- matrix(NA, nrow = length(usp), ncol = ncell*6)
#number of pairwise combos
N_edges <- rep(NA, length(usp))
#cell lats at each cell for each species
cell_lat2 <- matrix(NA, nrow = length(usp), ncol = ncell)
#start and end of cn_id for each species for spatial component
st <- rep(NA, length(usp))
end <- rep(NA, length(usp))
st[1] <- 1
for (k in 1:length(usp))
{
  #k <- 1
  t_idx <- which(mrg_f5$sp_id == usp[k])
  temp <- mrg_f5[t_idx,]
  t_ucell <- unique(temp$cell)
  t_ucid <- unique(temp[,c('cell', 'cn_id', 'cell_lat')])
  
  Csp[k,] <- c(t_ucid$cn_id, rep(-999, ncell - NROW(t_ucid)))
  tncell <- NROW(t_ucid)
  clen[k] <- tncell
  cell_lat2[k,] <- c(scale(t_ucid$cell_lat, scale = FALSE)[,1], 
                     rep(-999, ncell - NROW(t_ucid)))
  
  #cell centers
  tcellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_ucid$cell)
  
  #create adjacency matrix - 1 if adjacent to cell, 0 if not
  t_adjacency_matrix <- matrix(data = NA, nrow = tncell, ncol = tncell)
  for (i in 1:tncell)
  {
    #i <- 1
    for (j in i:tncell)
    {
      #j <- 2
      dists <- geosphere::distm(c(tcellcenters$lon_deg[i], tcellcenters$lat_deg[i]),
                                c(tcellcenters$lon_deg[j], tcellcenters$lat_deg[j]))
      t_adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
    }
  }
  
  #indices for 1s
  tninds <- which(t_adjacency_matrix == 1, arr.ind = TRUE)
  
  node1[k,] <- c(tninds[,1], rep(-999, NCOL(node1) - NROW(tninds)))
  node2[k,] <- c(tninds[,2], rep(-999, NCOL(node1) - NROW(tninds)))
  
  N_edges[k] <- NROW(tninds)
  
  #start and end of cn_id for each species for spatial component
  if (k > 1)
  {
    st[k] <- end[k-1] + 1
  }
  end[k] <- st[k] + tncell - 1 
}

#species index and cell_lat
u_cn_id <- sort(unique(mrg_f5$cn_id))
sp_id <- c()
cell_lat <- c()
for (i in 1:length(u_cn_id))
{
  #i <- 1
  s_idx <- which(mrg_f5$cn_id == u_cn_id[i])
  sn <- mrg_f5[s_idx,'sp_idx'][1]
  sp_id <- c(sp_id, sn)
  
  cl <- mrg_f5[s_idx,'cell_lat'][1]
  cell_lat <- c(cell_lat, cl)
}

tsp <- unique(mrg_f5[,c('species', 'cell', 'cell_lat', 'sp_idx')])

#unique species/cells
u5 <- unique(mrg_f5[,c('species', 'cell', 'cell_lat', 'cell_lng', 'sp_idx', 'cn_id')])

#no missing data (no need to impute values)
mrg_f6 <- dplyr::filter(mrg_f5, !is.na(arr_GAM_mean))

#response data
arr_obs <- mrg_f6$arr_IAR_mean
sd_arr <- mrg_f6$arr_IAR_sd


#PC1 and PC2
#inverse migration speed
f1 <- unique(mrg_f6[,c('species', 'beta_gamma_mean')])
#arrival mean
f2 <- aggregate(arr_IAR_mean ~ species, arr_master, mean)
f_mrg <- dplyr::left_join(f1, f2, by = 'species')

#study area grid
setwd(paste0(dir, 'Data/hex_grid_crop'))
hg_crop_comb_p <- rgdal::readOGR('hex_grid_crop_combine.shp')
#reproject
hg_crop_comb <- sp::spTransform(hg_crop_comb_p, sp::CRS("+proj=laea +lat_0=0 +lon_0=-70 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


#overwinter latitude
centers <- data.frame(species = rep(NA, length(sp_k2)),
                      nb_lat = NA, 
                      nb_lng = NA,
                      mb_lat = NA,
                      mb_lng = NA,
                      nb_lat_4326 = NA, 
                      nb_lng_4326 = NA,
                      mb_lat_4326 = NA,
                      mb_lng_4326 = NA,
                      dist = NA)

#reference key for species synonyms
setwd(paste0(dir, '../Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')
#get centroids
for (i in 1:length(sp_k2))
{
  #i <- 48
  print(paste('species ', i, ' of ', length(sp_k2)))
  
  args <- sp_k2[i]
  
  #change dir to shp files
  setwd(paste0(dir, '../Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
  
  #filter by breeding/migration cells
  #match species name to shp file name
  g_ind <- grep(args, sp_key$file_names_2016)
  
  #check for synonyms if there are no matches
  if (length(g_ind) == 0)
  {
    g_ind2 <- grep(args, sp_key$BL_Checklist_name)
  } else {
    g_ind2 <- g_ind
  }
  
  #get filename and read in
  fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
  sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
  
  #reproject
  sp_rng2 <- sp::spTransform(sp_rng, sp::CRS("+proj=laea +lat_0=0 +lon_0=-70 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  # #filter by non-breeding (3) range [reference: resident (1), breeding (2), passage (4)]
  nrng_nb <- sp_rng2[which(sp_rng2$SEASONAL == 3),]
  #buffer of 0 to deal with intersection issues
  nrng_nb2 <- rgeos::gBuffer(nrng_nb, byid = TRUE, width = 0)
  #get centroid
  center_nb <- rgeos::gCentroid(nrng_nb2)
  #back to 4326
  center_nb_4326 <- sp::spTransform(center_nb, sp::CRS('+init=epsg:4326'))
  
  #breeding and migration range
  nrng_mb <- sp_rng2[which(sp_rng2$SEASONAL == 2 | sp_rng2$SEASONAL == 4),]
  #buffer of 0 to deal with intersection issues
  nrng_mb2 <- rgeos::gBuffer(nrng_mb, byid = TRUE, width = 0)
  #clip to study area
  nrng_mb_cr <- rgeos::gIntersection(nrng_mb2, hg_crop_comb, byid = TRUE)
  #centroid
  center_mb <- rgeos::gCentroid(nrng_mb_cr)
  #back to 4326
  center_mb_4326 <- sp::spTransform(center_mb, sp::CRS('+init=epsg:4326'))
  
  nb_mb_dist <- geosphere::distm(center_nb_4326@coords, center_mb_4326@coords)
  
  #store in df
  centers$species[i] <- args
  centers$nb_lat[i] <- center_nb@coords[2]
  centers$nb_lng[i] <- center_nb@coords[1]
  centers$mb_lat[i] <- center_mb@coords[2]
  centers$mb_lng[i] <- center_mb@coords[1]
  centers$nb_lat_4326[i] <- center_nb_4326@coords[2]
  centers$nb_lng_4326[i] <- center_nb_4326@coords[1]
  centers$mb_lat_4326[i] <- center_mb_4326@coords[2]
  centers$mb_lng_4326[i] <- center_mb_4326@coords[1]
  centers$dist[i] <- nb_mb_dist / 1000
}

f_mrg2 <- dplyr::left_join(f_mrg, centers, by = 'species')

#multicolinearity
covs <- select(f_mrg2, beta_gamma_mean, arr_IAR_mean, nb_lat_4326)
colnames(covs) <- c('Inv mig speed', 'Mean arrival', 'Lat')
#cor(covs)

#PCA on variables
covs_pca <- prcomp(covs, center = TRUE, scale. = TRUE)


# Stan model --------------------------------------------------------------

#cell_lat2 species in rows (std lat), cell in cols -> pad with -999

DATA <- list(N = NROW(mrg_f6),
             NC = length(u_cn_id),
             Nsp = length(unique(sp_id)),
             arr_obs = arr_obs,
             sd_arr = sd_arr,
             sc_gr = mrg_f6$sc_gr,
             sp_id = sp_id,
             sp_id2 = mrg_f6$sp_idx,
             cn_id = mrg_f6$cn_id,
             cell_lat = cell_lat,
             cell_lat2 = cell_lat2[,1:max(clen)],
             N_edges = N_edges,
             MN_edges = max(N_edges),
             node1 = node1[,1:max(N_edges)],
             node2 = node2[,1:max(N_edges)],
             clen = clen,
             Sclen = sum(clen),
             Mclen = max(clen),
             Csp = Csp[,1:max(clen)],
             neg_nine = rep(-999, max(clen)),
             st = st,
             end = end,
             PC1 = covs_pca$x[,1],
             covs_pca = covs_pca,
             f_mrg2 = f_mrg2,
             u5 = u5,
             mrg_f6 = mrg_f6)

stanmodel1 <- "
data {
int<lower=0> N;                            // number of data points
int<lower=0> NC;                           // number of species/cells
int<lower=0> Nsp;                          // number of species
int<lower=0> clen[Nsp];                    // number of cell for each species
int<lower=0> Mclen;                        // max number of cell for any species
int<lower=0> Sclen;                        // number of cells used over all species
vector<lower=0>[N] arr_obs;
vector<lower=0>[N] sd_arr;
vector[N] sc_gr;
int<lower=0> sp_id[NC];                    // species ids by cell
int<lower=0> cn_id[N];                     // species/cell ids
vector<lower=0>[NC] cell_lat;
real cell_lat2[Nsp, Mclen];                // cell lat for every cell in model
int<lower = 0> N_edges[Nsp];               // number of edges in adjacency matrices
int<lower = 0> MN_edges;                   // max number of edges
int node1[Nsp, MN_edges];                  // node1[i] adjacent to node2[i]
int node2[Nsp, MN_edges];                  // and node1[i] < node2[i]
int Csp[Nsp, Mclen];
real neg_nine[Mclen];                      // array of -999 to fill phis not estimated
int st[Nsp];                               // start index for phi
int end[Nsp];                              // end index for phi
vector[Nsp] PC1;
}

parameters {
vector[N] mu_arr;
real<lower = 0> sigma_APG;
real mu_pi;
real mu_nu;
vector<lower = 0>[2] sigma_pn;
cholesky_factor_corr[2] L_Rho_pn;         // cholesky factor of corr matrix
matrix[2, Nsp] z_pn;
vector[NC] alpha_raw;
vector[NC] beta_raw;
real<lower = 0> sigma_alpha;
real<lower = 0> sigma_beta;
real<lower = 0> sigma_xi;
vector[Nsp] xi_raw;
real mu_gamma;
real<lower = 0> sigma_gamma;
vector<offset = mu_gamma, multiplier = sigma_gamma>[Nsp] gamma;
real phi_est[Sclen];                  // phis that are to be estimated
real lambda;
real<lower = 0> kappa;
vector[Nsp] sigma_phi_raw;
// vector<offset = lambda, multiplier = kappa><lower = 0>[Nsp] sigma_phi;
real alpha_xi;
real beta_xi;
}

transformed parameters {
vector[N] mu_APG;
matrix[Nsp, 2] pn;
matrix[2, 2] Rho_pn;                  // correlation matrix
vector[NC] mu_alpha;
vector[NC] mu_beta;
vector[NC] alpha;
vector[NC] beta;
vector[Nsp] pi;
vector[Nsp] nu;
real phi[Nsp, Mclen];                // spatial error component
vector<lower = 0>[Nsp] sigma_phi;
vector[Nsp] mu_xi;
vector[Nsp] xi;

// implies sigma_phi ~ lognormal(lambda, kappa)
sigma_phi = exp(sigma_phi_raw * kappa + lambda);

mu_xi = alpha_xi + beta_xi * PC1;

// implies xi ~ normal(mu_xi, sigma_xi)
xi = xi_raw * sigma_xi + mu_xi;

// cholesky factor of covariance matrix (i.e., diagonal matrix of scale times cholesky factor of correlation matrix) multiplied by z score
// implies pn ~ MVN(0, Sigma)
pn = (diag_pre_multiply(sigma_pn, L_Rho_pn) * z_pn)';
// implies Rho_pn = L_Rho_pn * L_Rho_pn';
Rho_pn = multiply_lower_tri_self_transpose(L_Rho_pn);
pi = mu_pi + pn[,1];
nu = mu_nu + pn[,2];

mu_alpha = pi[sp_id] + nu[sp_id] .* cell_lat;

// implies alpha[jk] ~ N(mu_alpha[jk], sigma_alpha)
alpha = alpha_raw * sigma_alpha + mu_alpha;

for (k in 1:Nsp)
{
  // estimate phi where cells are used in spatial component of model 
  phi[k, 1:clen[k]] = phi_est[st[k]:end[k]];
  // fill phi with -999 where cells are not used in spatial model (to prevent estimation)
  phi[k, (clen[k] + 1):Mclen] = neg_nine[1:(Mclen - clen[k])];
  
  // Csp = sp/cell numbers for that species (in order adj matrix was created with)
  // only use phi values where there are cells for that species
  mu_beta[Csp[k, 1:clen[k]]] = xi[k] + gamma[k] * to_vector(cell_lat2[k, 1:clen[k]]) + to_vector(phi[k, 1:clen[k]]) * sigma_phi[k];
}
beta = beta_raw * sigma_beta + mu_beta;

// cn_id is ids just for cells that have valid GAM HM
// alpha and beta are estimated for all cells in model (entire range of all species)
// response data are only available for some of these cells
mu_APG = alpha[cn_id] + beta[cn_id] .* sc_gr;

}

model {
// priors
sigma_APG ~ normal(0, 5);
sigma_alpha ~ normal(0, 10);
sigma_beta ~ std_normal();
mu_pi ~ normal(0, 40);
mu_nu ~ normal(0, 3);
sigma_pn[1] ~ normal(0, 50);
sigma_pn[2] ~ std_normal();
lambda ~ normal(-1, 2);
kappa ~ normal(0, 1);
// mu_xi ~ normal(0.5, 0.5); // based on phenological response of other taxa
sigma_xi ~ normal(0, 0.5);
mu_gamma ~ normal(0, 0.1);
sigma_gamma ~ normal(0, 0.1);
alpha_xi ~ normal(0, 2);
beta_xi ~ normal(0, 2);

// centered
mu_arr ~ normal(mu_APG, sigma_APG);

// non-centered
// xi ~ normal(mu_xi, sigma_xi);
gamma ~ normal(mu_gamma, sigma_gamma);
xi_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();
sigma_phi_raw ~ std_normal();
// sigma_phi ~ lognormal(lambda, kappa);

to_vector(z_pn) ~ std_normal();
L_Rho_pn ~ lkj_corr_cholesky(1);

for (k in 1:Nsp)
{
  int n1[N_edges[k]] = node1[k, 1:N_edges[k]];
  int n2[N_edges[k]] = node2[k, 1:N_edges[k]];

  // pairwise difference formulation
  // only include phi values where there are cells for that species
  target += -0.5 * dot_self(to_vector(phi[k, n1]) - to_vector(phi[k, n2]));
  // soft sum to 0 constraint
  sum(to_vector(phi[k, 1:clen[k]])) ~ normal(0, 0.001 * clen[k]);
}

// observation model for arrival
arr_obs ~ normal(mu_arr, sd_arr);

}

generated quantities {
real arr_rep[N];
vector[N] delta;

arr_rep = normal_rng(mu_arr, sd_arr);
delta = mu_arr - alpha[cn_id];
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.94
TREE_DEPTH <- 13
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 7000
WARM <- 3000

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date))

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   warmup = WARM,
                   cores = CHAINS,
                   #diagnostic_file = 'diagnostics.csv',
                   pars = c('alpha',
                            'beta',
                            'mu_alpha',
                            'sigma_alpha',
                            'mu_beta',
                            'sigma_beta',
                            'pi',
                            'nu',
                            'pn',
                            'mu_pi',
                            'mu_nu',
                            'sigma_pn',
                            'Rho_pn',
                            'xi',
                            'mu_xi',
                            'alpha_xi',
                            'beta_xi',
                            'sigma_xi',
                            'gamma',
                            'mu_gamma',
                            'sigma_gamma',
                            'phi',
                            'sigma_phi',
                            'lambda',
                            'kappa',
                            'sigma_APG',
                            'mu_arr',
                            'arr_rep',
                            'delta'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


# save object -------------------------------------------------------------

saveRDS(fit, file = paste0('arr-gr-SVC-sens-stan-output-', run_date, '.rds'))
saveRDS(DATA, file = paste0('arr-gr-SVC-sens-data-', run_date, '.rds'))
#fit <- readRDS(paste0('arr-gr-SVC-sens-stan-output-', run_date, '.rds'))
#DATA <- readRDS(paste0('arr-gr-SVC-sens-data-', run_date, '.rds'))


# Calc diagnostics ---------------------------------------------------

# library(shinystan)
# launch_shinystan(fit)

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mn_stepsize <- sapply(sampler_params, 
                      function(x) mean(x[, 'stepsize__']))
mn_treedepth <- sapply(sampler_params, 
                       function(x) mean(x[, 'treedepth__']))
accept_stat <- sapply(sampler_params, 
                      function(x) mean(x[, 'accept_stat__']))
num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))


# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, 
                                      n.eff = TRUE, 
                                      round = 3, 
                                      excl = 'arr_rep')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])


# PPC ---------------------------------------------------------------------

y_val <- DATA$arr_obs
y_rep <- MCMCvis::MCMCchains(fit, params = 'arr_rep')

bayesplot::ppc_dens_overlay(y_val, y_rep[1:100,])


# correlation -------------------------------------------------------------

#extract posteriors
delta <- MCMCvis::MCMCchains(fit, params = 'delta')
xi <- MCMCvis::MCMCchains(fit, params = 'xi')

#calculate sd of each species at each iteration
sp_sd <- matrix(NA, nrow = NROW(xi), ncol = dim(xi)[2])
for (k in 1:dim(xi)[2])
{
  #k <- 1
  t_idx <- which(DATA$sp_id2 == k)
  sp_sd[,k] <- apply(delta[,t_idx], 1, sd)
}

#mean interannual variation (sd) for each species
apply(sp_sd, 2, mean)

#mean interannual variation (sd) across species
#posterior mean of the standard deviation of \delta at each iteration.
mean(apply(delta, 1, sd))

#mean green-up variation (sd)
sd(DATA$sc_gr)

apply(xi, 2, mean)

#correlation between sd and sens at each iteration
cor_vec <- rep(NA, NROW(xi))
for (i in 1:NROW(xi))
{
  #i <- 1
  cor_vec[i] <- cor(sp_sd[i,], xi[i,])
}

mean(cor_vec)
quantile(cor_vec, probs = c(0.025, 0.975))

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date))
saveRDS(cor_vec, paste0('cor-vec-' , run_date, '.rds'))


# proportion of variance in arrival explained by greenup ------------------

alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')
mu_arr_ch <- MCMCvis::MCMCchains(fit, params = 'mu_arr')

uid <- unique(DATA$cn_id)
pv_mat <- matrix(NA, nrow = NROW(alpha_ch), ncol = length(uid))
mu <- matrix(NA, nrow = max(DATA$cn_id), ncol = 17)
#j = cell
#k = data point
for (i in 1:NROW(alpha_ch))
{
  #i <- 1
  print(i)
  for (j in 1:length(uid))
  {
    #j <- 455
    #which data points are associated with sp/cell id
    idx <- which(DATA$cn_id == uid[j])
    for (k in 1:length(idx))
    {
      #k <- 1
      mu[j,k] <- alpha_ch[i, uid[j]] + beta_ch[i, uid[j]] * DATA$sc_gr[idx[k]]
    }
    
    var_y <- var(mu_arr_ch[i,idx])
    var_mu <- var(mu[j,], na.rm = TRUE)
    pv_mat[i,j] <- var_mu / (var_y)
  }
}

med_pv <- apply(pv_mat, 2, function(x) median(x, na.rm = TRUE))
u5$pv <- NA
u5$pv[uid] <- med_pv

u6 <- dplyr::filter(u5, !is.na(pv))
mn_pv <- aggregate(pv ~ species, data = u6, mean)

range(mn_r2$pv)
mean(mn_r2$pv)

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date))
saveRDS(mn_pv, paste0('pv-arr-gr-' , run_date, '.rds'))


# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0('arr-gr-SVC-sens-stan-results-', run_date, '.txt'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output, na.rm = TRUE), ' \n'))
print(model_summary)
sink()


# PPO ---------------------------------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Results/arr-gr-SVC-sens-',
                          run_date, '/Figures')),
       dir.create(paste0(dir, 'Results/arr-gr-SVC-sens-',
                         run_date, '/Figures')),
       FALSE)

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date, '/Figures'))


#alpha_xi ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_xi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-alpha_xi-', run_date, '.pdf'))

#beta_xi ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_xi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-beta_xi-', run_date, '.pdf'))

#sigma_xi ~ HN(0, 0.5)
PR_p <- rnorm(10000, 0, 0.5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_xi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_xi-', run_date, '.pdf'))

#mu_gamma ~ N(0, 0.1)
PR <- rnorm(10000, 0, 0.1)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-mu_gamma-', run_date, '.pdf'))

#sigma_gamma ~ HN(0, 0.1)
PR_p <- rnorm(10000, 0, 0.1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_gamma-', run_date, '.pdf'))


#mu_pi ~ N(0, 40)
PR <- rnorm(10000, 0, 40)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_pi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-mu_pi-', run_date, '.pdf'))

#mu_nu ~ N(0, 3)
PR <- rnorm(10000, 0, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-mu_nu-', run_date, '.pdf'))

#sigma_pn[1] ~ HN(0, 50)
PR_p <- rnorm(10000, 0, 50)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_pn\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_pn[1]-', run_date, '.pdf'))

#sigma_pn[2] ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_pn\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_pn[2]-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_beta-', run_date, '.pdf'))

#sigma_APG ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_APG',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-sigma_APG-', run_date, '.pdf'))

#lambda ~ N(-1, 2)
PR <- rnorm(10000, -1, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'lambda',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-lambda-', run_date, '.pdf'))

#kappa ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'kappa',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('arr-gr-SVC-sens-trace-kappa-', run_date, '.pdf'))


# subset object -----------------------------------------------------------

y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary',
                                 func = mean)[[1]]
y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary',
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary',
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]
fit_ab <- MCMCvis::MCMCchains(fit, params = c('alpha', 'beta'), mcmc.list = TRUE)
fit_xgpn <- MCMCvis::MCMCchains(fit, params = c('xi', 'gamma', 'pi', 'nu'), mcmc.list = TRUE)
fit_subset <- MCMCvis::MCMCchains(fit, params = c('mu_pi', 'mu_nu', 'alpha_xi', 'mu_gamma', 'beta_xi'), 
                                  mcmc.list = TRUE)


#merge with names
beta_mean <- MCMCvis::MCMCpstr(fit_ab, params = 'beta', func = mean)[[1]]
beta_sd <- MCMCvis::MCMCpstr(fit_ab, params = 'beta', func = sd)[[1]]
u5$beta_mean <- beta_mean
u5$beta_sd <- beta_sd

xi_mean <- MCMCvis::MCMCpstr(fit_xgpn, params = 'xi', func = mean)[[1]]
xi_sd <- MCMCvis::MCMCpstr(fit_xgpn, params = 'xi', func = sd)[[1]]

gamma_mean <- MCMCvis::MCMCpstr(fit_xgpn, params = 'gamma', func = mean)[[1]]
gamma_sd <- MCMCvis::MCMCpstr(fit_xgpn, params = 'gamma', func = sd)[[1]]

pi_mean <- MCMCvis::MCMCpstr(fit_xgpn, params = 'pi', func = mean)[[1]]
pi_sd <- MCMCvis::MCMCpstr(fit_xgpn, params = 'pi', func = sd)[[1]]

nu_mean <- MCMCvis::MCMCpstr(fit_xgpn, params = 'nu', func = mean)[[1]]
nu_sd <- MCMCvis::MCMCpstr(fit_xgpn, params = 'nu', func = sd)[[1]]

ps_mrg <- data.frame(xi_mean, xi_sd, 
                     gamma_mean, gamma_sd, 
                     pi_mean, pi_sd, 
                     nu_mean, nu_sd, 
                     species = unique(u5$species))
psu_mrg <- dplyr::left_join(u5, ps_mrg, by = 'species')

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date))
saveRDS(psu_mrg, file = paste0('arr-gr-SVC-sens-psummary-', run_date, '.rds'))
saveRDS(fit_subset, file = paste0('arr-gr-SVC-sens-stan-output-subset-', run_date, '.rds'))

# rm(fit)
# gc()


# Map of trends ------------------------------------------

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date, '/Figures'))

#load maps
worldmap <- data.frame(maps::map("world", plot = FALSE)[c("x", "y")])

#min/max for plotting using output data
MIN <- min(u5$beta_mean)
MAX <- max(u5$beta_mean)


#plot
map_vis_fun <- function(SPECIES = 'Geothlypis_trichas', min = MIN, max = MAX)
{
  #SPECIES <- 'Catharus_guttatus'
  uj5t <- dplyr::filter(u5, species == SPECIES)
  
  if (NROW(uj5t) < 1)
  {
    stop(paste0('Species: ', SPECIES, ' not found!'))
  }
  
  #transform cells to grid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  cell_grid <- dggridR::dgcellstogrid(hexgrid6, uj5t$cell)
  cell_grid$cell <- as.numeric(cell_grid$cell)
  
  to_plt <- dplyr::inner_join(uj5t, cell_grid, by = 'cell')
  
  p_beta <- ggplot() +
    geom_path(data = worldmap,
              aes(x = x, y = y), color = 'black') +
    # geom_polygon(data = nrng.df,
    #           aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) +
    # geom_polygon(data = nrng_rm.df,
    #              aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.4) +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -65), ylim = c(23, 55)) +
    geom_polygon(data = to_plt, aes(x = long, y = lat,
                                    group = group, fill = beta_mean),
                 alpha = 0.5) +
    geom_path(data = to_plt, aes(x = long, y = lat, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradient2(low = 'indianred', high = 'royalblue', mid = 'lightgoldenrod',
                         limits = c(min, max), midpoint = 0) +
    labs(fill = 'Slope') +
    ggtitle(paste0(SPECIES, ' - Arrival ~ gr')) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')
  
  print(p_beta)
}

sps <- unique(u5$species)
#one species, all cells
for (i in 1:length(sps))
{
  #i <- 24
  tt <- dplyr::filter(u5, species == sps[i])
  pdf(paste0('arr-gr-SVC-sens-map-', sps[i], '-', run_date, '.pdf'))
  map_vis_fun(SPECIES = sps[i], 
              # min = min(tt$beta_mean), 
              # max = max(tt$beta_mean))
              min = MIN, 
              max = MAX)
  dev.off()
}


# param est --------------------------------------------------------------

setwd(paste0(dir, 'Results/arr-gr-SVC-sens-', run_date, '/Figures'))

sp_nms <- unique(u5[,c('species', 'sp_idx')])$species

pdf('xi.pdf', height = 12, width = 8)
MCMCvis::MCMCplot(fit_xgpn, params = 'xi',
                  labels = sp_nms,
                  guide_lines = TRUE,
                  sz_labels = 0.5,
                  #xlim = c(-0.2, 0.6),
                  rank = TRUE,
                  main = 'xi (change in days arr / change in day greenup)')
dev.off()

pdf('gamma.pdf', height = 12, width = 8)
MCMCvis::MCMCplot(fit_xgpn, params = 'gamma',
                  labels = sp_nms,
                  guide_lines = TRUE,
                  sz_labels = 0.5,
                  #xlim = c(-0.2, 0.6),
                  rank = TRUE,
                  main = 'gamma (change in xi / degree lat)')
dev.off()

