######################
# trends in greenup date - spatially varying coefficients model
#
# Model for 1) forest-filtered and 2) all land-cover greenup
######################

 
# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/pheno_trends/'

#Xanadu
#dir <- '/labs/Tingley/'


# model dir ------------------------------------------------------------

#run date
run_date <- '2020-08-20'
env_date <- '2020-08-06'
arr_date <- '2020-07-21'


# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)
library(dggridR)
library(geosphere)


# Filter data ------------------------------------------------------------------

setwd(paste0(dir, 'Data/environment/processed/', env_date))

#read in greenup
gr_forest <- readRDS(paste0('MidGreenup-', env_date, '-forest.rds'))
gr_all <- readRDS(paste0('MidGreenup-', env_date, '-all.rds'))

#read in arrival to get cells
setwd(paste0(dir, 'Data/arrival_master_', arr_date))
arr_master <- readRDS(paste0('arrival_master_', arr_date, '.rds'))
ucells <- unique(arr_master[,c('per_ovr', 'cell')])
ucells2 <- dplyr::filter(ucells, per_ovr >= 0.05)

#filter forest by cells
gr_forest2 <- dplyr::filter(gr_forest, cell %in% ucells2$cell, 
                            gr_ncell > 10000, year >= 2002, year <= 2017)
gr_all2 <- dplyr::filter(gr_all, cell %in% unique(gr_forest2$cell),
                         gr_ncell > 10000, year >= 2002, year <= 2017)


# data processing ----------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
gr_cell <- unique(gr_forest2[,c('cell', 'cell_lat')])
gr_cell2 <- gr_cell[order(gr_cell[,'cell']),]
ncell <- NROW(gr_cell2)

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, gr_cell2$cell)

#create adjacency matrix - 1 if adjacent to cell, 0 if not
adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)

for (i in 1:ncell)
{
  #i <- 1
  for (j in i:ncell)
  {
    #j <- 2
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)


# Stan model --------------------------------------------------------------

#forest
DATA1 <- list(N = NROW(gr_forest2),
             NC = ncell,
             y = gr_forest2$gr_mn,
             cn_id = as.numeric(factor(gr_forest2$cell)),
             year = as.numeric(factor(gr_forest2$year)),
             cell_lat = gr_cell2$cell_lat,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2])

#all
DATA2 <- list(N = NROW(gr_all2),
              NC = ncell,
              y = gr_all2$gr_mn,
              cn_id = as.numeric(factor(gr_all2$cell)),
              year = as.numeric(factor(gr_all2$year)),
              cell_lat = gr_cell2$cell_lat,
              N_edges = nrow(ninds), 
              node1 = ninds[,1],
              node2 = ninds[,2])


#env over time
stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
int<lower=0> NC;                      // number of cells
vector[N] y;                          // response
int<lower=0> cn_id[N];                // cell ids
vector<lower=0>[N] year;
vector<lower=0>[NC] cell_lat;
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = NC> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = NC> node2[N_edges];             // and node1[i] < node2[i]
}

parameters {
real<lower = 0> sigma;
real pi;
real nu;
real xi;
vector[NC] alpha_raw;
vector[NC] beta_raw;
real<lower = 0> sigma_alpha;
real<lower = 0> sigma_beta;
vector[NC] phi;                                     // spatial error componenet
real<lower = 0> sigma_phi;
}

transformed parameters {
vector[N] mu;
vector[NC] mu_alpha;
vector[NC] mu_beta;
vector[NC] alpha;
vector[NC] beta;

mu_alpha = pi + nu * cell_lat;
alpha = alpha_raw * sigma_alpha + mu_alpha;
mu_beta = xi + phi * sigma_phi;
beta = beta_raw * sigma_beta + mu_beta;

mu = alpha[cn_id] + beta[cn_id] .* year;
}

model {
// priors
sigma ~ normal(0, 10);
sigma_alpha ~ normal(0, 20);
sigma_beta ~ std_normal();
pi ~ normal(120, 50);
nu ~ normal(0, 5);
xi ~ std_normal();
sigma_phi ~ normal(0, 5);


// non-centered
alpha_raw ~ std_normal();
beta_raw ~ std_normal();

// pairwise difference formulation
target += -0.5 * dot_self(phi[node1] - phi[node2]);
// soft sum to 0 constraint
sum(phi) ~ normal(0, 0.001 * NC);

y ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu, sigma);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.85
TREE_DEPTH <- 15
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 8000


# fit gr ----------------------------------------------------------------

tt <- proc.time()
fit1 <- rstan::stan(model_code = stanmodel1,
                   data = DATA1,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'mu_alpha',
                            'sigma_alpha',
                            'mu_beta',
                            'sigma_beta',
                            'pi',
                            'nu',
                            'xi',
                            'phi',
                            'sigma_phi',
                            'sigma',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))

fit2 <- rstan::stan(model_code = stanmodel1,
                    data = DATA2,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'mu_alpha',
                             'sigma_alpha',
                             'mu_beta',
                             'sigma_beta',
                             'pi',
                             'nu',
                             'xi',
                             'phi',
                             'sigma_phi',
                             'sigma',
                             'y_rep'), 
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


# save object -------------------------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Results/gr-SVC-', run_date)),
       dir.create(paste0(dir, 'Results/gr-SVC-', run_date)),
       FALSE)

setwd(paste0(dir, 'Results/gr-SVC-', run_date))
saveRDS(fit1, file = paste0('gr-forest-SVC-stan-output-', run_date, '.rds'))
saveRDS(DATA1, file = paste0('gr-forest-SVC-data-', run_date, '.rds'))
saveRDS(fit2, file = paste0('gr-all-SVC-stan-output-', run_date, '.rds'))
saveRDS(DATA2, file = paste0('gr-all-SVC-data-', run_date, '.rds'))


# Analyze -----------------------------------------------------------------

# MCMCvis::MCMCsummary(fit, 
#                      params = c('pi', 'nu'), 
#                      round = 3)
# MCMCvis::MCMCsummary(fit, 
#                      params = c('mu_pi', 'mu_nu'), 
#                      round = 4)
# MCMCvis::MCMCsummary(fit, 
#                      params = c('mu_alpha', 'mu_beta'), 
#                      round = 4)
# MCMCvis::MCMCsummary(fit, 
#                      params = 'sigma',
#                      ISB = FALSE,
#                      round = 4)
# MCMCvis::MCMCplot(fit, params = 'beta', ref_ovl = TRUE)
# MCMCvis::MCMCplot(fit, params = 'mu_beta', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'pi')
# MCMCvis::MCMCplot(fit, params = 'nu')
# 
# beta_s <- MCMCvis::MCMCsummary(fit, params = 'beta')
# which(beta_s[,'97.5%'] < 0)


# Calc diagnostics ---------------------------------------------------

# library(shinystan)
# launch_shinystan(fit1)
# launch_shinystan(fit2)

sampler_params1 <- get_sampler_params(fit1, inc_warmup = FALSE)
mn_stepsize1 <- sapply(sampler_params1, 
                      function(x) mean(x[, 'stepsize__']))
mn_treedepth1 <- sapply(sampler_params1, 
                       function(x) mean(x[, 'treedepth__']))
accept_stat1 <- sapply(sampler_params1, 
                      function(x) mean(x[, 'accept_stat__']))
num_diverge1 <- rstan::get_num_divergent(fit1)
num_tree1 <- rstan::get_num_max_treedepth(fit1)
num_BFMI1 <- length(rstan::get_low_bfmi_chains(fit1))

sampler_params2 <- get_sampler_params(fit2, inc_warmup = FALSE)
mn_stepsize2 <- sapply(sampler_params2, 
                       function(x) mean(x[, 'stepsize__']))
mn_treedepth2 <- sapply(sampler_params2, 
                        function(x) mean(x[, 'treedepth__']))
accept_stat2 <- sapply(sampler_params2, 
                       function(x) mean(x[, 'accept_stat__']))
num_diverge2 <- rstan::get_num_divergent(fit2)
num_tree2 <- rstan::get_num_max_treedepth(fit2)
num_BFMI2 <- length(rstan::get_low_bfmi_chains(fit2))


# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary1 <- MCMCvis::MCMCsummary(fit1, Rhat = TRUE, 
                                      n.eff = TRUE, 
                                      round = 2, 
                                      excl = c('y_rep'))

model_summary2 <- MCMCvis::MCMCsummary(fit2, Rhat = TRUE, 
                                       n.eff = TRUE, 
                                       round = 2, 
                                       excl = c('y_rep'))

#extract Rhat and neff values
rhat_output1 <- as.vector(model_summary1[, grep('Rhat', colnames(model_summary1))])
neff_output1 <- as.vector(model_summary1[, grep('n.eff', colnames(model_summary1))])

rhat_output2 <- as.vector(model_summary2[, grep('Rhat', colnames(model_summary2))])
neff_output2 <- as.vector(model_summary2[, grep('n.eff', colnames(model_summary2))])


# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('gr-forest-SVC-stan-results-', run_date, '.txt'))
cat(paste0('Total minutes (both): ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge1, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree1, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI1, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize1), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth1), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat1), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output1, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output1, na.rm = TRUE), ' \n'))
print(model_summary1)
sink()

options(max.print = 50000)
sink(paste0('gr-all-SVC-stan-results-', run_date, '.txt'))
cat(paste0('Total minutes (both): ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge2, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree2, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI2, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize2), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth2), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat2), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output2, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output2, na.rm = TRUE), ' \n'))
print(model_summary2)
sink()


# PPC ---------------------------------------------------------------------

y_rep1 <- MCMCvis::MCMCchains(fit1, params = 'y_rep')
bayesplot::ppc_dens_overlay(DATA1$y, y_rep1)

y_rep2 <- MCMCvis::MCMCchains(fit2, params = 'y_rep')
bayesplot::ppc_dens_overlay(DATA2$y, y_rep2)


# PPO ---------------------------------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Results/gr-SVC-', 
                          run_date, '/Figures')),
       dir.create(paste0(dir, 'Results/gr-SVC-', 
                         run_date, '/Figures')),
       FALSE)

setwd(paste0(dir, 'Results/gr-SVC-', run_date, '/Figures'))


#pi ~ N(120, 50)
PR <- rnorm(10000, 120, 50)
MCMCvis::MCMCtrace(fit1,
                   params = 'pi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-pi-', run_date, '.pdf'))

#nu ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit1,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-nu-', run_date, '.pdf'))

#gamma ~ N(0, 1)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit1,
                   params = 'gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-gamma-', run_date, '.pdf'))

#sigma ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit1,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-sigma-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 20)
PR_p <- rnorm(10000, 0, 20)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit1,
                   params = 'sigma_alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit1,
                   params = 'sigma_beta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-sigma_beta-', run_date, '.pdf'))

#sigma_phi ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit1,
                   params = 'sigma_phi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-forest-SVC-trace-sigma_phi-', run_date, '.pdf'))

#pi ~ N(120, 50)
PR <- rnorm(10000, 120, 50)
MCMCvis::MCMCtrace(fit2,
                   params = 'pi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-pi-', run_date, '.pdf'))

#nu ~ N(0, 5)
PR <- rnorm(10000, 0, 5)
MCMCvis::MCMCtrace(fit2,
                   params = 'nu',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-nu-', run_date, '.pdf'))

#gamma ~ N(0, 1)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit2,
                   params = 'gamma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-gamma-', run_date, '.pdf'))

#sigma ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit2,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-sigma-', run_date, '.pdf'))

#sigma_alpha ~ HN(0, 20)
PR_p <- rnorm(10000, 0, 20)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit2,
                   params = 'sigma_alpha',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-sigma_alpha-', run_date, '.pdf'))

#sigma_beta ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit2,
                   params = 'sigma_beta',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-sigma_beta-', run_date, '.pdf'))

#sigma_phi ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit2,
                   params = 'sigma_phi',
                   priors = PR,
                   pdf = TRUE,
                   open_pdf = FALSE,
                   filename = paste0('gr-all-SVC-trace-sigma_phi-', run_date, '.pdf'))


# Plot trends ----------------------------------------------------------------

plt_fun <- function(DATA, GR, FIT, MAIN = 'forest')
{
  DATA_PLOT <- data.frame(y_obs = DATA$y,
                          x_obs = DATA$year,
                          cn_id = as.numeric(factor(GR$cell)),
                          cell_lat = GR$cell_lat)
  
  c_idx <- unique(DATA_PLOT$cn_id)
  
  a_ch <- MCMCvis::MCMCchains(FIT, ISB = FALSE, params = 'alpha')
  b_ch <- MCMCvis::MCMCchains(FIT, ISB = FALSE, params = 'beta')
  
  sim_x <- seq(min(DATA_PLOT$x_obs) - 1, max(DATA_PLOT$x_obs) + 1, 
               length = 100)
  
  FIT_PLOT <- data.frame()
  for (j in 1:length(c_idx))
  {
    #j <- 1
    mf <- matrix(nrow = NROW(a_ch), ncol = length(sim_x))
    for (i in 1:length(sim_x))
    {
      #i <- 1
      mf[,i] <- a_ch[,j] + b_ch[,j] * sim_x[i]
    }
    
    med_mf <- apply(mf, 2, median)
    LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
    UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
    tdf <- data.frame(cn_id = rep(c_idx[j], 100), sim_x, med_mf, LCI_mf, UCI_mf)
    FIT_PLOT <- rbind(FIT_PLOT, tdf)
  }
  
  p <- ggplot(data = DATA_PLOT, aes(x_obs, y_obs)) +
    #model fit
    # geom_ribbon(data = FIT_PLOT,
    #             aes(x = MN_X, ymin = LCI, ymax = UCI),
    #             fill = 'grey', alpha = 0.7,
    #             inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(sim_x, med_mf, color = factor(cn_id)),
              alpha = 0.7,
              inherit.aes = FALSE,
              size = 1.4) +
    geom_point(data = DATA_PLOT, aes(x_obs, y_obs, color = factor(cn_id)),
               inherit.aes = FALSE, size = 3, alpha = 0.3) +
    theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    theme(legend.position = 'none') +
    ylab('gr') +
    xlab('Year') +
    ggtitle(paste0(MAIN, ' greenup')) +
    theme(
      plot.title = element_text(size = 22),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
  
  ggsave(paste0('gr-', MAIN, '-SVC-ind-', run_date, '.pdf'), p)
}

plt_fun(DATA = DATA1, GR = gr_forest2, FIT = fit1, MAIN = 'forest')
plt_fun(DATA = DATA2, GR = gr_all2, FIT = fit2, MAIN = 'all')


# Map of trends ------------------------------------------

map_plt_fun <- function(DATA, GR, FIT, MAIN = 'forest', GR_CELL, RNG)
{

  #estimated slope in grey, sd in white
  #extract median and sd estimates for slopes
  GR_CELL$beta_med <- MCMCvis::MCMCpstr(FIT, params = 'beta', func = median)[[1]]
  GR_CELL$beta_sd <- MCMCvis::MCMCpstr(FIT, params = 'beta', func = sd)[[1]]
  
  #load maps
  worldmap <- ggplot2::map_data("world")
  
  #min/max for plotting using output data
  # MIN <- min(GR_CELL$beta_med)
  # MAX <- max(GR_CELL$beta_med)
  MIN <- RNG[1]
  MAX <- RNG[2]
  
  #plot
  
  #transform cells to grid
  cell_grid <- dggridR::dgcellstogrid(hexgrid6, GR_CELL$cell)
  cell_grid$cell <- as.numeric(cell_grid$cell)
  
  to_plt <- dplyr::inner_join(GR_CELL, cell_grid, by = 'cell')
  
  p_beta <- ggplot(data = worldmap, aes(x = long, y = lat, 
                                        group = group)) +
    geom_polygon(fill = alpha('black', 0.3), color = NA) +
    coord_map("ortho", orientation = c(35, -90, 0),
              xlim = c(-130, -45), ylim = c(18, 66)) +
    geom_polygon(data = to_plt, aes(x = long, y = lat, 
                                    group = group, fill = beta_med), 
                 alpha = 0.7) +
    geom_path(data = to_plt, aes(x = long, y = lat, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradient2(low = 'indianred', high = 'royalblue', mid = 'lightgoldenrod',
                         limits = c(MIN, MAX), midpoint = 0) +
    # scale_fill_gradientn(colors = c('orange', 'grey', 'light blue'),
    #                        #c(hcl(h = 240, c = 35, l = 35), hcl(h = 180, c = 15, l = 92)),
    #                      limits = c(MIN, MAX)) +
    # scale_fill_gradient2(low = 'indianred2', mid = 'grey60', high = 'lightblue2', 
    #                      limits = c(MIN, MAX), midpoint = 0) +
    labs(fill = 'Slope') +
    # annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5,
    #          label = round(to_plt2$med_beta, digits = 2), col = 'black', alpha = 0.2,
    #          size = 3) +
    # annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5,
    #          label = round(to_plt2$sd_beta, digits = 2), col = 'white', alpha = 1,
    #          size = 2) +
    ggtitle(paste0(MAIN)) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')
  
  
  ggsave(paste0('gr-', MAIN, '-SVC-map-', run_date, '.pdf'), p_beta)
}


f1_beta <- MCMCvis::MCMCpstr(fit1, params = 'beta', func = median)[[1]]
f2_beta <- MCMCvis::MCMCpstr(fit2, params = 'beta', func = median)[[1]]

map_plt_fun(DATA = DATA1, GR = gr_forest2, FIT = fit1, 
            MAIN = 'forest', GR_CELL = gr_cell2,
            RNG = range(c(f1_beta, f2_beta)))
map_plt_fun(DATA = DATA2, GR = gr_all2, FIT = fit2, 
            MAIN = 'all', GR_CELL = gr_cell2,
            RNG = range(c(f1_beta, f2_beta)))



# slope plots -------------------------------------------------------------

gr_cell2$cid <- as.numeric(factor(gr_cell2$cell))
srt_cg2 <- gr_cell2[order(gr_cell2$cell_lat),]
beta_srt <- paste0('beta\\[',srt_cg2$cid, '\\]')

fit1_beta <- MCMCvis::MCMCchains(fit1, params = 'beta', mcmc.list = TRUE)
fit2_beta <- MCMCvis::MCMCchains(fit2, params = 'beta', mcmc.list = TRUE)

#cells sorted from S to N
pdf('beta.pdf', height = 5, width = 15)
MCMCvis::MCMCplot(fit1_beta, fit2_beta, params = beta_srt, ISB = FALSE, 
                  horiz = FALSE, col = rgb(0,0,0,0.5), col2 = rgb(1,0,0,0.5),
                  sz_med = 1, sz_thick = 3, sz_thin = 1, offset = 0.15)
dev.off()
