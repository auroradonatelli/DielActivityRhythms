
# Load packages
library(rstan)
library(tidyverse)

# Upload dataset available at https://doi.org/10.5061/dryad.d51c5b0fd
df = read_csv("dataset.csv")

bear_lev=length(unique(df$id)) # nu. of individual bears
year_lev=length(unique(df$year)) # nu. of years nested in population

# counts of individual bears in each population
n_ap=length(unique(df$id[df$study==1])) # Apennine
n_cr=length(unique(df$id[df$study==2])) # Dinaric Pindos (Croatia)
n_se=length(unique(df$id[df$study==3])) # Dinaric Pindos (Serbia)
n_sw=length(unique(df$id[df$study==4])) # Scandinavian
n_fi=length(unique(df$id[df$study==5])) # Karelian
n_ye=length(unique(df$id[df$study==6])) # Yellowstone

# Specify model:
modelString = "
  data {
    int N;      ## number of rows in dataset 
    int n_ap;
    int n_cr;
    int n_se;
    int n_sw;
    int n_fi;
    int n_ye;
    int bear_lev;
    int year_lev;
    int year[N];   ## year nested in population
    real vel[N];   ## response variable, i.e., hourly movement rates
    int hour[N];   ## hour of the day
    int id[N];     ## bear id
    int study[N];  ## bear population
    int season[N]; ## season
    real index[N]; ## GHM index, i.e., human encroachment
    real ndvi[N];  ## NDVI
    real tmx[N];   ## maximum daily ambient temperature
    real dayle[N]; ## daylength  
    vector[24] onevec;
  }
  parameters {
    vector[24] h[6,3];      ## hour of day effect in interaction with season and population
    vector[bear_lev] i;     ## random effect of bear id
    vector[year_lev] y;     ## random effect of year nested in population
    vector[24] ind[6,3];    ## GHM effect in interaction with season and population
    vector[24] nd[6,3];     ## NDVI effect in interaction with season and population
    vector[24] tx[6,3];     ## maximum daily temperature effect in interaction with season and population
    vector[24] dl[6];       ## daylength effect in interaction with population
    real alpha;             ## intercept   
    real<lower=0> rigma;    ## model variance
    real<lower=0> sigma_ap; ## sigma_xx are variances associated with each population for random effects of bear id
    real<lower=0> sigma_cr;
    real<lower=0> sigma_se;
    real<lower=0> sigma_sw;
    real<lower=0> sigma_fi;
    real<lower=0> sigma_ye;
    real<lower=0> sigma_y;  ## variance associated with random effect of year
    real<lower=0.125,upper=3> nu_11;  ## parameters in covariance matrices
    real<lower=0.125,upper=3> nu_21;
    real<lower=0.125,upper=3> nu_31;
    real<lower=0.125,upper=3> nu_41;
    real<lower=0.125,upper=3> nu_51;
    real<lower=0.125,upper=3> nu_61;
    real<lower=0.125,upper=3> nu_12;
    real<lower=0.125,upper=3> nu_22;
    real<lower=0.125,upper=3> nu_32;
    real<lower=0.125,upper=3> nu_42;
    real<lower=0.125,upper=3> nu_52;
    real<lower=0.125,upper=3> nu_62;
    real<lower=0.125,upper=3> nu_13;
    real<lower=0.125,upper=3> nu_23;
    real<lower=0.125,upper=3> nu_33;
    real<lower=0.125,upper=3> nu_43;
    real<lower=0.125,upper=3> nu_53;
    real<lower=0.125,upper=3> nu_63;
    real<lower=0> psi_11;   ## parameters in covariance matrices
    real<lower=0> psi_21;
    real<lower=0> psi_31;
    real<lower=0> psi_41;
    real<lower=0> psi_51;
    real<lower=0> psi_61;
    real<lower=0> psi_12;
    real<lower=0> psi_22;
    real<lower=0> psi_32;
    real<lower=0> psi_42;
    real<lower=0> psi_52;
    real<lower=0> psi_62;
    real<lower=0> psi_13;
    real<lower=0> psi_23;
    real<lower=0> psi_33;
    real<lower=0> psi_43;
    real<lower=0> psi_53;
    real<lower=0> psi_63;
  }
    transformed parameters {
    cov_matrix[24] Sigma_11;  ## covariance matrices for the hour of day effect for each population and season
    cov_matrix[24] Sigma_21;
    cov_matrix[24] Sigma_31;
    cov_matrix[24] Sigma_41;
    cov_matrix[24] Sigma_51;
    cov_matrix[24] Sigma_61;
    cov_matrix[24] Sigma_12;
    cov_matrix[24] Sigma_22;
    cov_matrix[24] Sigma_32;
    cov_matrix[24] Sigma_42;
    cov_matrix[24] Sigma_52;
    cov_matrix[24] Sigma_62;
    cov_matrix[24] Sigma_13;
    cov_matrix[24] Sigma_23;
    cov_matrix[24] Sigma_33;
    cov_matrix[24] Sigma_43;
    cov_matrix[24] Sigma_53;
    cov_matrix[24] Sigma_63;

    
    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_11[k1,k2] = exp(-nu_11*min(k2-k1, 24+k1-k2))/psi_11;
    Sigma_11[k2,k1] = Sigma_11[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_11[k1,k1] = 1/psi_11;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_21[k1,k2] = exp(-nu_21*min(k2-k1, 24+k1-k2))/psi_21;
    Sigma_21[k2,k1] = Sigma_21[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_21[k1,k1] = 1/psi_21;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_31[k1,k2] = exp(-nu_31*min(k2-k1, 24+k1-k2))/psi_31;
    Sigma_31[k2,k1] = Sigma_31[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_31[k1,k1] = 1/psi_31;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_41[k1,k2] = exp(-nu_41*min(k2-k1, 24+k1-k2))/psi_41;
    Sigma_41[k2,k1] = Sigma_41[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_41[k1,k1] = 1/psi_41;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_51[k1,k2] = exp(-nu_51*min(k2-k1, 24+k1-k2))/psi_51;
    Sigma_51[k2,k1] = Sigma_51[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_51[k1,k1] = 1/psi_51;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_61[k1,k2] = exp(-nu_61*min(k2-k1, 24+k1-k2))/psi_61;
    Sigma_61[k2,k1] = Sigma_61[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_61[k1,k1] = 1/psi_61;
    }



    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_12[k1,k2] = exp(-nu_12*min(k2-k1, 24+k1-k2))/psi_12;
    Sigma_12[k2,k1] = Sigma_12[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_12[k1,k1] = 1/psi_12;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_22[k1,k2] = exp(-nu_22*min(k2-k1, 24+k1-k2))/psi_22;
    Sigma_22[k2,k1] = Sigma_22[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_22[k1,k1] = 1/psi_22;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_32[k1,k2] = exp(-nu_32*min(k2-k1, 24+k1-k2))/psi_32;
    Sigma_32[k2,k1] = Sigma_32[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_32[k1,k1] = 1/psi_32;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_42[k1,k2] = exp(-nu_42*min(k2-k1, 24+k1-k2))/psi_42;
    Sigma_42[k2,k1] = Sigma_42[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_42[k1,k1] = 1/psi_42;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_52[k1,k2] = exp(-nu_52*min(k2-k1, 24+k1-k2))/psi_52;
    Sigma_52[k2,k1] = Sigma_52[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_52[k1,k1] = 1/psi_52;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_62[k1,k2] = exp(-nu_62*min(k2-k1, 24+k1-k2))/psi_62;
    Sigma_62[k2,k1] = Sigma_62[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_62[k1,k1] = 1/psi_62;
    }



    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_13[k1,k2] = exp(-nu_13*min(k2-k1, 24+k1-k2))/psi_13;
    Sigma_13[k2,k1] = Sigma_13[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_13[k1,k1] = 1/psi_13;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_23[k1,k2] = exp(-nu_23*min(k2-k1, 24+k1-k2))/psi_23;
    Sigma_23[k2,k1] = Sigma_23[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_23[k1,k1] = 1/psi_23;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_33[k1,k2] = exp(-nu_33*min(k2-k1, 24+k1-k2))/psi_33;
    Sigma_33[k2,k1] = Sigma_33[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_33[k1,k1] = 1/psi_33;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_43[k1,k2] = exp(-nu_43*min(k2-k1, 24+k1-k2))/psi_43;
    Sigma_43[k2,k1] = Sigma_43[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_43[k1,k1] = 1/psi_43;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_53[k1,k2] = exp(-nu_53*min(k2-k1, 24+k1-k2))/psi_53;
    Sigma_53[k2,k1] = Sigma_53[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_53[k1,k1] = 1/psi_53;
    }


    for(k1 in 1:23){
    for(k2 in (k1+1):24){
    Sigma_63[k1,k2] = exp(-nu_63*min(k2-k1, 24+k1-k2))/psi_63;
    Sigma_63[k2,k1] = Sigma_63[k1,k2];
    }
    }
    for(k1 in 1:24){
    Sigma_63[k1,k1] = 1/psi_63;
    }


  }
  model {

    rigma ~ inv_gamma(1,1);
    
    nu_11 ~ uniform(0.125,3);
    nu_21 ~ uniform(0.125,3);
    nu_31 ~ uniform(0.125,3);
    nu_41 ~ uniform(0.125,3);
    nu_51 ~ uniform(0.125,3);
    nu_61 ~ uniform(0.125,3);
    nu_12 ~ uniform(0.125,3);
    nu_22 ~ uniform(0.125,3);
    nu_32 ~ uniform(0.125,3);
    nu_42 ~ uniform(0.125,3);
    nu_52 ~ uniform(0.125,3);
    nu_62 ~ uniform(0.125,3);
    nu_13 ~ uniform(0.125,3);
    nu_23 ~ uniform(0.125,3);
    nu_33 ~ uniform(0.125,3);
    nu_43 ~ uniform(0.125,3);
    nu_53 ~ uniform(0.125,3);
    nu_63 ~ uniform(0.125,3);

    psi_11 ~ gamma(1,1);
    psi_21 ~ gamma(1,1);
    psi_31 ~ gamma(1,1);
    psi_41 ~ gamma(1,1);
    psi_51 ~ gamma(1,1);
    psi_61 ~ gamma(1,1);
    psi_12 ~ gamma(1,1);
    psi_22 ~ gamma(1,1);
    psi_32 ~ gamma(1,1);
    psi_42 ~ gamma(1,1);
    psi_52 ~ gamma(1,1);
    psi_62 ~ gamma(1,1);
    psi_13 ~ gamma(1,1);
    psi_23 ~ gamma(1,1);
    psi_33 ~ gamma(1,1);
    psi_43 ~ gamma(1,1);
    psi_53 ~ gamma(1,1);
    psi_63 ~ gamma(1,1);

    alpha ~ normal(0, sqrt(1000));

      h[1,1,] ~ multi_normal(alpha*onevec, Sigma_11);
      h[2,1,] ~ multi_normal(alpha*onevec, Sigma_21);
      h[3,1,] ~ multi_normal(alpha*onevec, Sigma_31);
      h[4,1,] ~ multi_normal(alpha*onevec, Sigma_41);
      h[5,1,] ~ multi_normal(alpha*onevec, Sigma_51);
      h[6,1,] ~ multi_normal(alpha*onevec, Sigma_61);
      h[1,2,] ~ multi_normal(alpha*onevec, Sigma_12);
      h[2,2,] ~ multi_normal(alpha*onevec, Sigma_22);
      h[3,2,] ~ multi_normal(alpha*onevec, Sigma_32);
      h[4,2,] ~ multi_normal(alpha*onevec, Sigma_42);
      h[5,2,] ~ multi_normal(alpha*onevec, Sigma_52);
      h[6,2,] ~ multi_normal(alpha*onevec, Sigma_62);
      h[1,3,] ~ multi_normal(alpha*onevec, Sigma_13);
      h[2,3,] ~ multi_normal(alpha*onevec, Sigma_23);
      h[3,3,] ~ multi_normal(alpha*onevec, Sigma_33);
      h[4,3,] ~ multi_normal(alpha*onevec, Sigma_43);
      h[5,3,] ~ multi_normal(alpha*onevec, Sigma_53);
      h[6,3,] ~ multi_normal(alpha*onevec, Sigma_63);

    sigma_ap ~ inv_gamma(1,1);
    sigma_cr ~ inv_gamma(1,1);
    sigma_se ~ inv_gamma(1,1);
    sigma_sw ~ inv_gamma(1,1);
    sigma_fi ~ inv_gamma(1,1);
    sigma_ye ~ inv_gamma(1,1);

    i[1 : n_ap] ~ normal(0, sqrt(sigma_ap));
    i[(n_ap+1) : (n_ap+n_cr)] ~ normal(0, sqrt(sigma_cr));
    i[(n_ap+n_cr+1) : (n_ap+n_cr+n_se)] ~ normal(0, sqrt(sigma_se));
    i[(n_ap+n_cr+n_se+1) : (n_ap+n_cr+n_se+n_sw)] ~ normal(0, sqrt(sigma_sw));
    i[(n_ap+n_cr+n_se+n_sw+1) : (n_ap+n_cr+n_se+n_sw+n_fi)] ~ normal(0, sqrt(sigma_fi));
    i[(n_ap+n_cr+n_se+n_sw+n_fi+1) : (n_ap+n_cr+n_se+n_sw+n_fi+n_ye)] ~ normal(0, sqrt(sigma_ye));

    ind[1,1,1:24] ~ normal(0, sqrt(1000));
    ind[1,2,1:24] ~ normal(0, sqrt(1000));
    ind[1,3,1:24] ~ normal(0, sqrt(1000));
    ind[2,1,1:24] ~ normal(0, sqrt(1000));
    ind[2,2,1:24] ~ normal(0, sqrt(1000));
    ind[2,3,1:24] ~ normal(0, sqrt(1000));
    ind[3,1,1:24] ~ normal(0, sqrt(1000));
    ind[3,2,1:24] ~ normal(0, sqrt(1000));
    ind[3,3,1:24] ~ normal(0, sqrt(1000));
    ind[4,1,1:24] ~ normal(0, sqrt(1000));
    ind[4,2,1:24] ~ normal(0, sqrt(1000));
    ind[4,3,1:24] ~ normal(0, sqrt(1000));
    ind[5,1,1:24] ~ normal(0, sqrt(1000));
    ind[5,2,1:24] ~ normal(0, sqrt(1000));
    ind[5,3,1:24] ~ normal(0, sqrt(1000));
    ind[6,1,1:24] ~ normal(0, sqrt(1000));
    ind[6,2,1:24] ~ normal(0, sqrt(1000));
    ind[6,3,1:24] ~ normal(0, sqrt(1000));

    nd[1,1,1:24] ~ normal(0, sqrt(1000));
    nd[1,2,1:24] ~ normal(0, sqrt(1000));
    nd[1,3,1:24] ~ normal(0, sqrt(1000));
    nd[2,1,1:24] ~ normal(0, sqrt(1000));
    nd[2,2,1:24] ~ normal(0, sqrt(1000));
    nd[2,3,1:24] ~ normal(0, sqrt(1000));
    nd[3,1,1:24] ~ normal(0, sqrt(1000));
    nd[3,2,1:24] ~ normal(0, sqrt(1000));
    nd[3,3,1:24] ~ normal(0, sqrt(1000));
    nd[4,1,1:24] ~ normal(0, sqrt(1000));
    nd[4,2,1:24] ~ normal(0, sqrt(1000));
    nd[4,3,1:24] ~ normal(0, sqrt(1000));
    nd[5,1,1:24] ~ normal(0, sqrt(1000));
    nd[5,2,1:24] ~ normal(0, sqrt(1000));
    nd[5,3,1:24] ~ normal(0, sqrt(1000));
    nd[6,1,1:24] ~ normal(0, sqrt(1000));
    nd[6,2,1:24] ~ normal(0, sqrt(1000));
    nd[6,3,1:24] ~ normal(0, sqrt(1000));

    tx[1,1,1:24] ~ normal(0, sqrt(1000));
    tx[1,2,1:24] ~ normal(0, sqrt(1000));
    tx[1,3,1:24] ~ normal(0, sqrt(1000));
    tx[2,1,1:24] ~ normal(0, sqrt(1000));
    tx[2,2,1:24] ~ normal(0, sqrt(1000));
    tx[2,3,1:24] ~ normal(0, sqrt(1000));
    tx[3,1,1:24] ~ normal(0, sqrt(1000));
    tx[3,2,1:24] ~ normal(0, sqrt(1000));
    tx[3,3,1:24] ~ normal(0, sqrt(1000));
    tx[4,1,1:24] ~ normal(0, sqrt(1000));
    tx[4,2,1:24] ~ normal(0, sqrt(1000));
    tx[4,3,1:24] ~ normal(0, sqrt(1000));
    tx[5,1,1:24] ~ normal(0, sqrt(1000));
    tx[5,2,1:24] ~ normal(0, sqrt(1000));
    tx[5,3,1:24] ~ normal(0, sqrt(1000));
    tx[6,1,1:24] ~ normal(0, sqrt(1000));
    tx[6,2,1:24] ~ normal(0, sqrt(1000));
    tx[6,3,1:24] ~ normal(0, sqrt(1000));

    dl[1,1:24] ~ normal(0, sqrt(1000));
    dl[2,1:24] ~ normal(0, sqrt(1000));
    dl[3,1:24] ~ normal(0, sqrt(1000));
    dl[4,1:24] ~ normal(0, sqrt(1000));
    dl[5,1:24] ~ normal(0, sqrt(1000));
    dl[6,1:24] ~ normal(0, sqrt(1000));

    sigma_y ~ inv_gamma(1,1);
    y ~ normal(0, sqrt(sigma_y));

    for(t in 1:N){
      vel[t] ~ normal(h[study[t],season[t],hour[t]] + ind[study[t],season[t],hour[t]]*index[t] +
                nd[study[t],season[t],hour[t]]*ndvi[t] +
                tx[study[t],season[t],hour[t]]*tmx[t] +
                i[id[t]] + dl[study[t],hour[t]]*dayle[t] + y[year[t]],
                sqrt(rigma));
    }
  }
" # close quote for modelString

# Translate model to C++ and compile to DSO:
stanDso <- stan_model( model_code=modelString ) 

# Specify data:
N=nrow(df)
onevec = array(1, dim=c(24))
dataList = list(
  vel = df$movem_60,
  hour = df$cor_hour, 
  id=df$id,
  study=df$study,
  season=df$season,
  year=df$year,
  index=df$index,
  ndvi=df$ndvi,
  tmx=df$tmx,
  dayle = df$dayle,
  bear_lev=bear_lev, 
  year_lev = year_lev,
  N=N, 
  n_ap=n_ap,
  n_cr=n_cr,
  n_se=n_se,
  n_sw=n_sw,
  n_fi=n_fi,
  n_ye=n_ye,
  onevec=onevec 
)

# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 2 ,
                     cores = 2, 
                     iter = 5000 , 
                     warmup = 1000 ,
                     thin = 2 )

# save results
save.image("workspace.Rdata")
