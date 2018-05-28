import pymc3 as mc
import theano.tensor as T
import theano
from bird import BirdModel, get_pop
import numpy as np
import xarray as xr

num_years = 22
obs = xr.open_dataset('bird-pcout.nc')['pop'].values

def model_out(sup, scp, ssp):
    model = BirdModel(num_years=num_years, surv_prob=sup, scout_prob=scp, scout_surv_prob=ssp, seed=0, model_queries={'pop': get_pop})
    model.run_model()
    out = model.query_out.model_query_to_np()['pop']
    return out

with mc.Model() as model:
    # construct likelihood function
    # define priors/hyperpriors
    a = mc.Normal('a', mu=0, sd=10)
    b = mc.Lognormal('b', mu=np.log(0.9), sd=1)
#    rho = mc.Normal('rho', mu=0, sd=1)
    # model parameters
    sup = mc.Beta('surv_prob', alpha=5, beta=1)
    scp = mc.Beta('scout_prob', alpha=5, beta=3)
    ssp = mc.Beta('scout_surv_prob', alpha=8, beta=1)
    
    # evaluate model given inputs
        
#    tau = mc.Gamma('tau', alpha=10, beta=1)
#    err = mc.Normal('err', mu=rho, tau = tau)
#    sigma = mc.Gamma('sigma', alpha=2, beta=0.5)
    @theano.compile.ops.as_op(itypes=[T.dscalar, T.dscalar, T.dscalar, T.dscalar, T.dscalar], otypes=[T.dvector])
    def log_mu(a, b, sup, scp, ssp):
        return a + b * np.log(model_out(sup, scp, ssp))
    pop = mc.Poisson('pop', mu=np.exp(log_mu(a, b, sup, scp, ssp)), observed=obs)
    
with model:
    blocked_steps = [mc.Metropolis([a, b], blocked=False), mc.Metropolis([sup, scp, ssp], blocked=True)]
    samps = mc.sample(1000, tune=500, step=blocked_steps, chains=1)
#    samps = mc.step_methods.smc.sample_smc(samples=100, chains=100, homepath=test_folder)
    print(mc.summary(samps))