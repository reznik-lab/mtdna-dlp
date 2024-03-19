import pyro
import pyro.distributions as dist
from pyro.infer.autoguide import AutoDelta, AutoNormal
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate, infer_discrete, Trace_ELBO
import torch
from torch.distributions import constraints
from pyro.ops.indexing import Vindex
from pyro import poutine
from pyro.infer.autoguide.initialization import init_to_sample
import math
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from CloneAssignClone import CloneAssignClone

## Parse necessary arguments
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-s", "--sample_id",type=str, help="Sample ID")
parser.add_argument("-d", "--dataset",type=str, help="Dataset")
parser.add_argument("-c","--num_of_clones",type=int,help="Number of clones, default = 3",default = 3)

## read in arguments
args = parser.parse_args()
sample_id = args.sample_id
dataset=args.dataset
num_of_clones = args.num_of_clones

## load SNV input
# Alternate allele count matrix from scRNA (row represents a SNV, column represents a cell)
snv_allele = pd.read_csv(dirpath + sample_id + '_allele.csv')
cell_ids = snv_allele.columns.values

# Total count matrix at SNV from scRNA (row represents a SNV, column represents a cell)
snv_total = pd.read_csv(dirpath + sample_id + '_total.csv')

# Reformatting the data into tensors
cell_ids = snv_allele.columns.values
snv_allele = torch.tensor(snv_allele.values).T
snv_total = torch.tensor(snv_total.values).T

@config_enumerate
def cloneassign_pyro_model(num_of_clones, snv_allele, snv, temperature=0.5):
    '''
    clone assign model
    :param num_of_clones: integer
    :param snv_allele: torch.tensor
    :param snv: torch.tensor
    :param temperature: float
    :return: None
    '''
    
    num_of_cells, num_of_snps = snv_allele.shape
    assert snv_allele.shape == snv.shape
    
    snv_plate = pyro.plate('snv', num_of_snps, dim=-1)
    clone_plate = pyro.plate('clone', num_of_clones, dim=-2)
    
    with snv_plate:
        snv_het = pyro.sample('expose_snv_het', dist.Beta(torch.ones(1), torch.ones(1) * 10))
        
        with clone_plate:
            clone_het = pyro.sample('expose_clone_het', dist.Beta(torch.ones(1), torch.ones(1)))
        
        # draw snv_weights from uniform distribution.
        snv_weights_prob = pyro.sample('expose_snv_weights_prob', dist.Dirichlet(torch.ones(2) * 1))
        
        # the score reflects how much the SNV influences clonal assignment
        snv_weights = pyro.sample('expose_snv_weights', dist.RelaxedOneHotCategorical(
                                                              temperature=torch.tensor(temperature),
                                                              probs=snv_weights_prob))
        
    with pyro.plate('cell', num_of_cells):
        # Draw clone_assign_prob from Dirichlet
        clone_assign_prob = pyro.sample('expose_clone_assign_prob', dist.Dirichlet(torch.ones(num_of_clones) * 1))
        
        # Draw clone_assign from Categorical
        clone_assign = pyro.sample('clone_assign', dist.Categorical(clone_assign_prob))
        
        # Add the noise
        expected_clone_het = (snv_het * clone_het[clone_assign] * snv_weights[:,0] + snv_het * snv_weights[:,1])
        
        # Draw the observed variant alleles based on total snv count (observed) and heteroplasmy for this cell
        pyro.sample('het', dist.Binomial(snv, expected_clone_het).to_event(1), obs=snv_allele)

optim = pyro.optim.Adam({'lr': 0.01, 'betas': [0.8, 0.99]})
elbo = TraceEnum_ELBO(max_plate_nesting=2)

guide = AutoDelta(poutine.block(cloneassign_pyro_model, expose_fn=lambda msg: msg["name"].startswith("expose_")))
svi = SVI(cloneassign_pyro_model, guide, optim, loss=elbo)

losses = []
for i in range(1000):
    loss = svi.step(num_of_clones, snv_allele, snv_total)
    
    if i >= 300:
        loss_diff = abs((losses[-1] - loss) / losses[-1])
        if loss_diff < 5e-5:
            break
    
    losses.append(loss)
    print('.' if i % 200 else '\n', end='')
    
map_estimates = guide(num_of_clones, snv_allele, snv_total)

for k in pyro.get_param_store():
    print(k, pyro.get_param_store()[k].detach())

# Replay model
guide_trace = poutine.trace(guide).get_trace(num_of_clones, snv_allele, snv_total)
trained_model = poutine.replay(cloneassign_pyro_model, trace=guide_trace)

# Infer discrete sites and get model trace
inferred_model = infer_discrete(
    trained_model, temperature=0,
    first_available_dim=-3)
trace = poutine.trace(inferred_model).get_trace(num_of_clones, snv_allele, snv_total)

# Extract fitted parameters
clone_het = trace.nodes['expose_clone_het']['value'].detach().numpy()
clone_assign_prob = trace.nodes['expose_clone_assign_prob']['value'].detach().numpy()
clone_assign = trace.nodes['clone_assign']['value'].detach().numpy()
snv_weights = trace.nodes['expose_snv_weights']['value'].detach().numpy()

# Output the clone labels result
output = pd.DataFrame((cell_ids,clone_assign), index = ['cell_id','clone']).transpose()
output.to_csv(outputpath + sample_id + '_clonelabels.tsv', index = False, sep='\t')

# Output the SNV weights result
snv_allele = pd.read_csv(dirpath + sample_id + '_allele.csv')
output = pd.DataFrame(snv_weights)
output.index = snv_allele.index.values
output.to_csv(outputpath + sample_id + '_snv_weights.tsv', index = True, sep='\t')

