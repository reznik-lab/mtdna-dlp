#!/usr/bin/env python
# coding: utf-8

# In[6]:


"""
CloneAssign base class: initialize and clean data frame
"""

import torch
import numpy as np
import pandas as pd
from torch.nn import Softplus

import pyro
import pyro.optim
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.infer.autoguide.initialization import init_to_sample
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate
from pyro.ops.indexing import Vindex


def max_count(s):
    '''
    calculate the number of appearances of the most common item in pandas.Series
    :param s: (pandas.Series)
    :return: (int)
    '''
    return s.value_counts(dropna=False).max()


def inverse_softplus(x):
    '''
    inverse the softplus function
    :param x: number matrix (torch.tensor)
    :return: number matrix (torch.tensor)
    '''
    return x + torch.log(-torch.expm1(-x))


class CloneAssign():
    
    def construct_input(self, terminals):
        '''
        clean up mutations that are not informative
        :return: clone_df, snv_allele, snv
        '''
        
        snv = self.snv_df
        snv_allele = self.snv_allele_df
        num_of_clones = self.num_of_clones
        
        return num_of_clones, snv_allele, snv
    
    # inplace method
    def average_param_dict(self, param_dict):
        if param_dict is not None:
            for key in param_dict:
                param_dict[key] = sum(param_dict[key])/len(param_dict[key])
    
    # inplace method            
    def max_param_dict(self, param_dict):
        if param_dict is not None:
            for key in param_dict:
                param_dict[key] = max(param_dict[key])
    
    # inplace method
    def make_columns_consistent(self, *args):
        intersect_columns = None
        for arg in args:
            if hasattr(arg, 'columns'):
                if intersect_columns is None:
                    intersect_columns = arg.columns
                else:
                    intersect_columns = intersect_columns.intersection(arg.columns)
        for arg in args:
            if hasattr(arg, 'columns'):
                arg.drop(columns=[col for col in arg if col not in intersect_columns], inplace=True)
    
    def convert_df_to_torch(self, df):
        if df is not None:
            df_torch = torch.tensor(df.values, dtype=torch.float64)
            df_torch = torch.transpose(df_torch, 0, 1)
            return df_torch
        else:
            return None
    
    def generate_output(self):
        '''
        generate clone_assign_df and gene_type_score_df
        :return: cloneassign output (pandas.DataFrame)
        '''
        if self.clone_assign_df is None:
            self.clone_assign_df = pd.DataFrame.from_dict(self.clone_assign_dict.items())
            self.clone_assign_df.rename(columns={0: "cell_id", 1: "clone_id"}, inplace=True)

        return self.clone_assign_df

    def __init__(self, num_of_clones=None, snv_allele=None, snv=None, repeat=10,
                 min_clone_assign_prob=0.8, min_clone_assign_freq=0.7,min_consensus_snv_freq=0.6,
                 max_temp=1.0, min_temp=0.5, anneal_rate=0.01, learning_rate=0.1, max_iter=400, rel_tol=5e-5, 
                 record_input_output=False):
        '''
        initialize CloneAssign object
        :param repeat: num of times to run cloneassign to generate consensus results. (int)
        :param min_clone_assign_prob: assign cells to a clone if clone assignment prob reaches min_clone_assign_prob (float)
        :param min_clone_assign_freq: assign cells to a clone if a min proportion of runs generate the same results (float)
        :param max_temp: starting temperature in Gumbel-Softmax reparameterization. (float)
        :param min_temp: min temperature in Gumbel-Softmax reparameterization. (float)
        :param anneal_rate: annealing rate in Gumbel-Softmax reparameterization. (float)
        :param learning_rate: learning rate of Adam optimizer. (float)
        :param max_iter: max number of iterations of elbo optimization during inference. (int)
        :param rel_tol: when the relative change in elbo drops to rel_tol, stop inference. (float)
        :param record_input_output: record input output before and after cloneassign runs. (bool)
        '''
        
        self.map_estimates = None
        self.num_of_clones = num_of_clones
        self.snv_allele_df = snv_allele
        self.snv_df = snv
        
        self.min_consensus_snv_freq = min_consensus_snv_freq
        self.repeat = repeat

        self.min_clone_assign_prob = min_clone_assign_prob
        self.min_clone_assign_freq = min_clone_assign_freq

        self.max_temp = max_temp
        self.min_temp = min_temp
        self.anneal_rate = anneal_rate
        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.record_input_output = record_input_output
        
        self.rel_tol = rel_tol

        # output
        self.clone_assign_dict = dict()
        self.params_dict = dict()
        
        self.clone_assign_df = None

    @config_enumerate
    def cloneassign_pyro_model(num_of_clones, snv_allele, snv, temperature=0.5):
        '''
        original cloneassign model
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
            with clone_plate:
                clone_het = pyro.sample('expose_clone_het', dist.Beta(torch.ones(1), torch.ones(1)))

        with pyro.plate('cell', num_of_cells):

            # draw clone_assign_prob from Dir
            clone_assign_prob = pyro.sample('expose_clone_assign_prob', dist.Dirichlet(torch.ones(num_of_clones) * 1))

            # draw clone_assign from Cat
            clone_assign = pyro.sample('clone_assign', dist.Categorical(clone_assign_prob))

            # draw the observed variant alleles based on total snv count (observed) and heteroplasmy for this cell
            pyro.sample('het', dist.Binomial(snv, clone_het[clone_assign]).to_event(1), obs=snv_allele)
            
    def run_cloneassign_pyro(self, num_of_clones, snv_allele, snv, current_repeat):
        '''
        cloneassign inference
        :return: clone assignment, gene_type_score
        '''
        
        np_temp = self.max_temp
        losses = []

        optim = pyro.optim.Adam({'lr': self.learning_rate, 'betas': [0.8, 0.99]})
        elbo = TraceEnum_ELBO(max_plate_nesting=1)

        model = self.cloneassign_pyro_model

        def initialize(seed):
            global global_guide, svi
            pyro.set_rng_seed(seed)
            pyro.clear_param_store()
            global_guide = AutoDelta(poutine.block(model, expose_fn=lambda msg: msg["name"].startswith("expose_")), init_loc_fn=init_to_sample)
            svi = SVI(model, global_guide, optim, loss=elbo)
            return svi.loss(model, global_guide, snv_allele, snv, np_temp)

        loss, seed = min((initialize(seed), seed) for seed in range(current_repeat * 100, (current_repeat + 1) * 100))
        initialize(seed)
        print('seed = {}, initial_loss = {}'.format(seed, loss))

        # start inference
        print('Start Inference.')
        for i in range(self.max_iter):
            if i % 100 == 1:
                np_temp = np.maximum(self.max_temp * np.exp(-self.anneal_rate * i), self.min_temp)
            loss = svi.step(snv_allele, snv, np_temp)

            if i >= 1:
                loss_diff = abs((losses[-1] - loss) / losses[-1])
                if loss_diff < self.rel_tol:
                    print('ELBO converged at iteration ' + str(i))
                    break

            losses.append(loss)

        map_estimates = global_guide(snv_allele, snv)

        # also record inferred parameters
        self.map_estimates = map_estimates
        results = dict()
        for entry in map_estimates:
          entry_name = entry[3:]
          results[entry_name] = pd.DataFrame(map_estimates[entry].data.numpy())

        return results
    
    def run_cloneassign_pyro_repeat(self, num_of_clones, snv_allele_df, snv_df):
        '''
        call run_cloneassign_pyro() multiple times to generate consensus results
        :return: frequency of unassigned cells, clone assignment, gene_type_score
        '''
        
        torch.set_default_dtype(torch.float64)
        snv_allele = self.convert_df_to_torch(snv_allele_df)
        snv = self.convert_df_to_torch(snv_df)

        clone_assign_list = []
        other_params = dict()

        for i in range(self.repeat):
            current_results = self.run_cloneassign_pyro(num_of_clones, snv_allele, snv, i)

            current_clone_assign = current_results['clone_assign_prob']
            current_clone_assign_prob = current_clone_assign.max(1)
            current_clone_assign_discrete = current_clone_assign.idxmax(1)

            current_clone_assign_discrete[current_clone_assign_prob < self.min_clone_assign_prob] = np.nan
            clone_assign_list.append(current_clone_assign_discrete)

            for param_name in current_results:
                if param_name != 'clone_assign_prob':
                    if param_name not in other_params:
                        other_params[param_name] = []
                    other_params[param_name].append(current_results[param_name].iloc[:, 0])

        mean_params = dict()
        for param_name in other_params:
            other_params[param_name] = pd.concat(other_params[param_name], axis=1)
            mean_params['mean_' + param_name] = other_params[param_name].mean(1)

        other_params.update(mean_params)

        clone_assign = pd.concat(clone_assign_list, axis=1)
        clone_assign_max = clone_assign.mode(1, dropna=False)[0]
        clone_assign_freq = clone_assign.apply(max_count, axis=1) / self.repeat

        clone_assign_max[clone_assign_freq < self.min_clone_assign_freq] = np.nan

        # calculate NaN freq
        none_freq = clone_assign_max.isna().sum() / clone_assign_max.shape[0]

        return none_freq, clone_assign_max, clone_assign, other_params


    
    
    
    



