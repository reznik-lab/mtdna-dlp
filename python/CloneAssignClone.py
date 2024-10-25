#!/usr/bin/env python
# coding: utf-8

# In[6]:


"""
CloneAssignClone class
"""

import pandas as pd
import numpy as np
from CloneAssign import CloneAssign


class CloneAssignClone(CloneAssign):

    def __init__(self, num_of_clones=None, snv_allele=None, snv=None, repeat=10,
                 min_clone_assign_prob=0.8, min_clone_assign_freq=0.7, min_consensus_snv_freq=0.6,
                 max_temp=1.0, min_temp=0.5, anneal_rate=0.01, learning_rate=0.1, max_iter=400, rel_tol=5e-5, 
                 record_input_output=False):
        '''
        initialize CloneAssignClone object
        :param clone: groupings of cnv cells. (pandas.DataFrame)
        :param repeat: num of times to run cloneassign to generate consensus results. (int)
        :param min_clone_assign_prob: assign cells to a clone if clone assignment prob reaches min_clone_assign_prob (float)
        :param min_clone_assign_freq: assign cells to a clone if a min proportion of runs generate the same results (float)
        :param max_temp: starting temperature in Gumbel-Softmax reparameterization. (float)
        :param min_temp: min temperature in Gumbel-Softmax reparameterization. (float)
        :param anneal_rate: annealing rate in Gumbel-Softmax reparameterization. (float)
        :param learning_rate: learning rate of Adam optimizer. (float)
        :param max_iter: max number of iterations of elbo optimization during inference. (int)
        :param rel_tol: when the relative change in elbo drops to rel_tol, stop inference. (float)
        '''
        CloneAssign.__init__(self, num_of_clones, snv_allele, snv, repeat, min_clone_assign_prob,
                             min_clone_assign_freq, min_consensus_snv_freq, max_temp, min_temp, anneal_rate,
                             learning_rate, max_iter, rel_tol, record_input_output)
        
    def assign_cells_to_clones(self):
        '''
        assign cells from scRNA to clones based on clonal heteroplasmy
        :return: clone_assign_df (pandas.DataFrame) and gene_type_score_df (pandas.DataFrame)
        '''
        terminals = []
        if self.snv_df is not None:
            all_cells = self.snv_df.columns.values
        
        num_of_clones_input, snv_allele_input, snv_input = self.construct_input(terminals)
        
        # make columns consistent
        self.make_columns_consistent(snv_allele_input, snv_input)
        
        # run cloneassign
        snp_count = 0
        cell_count = 0
        
        snp_count = snv_input.shape[0]
        cell_count = snv_input.shape[1]
        print("snp count: " + str(snp_count))
        
        if snp_count == 0:
            raise ValueError('No valid snps exist in the matrix after filtering. Maybe loose the filtering criteria?')
        print("cell count: " + str(cell_count))
        
        # record input
        if self.record_input_output:
            self.params_dict = dict()
            self.params_dict['input'] = dict()
            self.params_dict['input']['num_of_clones'] = num_of_clones_input
            self.params_dict['input']['snv_allele'] = snv_allele_input
            self.params_dict['input']['snv'] = snv_input
        
        none_freq, clone_assign, clone_assign_df, params_dict = self.run_cloneassign_pyro_repeat(num_of_clones_input, snv_allele_input, snv_input)
        
        if self.record_input_output:
            self.params_dict['output'] = dict()
            self.params_dict['output']['none_freq'] = none_freq
            self.params_dict['output']['clone_assign'] = clone_assign
            self.params_dict['output']['clone_assign_df'] = clone_assign_df
            self.params_dict['output']['params_dict'] = params_dict
        
        cell_count = snv_input.shape[1]
        cell_names = snv_input.columns.values
        
        # record clone_assign
        for i in range(cell_count):
            if np.isnan(clone_assign.values[i]):
                self.clone_assign_dict[cell_names[i]] = np.nan
            else:
                self.clone_assign_dict[cell_names[i]] = clones[int(clone_assign.values[i])]
        
        # record clonal heteroplasmy level
        #self.allele_assign_prob_dict[hscn_input.index.values[i]] = [params_dict['mean_allele_assign_prob'][i]]

        return
