import abc
import datetime as dt
import random
import numpy as np
import pandas as pd

# define Model abstract class
class Model(abc.ABCMeta):
    
    """ Base class for models. """
    
    def __init__(self, seed=None, scheduler=None):
        """ create new model instance. """
        self.nagt = 0
        # initialize storage for agents
        self.agents = pd.DataFrame(columns=['uid', 'agent', 'active'])
        
        self.time = 0
        self.scheduler = scheduler
        
        # set seed for Python and numpy
        if seed is None:
            self.seed = dt.datetime.now()
        else:
            self.seed = seed
        random.seed(seed)
        np.random.seed(seed)

    @abc.abstractmethod
    def step(self):
        """ advance in time. """

# define scheduler functions
# basic scheduler which advances each active agent in turn
def basic_scheduler(model):
    if model.nagt == 0:
        raise IndexError('No agents in model!')
    active_agt = model.agents['agent'][model.agents['active']].tolist()
    for agent in active_agt:
        agent.step()
    model.time += 1
