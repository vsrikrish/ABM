import abc
import datetime as dt
import random
import numpy as np
import pandas as pd

# define Model abstract class
class Model(metaclass=abc.ABCMeta):
    
    """ Base class for models. """
    
    def __init__(self, seed=None):
        """ create new model instance. """
        super().__init__()
        self.nagt = 0
        # initialize storage for agents
        self.agents = pd.DataFrame(columns=['uid', 'agent', 'active'])
        
        self.time = 0
        
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
    """ step is the name of the agent method """
    if model.nagt == 0:
        raise IndexError('No agents in model!')
    active_agt = model.agents['agent'][model.agents['active']].tolist()
    for agent in active_agt:
        agent.step()

# random scheduler, which advances each active agent in random order
def random_scheduler(model):
    if model.nagt == 0:
        raise IndexError('No agents in model!')
    active_agt = model.agents.loc[model.agents['active'],'agent'].tolist()
    random.shuffle(active_agt)
    for agent in active_agt:
        agent.step()