from agent import Agent
from model import Model, random_scheduler
from space import Grid, wrap_tuple
from query import Query
import random
import numpy as np
from scipy.stats import genextreme
import operator
import itertools
import collections

# model query functions
def get_vac(model):
    """ gets the number of vacant lots when called """
    return sum([model.grid.is_empty((x, y))
                for ids, x, y in model.grid.iter_with_coords()
                ])
                
def get_states(model):
    """ gets the states (vacant or occupied) of each lot when called"""
    states = np.zeros((model.grid.width, model.grid.height), dtype=np.bool)
    for ids, x, y in model.grid.iter_with_coords():
        states[x][y] = not model.grid.is_empty((x, y))
    return states
                    
def get_floods(model):
    """ gets the number of floods over the previous 10 years for each lot """
    floods = np.zeros((model.grid.width, model.grid.height), dtype=np.int8)
    for ids, x, y in model.grid.iter_with_coords():
        floods[x][y] = model.grid.num_floods((x, y), 10)
    return floods

class HousingGrid(Grid):
    
    """ defines a grid of houses along a river """
    
    # rows are the number of housing rows
    # cols the number of housing columns
    # river_params is a dictionary with the GEV parameters for extreme water levels, with values "loc", "scale", and "shape"
    # kwargs should include at least one of "elev" (structure of elevations with dimension equal to (cols, rows)) or "slope" (dictionary of parameters to compute stochastic elevations, with entries "offset" and "mean" slope)
    def __init__(self, model, rows, cols, n_yr, **kwargs):
        
        super().__init__(width=cols, height=rows) # initialize superclass object
        self.model = model
        # if housing elevations are passed, save those
        elev = kwargs.get('elev', None)
        if elev is not None:
            elev = np.array(elev)
            if elev.shape == (cols, rows):
                self.elev = elev
            else:
                raise ValueError("Elevation array must have dimensions (cols, rows)!")
            
        else:
            slope = kwargs.get('slope', None)
            try:
                # calculate flood elevation heights for each house and initialize
                heights = []
                for row in range(rows):
                    base_height = row * slope['mean'] + slope['offset']
                    heights.append(base_height +
                                   np.random.normal(0, 1, cols))
                self.elev = np.stack(heights, axis=1)
            except TypeError:
                raise ValueError('Must pass valid elevations or slope!')
    
        # initialize flood history storage
        self.floods = np.empty((cols, rows, n_yr+50), dtype='bool')
        
    # get flood elevation for a particular parcel
    def get_elev(self, position):
        col, row = position
        return self.elev[col][row]

    # get number of floods for a particular pacel within the previous num_yr years
    def num_floods(self, position, num_yr):
        col, row = position
        eidx = self.model.time-self.model.start+51
        sidx = eidx-num_yr
        return sum(self.floods[col][row][sidx:eidx])
    
    # add to flood history for a parcel
    def update_history(self, position, flood):
        col, row = position
        self.floods[col][row].append(flood)
        
        
class Resident(Agent):
    
    """ defines an agent residing in a house """
    def __init__(self, uid, model, position,
                 move_params, mem_length, res_length):
        super().__init__(uid, model)
        self.location = position
        self.move_params = move_params
        self.mem_length = mem_length
        self.res_length = res_length
        
    def move(self):
        x, y = self.location
        # the number of years in an agent's memory should be no longer
        # than the agent's memory length, but no less than 5 years or
        # the amount of time the agent has lived in the house
        mem_years = min(max(self.res_length, 5), self.mem_length)
        n_flood = self.model.grid.num_floods(self.location, mem_years)
        x = self.move_params['int'] + \
            self.move_params['self_coef'] * n_flood / mem_years
        u = np.random.random()
        return u < 1/(1+np.exp(-x))
        
    def step(self):
        self.res_length += 1
        if self.move():
            self.model.grid.remove_agent(self)
            
class BasicFloodModel(Model):
    
    """ defines a simple model structure for flooding and migration with synthetic data history"""
    
    def __init__(self, grid_sz, move_params,
                 mem_length, seed, **kwargs):
        super().__init__(seed=seed)
        
        # set model time
        self.start = kwargs.get('start_year', 2018)
        self.end = kwargs.get('end_year', 2100)
        self.time = self.start

        # construct housing grid
        cols, rows = grid_sz    # unpack grid size tuple
        # get housing elevations or slope to compute them from kwargs
        elev = kwargs.get('elev', None)
        slope = kwargs.get('slope', None)
        if elev is not None:
            self.grid = HousingGrid(model=self, rows=rows, cols=cols,
                                    n_yr=self.end-self.time+1, elev=elev)
        elif slope is not None:
            self.grid = HousingGrid(model=self, rows=rows, cols=cols,
                                    n_yr=self.end-self.time+1, slope=slope)
        else:
            raise ValueError('Need to pass housing elevations or slope')
        
        # store agent parameters
        self.mem_length = mem_length
        self.move_params = move_params
        
        # get either river height elevation history or parameters for river height elevation samples
        self.river_params = kwargs.get('river_params', None)
        self.river_hist = kwargs.get('river_hist', None)
        if self.river_params is None and self.river_hist is None:
            raise ValueError('Need to pass river parameters or history!')
        if self.river_hist is not None:
            self.river_hist = np.array(self.river_hist)
        
        # initialize model queries
        model_queries = kwargs.get('query', {})
        self.query = Query(model_queries=model_queries)
        
        # set up 50 years of flood history
        try:
            water_elev = self.river_hist[0:50]
        except TypeError:
            if self.river_params is not None:
                water_elev = genextreme.rvs(c=self.river_params['shape'],
                                            loc=self.river_params['loc'],
                                            scale=self.river_params['scale'],       size=50)
            else:
                raise ValueError('No valid way to get water elevations: need river parameters or history!')
        floods = [self.grid.elev < hgt for hgt in water_elev]
        self.grid.floods[:, :, 0:50] = np.transpose(np.array(floods),
                                                    (1, 2, 0))
                
        # initialize agents. We assume they have fixed sensitivity to floods for now and that there is a fixed probability of each house being occupied
        p_occ = kwargs.get('p_occ', 0.99)
        for row in range(rows):
            for col in range(cols):
                if np.random.random() < p_occ:
                    self.create_agent((col, row), np.random.randint(5, 50))
                
    def create_agent(self, position, res_length):
        uid = random.choice(list(set(range(0, 10001))-self.agent_ids))
        agent = Resident(model=self, uid=uid, position=position,
                         move_params={k:v
                                      for k, v in self.move_params.items()
                                      if k not in ['fill_prob']
                                     },
                         mem_length=self.mem_length, res_length=res_length)
        self.grid.place_agent(agent, position)
        self.agents[uid] = agent
        self.agent_ids.add(uid)
        self.nagt += 1
        
    def step(self):
        """ advances time step by 1 """
        # generate extreme water level
        try:
            water_elev = self.river_hist[self.time-self.start+50]
        except TypeError:
            if self.river_params is not None:
                water_elev = genextreme.rvs(c=self.river_params['shape'],
                                            loc=self.river_params['loc'],
                                            scale=self.river_params['scale'],
                                            size=1)
            else:
                raise ValueError('No valid way to get water elevations: need river parameters or history!')

        # see if any houses flood
        floods = self.grid.elev < water_elev
        # update flood history
        self.grid.floods[:, :, self.time-self.start+50] = floods
        # see if any agents move in or out
        for ids, col, row in self.grid.iter_with_coords():
            if len(ids) == 0:
                if np.random.random() < self.move_params['fill_prob']:
                    self.create_agent((col, row), 0)
            else:
                id_list = list(ids)
                for id in id_list:
                    agent = self.agents[id]
                    agent.step()
        # run model queries
        self.query.query_model(self)
        # advance time
        self.time += 1

    def run(self):
        while self.time <= self.end:
            self.step()
 
