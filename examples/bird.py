import sys

sys.path.append('/storage/work/vxs914/ABM/scripts')
from agent import Agent
from model import Model, random_scheduler
from space import Grid, wrap_tuple
from query import Query
import random
import numpy as np
import operator
import itertools

# model query functions
def get_pop(model):
    """ gets the bird population of the model when called """
    return len(set.union(*(ids for ids, x, y in model.grid.iter_with_coords())))   
    
def get_vacancies(model):
    """ gets the percentage of territories missing one or both alphas """
    m_alpha_list = []
    f_alpha_list = []
    for ids, x, y in model.grid.iter_with_coords():
        m_alpha_list.append(model.grid.contains_alpha((x,y), 'M'))
        f_alpha_list.append(model.grid.contains_alpha((x,y), 'F'))
    alpha_vacancy = [not (m and f) for m, f in zip(m_alpha_list, f_alpha_list)]
    return sum(alpha_vacancy) / (model.grid.height * model.grid.width)

class TerritoryGrid(Grid):

    """ defines a 1d grid of territories, arranged in a ring """

    def __init__(self,
                 model,
                 terr_count  = 25,
                 scout_dist = 5):

        super().__init__(width=terr_count, height=1, wrap=True)
        self.model = model
        self.scout_dist = scout_dist

    def contains_alpha(self, position, sex):
        """ determines if the territory contains an alpha of the
            desired sex """

        return any([b.is_alpha for b in self.iter_bird(position) if b.sex == sex])
        
    def contains_alpha_list(self, cell_iter, sex):
        """ determines if territories contain alphas of the given sex """
        return [self.contains_alpha(cell, sex) for cell in cell_iter]

    def scout_iter(self, position, direction):
        """ returns an iterator of territories for a bird to scout from its given location """
        x, y = position
        for dx in range(1, self.scout_dist + 1):
            if self.out_of_bounds((x+dx, 0)):
                coords = self.wrap_coords((x+dx, 0))
            else:
                coords = (x+dx, 0)
            yield coords
            
    def iter_bird(self, position):
        """ returns an iterator over the birds in the territory """
        x,y = position
        return (self.model.agents[id] for id in self.grid[x][y])
        
class Bird(Agent):

    """ Birds are the model agents. They can be male or female, and alpha or non-alpha.
        Non-alpha birds will look for new territories via scouting where they can be
        alphas. """

    def __init__(self,
                 model,                     # model instance
                 grid,                      # grid instance
                 uid,                       # unique ID
                 location,                  # location of cell bird spawns in
                 sex = None,                # sex of bird; can be 'M' or 'F' or None (will be randomly assigned)
                 age = 0,                   # age in months of the bird at initialization
                 scout_surv_prob = 0.8,     # probability a bird survives a scouting trip
                 surv_prob = 0.99,          # probably of surviving each month
                 scout_prob = 0.5):         # probability of going on a scouting trip

        super().__init__(uid, model)
        self.grid = grid
        self.location = location
        if sex is not None:
            self.sex = sex
        else:
            self.sex = random.sample(('M', 'F'), 1)
        self.age = age
        self.scout_surv_prob = scout_surv_prob
        self.surv_prob = surv_prob
        self.scout_prob = 0.5
        self.is_alpha = False
        
    def scout_trip(self):
        # if the bird survives the scouting trip,
        # scout. Otherwise, remove from the grid.
        if self.scout_surv_prob >= random.uniform(0,1):
            scout_dir = random.sample([-1, 1], 1)
            scout_iter = self.grid.scout_iter(self.location, scout_dir)
            for position in scout_iter:
                # if any of the scouted territories does not have an alpha
                # of the same sex as the bird, the bird stays in the same position.
                # Otherwise, the bird moves to the first territory found with 
                if not self.grid.contains_alpha(position, self.sex):
                    self.grid.move_agent(self, position)
                    break
        else:
            self.grid.remove_agent(self)
                                
    def step(self):
        """ advances the agent one step, which involves scouting """
        x, y = self.location
        if not self.is_alpha:
            # if the bird is the oldest subordinate bird
            try: 
                oldest_sub = max([self.model.agents[id] for id in self.grid[x][y] if (not self.model.agents[id].is_alpha) and (self.model.agents[id].sex == self.sex)], key=operator.attrgetter('age'))
            except ValueError:
                oldest_sub = None
            if oldest_sub != self:
                if random.uniform(0, 1) <= self.scout_prob:
                    self.scout_trip()            

class BirdModel(Model):

    """ Bird ecology model from Thiele et al (2014) """

    def __init__(self,
                 num_years = 22,                    # number of years for simulation
                 num_territory = 25,                # number of territories in grid
                 fecundity = 2,                     # number of offspring for each alpha pair
                 scout_dist = 5,                    # distance a bird can scout
                 scout_surv_prob = 0.8,             # probability of surviving a scouting trip
                 surv_prob = 0.99,                  # probability of surviving each month
                 scout_prob = 0.5,                  # probability of going on a scouting trip
                 num_initial_agents = 2,            # initial number of birds of each sex in each territory
                 seed = None,                       # RNG seed
                 scheduler = random_scheduler,      # scheduling function
                 model_queries = {'pop': get_pop,
                                  'vacancy': get_vacancies}):   # model reporting functions

        super().__init__(seed)
        self.t_end = num_years * 12 # convert number of years to monthly timesteps
        
        self.fecundity = fecundity
        self.scout_surv_prob = scout_surv_prob
        self.scout_prob = scout_prob
        self.surv_prob = surv_prob
        self.scout_dist = scout_dist
        self.scheduler = scheduler
        self.grid = TerritoryGrid(model=self, terr_count=num_territory, 
                                    scout_dist=self.scout_dist)
                                    
        self.query_out = Query(model_queries=model_queries)
        
        # initialize bird population
        for x in range(num_territory):
            for sex in ['M', 'F']:
                new_agents = self.create_agents(n_agents=2, position=(x, 0),
                                   sex=sex, age=[random.randint(0, 24) for i in range(2)])
                # the oldest new agent of each sex becomes the alpha
                alpha = max(new_agents, key=operator.attrgetter('age'))
                alpha.is_alpha = True
                
            
    def create_agents(self, n_agents, position, sex=None, age=0):
        # increment number of total agents
        new_agents = [None] * n_agents
        if sex is None:
            sex = [random.choice(['M', 'F']) for i in range(n_agents)]
        elif isinstance(sex, str):
            sex = [sex] * n_agents
        if not isinstance(age, list):
            age = [age] * n_agents
        for i in range(n_agents):
            uid = random.choice(list(set(range(0, 10001))-self.agent_ids))
            agent = Bird(model=self, grid=self.grid, uid=uid, location=position, sex=sex[i], age=age[i], 
                    scout_surv_prob=self.scout_surv_prob, scout_prob=self.scout_prob, 
                    surv_prob=self.surv_prob)
            new_agents[i] = agent
            self.grid.place_agent(agent, position)
            self.agents[uid] = agent
            self.agent_ids.add(uid)
            self.nagt += 1
        return new_agents
                    
    def reproduce(self, position):
        """ if a territory contains alphas of both sexes, they reproduce """
        if self.grid.contains_alpha(position, 'M') and self.grid.contains_alpha(position, 'F'):
            sex = [random.choice(['M', 'F']) for i in range(self.fecundity)]
            self.create_agents(self.fecundity, position, sex=sex, age=0)
            
    def step(self):
        """ advances time step by one """
        self.time += 1  # update time
        # update bird ages and fill alphas for each territory
        for ids, x, y in self.grid.iter_with_coords():
            for agent in self.grid.iter_bird((x,y)):
                agent.age += 1
            for sex in['M', 'F']:
                if not self.grid.contains_alpha((x,y), sex):
                    try:
                        new_alpha = max([bird for bird in self.grid.iter_bird((x, y)) if bird.age >= 12 and bird.sex == sex], key=operator.attrgetter('age'))
                        new_alpha.is_alpha = True
                    except ValueError:
                        continue
        
        # step each bird (scouting and movement) according to the scheduler
        self.scheduler(self)
        
        # make observations in November of each year
        if self.time % 12 == 11:
            self.query_out.query_model(self)
        
        # in last month of each year, alpha females reproduce if there is an alpha male
        if self.time % 12 == 0:
            for ids, x, y in self.grid.iter_with_coords():
                if self.grid.contains_alpha((x,y), 'M') and self.grid.contains_alpha((x,y), 'F'):
                    self.reproduce((x,y))
        # birds experience mortality
        for id, agent in self.agents.items():
            if (agent.location is not None) and (self.surv_prob < random.uniform(0, 1)):
                self.grid.remove_agent(agent)
                
    def run_model(self):
        """ model to invoke to run model from start to finish """
        for t in range(self.t_end):
            self.step()

if __name__ == "main":
    import pandas as pd
    import xarray as xr
    
    m = bird.BirdModel(seed=0)
    m.run_model()
    d = m.query_out['pop']
    print(d)
#    ds = xr.Dataset.from_dataframe(df)
#    ds.to_netcdf('bird-pcout.nc')
