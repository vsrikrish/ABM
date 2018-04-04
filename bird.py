from agent import Agent
from model import Model, random_scheduler
from space import Grid, wrap_tuple
import random
import pandas
import numpy as np
import operator


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

        x, y = position
        return any([(row['agent'].is_alpha) and (row['agent'].sex == sex) for index, row in self.model.agents[self.model.agents['uid'].isin(self[x][y])].iterrows()])

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
        return (row['agent'] for index, row in self.model.agents[self.model.agents['uid'].isin(self[x][y])].iterrows())

class Bird(Agent):

    """ Birds are the model agents. They can be male or female, and alpha or non-alpha.
        Non-alpha birds will look for new territories via scouting where they can be
        alphas. """

    def __init__(self,
                 model,                     # model instance
                 uid,                       # unique ID
                 location,                  # location of cell bird spawns in
                 sex = None,                # sex of bird; can be 'M' or 'F' or None (will be randomly assigned)
                 age = 0,                   # age in months of the bird at initialization
                 scout_surv_prob = 0.8,     # probability a bird survives a scouting trip
                 surv_prob = 0.99,          # probably of surviving each month
                 scout_prob = 0.5):         # probability of going on a scouting trip

        super().__init__(uid, model)
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

    def become_alpha(self):
        self.is_alpha = True
        
    def scout_trip(self):
        # if the bird survives the scouting trip,
        # scout. Otherwise, remove from the grid.
        if self.scout_surv_prob >= random.uniform(0,1):
            scout_dir = random.sample([-1, 1], 1)
            scout_iter = self.model.grid.scout_iter(self.location, scout_dir)
            for position in scout_iter:
                # if any of the scouted territories does not have an alpha
                # of the same sex as the bird, the bird stays in the same position.
                # Otherwise, the bird moves to the first territory found with 
                if not self.model.grid.contains_alpha(position, self.sex):
                    self.model.grid.move_agent(self, position)
                    break
        else:
            self.model.grid.remove_agent(self)
            self.model.agents.loc[self.model.agents['uid'] == self.get_id(),'active'] = False
            
    def step_age(self):
        self.age += 1
        
    def step(self):
        """ advances the agent one step, which involves scouting """
        if not self.is_alpha:
            ages = [(bird.get_id(), bird.age, bird.is_alpha) for bird in self.model.grid.iter_bird(self.location) if bird.sex == self.sex]
            oldest_bird = max([a for a in ages if not a[2]], key=operator.itemgetter(1))[0]
            if oldest_bird != self.get_id():
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
                 scheduler = random_scheduler):    # scheduling function

        super().__init__(seed)
        self.t_end = num_years * 12 # convert number of years to monthly timesteps
        
        self.fecundity = fecundity
        self.scout_surv_prob = scout_surv_prob
        self.scout_prob = scout_prob
        self.surv_prob = surv_prob
        self.scout_dist = scout_dist
        self.agents = self.agents.reindex(
            columns=['uid', 'agent', 'sex', 'age', 'active'])
        self.scheduler = scheduler
        self.grid = TerritoryGrid(model=self, terr_count=num_territory, 
                                    scout_dist=self.scout_dist)
        # initialize bird population
        for x in range(num_territory):
            for sex in ['M', 'F']:
                self.create_agents(n_agents=2, position=(x, 0),
                                   sex=sex, age=[random.randint(0, 24) for i in range(2)])
                # the oldest new agent of each sex becomes the alpha
                new_agts = self.agents[(self.agents['uid'].isin(self.grid[x][0])) & (self.agents['sex'] == sex)]
                alpha = new_agts.loc[new_agts['age'].idxmax()]
                alpha.agent.become_alpha()
                
            
    def create_agents(self, n_agents, position, sex=None, age=0):
        # increment number of total agents
        if sex is None:
            sex = random.choices(['M', 'F'], n_agents)
        elif isinstance(sex, str):
            sex = [sex] * n_agents
        if not isinstance(age, list):
            age = [age] * n_agents
        for i in range(n_agents):
            agent = Bird(model=self, uid=self.nagt, location=position, sex=sex[i], age=age[i], 
                    scout_surv_prob=self.scout_surv_prob, scout_prob=self.scout_prob, surv_prob=self.surv_prob)
            self.grid.place_agent(agent, position)
            self.agents = self.agents.append({'uid': agent.get_id(), 'agent': agent, 'sex': agent.sex, 'age': agent.age, 'active': True}, ignore_index=True)
            self.nagt += 1
                    
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
            self.agents.loc[self.agents['uid'].isin(ids),'age'] += 1
            for agent in self.agents.loc[self.agents['uid'].isin(ids), 'agent'].tolist():
                agent.step_age()
            for sex in['M', 'F']:
                if not self.grid.contains_alpha((x,y), sex):
                    cands = self.agents.loc[(self.agents['uid'].isin(ids)) & (self.agents['sex'] == sex),'age']
                    if cands.empty:
                        continue
                    else:
                        new_alpha = cands.idxmax()
                        if self.agents.at[new_alpha, 'age'] >= 12:
                            self.agents.at[new_alpha, 'agent'].become_alpha()
        # step each bird (scouting and movement) according to the scheduler
        self.scheduler(self)
        # in last month of each year, alpha females reproduce if there is an alpha male
        if self.time % 12 == 0:
            for ids, x, y in self.grid.iter_with_coords():
                if self.grid.contains_alpha((x,y), 'M') and self.grid.contains_alpha((x,y), 'F'):
                    self.reproduce((x,y))
        # birds experience mortality
        for index, row in self.agents[self.agents['active']].iterrows():
            if self.surv_prob <= random.uniform(0, 1):
                self.grid.remove_agent(self.agents.at[index, 'agent'])
                self.agents.at[index, 'active'] = False
                
    def run_model(self):
        """ model to invoke to run model from start to finish """
        for t in range(self.t_end):
            self.step()
