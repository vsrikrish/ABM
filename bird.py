from agent import Agent
from model import Model, basic_scheduler
from space import Grid, wrap_tuple
import random
import pandas
import numpy

class TerritoryGrid(Grid):

    """ defines a 1d grid of territories, arranged in a ring """

    def __init__(self,
                 terr_count  = 25,
                 scout_dist = 5,
                 model):

        super().__init(self, width=terr_count, height=1, wrap=True)
        self.model = model
        self.scout_dist = scout_dist

    def contains_alpha(self, position, sex):
        """ determines if the territory contains an alpha of the
            desired sex """

        x, y = position
        return any((row['agent'].is_alpha()) and (row['agent'].sex == sex)
                   for index, row in self.model.agents[self.model.agents['uid'].isin(self[x][y])].iterrows())

    def contains_alpha_list(self, cell_iter, sex):
        """ determines if territories contain alphas of the given sex """
        return [self.contains_alpha(cell, sex) for cell in cell_iter]

    def scout_iter(self, position, direction):
        """ returns an iterator of territories for a bird to scout from its given location """
        yield

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
        self.scout_dist = scout_dist
        self.scout_surv_prob = scout_surv_prob
        self.surv_prob = surv_prob
        self.scout_prob = 0.5
        self.is_alpha = False

    def become_alpha(self):
        self.is_alpha = True

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
                 scheduler = basic_scheduler()):    # scheduling function


