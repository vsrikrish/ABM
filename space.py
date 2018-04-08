import itertools
import numpy as np
import operator
import warnings

# define decorator to wrap single-tuple arguments in a list
def wrap_tuple(func):
    def wrapper(*args):
        if isinstance(args[1], tuple) and len(args[1]) == 2:
            return func(args[0], [args[1]])
        else:
            return func(*args)
            
    return wrapper
    

class Grid(object):
    """ base class for spatial grids """
    
    def __init__(self, width, height, wrap=False):
        super().__init__()
        self.width = width
        self.height = height
        self.wrap = wrap
        
        self.grid = []
        
        for x in range(self.width):
            col = []
            for y in range(self.height):
                col.append(self.default_pop())
            self.grid.append(col)

    @staticmethod
    def default_pop():
        """ Default population for cell elements """
        return set()
        
    def __getitem__(self, position):
        """ Get grid element by position """
        return self.grid[position]   
    
    def __iter__(self):
        """ Define iterator over grid cells """
        return itertools.chain(*self.grid)
                
    def wrap_coords(self, position):
        """ Handle wrapping around the boundary """
        if not self.out_of_bounds(position):
            return position
        elif not self.wrap:
            raise IndexError('Desired cell out of bounds')
        else:
            return position[0] % self.width, position[1] % self.height
            
    def is_empty(self, position):
        """ Return True if cell has no contents """
        x, y = position
        return True if self.grid[x][y] == self.default_pop() else False

    def out_of_bounds(self, position):
        """ Test if cell position is off of the grid """
        x, y = position
        return (not ((-1 < x < self.width) and (-1 < y < self.height)))
    
    def _place_agent(self, agent_id, position):
        """ Place an agent on the grid at position """
        x, y = position
        self.grid[x][y].add(agent_id)
        
    def place_agent(self, agent, position):
        """ Place an agent on the grid at position """
        self._place_agent(agent.get_id(), position)
        agent.location = position
        
    def _remove_agent(self, agent_id, position):
        """ Remove an agent from grid position """
        x, y = position
        try: 
            self.grid[x][y].remove(agent_id)
        except ValueError:
            warnings.warn('Agent {} not in cell {}!'.format(agent_id, position))
        
    def remove_agent(self, agent):
        """ Remove an agent from the grid """
        position = agent.location
        self._remove_agent(agent.get_id(), position)
        agent.location = None
        
    def move_agent(self, agent, position):
        """ Move agent from its current location to position """
        position = self.wrap_coords(position)
        self._remove_agent(agent.get_id(), agent.location)
        self._place_agent(agent.get_id(), position)
        agent.location = position

    def iter_with_coords(self):
        """ Return iterator with cell contents and coordinates """
        for x in range(self.width):
            for y in range(self.height):
                yield self.grid[x][y], x, y
                
    @wrap_tuple
    def iter_cell_list_contents(self, cell_lst):
        """ Return an iterator of contents of cells in cell_lst """
        return (self[x][y] for x, y in cell_lst if not self.is_empty((x, y)))

    @wrap_tuple
    def cell_list_contents(self, cell_lst):
        """ Return a list of contents of cells in cell_lst """
        return list(self.iter_cell_list_contents(cell_lst))
            
    def _iter_nghd(self, position, r=1, center=False):
        """ Return iterator over cells within an r-neighborhood
            of the cell referenced by position. Cells
            are referenced by their coordinate tuples. """
        
        x, y = position
        nghd = set()
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                # skip center cell if not desired
                if (dx == 0) and (dy == 0) and (not center):
                    continue
                # skip cells exceeding Manhattan distance
                if abs(dx) + abs(dy) > r:
                    continue
                # if grid does not wrap, skip when boundary is crossed
                if (not self.wrap) and self.out_of_bounds((x + dx, y + dy)):
                    continue

                coords = self.wrap_coords((x + dx, y + dy))
                
                if coords not in nghd:
                    nghd.add(coords)
                    yield coords
        
    def nghd_list(self, position, r=1, center=False):
        """ Return list of neighboring cells """
        return list(self._iter_nghd(position, r, center))
            
    def iter_nghd_cont(self, position):
        """ Iterate over neighboring cell contents """
        nghd = self._iter_nghd(position)
        return self.iter_cell_list_contents(nghd)

def ContinuousDomain(object):
    """ Continuous domain, where agents can have arbitrary
        positions instead of grid coordinates. The domain 
        object stores agent positions as a numpy array """
        
    def __init__(self, x=(0,1), y=(0,1), wrap=False):
        """ x and y are tuples giving lower and upper bounds
            on those coordinates. Wrap is a boolean for whether
            the grid wraps around as a torus. """
            
        super().__init__()
        self.xlim = x
        self.ylim = y
        
        self.width = x[1] - x[0]
        self.height = y[1] - y[0]
        self.wrap = wrap
        
        # initialize agent storage array
        self._agent_pos = np.zeros(shape=(0,2))
        
        # initialize list mapping position indices and agent ids 
        self._idx_to_agent = []
        
    def out_of_bounds(self, position):
        """ determine if position is out of bounds """
                x, y = position
        return (not ((self.xlim[0] <= x <= self.xlim[1]) and (self.ylim[0] < y < self.ylim[1])))
        
    def wrap_coords(self, position):
        if not self.out_of_bounds(position):
            return position
        elif not self.wrap:
            raise IndexError('Position out of bounds!')
        else:
            new_x = self.xlim[0] + (position[0] - self.xlim[0]) % self.width
            new_y = self.ylim[0] + (position[1] - self.ylim[0]) % self.height
            if isinstance(position, tuple):
                return (new_x, new_y)
            elif isinstance(position, np.ndarray):
                return np.array((new_x, new_y))
                
    def _place_agent(self, agent_id, position):
        """ private method to assign a location to the agent id """
        position = wrap_coords(position)
        # add position to array
        self._agent_pos = np.vstack((self._agent_pos, np.array(position)))
        # store index-agent mapping
        self._idx_to_agent.append((self_agent_pos.shape[0] - 1, agent_id))
        
    def place_agent(self, agent, position):
        """ place an agent in the model """
        # assign position to agent id
        self._place_agent(agent.get_id(), position)
        # set agent location
        agent.location = position
        
    def _remove_agent(self, agent_id, position):
        """ private method to remove an agent id from the model """
        # get index in position array 
        try:
            lidx, pidx = zip(*[i, v[0] for i, v in enumerate(self._idx_to_agent) if v[1] == agent_id])
        except ValueError:
            raise ValueError('Agent {} is not placed in the model!'.format(agent_id))
        if len(lidx) > 1:
            raise ValueError('Agent {} is placed multiple times!'.format(agent_id))
        else:
            # remove agent from position array
            self._agent_pos = np.delete(self._agent_pos, pidx, axis=0)
            # remove agent from index-agent tuple list
            del self._idx_to_agent[lidx]
    
    def remove_agent(self, agent, position):
        """ remove an agent from the model """
        self._remove_agent(agent.get_id(), position)
        agent.location = None
    
    def _move_agent(self, agent_id, position):
        """ private method to move agent """
        idx = [v[0] for v in self._idx_to_agent if v[1] == agent_id]
        self._agent_pos[idx] = position
    
    def move_agent(self, agent, position):
        """ move agent from its current position
            to the specified position """
        position = self.wrap_coords(position)
        self._move_agent(agent.get_id(), position)
        agent.location = position
        
    def iter_neighbors(self, position, r, center = True):
        """ returns an iterator of the
            agent ids within radius r from
            position. if center is True, this will
            include any agents at the coordinates
            specified by position. """
        position = wrap_coords(position)
        dists = np.linalg.norm(self._agent_pos - position)
        idx, = np.where(dists <= r ** 2)
        return (v[1] for v in self._idx_to_agent 
                if (v[0] is in idx) and (center or dists[v[0]] > 0))
            