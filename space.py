import itertools
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
        self.model.agents.loc[self.model.agents['uid'] == agent.get_id(), 'active'] = False
        
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

