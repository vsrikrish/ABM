import itertools
import numpy as np
import warnings

# define decorator to wrap single-tuple arguments in a list
def wrap_tuple(decorated_function):
    def wrapper(*args):
        if isinstance(args[1], tuple) and len(args[1]) == 2:
            return decorated_function(args[0], [args[1]])
        else:
            return decorated_function(*args)
            
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
            
    def default_pop(self):
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
    
    def _place_agent(self, agent, position):
        """ Place an agent on the grid at position """
        x, y = position
        self.grid[x][y].add(agent.get_id())
        
    def place_agent(self, agent, position):
        """ Place an agent on the grid at position """
        self._place_agent(agent, position)
        agent.location = position
        
    def _remove_agent(self, agent_id, position):
        """ Remove an agent from grid position """
        x, y = position
        try: 
            self.grid[x][y].remove(agent_id)
        except ValueError:
            warnings.warn(f'Agent {agent_id} not in cell {position}!')
        
    def remove_agent(self, agent):
        """ Remove an agent from the grid """
        position = agent.location
        self._remove_agent(agent.get_id(), position)
        agent.location = None
        
    def move_agent(self, agent, position):
        """ Move agent from its current location to position """
        position = self.wrap_coords(position)
        self._remove_agent(agent.get_id(), agent.location)
        
        
        
    def iter_with_coords(self):
        """ Return iterator with cell contents and coordinates """
        for x in range(self.width):
            for y in range(self.height):
                yield self.grid[x][y], x, y
                
    @wrap_tuple
    def iter_cell_list_cont(self, lst):
        """ Return an iterator of contents of cells in lst """
        return (self[x][y] for x, y in lst if not self.is_empty((x, y)))
            
    def _iter_nghd(self, position, r=1, center=False):
        """ Return iterator over cells within an r-neighborhood
            of the cell referenced by position. Cells
            are referenced by their coordinate tuples. """
        
        x, y = position
        nghd = set()
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                # skip center cell if not desired
                if dx == 0 and dy == 0 and not center:
                    continue
                # skip cells exceeding Manhattan distance
                if abs(dx) + abs(dy) > r:
                    continue
                # if grid does not wrap, skip when boundary is crossed
                if not self.wrap and (not (0 <= dx + x < self.width) or not ((not (0 <= dy + y < self.height)):
                    continue
                    
                coords = self.wrap_coords((x + dx, y + dy))
                
                if coords not in nghd:
                    nghd.add(coords)
                    yield coords
        
    def nghd_list(self, position, r=1, center=False):
        """ Return list of neighboring cells """
        return list(self._iter_nghd(position, r, center))
            
    def nghd_cont_iter(self, position):
        """ Iterate over neighboring cell contents """
        nghd = self.iter_nghd(position)
        return self.iter_cell_list_cont(nghd)
        
        
    