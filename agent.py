import abc

class Agent(metaclass=abc.ABCMeta):
 """ base class for agent objects """
    
    super().__init__()
    def __init__(self, uid):
        self._uid = uid
        self.location = None
        
    def get_id(self):
        return self._uid
        
    def get_location(self):
        return self.location
        
    @abc.abstractmethod
    def step(self):
        pass