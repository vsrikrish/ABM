import abc

class Agent(abc.ABCMeta):
 """ base class for agent objects """
    
    def __init__(self, uid):
        self._uid = uid
        
    def get_id(self):
        return self.uid
        
    @abc.abstractmethod
    def step(self):
        pass