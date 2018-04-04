import abc

class Agent(metaclass=abc.ABCMeta):

    """ base class for agent objects """

    def __init__(self, uid, model):
        super().__init__()
        self._uid = uid
        self.location = None
        self.model = model
        
    def get_id(self):
        return self._uid
        
    def set_id(self, uid):
        self._uid = uid
        
    def get_location(self):
        return self.location
        
    @abc.abstractmethod
    def step(self):
        pass