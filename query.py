import pandas as pd
import warnings
from collections import defaultdict

class Query(object):
    """ class for querying agents and models """
    
    def __init__(self, 
                 model_queries = None, 
                 agent_queries = None):
                 
        """ queries are model or agent functions
            which run to gather data or output. """
                 
        super().__init__()
        self.model_queries = {}
        self.agent_queries = {}
        
        self.model_vals = {}
        self.agent_vals = {}
        
        self.model_qtimes = []
        self.agent_qtimes = []
        
        # if queries are passed, read them in
        if model_queries is not None:
            for name, query in model_queries.items():
                self._new_model_query(name, query)
        if agent_queries is not None:
            for name, query in agent_queries.items():
                self._new_agent_query(name, query)
                
    @staticmethod
    def _query_attribute(attribute):
        """ returns a function which queries the
            named attribute. """
            
        def attribute_query(obj):
            return getattr(obj, attribute)
            
        return attribute_query
                
    def _new_model_query(self, name, query):
        """ add new model query """
        if isinstance(query, str):
            query = _query_attribute_(query)
        self.model_queries[name] = query
        self.model_vals[name] = []
        
    def _new_agent_query(name, query):
        """ add new agent query """
        if isinstance(query, str):
            query = _query_attribute(query)
        self.agent_queries[name] = query
        self.agent_vals[name] = []
            
    def query_model(self, model):
        """ when called, runs the model-level queries """
        if self.model_queries:
            for name, query in self.model_queries.items():
                self.model_vals[name].append(query(model))
            self.model_qtimes.append(model.time)
        else:
            warnings.warn("No model queries specified!")
                
    def query_agents(self, model):
        """ when called, runs the agent-level queries """
        if self.agent_queries:
            for name, query in self.agent_queries.items():
                for index, row in model.agents.loc[model.agents['active']].iterrows():
                    self.agent_vals[name].append((model.agents.at[index, 'uid'],
                        query(model.agents.at[index, 'agent'])))
            self.agent_qtimes.append(model.time)
        else:
            warnings.warn("No agent queries specified!")    

    def model_query_to_df(self):
        """ convert the dictionary of model queries to a pandas dataframe.
            The dataframe has one column for each variable, and another for the
            model times of query. """
        df = pd.DataFrame(self.model_vals)
        df['time'] = pd.Series(model.qtimes)
        return df
        
    def agent_query_to_df(self):
        """ convert the dictionary of agent queries to a pandas dataframe.
            The dataframe has one column for each variable, and is indexed by
            agent unique id and model times of query. """
        data = defaultdict(dict)
        for name, vals in self.agent_vals.items():
            for time, entry in enumerate(vals):
                agent_id = entry[0]
                val = entry[1]
                data[(time, agent_id)][name] = val
        df = pd.DataFrame_from_dict(data, orient='index')
        df.index.names = ["Time", "UID"]
        return df
