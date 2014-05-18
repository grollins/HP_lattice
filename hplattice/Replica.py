from string import zfill

from .Chain import Chain
from .Monty import Monty


class Replica:
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self, config, repnum):
        """Initialize the Replica() object."""
        temp = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.repname = 'rep' + zfill( str(repnum), 2 )
        self.repdir = config.DATADIR + self.repname
        self.chain = Chain(config)
        self.mc = Monty(config, temp, self.chain)    

