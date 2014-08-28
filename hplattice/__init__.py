from .Chain import Chain
from .Config import Config
from .Trajectory import Trajectory
from .Monty import Monty
from .Replica import Replica


class LatticeFactory(object):
    """docstring for LatticeFactory"""
    def __init__(self, chain_cls=Chain, replica_cls=Replica, monty_cls=Monty,
                 trajectory_cls=Trajectory, conf_cls=Config):
        self.chain_cls = chain_cls
        self.replica_cls = replica_cls
        self.monty_cls = monty_cls
        self.trajectory_cls = trajectory_cls
        self.conf_cls = conf_cls

    def make_chain(self, *args, **kwargs):
        return self.chain_cls(*args, **kwargs)

    def make_replica(self, *args, **kwargs):
        return self.replica_cls(*args, **kwargs)

    def make_monty(self, *args, **kwargs):
        return self.monty_cls(*args, **kwargs)

    def make_trajectory(self, *args, **kwargs):
        return self.trajectory_cls(*args, **kwargs)

    def make_configuration(self, *args, **kwargs):
        return self.conf_cls(*args, **kwargs)
