from .Chain import Chain
from .Config import Config
from .Trajectory import Trajectory
from .Monty import Monty
from .Replica import Replica


class LatticeFactory(object):
    """
    *LatticeFactory* objects orchestrate the creation of other objects
    in the HP model, namely chains, replicas, monte carlo samplers,
    trajectories, and configurations. 

    :param callable chain_cls: optional, chain class
    :param callable replica_cls: optional, replica class
    :param callable monty_cls: optional, monte carlo sampler class
    :param callable trajectory_cls: optional, trajectory class
    :param callable conf_cls: optional, HP configuration class

    """
    def __init__(self, chain_cls=Chain, replica_cls=Replica, monty_cls=Monty,
                 trajectory_cls=Trajectory, conf_cls=Config):
        self.chain_cls = chain_cls
        self.replica_cls = replica_cls
        self.monty_cls = monty_cls
        self.trajectory_cls = trajectory_cls
        self.conf_cls = conf_cls

    def make_chain(self, *args, **kwargs):
        """
        Make a chain object.
        """
        return self.chain_cls(*args, **kwargs)

    def make_replica(self, *args, **kwargs):
        """
        Make a replica object.
        """
        return self.replica_cls(*args, **kwargs)

    def make_monty(self, *args, **kwargs):
        """
        Make a monte carlo sampler object.
        """
        return self.monty_cls(*args, **kwargs)

    def make_trajectory(self, *args, **kwargs):
        """
        Make a trajectory object.
        """
        return self.trajectory_cls(*args, **kwargs)

    def make_configuration(self, *args, **kwargs):
        """
        Make a configuration object.
        """
        return self.conf_cls(*args, **kwargs)
