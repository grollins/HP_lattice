import pytest
from .. import LatticeFactory
from ..Replica import _compute_boltz_factor


@pytest.fixture
def replica():
    lattice_factory = LatticeFactory()
    conf = lattice_factory.make_configuration()
    conf.HPSTRING = 'HHPPHH'
    conf.INITIALVEC = [1,0,1,2,1]
    conf.NREPLICAS = 2
    conf.REPLICATEMPS = [275.0, 375.0]
    r = lattice_factory.make_replica(lattice_factory, conf, 0,
                                     nativeclist=[(0, 5), (1, 4)])
    r.init_mc_stats()
    return r

@pytest.fixture
def replica2():
    lattice_factory = LatticeFactory()
    conf = lattice_factory.make_configuration()
    conf.HPSTRING = 'HHPPHH'
    conf.INITIALVEC = [1,1,2,3,3]
    conf.NREPLICAS = 2
    conf.REPLICATEMPS = [275.0, 375.0]
    r = lattice_factory.make_replica(lattice_factory, conf, 1,
                                     nativeclist=[(0, 5), (1, 4)])
    r.init_mc_stats()
    return r

def test_mc_stats_init_to_zero(replica):
    assert replica.steps == 0
    assert replica.viablesteps == 0
    assert replica.acceptedsteps == 0
    assert replica.move_viability == 0
    assert replica.acceptance == 0

def test_record_stats_increments_counters_for_viable_moves(replica):
    replica.record_stats(move_is_viable=False, move_is_accepted=False)
    assert replica.viablesteps == 0
    assert replica.acceptedsteps == 0
    replica.record_stats(move_is_viable=True, move_is_accepted=True)
    assert replica.viablesteps == 1
    assert replica.acceptedsteps == 1

def test_compute_mc_ratios_correctly(replica):
    replica.steps = 100
    replica.viablesteps = 10
    replica.acceptedsteps = 5
    replica.compute_mc_acceptance()
    assert replica.move_viability == 0.1
    assert replica.acceptance == 0.5

def test_is_native_returns_true_when_contacts_are_native(replica, replica2):
    assert replica.is_native() is False
    assert replica2.is_native() is True

def test_compute_boltz_factor(replica, replica2):
    '''
    Should always swap if doing so would drop the lower energy structure
    to a lower temperature. In this test, replica2 is a lower energy structure
    at a higher temperature, so it should swap with replica. A boltz factor of
    less than one would mean that there's a chance they won't swap.
    '''
    bf = _compute_boltz_factor(replica, replica2)
    assert bf >= 1
