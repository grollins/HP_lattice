import pytest
from mock import Mock
from .. import LatticeFactory


@pytest.fixture
def replica():
    lattice_factory = LatticeFactory()
    conf = lattice_factory.make_configuration()
    conf.HPSTRING = 'HHPPHH'
    conf.INITIALVEC = [0,0,1,2,1]
    conf.NREPLICAS = 1
    conf.REPLICATEMPS = [275.0]
    r = lattice_factory.make_replica(lattice_factory, conf, 0,
                                     nativeclist=[(0, 5), (1, 4)])
    return r

@pytest.fixture
def mock_chain():
    m = Mock()
    m.get_vec_length.return_value = 5
    mock_nextvec = Mock()
    mock_nextvec.__len__ = Mock()
    mock_nextvec.__len__.return_value = 5
    m.nextvec = mock_nextvec
    return m

def test_move1_does_rigid_rotation_of_N_terminus(replica, mock_chain):
    vecindex = 0
    direction = -1
    replica.mc.move1(mock_chain, vecindex=vecindex, direction=direction)
    mock_chain.do_rigid_rot.assert_called_once_with(vecindex, direction)

def test_move1_does_3bead_flip_of_nontermini(replica, mock_chain):
    vecindex = 2
    direction = -1
    replica.mc.move1(mock_chain, vecindex=vecindex, direction=direction)
    mock_chain.do_three_bead_flip.assert_called_once_with(vecindex)

def test_move1_does_rigid_rotation_of_C_terminus(replica, mock_chain):
    vecindex = 4
    direction = 1
    replica.mc.move1(mock_chain, vecindex=vecindex, direction=direction)
    mock_chain.do_rigid_rot.assert_called_once_with(vecindex, direction)

def test_move2_does_3bead_flip(replica, mock_chain):
    vecindex = 2
    direction = 1
    replica.mc.move2(mock_chain, vecindex=vecindex, direction=direction,
                     moveseed=0.1)
    mock_chain.do_three_bead_flip.assert_called_once_with(vecindex)

def test_move2_does_crankshaft(replica, mock_chain):
    vecindex = 2
    direction = 1
    replica.mc.move2(mock_chain, vecindex=vecindex, direction=direction,
                     moveseed=0.5)
    mock_chain.do_crankshaft.assert_called_once_with(vecindex)

def test_move2_does_rigid_rotation(replica, mock_chain):
    vecindex = 2
    direction = 1
    replica.mc.move2(mock_chain, vecindex=vecindex, direction=direction,
                     moveseed=0.8)
    mock_chain.do_rigid_rot.assert_called_once_with(vecindex, direction)

def test_move3_does_rigid_rotation(replica, mock_chain):
    vecindex = 0
    direction = -1
    replica.mc.move3(mock_chain, vecindex=vecindex, direction=direction)
    mock_chain.do_rigid_rot.assert_called_once_with(vecindex, direction)
