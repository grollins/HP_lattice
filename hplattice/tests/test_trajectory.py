import pytest
from mock import Mock

from ..Chain import Chain
from ..Trajectory import Trajectory


@pytest.fixture
def chain():
    hpstring = 'PHPPHPHPPHH'
    initial_vec = [0, 0, 3, 0, 0, 0, 0, 0, 0, 0]
    return Chain(hpstring, initial_vec)

def test_write_coords_to_trajectory(chain):
    mock_stream_factory = Mock()
    mock_stream = Mock()
    mock_stream_factory.return_value = mock_stream
    traj = Trajectory(save_trajectory=True, trajectory_filename='temp.txt',
                   open_stream_fcn=mock_stream_factory)
    traj.snapshot(chain)
    assert mock_stream.write.call_count == 2
    traj.finalize()
    assert mock_stream.close.called
