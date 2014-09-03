import pytest
from math import sqrt

from ..Chain import Chain


@pytest.fixture
def coords():
    return Chain.Coords(num_monomers=5)

@pytest.fixture
def bad_coords():
    c = Chain.Coords(num_monomers=5)
    c.set((0,0), idx=4)
    return c

@pytest.fixture
def vec():
    return Chain.Vectors([0,0,0,1])

def test_coords_are_correct_length(coords):
    assert len(coords) == 5

def test_grow_adds_monomer_in_y_direction(coords):
    n = len(coords)
    xy = coords.get(-1)
    coords.grow()
    new_n = len(coords)
    new_xy = coords.get(-1)
    assert new_n == (n + 1)
    assert new_xy[0] == xy[0]
    assert new_xy[1] == (xy[1] + 1)

def test_pop_removes_last_monomer(coords):
    n = len(coords)
    xy = coords.get(-2)
    coords.pop()
    new_n = len(coords)
    new_xy = coords.get(-1)
    assert new_n == (n - 1)
    assert new_xy[0] == xy[0]
    assert new_xy[1] == xy[1]

def test_rotate_up_to_right(coords):
    xy = coords.get(-1)
    x, y = xy[0], xy[1]
    coords.rotate_up_to_right(-1)
    new_xy = coords.get(-1)
    assert new_xy[0] == (x + 1)
    assert new_xy[1] == (y - 1)

def test_rotate_right_to_down(coords):
    xy = coords.get(-1)
    x, y = xy[0], xy[1]
    coords.rotate_right_to_down(-1)
    new_xy = coords.get(-1)
    assert new_xy[0] == (x - 1)
    assert new_xy[1] == (y - 1)

def test_rotate_down_to_left(coords):
    xy = coords.get(-1)
    x, y = xy[0], xy[1]
    coords.rotate_down_to_left(-1)
    new_xy = coords.get(-1)
    assert new_xy[0] == (x - 1)
    assert new_xy[1] == (y + 1)

def test_convert_vec_to_coords(coords, vec):
    coords.vec2coords(vec)
    assert coords.get(0).tolist() == [0,0]
    assert coords.get(1).tolist() == [0,1]
    assert coords.get(2).tolist() == [0,2]
    assert coords.get(3).tolist() == [0,3]
    assert coords.get(4).tolist() == [1,3]

def test_overlapping_chain_is_not_viable(coords, vec, bad_coords):
    coords.vec2coords(vec)
    assert coords.is_viable() == True
    assert bad_coords.is_viable() == False

def test_compute_correct_distance_between_monomers(coords, vec):
    coords.vec2coords(vec)
    assert coords.distance_between_pts(0,4) == sqrt(10)

def test_modify_copy_does_not_modify_original(coords):
    other_coords = coords.copy()
    other_coords.set((2,2), 4)
    assert other_coords.get(4).tolist() == [2,2] 
    assert coords.get(4).tolist() == [0,0]
