import pytest

from ..Chain import Chain


@pytest.fixture
def vec():
    return Chain.Vectors([0,0,0,1])

def test_modify_copy_does_not_modify_original(vec):
    other_vec = vec.copy()
    other_vec.set_idx(idx=3, value=2)
    assert other_vec.get(3) == 2
    assert vec.get(3) == 1

def test_grow_adds_monomer_in_y_direction(vec):
    n = len(vec)
    vec.grow()
    new_n = len(vec)
    new_v = vec.get(-1)
    assert new_n == (n + 1)
    assert new_v == 0

def test_pop_removes_last_monomer(vec):
    n = len(vec)
    v = vec.get(-2)
    vec.pop()
    new_n = len(vec)
    new_v = vec.get(-1)
    assert new_n == (n - 1)
    assert new_v == v
