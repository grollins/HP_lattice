import pytest

from ..Chain import Chain


@pytest.fixture
def chain1():
    hpstring = 'PHPPHPHPPHH'
    initial_vec = [0, 0, 3, 0, 0, 0, 0, 0, 0, 0]
    return Chain(hpstring, initial_vec)

@pytest.fixture
def chain2():
    hpstring = 'HPPH'
    initial_vec = [1, 2, 3]
    return Chain(hpstring, initial_vec)

@pytest.fixture
def chain3():
    hpstring = 'HPPH'
    initial_vec = [0, 1, 2]
    return Chain(hpstring, initial_vec)

def test_chain_is_correct_length(chain1):
    assert len(chain1) == 11
    assert chain1.get_vec_length() == 10

def test_convert_hpstring_to_binary(chain1):
    binstr = chain1.hpstr2bin()
    assert ''.join([str(b) for b in binstr]) == '01001010011'

def test_correctly_identify_HH_contacts(chain2):
    contacts = chain2.contactstate()
    assert len(contacts) == 1
    assert contacts[0] == (0, 3)

def test_compute_correct_energy(chain2):
    E, contacts = chain2.energy(epsilon=-2.0)
    assert E == -2.0
    assert len(contacts) == 1
    assert contacts[0] == (0, 3)

def test_rigid_rotation(chain1):
    vecindex = 0
    assert chain1.vec.get(vecindex) == 0
    for i in [3, 2, 1, 0]:
        chain1.do_rigid_rot(vecindex=vecindex, direction=-1)
        chain1.update_chain()
        assert chain1.vec.get(vecindex) == i
    for i in [1, 2, 3, 0]:
        chain1.do_rigid_rot(vecindex=vecindex, direction=1)
        chain1.update_chain()
        assert chain1.vec.get(vecindex) == i

def test_three_bead_flip(chain2):
    vecindex = 1
    assert chain2.vec.get(vecindex) == 2
    assert chain2.vec.get(vecindex+1) == 3
    chain2.do_three_bead_flip(vecindex)
    chain2.update_chain()
    assert chain2.vec.get(vecindex) == 3
    assert chain2.vec.get(vecindex+1) == 2

def test_crankshaft(chain2):
    vecindex = 0
    assert chain2.vec.get(vecindex) == 1
    assert chain2.vec.get(vecindex+2) == 3
    chain2.do_crankshaft(vecindex)
    chain2.update_chain()
    assert chain2.vec.get(vecindex) == 3
    assert chain2.vec.get(vecindex+2) == 1

def test_shift_removes_trailing3_and_increments_last_vec(chain2):
    chain2.shift()
    assert len(chain2.vec) == 2
    assert chain2.vec.get(0) == 1
    assert chain2.vec.get(1) == 3

def test_leading0_first_turn1_is_nonsymmetric(chain1, chain2, chain3):
    assert not chain1.nonsym()
    assert not chain2.nonsym()
    assert chain3.nonsym()
