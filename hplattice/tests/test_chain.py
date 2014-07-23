import pytest

from ..Chain import Chain


@pytest.fixture
def chain():
    hpstring = 'PHPPHPHPPHH'
    initial_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    return Chain(hpstring, initial_vec)

@pytest.fixture
def chain2():
    hpstring = 'HPPH'
    initial_vec = [1, 2, 3]
    return Chain(hpstring, initial_vec)

def test_chain_is_correct_length(chain):
    assert len(chain) == 11
    assert chain.get_vec_length() == 10

def test_convert_hpstring_to_binary(chain):
    binstr = chain.hpstr2bin()
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
