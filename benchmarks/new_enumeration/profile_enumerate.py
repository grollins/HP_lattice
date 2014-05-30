import cProfile
import pstats
from os.path import join

from hplattice.Config import Config
from hplattice.Enumerator import Enumerator


HP_STRING_SET = (10, 'HPPHPPHPPH')


def load_configuration(hp_string, initial_vec):
    configfile = 'enumerate.conf'
    config = Config( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    return config

def main():
    N, hp_string = HP_STRING_SET
    initial_vec = [0] * (N - 1)
    config = load_configuration(hp_string, initial_vec)
    en = Enumerator(config, verbose=False)
    en.enumerate_states()

if __name__ == '__main__':
    filename = './profile_stats.stats'
    cProfile.run('main()', filename)
    stats = pstats.Stats(filename)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats(20)
