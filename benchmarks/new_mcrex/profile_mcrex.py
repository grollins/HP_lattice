import cProfile
import pstats
from os.path import join

from hplattice.Config import Config
from hplattice.MCSampler import MCSampler


CLIST_PATH = '/Users/grollins/work_code/HP_lattice/sequences/clist'
HP_STRING_SET = (10, 'HHPPHPPHPH', join(CLIST_PATH, 'hp10'))


def load_configuration(hp_string, initial_vec, clist_dir):
    # load in config file
    configfile = 'mcrex.conf'
    config = Config( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    config.NATIVEDIR = clist_dir
    return config

def main():
    N, hp_string, clist_dir = HP_STRING_SET
    print N, hp_string, clist_dir
    initial_vec = [0] * (N - 1)
    config = load_configuration(hp_string, initial_vec, clist_dir)
    s = MCSampler(config)
    s.do_mc_sampling()

if __name__ == '__main__':
    filename = './profile_stats.stats'
    cProfile.run('main()', filename)
    stats = pstats.Stats(filename)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats(20)
