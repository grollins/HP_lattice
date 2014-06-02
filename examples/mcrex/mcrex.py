#! /usr/bin/python

from os.path import join

from hplattice.Config import Config
from hplattice.MCSampler import MCSampler


CLIST_PATH = '../../sequences/clist'
HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH'], join(CLIST_PATH, 'hp04')),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH'], join(CLIST_PATH, 'hp06')),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH'], join(CLIST_PATH, 'hp10')),
     # (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH'], join(CLIST_PATH, 'hp12')),
     # (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH'], join(CLIST_PATH, 'hp14')),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)
VERBOSE = False


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
    for N, hp_string_list, clist_dir in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string, clist_dir
            initial_vec = [0] * (N - 1)
            config = load_configuration(this_hp_string, initial_vec,
                                        clist_dir)
            s = MCSampler(config, VERBOSE)
            s.do_mc_sampling()

if __name__ == '__main__':
    main()
