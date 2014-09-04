#! /usr/bin/python

from os.path import join

from hplattice import LatticeFactory
from hplattice.MCSampler import MCSampler


CLIST_PATH = '../../HP-sequences/sequences/clist'
HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH'], join(CLIST_PATH, 'hp04')),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH'], join(CLIST_PATH, 'hp06')),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH'], join(CLIST_PATH, 'hp10')),
     # (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH'], join(CLIST_PATH, 'hp12')),
     # (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH'], join(CLIST_PATH, 'hp14')),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)


def load_configuration(lattice_factory, hp_string, initial_vec, clist_dir):
    # load in config file
    configfile = 'mcrex.conf'
    config = lattice_factory.make_configuration( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    config.NATIVEDIR = clist_dir
    return config

def main():
    lattice_factory = LatticeFactory()
    for N, hp_string_list, clist_dir in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string, clist_dir
            initial_vec = [0] * (N - 1)
            config = load_configuration(lattice_factory, this_hp_string,
                                        initial_vec, clist_dir)
            s = MCSampler(lattice_factory, config)
            s.do_mc_sampling(save_trajectory=True,
                             trajectory_filename='%02d_traj.xyz' % N)

if __name__ == '__main__':
    main()
