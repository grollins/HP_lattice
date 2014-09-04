from hplattice import LatticeFactory
from hplattice.Enumerator import Enumerator


HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH']),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH']),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH']),
     # (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH']),
     # (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH']),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)


def load_configuration(lattice_factory, hp_string, initial_vec):
    configfile = 'enumerate.conf'
    config = lattice_factory.make_configuration( filename=configfile )
    config.HPSTRING = hp_string.strip()
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    return config

def main():
    lattice_factory = LatticeFactory()
    for N, hp_string_list in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string
            initial_vec = [0] * (N - 1)
            config = load_configuration(lattice_factory, this_hp_string,
                                        initial_vec)
            en = Enumerator(config)
            en.enumerate_states(save_trajectory=True,
                                trajectory_filename='%d_%s_traj.xyz' % \
                                (N, this_hp_string))

if __name__ == '__main__':
    main()
