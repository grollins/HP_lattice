import sys
sys.path.append('../')

from timer import Timer

from hplattice.Config import Config
from hplattice.Enumerator import Enumerator


HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH']),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH']),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH']),
     (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH']),
     (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH']),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)


def load_configuration(hp_string, initial_vec):
    configfile = 'enumerate.conf'
    config = Config( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    return config

def main():
    timing_data = []
    for N, hp_string_list in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string
            initial_vec = [0] * (N - 1)
            with Timer() as t:
                config = load_configuration(this_hp_string, initial_vec)
                en = Enumerator(config)
                en.enumerate_states()
            time_elapsed = t.interval
            this_data = (N, time_elapsed, this_hp_string)
            timing_data.append(this_data)

    for td in timing_data:
        print "%d  %.2e  %s" % td

if __name__ == '__main__':
    main()