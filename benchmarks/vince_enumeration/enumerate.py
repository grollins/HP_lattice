import sys
sys.path.append('../')

from timer import Timer

from hplattice.Config import Config
from hplattice.Replica import Replica

HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH']),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH']),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH']),
     (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH']),
     (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH']),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)


def enumerate_states(hp_string, initial_vec):
    configfile = 'enumerate.conf'
    config = Config( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []

    # create a single Replica
    replicas = [ Replica(config, 0) ]
    
    nconfs = 0
    # dictionary of {repr{contact state}: number of conformations}
    contact_states = {}
    # dictionary of {number of contacts: number of conformations}
    contacts = {}

    #################
    #
    # This is a useful subroutine for enumerating all conformations
    # of an HP chain
    #
    # NOTE: in order for this to work correctly, the initial starting
    # vector must be [0,0,0,....,0]
    # 
    done = 0
    while not(done):
        if len(replicas[0].chain.vec) == replicas[0].chain.n - 1:    
            if replicas[0].chain.viable:                
                if replicas[0].chain.nonsym():
                    # tally the number of contacts
                    state = replicas[0].chain.contactstate()
                    ncontacts = len(state)
                    if contacts.has_key(ncontacts) == False:
                        contacts[ncontacts] = 1
                    else:
                        contacts[ncontacts] = contacts[ncontacts] + 1

                    # tally the contact state
                    this_state_repr = repr(state)
                    if contact_states.has_key(this_state_repr) == False:
                        contact_states[this_state_repr] = 1
                    else:
                        contact_states[this_state_repr] = \
                            contact_states[this_state_repr] + 1

                    # tally the number of conformations
                    nconfs = nconfs + 1

                done = replicas[0].chain.shift()
                    
            else:
                done = replicas[0].chain.shift()

        else:
            if replicas[0].chain.viable:
                replicas[0].chain.grow()
            else:
                done = replicas[0].chain.shift()

        if replicas[0].chain.vec[0] == 1:
            # skip the other symmetries
            break        
    #
    #
    #################    

    # print out the density of contact states
    print
    print 'DENSITY of CONTACT STATES:'
    print '%-40s %s' % ('contact state','number of conformations')
    for state, num_confs in contact_states.items():
        print '%-40s %d' % (state, num_confs)

    # print out the density of states (energies)
    print 
    print 'DENSITY of STATES (in energies/contacts):'
    print '%-20s %-20s %s' % \
            ('number of contacts', 'energy (kT)',
             'number of conformations')
    for num_contacts, num_confs in contacts.items():
        print '%-20d %-20d %d' % \
            (num_contacts, config.eps * num_contacts, num_confs)
    print
    print 'at T = %4.1f K' % config.T

def main():
    timing_data = []
    for N, hp_string_list in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string
            initial_vec = [0] * (N - 1)
            with Timer() as t:
                enumerate_states(this_hp_string, initial_vec)
            time_elapsed = t.interval
            this_data = (N, time_elapsed, this_hp_string)
            timing_data.append(this_data)

    for td in timing_data:
        print "%d  %.2e  %s" % td

if __name__ == '__main__':
    main()