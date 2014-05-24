from hplattice.Chain import Chain

class Enumerator(object):
    """docstring for Enumerator"""
    def __init__(self, config, verbose):
        self.config = config
        self.verbose = verbose
        self.chain = Chain(self.config)

    def enumerate_states(self):
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
            if len(self.chain.vec) == self.chain.n - 1:    
                if self.chain.viable:                
                    if self.chain.nonsym():
                        # tally the number of contacts
                        state = self.chain.contactstate()
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

                    done = self.chain.shift()
                        
                else:
                    done = self.chain.shift()

            else:
                if self.chain.viable:
                    self.chain.grow()
                else:
                    done = self.chain.shift()

            if self.chain.vec[0] == 1:
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
                (num_contacts, self.config.eps * num_contacts, num_confs)
        print
        print 'at T = %4.1f K' % self.config.T
