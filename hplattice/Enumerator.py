
class Enumerator(object):
    """
    *Enumerator* objects are used to enumerate all conformations of
    an HP chain. The HP chain is defined in a configuration file, specified
    by the *config* parameter.

    :param lattice_factory: factory object that knows how to create chains and
                            trajectories
    :type lattice_factory: :class:`hplattice.LatticeFactory`
    :param str config: path to configuration file
    """
    def __init__(self, lattice_factory, config):
        self.lattice_factory = lattice_factory
        self.config = config
        self.chain = lattice_factory.make_chain(self.config.HPSTRING,
                                                self.config.INITIALVEC)

    def enumerate_states(self, save_trajectory=False,
                         trajectory_filename='traj.xyz'):
        """
        Enumerate all conformations of an HP chain.
        Prints density of contact states to stdout.

        :param bool save_trajectory: Generate an xyz coordinate trajectory
                                     when ``True``.
        :param str trajectory_filename: optional, save trajectory to this path
        """
        nconfs = 0
        # dictionary of {repr{contact state}: number of conformations}
        contact_states = {}
        # dictionary of {number of contacts: number of conformations}
        contacts = {}

        traj = self.lattice_factory.make_trajectory(save_trajectory,
                                                    trajectory_filename)

        #################
        #
        # This is a useful subroutine for enumerating all conformations
        # of an HP chain
        #
        # NOTE: in order for this to work correctly, the initial starting
        # vector must be [0,0,0,....,0]
        # 
        done = False
        while not done:
            # print self.chain
            if len(self.chain.vec) == (self.chain.n - 1):
                if self.chain.is_viable():
                    if self.chain.nonsym():
                        # tally the number of contacts
                        # state = self.chain.contactstate()
                        E, state = self.chain.energy()
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
                        # save configuration
                        traj.snapshot(self.chain)

                    done = self.chain.shift()

                else:
                    done = self.chain.shift()

            else:
                if self.chain.is_viable():
                    self.chain.grow()
                else:
                    done = self.chain.shift()

            if self.chain.is_first_vec_one():
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

        traj.finalize()
