from numpy import zeros, nonzero
from .Replica import attemptswap


class MCSampler(object):
    """
    *MCSampler* objects are used to run replica exchange monte carlo simulations
    of an HP chain. The HP chain is defined in a configuration file, specified
    by the *config* parameter. The configuration file also specifies the
    replica exchange parameters, such as the temperature of each replica.

    :param lattice_factory: factory object that knows how to create replicas and
                            trajectories.
    :type lattice_factory: :class:`hplattice.LatticeFactory`
    :param str config: path to configuration file
    """
    def __init__(self, lattice_factory, config):
        self.lattice_factory = lattice_factory
        self.config = config
        self.native_contacts = self._load_native_contacts()
        self.replicas = []
        for i in range(0, self.config.NREPLICAS):
            r = lattice_factory.make_replica(lattice_factory, self.config, i,
                                             self.native_contacts)
            self.replicas.append(r)

    def _load_native_contacts(self):
        if self.config.STOPATNATIVE == 1:
            nativeclistfile = self.config.NATIVEDIR + '/' + self.config.HPSTRING + '.clist'
            with open(nativeclistfile,'r') as fnative:
                nativeclist_str = fnative.readline()
                nativeclist = eval(nativeclist_str)
        else:
            nativeclist = None
        return nativeclist

    def _init_mc_stats(self):
        ### Initialize replica exchange stats 
        # the number of attemped replica swaps
        self.swaps = zeros(len(self.replicas))
        # the number of viable swaps
        self.viable_swaps = zeros(len(self.replicas))
        # the fraction of swaps that were accepted
        self.swap_acceptance = zeros(len(self.replicas))
        # initialize replica-level stats
        for r in self.replicas:
            r.init_mc_stats()
        # T-based stats
        self.accepted_steps_at_T = {}
        for T in self.config.REPLICATEMPS:
            self.accepted_steps_at_T[T] = 0

    def _update_swap_stats(self, i, j, swap_sucess):
        ### increment swap stats for replica i and replica j
        self.swaps[i] += 1
        self.swaps[j] += 1
        if swap_sucess:
            self.viable_swaps[i] += 1
            self.viable_swaps[j] += 1

    def _compute_swap_acceptance(self):
        ### compute the fraction of swaps that have been viable
        inds = nonzero(self.swaps)
        self.swap_acceptance[inds] = \
            (1.*self.viable_swaps[inds]) / self.swaps[inds]

    def _output_stats(self, prodstep):
        ### Output the status of the simulation
        print prodstep, 'production steps'
        print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s ' % \
              ('replica','viablesteps','steps','MCaccept','viableswaps',
               'swaps','SWAPaccept')
        for idx, rep in enumerate(self.replicas):
            print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s %d %d' % \
                (idx, rep.viablesteps, rep.steps,
                 '%1.3f' % rep.acceptance, self.swaps[idx],
                 self.viable_swaps[idx], '%1.3f' % self.swap_acceptance[idx],
                 rep.mc.tempfromrep, rep.mc.temp)
        if self.config.STOPATNATIVE == 1:
            print 'NATIVE CLIST:', self.native_contacts
        print '%-8s %-12s %-12s' % \
            ('replica', 'foundnative', 'contact state') 
        for rep in self.replicas:
            print '%-8d %-12d %s %s' % \
                (rep.repnum, rep.is_native(), rep.contactstate(), rep.get_vec())
        print self.accepted_steps_at_T

    def do_mc_sampling(self, save_trajectory=False, trajectory_filename='traj.xyz'):
        """
        Run replica exchange monte carlo of the HP chain.

        :param bool save_trajectory: Generate xyz coordinate trajectories
                                     when ``True``. There will be separate
                                     trajectory for each replica.
        :param str trajectory_filename: optional, save trajectory to this path.
                                        Replica numbers will be prepended to the
                                        name specified here.
        """
        traj_dict = {}
        for i, r in enumerate(self.replicas):
            traj = self.lattice_factory.make_trajectory(save_trajectory,
                    '%03d_%s' % (i, trajectory_filename))
            traj_dict[r] = traj

        self._init_mc_stats()
        found_native = False
        
        for prodstep in xrange(self.config.MCSTEPS):
            # Run the replicas for a production cycle...
            for r in self.replicas:
                move_is_viable = r.propose_move()
                if move_is_viable:
                    move_is_accepted = r.metropolis_accept_move()
                    if move_is_accepted:
                        self.accepted_steps_at_T[r.get_T()] += 1
                else:
                    move_is_accepted = False
                r.record_stats(move_is_viable, move_is_accepted)
                if self.config.STOPATNATIVE == 1 and r.is_native():
                    found_native = True

            if self.config.STOPATNATIVE == 1 and found_native:
                break

            # After the production cycle,      
            if (prodstep % self.config.SWAPEVERY) == 0:
                ### ...after every production run, attempt a SWAP
                swap_results = \
                    attemptswap(self.config.SWAPMETHOD, self.replicas)
                self._update_swap_stats(*swap_results)

            # Print status
            if (prodstep % self.config.PRINTEVERY) == 0:
                # calc MC move acceptance
                for r in self.replicas:
                    r.compute_mc_acceptance()       
                # calc replica swap acceptance
                self._compute_swap_acceptance()
                # self._output_stats(prodstep)
                for rep, traj in traj_dict.iteritems():
                    traj.snapshot(rep.chain)

        self._output_stats(prodstep)
        for traj in traj_dict.itervalues():
            traj.finalize()
