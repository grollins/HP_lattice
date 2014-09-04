.. _article2:

================
 HP Simulations
================

Replica Exchange Monte Carlo
============================

The basic steps involved in running a monte carlo simulation are:

1. Create a *LatticeFactory*.
2. Use the *LatticeFactory* to load a *Configuration*.
3. Create an *MCSampler*, passing the *LatticeFactory* and *Configuration* as
   arguments.
4. Call the ``do_mc_sampling`` method of the *MCSampler*.

.. code-block:: python

    from hplattice import LatticeFactory
    from hplattice.MCSampler import MCSampler 

    lattice_factory = LatticeFactory()
    config = lattice_factory.make_configuration(filename='mcrex.conf')

    mc = MCSampler(lattice_factory, config)
    mc.do_mc_sampling(save_trajectory=True, trajectory_filename='traj.xyz')

Conformational State Enumeration
================================

The basic steps involved in running an enumeration simulation are:

1. Create a *LatticeFactory*.
2. Use the *LatticeFactory* to load a *Configuration*.
3. Create an *Enumerator*, passing the *LatticeFactory* and *Configuration* as
   arguments.
4. Call the ``enumerate_states`` method of the *Enumerator*.

.. code-block:: python

    from hplattice import LatticeFactory
    from hplattice.Enumerator import Enumerator

    lattice_factory = LatticeFactory()
    config = lattice_factory.make_configuration(filename='enumerate.conf')

    en = Enumerator(lattice_factory, config)
    en.enumerate_states(save_trajectory=True, trajectory_filename='traj.xyz')
