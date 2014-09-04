.. _article1:

===========================
 Configuration File Format
===========================

*hplattice* simulations require a configuration file. The file should have
formatted rows consisting of two fields, separated by white-space (or any non-printing characters, like tabs)::

    HPSTRING                PHPPHPHPPHH
    INITIALVEC              [1, 0, 1, 2, 1, 2, 1, 2, 3, 3]
    EPS                     -5.0
    RESTRAINED_STATE        [(1, 4), (6, 9)]
    KSPRING                 0.0
    NREPLICAS               9
    REPLICATEMPS            [275.0, 300.0, 325.0, 350.0, 400.0, 450.0, 500.0, 600.0,    1000.0]
    MCSTEPS                 1000
    SWAPEVERY               50
    SWAPMETHOD              random pair
    MOVESET                 MS2
    PRINTEVERY              1
    NATIVEDIR               ../../HP-sequences/sequences/clist/hp11
    STOPATNATIVE            False

Here is the full list of the parameters and what each one represents. If a
parameter is not specified in the file, it will be set to a default value.

HPSTRING
    The *HPSTRING* specifies the chemical nature of each monomer. The only
    supported options are Hydrophobic (H) or Polar (P).

INITIALVEC
    The *INITIALVEC* is a list of integers that specifies the direction of the
    chain between two neighboring monomers: 0 (up), 1 (left), 2 (down) or 
    3 (right). For example, ``[0,0,0,0]`` would correspond to a 5-mer that points straight up from the origin (in the positive-y direction).

EPS
    The energy of a hydrophobic contact (two H's in adjacent lattice spaces that
    are not :math:`i+1` or :math:`i+2` neighbors along the chain).

NREPLICAS
    The number of replicas for replica exchange simulations. This is usually
    set to ``1`` for enumeration simulations because those simulations do not
    involve random moves in conformational space.

REPLICATEMPS
    A list of floats that specifies the temperature (K) of each replica. The 
    length of the *REPLICATEMPS* list should be equal to *NREPLICAS*.

MCSTEPS
    The number of monte carlo steps to run. Each replica will run for this
    number of steps, and they will periodically attempt to swap temperatures.

SWAPEVERY
    The number of steps between replica swap attempts.

SWAPMETHOD
    How to swap replicas. ``random pair`` to randomly choose two replicas
    to swap; ``neighbors`` to randomly choose one replica ``i`` and swap it
    with its ``i+1`` neighbor.

MOVESET
    Select which type of monte carlo moves will be used to sample conformational
    space: ``MS1`` for three-bead flips and rigid rotations; ``MS2`` for
    three-bead flips, crankshaft moves, and rigid rotations; and ``MS3`` for
    rigid rotations only.

RESTRAINED_STATE
    A list of tuples that specifies contacts that should be harmonically
    restrained. Each tuple in the list should contain a pair of integers that
    correspond to the indices of the monomers that should be restrained. An
    example would be ``[(1, 4), (6, 9)]`` which would add restraints to the
    monomer1-monomer4 contact and the monomer6-monomer9 contact. Note that the
    indices are 0-indexed, so monomer0 is the first monomer in the chain.

KSPRING
    The force constant of the harmonic restraints specified in *RESTRAINED_STATE*.

PRINTEVERY
    In a monte carlo simulation, save coordinates to trajectory after this
    number of steps.

NATIVEDIR
    The path to the file that specifies what the native contacts are for the
    chain specified by *HPSTRING*.

STOPATNATIVE
    If the monte carlo simulation finds the native conformation of the chain
    (as defined by the contacts in *NATIVEDIR*), then halt the simulation if
    *STOPATNATIVE* is ``True``.
