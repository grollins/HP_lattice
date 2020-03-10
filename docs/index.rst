===========
 HPlattice
===========

HPlattice is a Python library for the HP lattice model of Dill and Chan. It is ideally used as a teaching tool, or as a way to quickly prototype 2D lattice simulation ideas with easy-to-use extensible code. HPlattice can either 1) enumerate conformations, or 2) perform replica exchange monte carlo "dynamics" for 2-dimensional, square-lattice "bead-on-a-string" type chains.

The easiest way to get started with HPlattice is to download and install the Anaconda Python distribution. After installing Anaconda or another Python distribution of your choice, here are the steps for installing HPlattice:

1. `Download Bento <https://github.com/cournape/Bento/archive/master.zip>`_
2. Install Bento

.. code-block:: bash

    unzip Bento-master.zip
    cd Bento-master
    python bootstrap.py
    ./bentomaker configure
    ./bentomaker build
    ./bentomaker install

3. `Download HPlattice <https://www.dropbox.com/s/dzxqoe7ye1rb2kk/hplattice-1.0.tar.gz?dl=1>` (or from `GitHub <https://github.com/grollins/HP_lattice>`)
4. Install HPlattice

.. code-block:: bash

    tar xzf hp-lattice-1.0.tar.gz
    cd hp-lattice-1.0
    tar xzf HP-sequences.tgz # this step may take 30min to complete
    cd hplattice/util
    python setup.py build_ext --inplace
    cd ../..
    bentomaker install

5. Run unit tests (optional)

.. code-block:: bash

    conda install pytest # (or pip install pytest if Anaconda not installed)
    py.test hplattice

6. Try examples

.. code-block:: bash

    cd examples/enumerate
    python enumerate.py
    cd ../mcrex
    python mcrex.py

Contents
========

.. toctree::
    :maxdepth: 2

    userguide/index
    reference/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

