Quickstart guid for cms-l1t-analysis
============================================
Software package to analyse L1TNtuples

Quick-Start Guide for lxplus
--------------------------------------------
 1. Read the README etc if you want to dive deeper!
 2. To get cracking, follow the instructions below:

On Scientific Linux 6 with CVMFS available, e.g. lxplus, Soolin etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This includes nodes lxplus.cern.ch & private clusters

For now let's use this repository, you can fork your own of course!

.. code-block:: bash

  git clone -b weekly-checks https://github.com/professor-calculus/cms-l1t-analysis.git
  cd cms-l1t-analysis
  source setup.sh
  # you might need your grid cert -- but you shouldn't need it on lxplus with NTuples on EOS
  voms-proxy-init --voms cms
  make setup

make setup should work, but in the case that it throws Error 54 ignore it for now. This just means it cannot find
the test NTuples, but everything else is compiled ok and ready to rock 'n' roll.

running tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Tests can be run either on an SL 6 machine or in the Vagrant box:

.. code-block:: bash

  make test
  # if a grid proxy is provided (e.g. via voms-proxy-init --voms cms)
  # you can also run tests that require grid access:
  make test-all


This should fail if make setup threw Error 54.


Implementing and running an analysis script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"Analyzers" are the parts of the code that receive events from the input tuples, extracts the relevant data and puts this into the histograms.

You can see an example of an analyzer at: `cmsl1t/analyzers/demo_analyzer.py`.

To implement your own analyzer, all you need to do is make a new class in a file under `cmsl1t/analyzers/` which inherits from `cmsl1t.analyzers.BaseAnalyzer.BaseAnalyzer`.
You then need to implement a few methods:
 - `prepare_for_event`
 - `fill_histograms`
 - `write_histograms`
 - `make_plots`.

 See the BaseAnalyzer class and the demo_analyzer for examples and documentation of these methods.

Once you have implemented an analyzer and written a simple configuration for it, you can run it with `cmsl1t` command:

.. code-block:: bash

  cmsl1t -c config/demo.yaml -n 1000


Get help on the command line options by doing:

.. code-block:: bash

  cmsl1t --help


Running the weekly checks:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the weekly checks of jets & sums efficiencies, 2D plots etc, we have an analyzer:

Currently this is `cmsl1t/analyzers/jetMet_analyzer.py` but in the future it will likely get rolled into weekly_analyzer.py.

To run this, we want the config file `configs/offline_met_studies.yaml`.

Lines in the yaml give the paths for the NTuples, starting with 'root://' and ending with '/L1NTuple_*.root' -- You should edit these to point to the NTuples you're running on.

Further down the yaml file you can specify relative/absolute paths for the output.

Once you're ready, run as so:

.. code-block:: bash

  cmsl1t -c config/offline_met_studies.yaml -n <no_of_entries>


Or, for large running you can run on lxbatch with the NTuple list broken up into many jobs. Files per job specifiable with -f <number>.
Try to keep max ~1000 jobs else combining later is a nightmare.

Also more than about 20 files per job will lead to some going over the walltime...
Edit bin/cmsl1t_dirty_batch to change queue from 1nh to something else if this happens a lot.

For small runs I reccommend changing queue to 8nm and running ~4-6 files per job.

Once this submits the jobs it will tell you how to combine them together later :-)

.. code-block:: bash

  cmsl1t_dirty_batch -c config/offline_met_studies.yaml -f <ntuple_root_files_per_job>


HW vs Emu at Constant Rate:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code now has facility to calculate for a given HW L1 threshold the threshold at which the Emulated quantity gives the same rate.

Then, we may plot the turnons etc for HW and Emu at their respective rates together, using the new analysers.

To run this, first we want the config file `configs/HW_Emu_jet_rates.yaml`.

Lines in the yaml give the paths for the NTuples, starting with 'root://' and ending with '/L1NTuple_*.root' -- You should edit these to point to the ZeroBias NTuples you're running on.

Further down the yaml file you can specify relative/absolute paths for the output.

You'll notice in the analyser section we have something like
.. code-block:: bash

  thresholds:
    HTT: [val1, val2, ...]
    etc...

Change these to the thresholds you want for HW quantities.

Once you're ready, run as so:

.. code-block:: bash

  cmsl1t -c config/HW_Emu_jet_rates.yaml -n <no_of_entries>


Or, for large running you can run on lxbatch with the NTuple list broken up into many jobs using the command cmsl1t_dirty_batch. Files per job specifiable with -f <number>.
Try to keep max ~1000 jobs else combining later is a nightmare, but also if <number> is greater than baout 4 you might struggle for walltime per job.

For example, you could do:

.. code-block:: bash

  cmsl1t_dirty_batch -c config/HW_Emu_jet_rates.yaml -f 4

Then one can combine the output with the command which will be given to you at this stage.

Finally, the output of the code contains something like:

.. code-block:: bash

  thresholds:
    HTT: [val1, val2, ...]
    HTT_Emu: [val1_, val2_, ...]
    etc...

Copy this, being careful to keep the formatting...

Now we take a look at config/HW_Emu_constant_rate_turnons.yaml:

There are similar lines in here. Set the input ntuples to the SingleMu you wish to run on.

Now again in the analyser section we have the thresholds listed in the same format, but now with emulated quantities too.

Replace this with the stuff you just copied -- this format is interpreted as a python dictionary, so whitespace matters. thresholds: should be indented 2 spaces wrt lines above,
and HTT etc indented 2 spaces wrt thresholds.

Finally, we may run this like we did the previous step:

.. code-block:: bash

  cmsl1t_dirty_batch -c config/HW_Emu_constant_rate_turnons.yaml -f 4

Or of course

.. code-block:: bash

  cmsl1t -c config/HW_Emu_jet_rates.yaml -n <number_of_events>

if one wants a quick test job.
