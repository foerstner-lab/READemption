Troubleshooting
===============

Pickling Error
--------------

Problem
~~~~~~~

When running reademption (e.g. with Python 3.2) you get an error like this:
::
   Traceback (most recent call last):
     File "/usr/lib/python3.2/multiprocessing/queues.py", line 272, in _feed
      send(obj)
   _pickle.PicklingError: Can't pickle <class 'method'>: attribute lookup builtins.method failed

Solution
~~~~~~~~

Switch to Python 3.3 or higher. The newer versions can handle this.
