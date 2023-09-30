.. code-block:: console

    $ scrub_deveny_pickup -h
    usage: scrub_deveny_pickup [-h] [--debug_plots] [--overwrite_raw]
                               file [file ...]
    
    Clean RF pickup noise from DeVeny raw frames
    
    positional arguments:
      file             File(s) to clean
    
    options:
      -h, --help       show this help message and exit
      --debug_plots    Create plots during the analysis for debugging purposes
                       (default: True)
      --overwrite_raw  Overwrite the raw file rather than create a new file with the
                       '_scrub' suffix (default: False)
    