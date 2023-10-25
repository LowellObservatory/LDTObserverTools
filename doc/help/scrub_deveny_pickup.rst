.. code-block:: console

    $ scrub_deveny_pickup -h
    usage: scrub_deveny_pickup [-h] [--proc_dir PROC_DIR] [--overwrite_raw] [-d]
                               [-n]
                               file [file ...]
    
    Clean RF pickup noise from DeVeny raw frames
    
    positional arguments:
      file                 File(s) to clean
    
    options:
      -h, --help           show this help message and exit
      --proc_dir PROC_DIR  Path to the directory containing the .pypeit file used to
                           process `file` (ldt_deveny_?) (default: None)
      --overwrite_raw      Overwrite the raw file rather than create a new file with
                           the '_scrub' suffix (default: False)
      -d, --diagnostics    Output additional information and plots during the
                           analysis for debugging purposes (default: False)
      -n, --no_refit       Force no refit of 'bad' RMS values (default: False)
    