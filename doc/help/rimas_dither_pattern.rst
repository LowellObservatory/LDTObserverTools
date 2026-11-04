.. code-block:: console

    $ rimas_dither_pattern -h
    usage: rimas_dither_pattern [-h] -s FRAME [FRAME ...] --date DATE
                                [--arms {YJ,HK} [{YJ,HK} ...]] [-o OUTPUT]
                                [--source SOURCE] [--overwrite]
                                {ABBA,ONOFF}
    
    Create an ECSV sidecar correcting RIMAS dither metadata.
    
    positional arguments:
      {ABBA,ONOFF}
    
    options:
      -h, --help            show this help message and exit
      -s, --sequence FRAME [FRAME ...]
                            Full frame list for one commanded sequence, including
                            frames not written to disk. Repeat for each sequence;
                            Bash brace expansion such as -s {86..89} is supported.
      --date DATE           Observing date (YYYYMMDD)
      --arms {YJ,HK} [{YJ,HK} ...]
      -o, --output OUTPUT
      --source SOURCE
      --overwrite
    
    LDT Observer Tools (obstools) version 0.7.5
