.. code-block:: console

    $ fix_ldt_header -h
    usage: fix_ldt_header [-h] [--rimas_exp] file [file ...] keyword new_value
    
    Fix a keyword in LDT FITS headers
    
    positional arguments:
      file         File(s) on which to operate
      keyword      FITS keyword to change
      new_value    New header keyword value to insert
    
    options:
      -h, --help   show this help message and exit
      --rimas_exp  Fix older RIMAS exposure time keywords (default: False)
    
    LDT Observer Tools (obstools) version 0.7.5
