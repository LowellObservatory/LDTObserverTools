# PyDeVeny
Port of existing IDL programs for the DeVeny spectrograph on the LDT (dct-obs1) to python.  

*[Existing IDL code in this repository was exported from dct-obs1 on 31 August 2020.]*

The two major programs contained here are:

   - **`deveny_grangle`**: Compute the desired grating angle based on selected grating and desired central wavelength.  I would like this to be constructed as a GUI (probably with Tkinter rather than PyQt), with a dropdown containing the available gratings.

   - **`dfocus`**: Compute the needed collimator focus based on a series of arc line frames taken at various collimator settings.  Read in the arc lamp frames in the current night's focus directory, find the appropriate spectral lines in each frame, compute the FHWM (or other measure) of those lines, plot the FHWM as a function of collimator position and suggest the optimal focus position.  The IDL version runs in the command line and displays 3 plot windows containing the various information.  Maybe start with a Tkinter GUI for selecting the correct focus directory, and a "Compute Focus" button... the other plots can be displayed in 1 or more Tkinter windows.
    
