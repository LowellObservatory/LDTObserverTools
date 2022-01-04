# LDTObserverTools
## Collection of Observer Tools for the Lowell Discovery Telescope

Major programs contained here:

   - **`deveny_grangle`**: Compute the desired grating angle based on selected grating and desired central wavelength.  This routine comes with two interfaces.  The default GUI features a dropdown menu for grating selection and contains error checking on the input for central wavelength.  There is a ``MAX_GUI`` option for computing central wavelength given the grating angle in addition to the standard GUI features.  Also included is a command-line interface, identical to the old IDL function.  Online help is available with the ``-h`` option.  [Completed: 2021-01-26]

   - **`dfocus`**: Compute the needed collimator focus based on a series of arc line frames taken at various collimator settings.  Read in the arc lamp frames in the current night's focus directory, find the appropriate spectral lines in each frame, compute the FHWM (or other measure) of those lines, plot the FHWM as a function of collimator position and suggest the optimal focus position.  This program is executed identically to the old IDL version.  The python version uses `scipy.signal` processing routines for identifying line peaks and widths, resulting in more accurate and consistent estimates of the correct collimator focus position.  Rather than separately producing plots to the screen and disk, this version writes all plots to a PDF, then launches `Preview.app` to display the plots on the screen.  Newly added is a readout of the mount temperature so the user can determine when/if the collimator needs to be refocused during the night.  Online help is available with the ``-h`` option.  [Completed: 2021-11-19]

Future programs:

 - **`lmi_etc`**: Python GUI version of the LMI exposure time calculator (http://www2.lowell.edu/users/massey/LMI/etc_calc.php).

 - **`deveny_collfoc_range`**: Use the specified grating angle and mount temperature to suggest a range for use with the DeVeny LOUI collimator focus sequence function.  This is important because, unlike all other focus routines at LDT, this function takes the _**starting point**, step, and number of exoposures_ rather than the _**expected focus value**, step, and number of exposures_.  Having a routine to compute the desired values would make this step easier and less prone to error (_i.e._, only searching on one side of the expected focus).
