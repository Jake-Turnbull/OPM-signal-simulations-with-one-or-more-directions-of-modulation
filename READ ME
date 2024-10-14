Instructions for use:
 
For the comparison code, the entire folder should be downloaded, and the code opened in spyder (python 3.7.9) to guarantee the correct outputs. 
Ensure they are all saved in the same folder and that spyder is using that folder as its path so that it can initialise the arrays from the saved data.
The arrays are labelled in the following structure: DATATYPE and direction_DIRECTION OF OSCILLATION. For example, pxxval_xyz are the polarisation values in the x-direction when there is modulation in x, y and z. Another example being valzc_xy, these are the cosine modulated components of the polarisation in z when there is modulation in x and y. 
The two scripts in this folder, ‘BY_BZ_contours’ and ‘Comparisson_code_18000’ serve slightly different purposes. The latter shows the polarisations in x, y and z as well as the frequency components taking Bz as our slice when plotting the contours. The former only shows the contours of the frequency components when is our axis to be sliced.
 
This folder contains 3 pieces of code, the first ‘zero_crossing_zhang_comparisson’ simply has the recorded values of the zero crossings when both By and Bz are non-zero with the same magnitude over the range shown in the report compared to the predicted analytical form.
The second, ‘Zhang_comparisson_both_x_and_xy_modulation’, provides two figures and real time simulation using the same method as the code used for our analysis. The figures respectively show the shape of the sine modulated frequency response as its predicted analytical form and our simulated form for both x modulation and y modulation. It also shows the amplitude difference between the current value set by the slider controlling magnetic field with that of the case where Bz=By=0.
The third piece of code ‘Zhang_single_mod_comparisson’ shows the same as the prior but only for single axis modulation along the x-axis and was used to verify out method.
 
The third folder contains out optimised code: the scripts beginning with the ‘varying’ take along time to run but were the scripts used to calculate our values for the analysis. Rather than running these scripts it is probably better to use the already obtained data in the first folder. 
However, if you want to get your own simulated data at different field ranges, the names of the code represent the amount of modulation being applied.
The ‘contour_plot_test’ code does the same as the simulation code, but for far less values and therefore is easier to change and observe y.
Finally, the script titles ‘Optimised_with_freq_plus_linearity_calc’ calculates the linear equations of the frequency response for a set range of values of Bz and Bx, when By is zero. This can be easily m
