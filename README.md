# Feature map statistics
A collection of MATLAB functions for generating and analysing cortical feature
maps. All functions were developed and tested with MATLAB R2014b. 
 
This code was used to perform the analyses in the following paper:
Cloherty, SL, Hughes, NJ, Hietanen, MA, Bhagavatula, PS, Goodhill, GJ,
Ibbotson, MR (2016) "Sensory experience modifies feature map relationships in
visual cortex."

`tests.m` contains a collection of tests and examples of the provided functions.

## Functions

### Map pre-processing and generation
- `esd`: Performs extended spatial decorrelation. Follows the method described
         in section 4.5.1 of http://webdoc.sub.gwdg.de/ebook/diss/2003/tu-berlin/diss/2001/schiessl_ingo.pdf
- `align_images`: Aligns a collection of images to a reference image using
                 spatial translation by maximising linear correlation.

### Map statistics and analysis
- `orientation_hist`: Generates a histogram of orientation preference from an
                      OP map.
- `crossing_angle_dist`: Calculates a distribution of the crossing-angles
                         between the contours of an OP and an OD map.
- `od_op_crossing`: Calculates the angle of intersection between all crossings
                    of the contours of an OP and an OD map. Used by
                    `crossing_angle_dist`.
- `pinwod`: Calculates a histogram of pinwheel locations in an OP map relative
            to the accompanying OD map.
- `binned_selectivity`: Calculates a histogram of orientation selectivity
                        in the same OD bins as used by `pinwod`.
- `pinw_density`: Calculates the pinwheel density of an OP map.
- `locate_pinwheels`: Provides the location of all pinwheels in an OP map. Used
                      by `pinw_density`.
- `fourier_wavelength`: Calculates the average wavelength of an OP map using
                        its Fourier spectrum. Used by `pinw_density`.

### Plotting 
- `op_contours` - Calculates the contours of an OP map for pretty plotting.

### Other
- `corr_matrix`: Calculates a shifted cross-correlation matrix. Used by `esd`.
- `lattice_gradient_descent`: Performs gradient descent minimisation on a
                              lattice. Used by `align_images`.
