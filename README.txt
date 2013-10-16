This is Muse, which stands for something I forget.  Anyway, it does
sound-source localization of ultrasonic mouse vocalizations.

generate_r_est_raw_for_single_mouse_data.m: Estimates the position of
  each vocalization in the single-mouse data set, and ancillary things
  like the MSE at the optimal position, and others.  Generates a file
  called r_est_raw_for_single_mouse_data.mat, in which all of this
  stored.  r_est_raw_for_single_mouse_data.mat contains no information
  about P-values or confidence regions.  That's what makes it "raw".

generate_cdf_dJ_emp_mat_file.m: Reads the file
  r_est_raw_for_single_mouse_data.mat, outputs a file called
  cdf_dJ_emp_unique.mat, which is needed for computation of P values
  and confidence regions.  It contructs an empirical CDF of the
  delta-J values at the mouse head for each vocalization and fits a
  scaled chi-squared CDF to this data.  This fit is saved to
  cdf_dJ_emp_unique.mat, and is used to map delta-J values to P
  values.  This enables both the calculation of confidence regions and
  P-values for mouse positions known from video.

cross_validate_crs_fit_cdf.m: Do cross-validation of the coverage of
  the confidence regions.  N-fold Cross-validation is done over the N
  trials.

generate_pdfs_for_single_mouse_data.m: Generates a directory full of
  PDFs with one page per vocalization from the single-mouse data.  Each
  page shows the estimate and the map, plus some other information.

collect_pdfs_for_single_mouse_data.m: Function to collect all the
  single-page PDFs for the single-mouse data into a single PDF, named
  single_mouse_data.pdf.

generate_pdfs_for_08212012_B_001_to_500.m: Generates a directory full of
  PDFs with one page per vocalization from the 08-12-2012-B data.
  Only does the first 500 vocalizations.

collect_pdfs_for_08212012_B_001_to_500.m: Function to collect all the
  single-page PDFs for the 08-12-2012-B data data into a single PDF, named
  08212012_B_001_to_500.pdf.

generate_pdfs_for_10072012_B.m: Generates a directory full of
  PDFs with one page per vocalization from the 10-07-2012-B data.

collect_pdfs_for_10072012_B.m: Function to collect all the
  single-page PDFs for the 10-07-2012-B data data into a single PDF, named
  10072012_B.pdf.

This version of Muse depends upon the Taylor Matlab Toolbox, release 1.14.

All code in Muse, except that in toolbox/snippeter, is copyright Adam
L. Taylor, 2013.  It is licensed under the BSD 2-clause license.  (See
below.)

Copyright (c) 2013, Adam L. Taylor
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of Adam L. Taylor or the Howard Hughes Medical Institute.

