/***************************************************************************************
        What: Gaussian Subroutine
        Author: Trevor Vincent (vincenttrevor@gmail.com)
        Description:

        Creates data which can then be used to plot a gaussian
 ****************************************************************************************/

void gaussian(FILE *fh, double displacement[], int num_of_particles,
              int increments, double D_extra, double D_intra, double timestep,
              double number_of_steps, double mss_factor) {

    bool D_extra_is_bigger = false;

    if (D_extra > D_intra) {
        D_extra_is_bigger = true;
    }

    double max_step_size; // absolute max DISPLACEMENT the particle can travel
                          // after the diffusion time. Note: since this is a
                          // random walk, it is more than likely that this will
                          // never be achieved.

    if (D_extra_is_bigger == true) {
        max_step_size =
            timestep * number_of_steps * sqrt(6 * D_extra / timestep);
    } else {
        max_step_size =
            timestep * number_of_steps * sqrt(6 * D_intra / timestep);
    }

    int(*num_part_wdisp) =
        new int[increments]; // the number of particles within a certain
                             // increment of displacement, i.e. the number of
                             // particles with displacement between 0 and 1 mm

    max_step_size = max_step_size * mss_factor;

    for (int ii = 0; ii < increments; ii++) {
        num_part_wdisp[ii] = 0;
    }

    for (int i = 0; i < num_of_particles; i++) {

        for (int j = 0; j < increments; j++) {

            // if the particles displacement is between a certain increment, the
            // 2 is there because the total width (at most) of the gaussian
            // curve is TWICE the max_step_size
            if (displacement[i] <
                    -max_step_size + (j + 1) * 2 * max_step_size / increments &&
                displacement[i] >=
                    -max_step_size + (j)*2 * max_step_size / increments) {
                num_part_wdisp[j] = num_part_wdisp[j] + 1;
            }
        }
    }

    /* print to file */
    for (int jj = 0; jj < increments; jj++) {

        if (num_part_wdisp[jj] != 0) {
            fprintf(fh, "%.14f %d\n",
                    -max_step_size + (jj)*2 * max_step_size / increments,
                    num_part_wdisp[jj]);
        }
    }

    cout << endl;
    cout << "Gaussian Plot Info " << endl;
    cout << "Max Step Size = " << max_step_size << endl;
    cout << "Increments = " << increments << endl;
    cout << "Increment Size " << max_step_size / increments << endl;
    cout << "Max Step Size factor " << mss_factor << endl;

    fclose(fh);
}
