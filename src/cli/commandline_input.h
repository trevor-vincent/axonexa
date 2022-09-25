//#include <string.h>
//#include "getopt.h"
//#include <iostream>
//#include <cmath>
//#include <ctype.h>

// using namespace std;

/*

We define here a stucture called VARIABLES which holds all the
parameters given at command line plus the ones that weren't

*/

typedef struct {
    int geom;
    int num;
    int increments; // number of times ADC will be computed during a simulation
    double de;      // extracellular diffusion coefficient
    double di;      // intracellular diffusion coefficient
    double ds; // diffusion coefficient of the sheath, atm this is NOT BEING
               // USED by any of the code
    double dc; // diffusion coefficient of the core, atm this is NOT BEING USED
               // by any of the code
    double rs; // radius of the sheath, atm this is NOT BEING USED by any of the
               // code
    double
        rc; // radius of the core, atm this is NOT BEING USED by any of the code
    double T2e; // T2 extracellular
    double T2i; // T2 intracellular
    double
        T2s; // T2 for the sheath, atm this is NOT BEING USED by any of the code
    double T2c; // T2 for the core, atm this is NOT BEING USED by any of the
                // code
    double gt;  // gradient duration
    double gs;  // gradient spacing
    double r;   // radius of cylinder
    double csx; // x dimension cell size for szafer box cells
    double csy; // y dimension cell size for szafer box cells
    double csz; // z dimension cell size for szafer box cells (NOT BEING USED BY
                // ANYTHING). But I would still recommend setting csz = ls, just
                // in case it is being used
    double p;   // permeability
    double dt;  // timestep
    double ls;  // lattice size
    double g;   // initial gradient strength
    double hex_sidelength;  // hexagon sidelength d (see lattice.h for more info
                            // about this)
    bool printparticleinfo; // prints particle info, NOT SUPPORTED BY THE CODE
                            // ANYMORE
    bool printparticleposition; // prints the trajectories of the particles to
                                // fiel
    bool printsignalinfo;
    bool printinitialvars;
    bool printtofile; // prints ADCx and ADCz to file
    bool gaussian;    // prints gaussian data to file
    bool use_T2;      // flag which determines if T2 is used or not
    int graddir;      // direction of the gradient (0 = x-axis, 1 = y-axis, 2=
                      // z-axis)
    int num_of_trial; // number of different b values used for linear regression
                      // calculation
    int num_of_repeat; // number of times to repeat the simulation (back to back
                       // simulations with the same initial parameters)
    FILE *fh;          // main file handle for printing ADCx and ADCz
    FILE *fh_g;        // file handle for gaussian
    FILE *fh_pp;       // defunct
    /* file handles for incorrect linear regression investiagetion */
    FILE *fh_incorrectlinreg_ADCx;
    FILE *fh_incorrectlinreg_ADCy;
    FILE *fh_incorrectlinreg_FA;
    FILE *fh_incorrectlinreg_ADCz;
    int num_of_cyl;    // num_of_cyl in random lattice main voxel
    double max_radius; // max radius
    double min_radius; // min radius
    int cyl_seed; // seed for random number generator used to generate cylinder
                  // positions and radii
    double rand_mean;
    double rand_stddev;

} VARIABLES;

VARIABLES commandline_input(int argc, char *argv[]) {

    VARIABLES vars;

    // default variables
    vars.de = 2.5 * pow(10.0, -6.0);
    vars.di = 1.0 * pow(10.0, -6.0);
    vars.ds = 1.5 * vars.di;
    vars.dc = vars.di;
    vars.rs = 0;
    vars.rc = 0;
    vars.T2e = 200;
    vars.T2i = 200;
    vars.T2s = vars.T2i * 1.5;
    vars.T2c = vars.T2i;
    vars.gt = .001;
    vars.gs = 5;
    vars.ls = .003;
    vars.csx = 0.0026833;
    vars.csy = 0.0026833;
    vars.csz = vars.ls;
    vars.p = .005;
    vars.dt = .1;
    vars.g = 0.0001609114;
    vars.hex_sidelength = .5;
    vars.printparticleinfo = false;
    vars.printparticleposition = false;
    vars.printsignalinfo = false;
    vars.printinitialvars = false;
    vars.printtofile = false;
    vars.gaussian = false;
    vars.graddir = 0;
    vars.increments = 1;
    vars.r = 0;
    vars.use_T2 = true;
    vars.num_of_trial = 2;
    vars.num_of_repeat = 1;
    vars.num_of_cyl = 1;
    vars.max_radius = vars.ls / 2;
    vars.min_radius = vars.ls / 20;
    vars.cyl_seed = 12345678;
    vars.rand_mean = .00535 / 2.0;
    vars.rand_stddev = .00252 / 2.0;

    const char *optv = "a:b:cdefghijklmnopqrstuvwxy";

    /* Here we connect the parameter (e.g. geom) with a shortcut character (e.g.
     * 'a') */
    static const struct option longOpts[] = {
        {"geom", required_argument, NULL, 'a'},
        {"num", required_argument, NULL, 'b'},
        {"de", required_argument, NULL, 'c'},
        {"di", required_argument, NULL, 'd'},
        {"ds", required_argument, NULL, 'e'},
        {"dc", required_argument, NULL, 'f'},
        {"rs", required_argument, NULL, 'g'},
        {"rc", required_argument, NULL, 'h'},
        {"T2e", required_argument, NULL, 'i'},
        {"T2i", required_argument, NULL, 'j'},
        {"T2s", required_argument, NULL, 'k'},
        {"T2c", required_argument, NULL, 'l'},
        {"gt", required_argument, NULL, 'm'},
        {"gs", required_argument, NULL, 'n'},
        {"csx", required_argument, NULL, 'p'},
        {"csy", required_argument, NULL, 'q'},
        {"csz", required_argument, NULL, 'r'},
        {"perm", required_argument, NULL, 's'},
        {"dt", required_argument, NULL, 't'},
        {"ls", required_argument, NULL, 'u'},
        {"grad", required_argument, NULL, 'v'},
        {"hex", required_argument, NULL, 'w'},
        {"printparticleinfo", no_argument, NULL, 0},
        {"printsignalinfo", no_argument, NULL, 2},
        {"printinitialvars", no_argument, NULL, 3},
        {"printtofile", required_argument, NULL, 4},
        {"graddir", required_argument, NULL, 5},
        {"gaussian", required_argument, NULL, 6},
        {"increments", required_argument, NULL, 7},
        {"radius", required_argument, NULL, 8},
        {"printparticleposition", no_argument, NULL, 9},
        {"incorrectlinear_ADCx", required_argument, NULL, 10},
        {"dontuseT2", no_argument, NULL, 11},
        {"incorrectlinear_ADCz", required_argument, NULL, 12},
        {"incorrectlinear_FA", required_argument, NULL, 13},
        {"num_of_trial", required_argument, NULL, 14},
        {"num_of_repeat", required_argument, NULL, 15},
        {"num_of_cyl", required_argument, NULL, 16},
        {"max_radius", required_argument, NULL, 17},
        {"min_radius", required_argument, NULL, 18},
        {"cyl_seed", required_argument, NULL, 19},
        {"rand_mean", required_argument, NULL, 20},
        {"rand_stddev", required_argument, NULL, 21},
        {"incorrectlinear_ADCx", required_argument, NULL, 22},
    };

    /* This is the menu. It comes up if all you type in the command line is the
     * "executable name" (e.g. test6) */
    if (argc ==
        1) { // if arguments == 1, i.e. only the executable name has been typed.
        std::cout
            << "Monte Carlo Simulation of Diffusion with restricting geometries"
            << std::endl
            << std::endl;
        std::cout << "Usage: test [options] -geom <geometry packing> -num "
                     "<number of particles> "
                  << std::endl;
        std::cout << "       test -geom <geometry packing> -num <number of "
                     "particles> [options] "
                  << std::endl
                  << std::endl
                  << std::endl;
        std::cout << "Available Options are: " << std::endl << std::endl;
        std::cout
            << "       -geom <geometry packing>                         (0 = "
               "Szafer Box Lattice, 1 = Cylinder Box Lattice, no default)"
            << std::endl;
        std::cout << "       -num  <number of particles>                      "
                     "(no default)"
                  << std::endl;
        std::cout << "       -de <Extracellular Diffusion Coefficient>        "
                     "(geometry 0,1 default = 2.5e-6 mm^2/ms)"
                  << std::endl;
        std::cout << "       -di <Intracellular Diffusion Coefficient>        "
                     "(geometry 0 default = 1.0e-6 mm^2/ms)"
                  << std::endl;
        std::cout << "       -T2e <Extracellular T2>                          "
                     "(geometry 0,1 default = 200 ms)"
                  << std::endl;
        std::cout << "       -T2i <Intracellular T2>                          "
                     "(geometry 0 default = 70 ms)"
                  << std::endl;
        std::cout << "       -radius <Cell Cylinder Radius>                   "
                     "(geometry 1 default = 0 mm)"
                  << std::endl;
        std::cout << "       -gt <gradient time>                              "
                     "(geometry 0,1 default = 1 ms)"
                  << std::endl;
        std::cout << "       -gs <gradient spacing>                           "
                     "(if dtg = 0, this is the diffusion time. geometry 0,1 "
                     "default = 10 ms)"
                  << std::endl;
        std::cout << "       -csx <Cell Size in x dimension>                  "
                     "(geometry 0 default = 0)"
                  << std::endl;
        std::cout << "       -csy <Cell Size in y dimension>                  "
                     "(geometry 0 default = 0)"
                  << std::endl;
        std::cout << "       -perm <permeability>                             "
                     "(geometry 0,1 default =  .005)"
                  << std::endl;
        std::cout << "       -dt <timestep>                                   "
                     "(geometry 0,1 default =  .1 ms)"
                  << std::endl;
        std::cout << "       -ls <lattice size>                               "
                     "(geometry 0,1 default = 1 mm) "
                  << std::endl;
        std::cout << "       -grad <initial gradient strength>                "
                     "(geometry 0,1 default = 0.0001609114 T/mm)"
                  << std::endl;
        std::cout << "       -help                                            "
                     "(not available atm)                            "
                  << std::endl;
        std::cout
            << "       -dim <dimension>                                 "
               "(default = 3, no other dimensions are supported at this time)"
            << std::endl;
        std::cout << "       -printparticleposition                           "
                     "(default = false)"
                  << std::endl;
        std::cout << "       -printsignalinfo                                 "
                     "(default = false)"
                  << std::endl;
        std::cout << "       -printinitialvars                                "
                     "(default = false)"
                  << std::endl;
        std::cout << "       -printtofile <text file>                         "
                     "(default = false, prints (timestep, ADC) for graddir to "
                     "text file)"
                  << std::endl;
        std::cout << "       -graddir <gradient direction>                    "
                     "(default = 0(x axes), 1(y axes), 2(z axes)"
                  << std::endl;
        std::cout << "       -gaussian <text file>                            "
                     "(default = false) "
                  << std::endl;
        std::cout << "       -increments <number of times ADC is computed>    "
                     "(default = 1) "
                  << std::endl;
        std::cout << "       -hex <hexagonal sidelength>                      "
                     "(default = .5) "
                  << std::endl;
        std::cout << "       -dontuseT2                                       "
                     "(default = true) "
                  << std::endl;
        std::cout << "       -num_of_trial                                    "
                     "(default = true) "
                  << std::endl;
        std::cout << "       -num_of_repeat                                   "
                     "(default = 1) "
                  << std::endl;
        std::cout << "       -num_of_cyl                                      "
                     "(geom 3 default = 1) "
                  << std::endl;
        std::cout << "       -max_radius                                      "
                     "(geom 3 default = no default) "
                  << std::endl;
        std::cout << "       -min_radius                                      "
                     "(geom 3 default = no default) "
                  << std::endl;
        std::cout << "       -cyl_seed                                        "
                     "(geom 3 default = 12345678.12345678) "
                  << std::endl;
        exit(1);
    }

    /*
            if more than just the executable name is typed, this else block is
       executed it contains the instructions to be executed if one of the
       parameters is selected Usually the instructions just place the command
       line user-entered parameter into the variables struct under the
       corresponding name.
     */
    else {

        int errors = 0, count = 0;
        int optchar, longindex;

        while ((optchar = getopt_long_only(argc, argv, optv, longOpts,
                                           &longindex)) != -1) {

            switch (optchar) {

            case 'a':
                if (optarg == NULL)
                    ++errors;
                else
                    vars.geom = (int)atol(optarg);
                count++;

                break;
            case 'b':
                if (optarg == NULL)
                    ++errors;
                else
                    vars.num = (int)atol(optarg);
                count++;
                break;
            case 'c':
                if (optarg != NULL)
                    vars.de = atof(optarg);
                break;
            case 'd':
                if (optarg != NULL)
                    vars.di = atof(optarg);
                break;
            case 'e':
                if (optarg != NULL)
                    vars.ds = atof(optarg);
                break;
            case 'f':
                if (optarg != NULL)
                    vars.dc = atof(optarg);
                break;
            case 'g':
                if (optarg != NULL)
                    vars.rs = atof(optarg);
                break;
            case 'h':
                if (optarg != NULL)
                    vars.rc = atof(optarg);
                break;
            case 'i':
                if (optarg != NULL)
                    vars.T2e = atof(optarg);
                break;
            case 'j':
                if (optarg != NULL)
                    vars.T2i = atof(optarg);
                break;
            case 'k':
                if (optarg != NULL)
                    vars.T2s = atof(optarg);
                break;
            case 'l':
                if (optarg != NULL)
                    vars.T2c = atof(optarg);
                break;
            case 'm':
                if (optarg != NULL)
                    vars.gt = atof(optarg);
                break;
            case 'n':
                if (optarg != NULL)
                    vars.gs = atof(optarg);
                break;
            case 'p':
                if (optarg != NULL)
                    vars.csx = atof(optarg);
                break;
            case 'q':
                if (optarg != NULL)
                    vars.csy = atof(optarg);
                break;
            case 'r':
                if (optarg != NULL)
                    vars.csz = atof(optarg);
                break;

            case 's':
                if (optarg != NULL)
                    vars.p = atof(optarg);
                break;

            case 't':
                if (optarg != NULL)
                    vars.dt = atof(optarg);
                break;

            case 'u':
                if (optarg != NULL)
                    vars.ls = atof(optarg);
                break;

            case 'v':
                if (optarg != NULL)
                    vars.g = atof(optarg);
                break;
            case 'w':
                if (optarg != NULL)
                    vars.hex_sidelength = atof(optarg);
                break;
            case 0:
                vars.printparticleinfo = true;
                break;
            case 2:
                vars.printsignalinfo = true;
                break;
            case 3:
                vars.printinitialvars = true;
                break;
            case 4:
                vars.printtofile = true;
                vars.fh = fopen(optarg, "w");
                break;
            case 5:
                if (optarg != NULL)
                    vars.graddir = (int)atol(optarg);
                break;
            case 6:
                vars.gaussian = true;
                vars.fh_g = fopen(optarg, "w");
                break;
            case 7:
                if (optarg != NULL)
                    vars.increments = (int)atol(optarg);
                break;
            case 8:
                if (optarg != NULL)
                    vars.r = atof(optarg);
                break;
            case 9:
                vars.printparticleposition = true;
                vars.fh_pp = fopen(optarg, "w");
                break;
            case 10:
                vars.fh_incorrectlinreg_ADCx = fopen(optarg, "w");
                break;
            case 11:
                vars.use_T2 = false;
                break;
            case 12:
                vars.fh_incorrectlinreg_ADCz = fopen(optarg, "w");
                break;
            case 13:
                vars.fh_incorrectlinreg_FA = fopen(optarg, "w");
                break;
            case 14:
                vars.num_of_trial = (int)atol(optarg);
                break;
            case 15:
                vars.num_of_repeat = (int)atol(optarg);
                break;
            case 16:
                vars.num_of_cyl = (int)atol(optarg);
                break;
            case 17:
                vars.max_radius = atof(optarg);
                break;
            case 18:
                vars.min_radius = atof(optarg);
                break;
            case 19:
                vars.cyl_seed = (int)atol(optarg);
                break;
            case 20:
                vars.rand_mean = atof(optarg);
                break;
            case 21:
                vars.rand_stddev = atof(optarg);
            case 22:
                vars.fh_incorrectlinreg_ADCy = fopen(optarg, "w");
                break;

            case '?':
                if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr, "Unknown option character `\\x%x'.\n",
                            optopt);
                exit(1);
            }
        }
        if (count != 2) {
            std::cout << "ERROR: both -geom and -num parameters are required"
                      << std::endl;
            exit(1);
        }
    }

    return vars;
}

/* Prints all of the variables to the screen*/
void printinitvariables(VARIABLES vars) {
    std::cout << std::endl << "INITIAL VARIABLES" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "Dextra " << vars.de << std::endl;
    std::cout << "Dintra" << vars.di << std::endl;
    std::cout << "T2extra " << vars.T2e << std::endl;
    std::cout << "T2intra " << vars.T2i << std::endl;
    std::cout << "radius of cylinder" << vars.r << std::endl;
    std::cout << "Gradient time " << vars.gt << std::endl;
    std::cout << "Gradient Spacing " << vars.gs << std::endl;
    std::cout << "Number of Trials " << vars.num_of_trial << std::endl;
    std::cout << "Number of Repeats " << vars.num_of_repeat << std::endl;
    std::cout << "Cell size x " << vars.csx << std::endl;
    std::cout << "Cell size y " << vars.csy << std::endl;
    std::cout << "Cell size z " << vars.csz << std::endl;
    std::cout << "Permeability " << vars.p << std::endl;
    std::cout << "Timestep " << vars.dt << std::endl;
    std::cout << "Lattice Size " << vars.ls << std::endl;
    std::cout << "Initial Gradient Strength" << vars.g << std::endl;
    std::cout << "Hexagonal Sidelength " << vars.hex_sidelength << std::endl;
}
