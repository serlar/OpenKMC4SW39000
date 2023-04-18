/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_APP_LATTICE_H
#define SPK_APP_LATTICE_H

#include "stdio.h"
#include "app.h"

#include <stack>
#include <map> //baodong

#include "cvector.h" // xulei

extern "C"
{
  void call_CRTS_sync_master_array();
}

namespace SPPARKS_NS
{

  class AppLattice : public App
  {
    friend class CommLattice;

  public:
    int delpropensity; // # of sites away needed to compute propensity
    int delevent;      // # of sites away affected by an event

    int nmax;       // max # of sites per-site arrays can store
    int maxneigh;   // max neighbors of any site in entire system
    int *numneigh;  // # of neighbors of each site
    int **neighbor; // local indices of neighbors of each site

    class CommLattice *comm;
    int nset; // # of subsets of lattice sites

    AppLattice(class SPPARKS *, int, char **);
    virtual ~AppLattice();
    void input(char *, int, char **);
    void init();
    void setup();
    void iterate();
    // baodong

    void grow(int);
    void add_site(tagint, double, double, double);
    void add_ghost(tagint, double, double, double);
    void add_neighbors(int, int, char **);
    void add_values(int, char **);
    void add_value(int, int, int, char *);
    void print_connectivity();

    // pure virtual functions, must be defined in child class

    virtual void grow_app() = 0;
    virtual double site_energy(int) = 0;
    // baodong
    //  virtual void site_event_rejection(int, class RandomPark *) = 0;
    virtual double site_propensity(int, double *, int *, int *) = 0;
    virtual void site_event(int, int *, class RandomPark *) = 0;

    // virtual functions, may be overridden by child class

    virtual void input_app(char *, int, char **);
    virtual void init_app() {}
    virtual void setup_app() {}
    virtual void setup_end_app() {}
    virtual void *extract_app(char *) { return NULL; }

    virtual void push_new_site(int, int *, int, std::stack<int> *);
    virtual void push_connected_neighbors(int, int *, int, std::stack<int> *);
    virtual void connected_ghosts(int, int *, class Cluster *, int);

    // baodong
    virtual void compute_e(int) {}
    virtual int get_cu_alone(int) {}
    virtual void get_e(int, double *, double *, double *) {}
    virtual void compute_e(int, int) {}
    virtual void count_vac(int, int *, cvector_vector_type(int) *) {}
    virtual int type(int) {}
    virtual void communicate() {}
    virtual void communicate_priority(int) = 0;
    virtual void communicate_perform_set(int) = 0;
    virtual void get_neighbor(int, int *) {}
    virtual void user_update(double) {}
    virtual void create_comm_all() {}
    virtual void perform_comm_all() {}
    virtual void set_it(int, int, int, int) {}
    virtual void transport() {}

    // xulei
    virtual int is_in_set(int) {}
    virtual void update_site_propensity(cvector_vector_type(int), double *, double **,
                                        int **, int *) {}
    virtual void call_athread_join() = 0;

  protected:
    int me, nprocs;

    bigint naccept, nattempt;      // number of accepted and attempted events
    bigint nsweeps;                // number of sweeps performed
    double temperature, t_inverse; // temperature settings
    double dt_sweep;               // rKMC time for nglobal attemped events
    double dt_rkmc;                // rKMC time for one pass thru all sectors
    double dt_kmc;                 // KMC time for one pass thru all sectors
    double dt_step;                // KMC time for single global KMC step

    int allow_kmc;       // 1 if app supports KMC
    int allow_rejection; // 1 if app supports rejection KMC
    int allow_masking;   // 1 if app supports rKMC masking
    int allow_update;    // 1 if app provides user_update()
    int numrandom;       // # of RN used by rejection routine

    int sweepflag;    // set if rejection KMC solver
    int sectorflag;   // 1 if partition my domain into sectors
    int nsector;      // 1,2,4,8 = # of sectors
    int nsector_user; // 0 if default, else 2,4,8
    int ncolors;      // # of colors, depends on lattice
    int bothflag;     // 1 if both sectors and colors
    // baodong
    int update_only; // 1 if skip other iteration techniques

    class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
    class RandomPark *ranstrict; // RN generator for per-site strict rKMC
    int *siteseeds;              // per-site seeds for ransite
    int *sitelist;               // randomized list of site indices

    // baodong
    // bool Lmask;                  // masking on/off
    // char *mask;                  // size of nlocal + nghost sites

    double Ladapt; // adaptive sector time increments for KMC
    double tstop;  // requested time increment in sector
    double nstop;  // requested events per site in sector

    // arrays for owned + ghost sites
    // int *owner;                  // proc who owns the site, delte by smz 170210
    // int *index;                  // index of site on owning proc

    double *propensity;      // probabilities for each owned site
    double **propensity_dir; // probabilities for each owned site
    int **jumpsite_dir;      // probabilities for each owned site
    int *num_dir;

    cvector_vector_type(int) *vacant = NULL;
    // std::map<int,int> *vacant;

    int current_set;
    double xmid;
    double ymid;
    double zmid;

    struct Set
    {              // subset of lattice sites I own
      int nlocal;  // # of owned sites in set
      int nselect; // # of selections from set for rKMC
      int nloop;   // # of loops over set for rKMC
      int nborder; // # of sites with non-set site as neighbor
      int *border; // lattice index for each border site
      // int *bsites;               // list of border sites to pass to solver
      class Solve *solve;      // KMC solver
      double *propensity;      // propensities for set sites
      double **propensity_dir; // propensities for set sites
      int **jumpsite_dir;      // probabilities for each owned site
      int *num_dir;
      int *site2i; // map from set sites to lattice index
      // int *i2site;               // map from lattice index to set sites
      // std::map<int,int> vacant;
      cvector_vector_type(int) vacant = NULL;
    };
    Set *set; // list of subsets

    struct Site
    {
      int id, proc, index;
      double x, y, z;
    };

    void iterate_kmc_global(double);
    void iterate_kmc_sector(double);

    void ghosts_from_connectivity();
    void connectivity_within_cutoff();

    void create_set(int, int, class Solve *);
    class Solve *free_set(int);

    void stats(char *);
    void stats_header(char *);

    // baodong
    void set_sector(int, char **);
    void set_temperature(int, char **);

    void bounds(char *, int, int, int &, int &);
  };

}

#endif

/* ERROR/WARNING messages:

E: App needs a KMC or rejection KMC solver

You must define either a solver or sweep option.

E: App cannot use both a KMC and rejection KMC solver

You cannot define both a solver and sweep option.

E: KMC events are not implemented in app

Not every application supports KMC solvers.

E: Rejection events are not implemented in app

Self-explanatory.

E: Mask logic not implemented in app

Not every application supports masking.

E: Cannot use KMC solver in parallel with no sectors

Self-explanatory.

E: Cannot use random rejection KMC in parallel with no sectors

Self-explanatory.

E: Cannot use raster rejection KMC in parallel with no sectors

Self-explanatory.

E: Cannot use color/strict rejection KMC with sectors

Self-explanatory.

E: App did not set dt_sweep

Internal SPPARKS error.

E: Invalid number of sectors

Self-explanatory.

E: Cannot color without a lattice definition of sites

UNDOCUMENTED

E: Cannot color without contiguous site IDs

UNDOCUMENTED

E: Cannot color this combination of lattice and app

Coloring is not supported on this lattice for the neighbor
dependencies of this application.

E: Choice of sector stop led to no rKMC events

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: App does not permit user_update yes

UNDOCUMENTED

E: Per-processor system is too big

UNDOCUMENTED

*/
