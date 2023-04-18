/* ----------------------------------------------------------------------
   Vacancy tries to build a model for vacancy diffusion. The calculation are based on the
   diffusion model and LAKIMOCA.
   By SMZ, 2015-12.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(vacancy, AppVacancy)

#else

#ifndef SPK_APP_VACANCY_H
#define SPK_APP_VACANCY_H

#include "app_lattice.h"
#include "swarg.h"

#define NPT_EAM 6500   /* max number of point of the pair or density parts of the potential */
#define NPT_EMBED 6200 /* max number of point of the embeded part of the EAM potential */
#define NSPECIES 6

extern "C"
{
void kernel_site_propensity(struct _sw_site_propensity_arg * arg);

  void kernel_func(struct _swarg *arg);

  void c_call_athread_join();
}

namespace SPPARKS_NS
{

  class AppVacancy : public AppLattice
  {
    // friend class DiagVacancy;

  public:
    AppVacancy(class SPPARKS *, int, char **);
    ~AppVacancy();
    void input_app(char *, int, char **);
    void grow_app();
    void init_app();
    void setup_app();

    double site_energy(int);
    double site_propensity(int, double *, int *, int *);
    void site_event(int, int *, class RandomPark *);

    void count_vac(int, int *, cvector_vector_type(int) *);
    int type(int);
    void transport();
    void box2sub(int *, int *, int *);
    void set_it(int, int, int, int);
    void set_it_e(int, int, int, int, double, double, double);
    int set_it_e_set(int, int, int, int, double, double, double, int &);
    void get_it_e(int, int, int, int *, double *, double *, double *);

  private:
    int *lattice;
    int allocated;

    struct Event
    {                    // one event for an owned site
      double propensity; // propensity of this event
      int destination;   // local ID of destination site
      int next;          // index of next event for this site
    };

    Event *events;   // list of events for all owned sites
    int nevents;     // # of events for all owned sites
    int maxevent;    // max # of events list can hold
    int *firstevent; // index of 1st event for each owned site
    int freeevent;   // index of 1st unused event in list

    double *E_V;
    double *E_R;
    double *E_S;

  public:
    double unitdis;
    int boxx, boxy, boxz, boxxlo, boxylo, boxzlo, boxxhi, boxyhi, boxzhi;
    int subx, suby, subz, subxlo, subylo, subzlo, subxhi, subyhi, subzhi, subxlen, subylen, subzlen;
    int xghost, yghost, zghost;

    int ***POS_ID;
    double *barrier; // barrier for fe, ni, si, mn, cu
    double a;
    int tag_COHESIVE_MODEL;
    std::map<int, int> wait_send;
    MPI_Datatype ctype;

    double vij1[NSPECIES][NSPECIES]; /* vij1 larger for several solute type Cu,Ni,Mn,Si */
    double rho1[NSPECIES][NSPECIES];
    double srho1[NSPECIES][NSPECIES];

    double vij2[NSPECIES][NSPECIES], rho2[NSPECIES][NSPECIES];
    double vij3[NSPECIES][NSPECIES], rho3[NSPECIES][NSPECIES], vij4[NSPECIES][NSPECIES], rho4[NSPECIES][NSPECIES];
    double vij[5][NSPECIES][NSPECIES], rho[5][NSPECIES][NSPECIES];
    double srho2[NSPECIES][NSPECIES], srho3[NSPECIES][NSPECIES], srho4[NSPECIES][NSPECIES];
    double srho[5][NSPECIES][NSPECIES];

    /* specific for tabulated EAM pot */
    double FEMBED[3][NPT_EMBED]; /* fct entourage */
    double FEMBED_dr[3];
    double FEMBED_rmin[3], FEMBED_rmax[3];

    /* specific for modified EAM pot */
    double SEMBED[3][NPT_EMBED]; /* fct entourage */
    double SEMBED_dr[3];
    double SEMBED_rmin[3], SEMBED_rmax[3];

    struct _sw_site_propensity_arg argOfSite_propensity;
    struct  _swarg argOfKernel;
    

    int get_POS_ID(int, int, int);
    int get_POS_ID_v2(int, int, int, int &);

    void read_pot_file(char *);
    double interpol(int, double *, double *, double);

    void init_potential_pair_file(char *);

    void find_border();
    void get_neighbor(int, int *);

    void compute_energy();
    double single_energy(int, double, double, double);
    void eatom(int, int *, int *, int *, double *, double *, double *);
    void set_e(int, double, double, double);
    void compute_e(int);
    int get_cu_alone(int);

    double delete1(int, double *, double *, double *);
    double calcul_de(int, int, int, int, int *);
    double calcul_dene_v1(int, int, int, int, int, double);
    double calcul_dene_v2(int, int, int, int, int, int);
    void do_jump(int, int);
    void do_jump_v1(int, int, int, int, int, double);
    void do_jump_v2(int, int, int, int, int, int);

    double eval_embed(int, double);
    double eval_sembed(int, double);

    void clear_events(int);
    void add_event(int, int, double);

    void allocate_data();

    double compute_linedis(int, int);
    void inc_px(int, int, int *);
    void inc_py(int, int, int *);
    void inc_pz(int, int, int *);
    int is_in_set(int);
    int in_which_set(int);

    void add_to_wait(int, int, int);
    void add_to_wait4TypeChange(int, int, int);
    void add_to_send();
    void communicate();
    void communicate_priority(int);
    void communicate_perform_set(int);

    // not use
    void get_e(int, double *, double *, double *);

    void update_site_propensity(cvector_vector_type(int), double *, double **,
                                int **, int *);
    void call_athread_join();
    void set_arg_of_kernel_func(int *, int *, int, int, int *, int);
    void update_wait_send_map(int, int);
  };

}

#endif
#endif

    /* ERROR/WARNING messages:

    E: Illegal ... command

    Self-explanatory.  Check the input script syntax and compare to the
    documentation for the command.  You can use -echo screen as a
    command-line option when running SPPARKS to see the offending
    line.

    E: Cannot use %s command until sites exist

    This command requires sites exist before using it in an input script.

    E: Can only use ecoord command with app_style diffusion nonlinear

    Self-explanatory.

    E: Cannot define Schwoebel barrier without Schwoebel model

    Self-explanatory.

    E: Unrecognized command

    The command is assumed to be application specific, but is not
    known to SPPARKS.  Check the input script.

    E: Cannot perform deposition in parallel

    UNDOCUMENTED

    E: Cannot perform deposition with multiple sectors

    UNDOCUMENTED

    E: One or more sites have invalid values

    The application only allows sites to be initialized with specific
    values.

    E: Did not reach event propensity threshhold

    UNDOCUMENTED

    E: BAD DONE

    UNDOCUMENTED

    */
