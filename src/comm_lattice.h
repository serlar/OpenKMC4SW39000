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

#ifndef SPK_COMM_LATTICE_H
#define SPK_COMM_LATTICE_H

#include "mpi.h"
#include "pointers.h"
#include <map>
#include "cvector.h"

namespace SPPARKS_NS
{

  class CommLattice : protected Pointers
  {
  public:
    // add by xulei for
    cvector_vector_type(int) listVacanciesUpdated[8];

    CommLattice(class SPPARKS *);
    ~CommLattice();
    void init(int, int, int, int *);
    // void all();
    // void all_reverse();
    // void sector(int);
    // void reverse_sector(int);

    void create_set();
    void create_all();
    void set_iset(int);
    void add_to_send(int, int, int, int, int);
    void perform_set();
    void perform_all();

    void perform_priority(int);
    void set_PriorityTargetList(int *);
    void set_PrioritySourceList(int *);
    void get_Priority_target(int, int *, int **);
    void get_Priority_source(int, int *, int **);

  private:
    // MPI_Comm comm_3d;
    int me, nprocs;
    int delghost, delreverse;

    int PriorityTargetList[24];
    int PrioritySourceList[24];

    struct Particle
    {
      int x, y, z, it, change_type; // 1 means type_change, 0 means energy_change.
      double ev, er, es;
    };

    int ninteger, ndouble;

    int **iarray;
    double **darray;

    // add
    int nset;
    std::map<int, int *> needer;
    struct Swap
    {
      int nsend, nrecv; // number of messages to send/recv
      int *sproc;       // proc for each send message
      int *scount;      // size of each send message in sites
      int *smax;
      std::map<int, int> shash;
      Particle **spart;
      int *rproc;  // proc for each recv message
      int *rcount; // size of each recv message in sites
      std::map<int, int> rhash;
      Particle **rpart;
      MPI_Request *request; // MPI datums for each recv message
      MPI_Status *status;
    };
    Swap *setswap;
    Swap *cur_setswap;
    Swap allswap;
    MPI_Datatype ctype;

    // Swap *create_swap_all();
    // void free_swap(Swap *);

    // void create_send_from_list(int, Site *, Swap *);
    // void create_send_from_recv(tagint, tagint, Site *, Swap *);
    // void create_recv_from_send(int, int, Site *, Swap *);
    // void create_recv_from_list(tagint, Site *, Swap *);

    // void perform_swap_site(Swap *);
    // void perform_swap_int(Swap *);
    // void perform_swap_double(Swap *);
    // void perform_swap_general(Swap *);
    // add
    void new_type_particle(MPI_Datatype *);
    void destroy_map(std::map<int, int *> *);
    void free_swap(Swap *);
  };

}

#endif

/* ERROR/WARNING messages:

E: Site-site interaction was not found

Internal SPPARKS error.

*/
