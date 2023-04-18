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

#include "spktype.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice.h"
#include "app.h"
#include "app_lattice.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "math.h"
#include "timer.h"
#include "app_vacancy.h"
#include "lattice.h"

#include <map>
#include <set>

using namespace SPPARKS_NS;

#define DELTA 16384
#define PARTNUM 100

/* ---------------------------------------------------------------------- */

CommLattice::CommLattice(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  for (int i = 0; i < 8; i++)
  {
    listVacanciesUpdated[i] = NULL;
  }
  nset = 0;
  setswap = NULL;
}

/* ---------------------------------------------------------------------- */

CommLattice::~CommLattice()
{
  // printf("me= %d Numscount= %d Zeroscount= %d Numrcount= %d Zerorcount= %d\n", me, Numscount, Zeroscount, Numrcount, Zerorcount); //wubd 2016-12-27

  // if (allswap) free_swap(allswap);
  // add
  MPI_Type_free(&ctype);

  if (setswap)
  {
    for (int i = 0; i < nset; i++)
    {
      free_swap(&setswap[i]);
    }
    delete[] setswap;
  }
  free_swap(&allswap);
  destroy_map(&needer);

  for (int i = 0; i < 8; i++)
  {
    cvector_free(listVacanciesUpdated[i]);
  }
}

/* ----------------------------------------------------------------------
   setup comm pattern
   3 kinds of patterns: all ghosts, sector ghosts, sector reverse
   if nsector_request = 1, just ghosts for entire proc domain
   if nsector_request > 1, do all and sector ghosts, reverse if needed
   array = NULL = communicate iarray/darray from app
   array = non-NULL = communicate passed-in array (from diagnostic)
------------------------------------------------------------------------- */

void CommLattice::init(int nsector_request, int delpropensity, int delevent,
                       int *array)
{
  // comm_3d = domain->comm_3d;
  delghost = delpropensity;
  delreverse = delevent;

  AppLattice *applattice = (AppLattice *)app;

  ninteger = app->ninteger;
  ndouble = app->ndouble;
  iarray = app->iarray;
  darray = app->darray;

  nset = nsector_request;
  setswap = NULL;

  new_type_particle(&ctype);
  // create new swaps as requested
  create_all();
  create_set();
}

/* ----------------------------------------------------------------------
   acquire ghost values for entire proc sub-domain
------------------------------------------------------------------------- */
// delete by smz 170302
// void CommLattice::all()

/* ----------------------------------------------------------------------
   reverse communicate changed border values for entire proc sub-domain
------------------------------------------------------------------------- */
// deleted by bdwu 170228
// void CommLattice::all_reverse()

/* ----------------------------------------------------------------------
   acquire ghost values for one sector
------------------------------------------------------------------------- */
// delete by smz 161017
// void CommLattice::sector(int isector)

/* ----------------------------------------------------------------------
   reverse communicate changed border values for one sector
------------------------------------------------------------------------- */
// delete by smz 161017
// void CommLattice::reverse_sector(int isector)

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost sites for my entire subdomain
------------------------------------------------------------------------- */
// delete by smz 170302
// CommLattice::Swap *CommLattice::create_swap_all()

/* ----------------------------------------------------------------------
   create send portion of a Swap communication pattern
   start with list of sites I need to recv
   circulate list to all procs
   when comes back to me, will also be flagged with who I recv sites from
------------------------------------------------------------------------- */
// delete by smz 170302
// void CommLattice::create_send_from_recv(tagint nsite, tagint maxsite,
// 					Site *buf, Swap *swap)

/* ----------------------------------------------------------------------
   create recv portion of a Swap communication pattern
   create from list of sites I need to recv and procs who will send them
   no communication required
------------------------------------------------------------------------- */
// delete by smz 170302
// void CommLattice::create_recv_from_list(tagint nsite, Site *buf, Swap *swap)

/* ----------------------------------------------------------------------
   free a Swap object
------------------------------------------------------------------------- */
// delete by smz 170302
// void CommLattice::free_swap(Swap *swap)

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use iarray and darray as source/destination, mixed integer/double data
------------------------------------------------------------------------- */
// delete by smz 170302
// void CommLattice::perform_swap_general(Swap *swap)

void CommLattice::create_all()
{
  AppLattice *applattice = (AppLattice *)app;
  int nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  double **xyz = applattice->xyz;

  int XN, XP, YN, YP, ZN, ZP;
  int XNYN, XNYP, XPYN, XPYP, XNZN, XNZP, XPZN, XPZP, YNZN, YNZP, YPZN, YPZP;
  int XNYNZN, XNYNZP, XNYPZN, XNYPZP, XPYNZN, XPYNZP, XPYPZN, XPYPZP;
  int XDELTA, YDELTA, ZDELTA, _XDELTA, _YDELTA, _ZDELTA;

  XN = domain->procneigh[0][0];
  XP = domain->procneigh[0][1];
  YN = domain->procneigh[1][0];
  YP = domain->procneigh[1][1];
  ZN = domain->procneigh[2][0];
  ZP = domain->procneigh[2][1];

  XDELTA = XP - me;
  YDELTA = YP - me;
  ZDELTA = ZP - me;
  _XDELTA = XN - me;
  _YDELTA = YN - me;
  _ZDELTA = ZN - me;

  XPYP = YP + XDELTA;
  XPYN = YN + XDELTA;
  XNYP = YP + _XDELTA;
  XNYN = YN + _XDELTA;

  XPZP = ZP + XDELTA;
  XPZN = ZN + XDELTA;
  XNZP = ZP + _XDELTA;
  XNZN = ZN + _XDELTA;

  YPZP = ZP + YDELTA;
  YPZN = ZN + YDELTA;
  YNZP = ZP + _YDELTA;
  YNZN = ZN + _YDELTA;

  XPYPZP = XPYP + ZDELTA;
  XPYPZN = XPYP + _ZDELTA;
  XPYNZP = XPYN + ZDELTA;
  XPYNZN = XPYN + _ZDELTA;
  XNYPZP = XNYP + ZDELTA;
  XNYPZN = XNYP + _ZDELTA;
  XNYNZP = XNYN + ZDELTA;
  XNYNZN = XNYN + _ZDELTA;

  int nsend = 0;
  int nrecv = 0;
  int procs[26] = {XN, XP, YN, YP, ZN, ZP, XNYN, XNYP, XPYN, XPYP, XNZN, XNZP, XPZN, XPZP, YNZN, YNZP, YPZN, YPZP, XNYNZN, XNYNZP, XNYPZN, XNYPZP, XPYNZN, XPYNZP, XPYPZN, XPYPZP};

  int PrioritySequence[24] = {XN, XP, YN, XPYN, XN, XP, YP, ZN, XPYP, XPZN, YPZN, XPYPZN, XN, XP, YN, XPYN, XN, XP, YP, ZP, XPYP, XPZP, YPZP, XPYPZP};
  set_PriorityTargetList(&PrioritySequence[0]);
  int PrioritySequence2[24] = {XP, XN, YP, XNYP, XP, XN, YN, ZP, XNYN, XNZP, YNZP, XNYNZP, XP, XN, YP, XNYP, XP, XN, YN, ZN, XNYN, XNZN, YNZN, XNYNZN};
  set_PrioritySourceList(&PrioritySequence2[0]);

  for (int i = 0; i < 26; i++)
  {
    if (procs[i] != me && allswap.shash.find(procs[i]) == allswap.shash.end())
    {
      allswap.shash.insert(std::pair<int, int>(procs[i], nsend++));
    }
    if (procs[i] != me && allswap.rhash.find(procs[i]) == allswap.rhash.end())
    {
      allswap.rhash.insert(std::pair<int, int>(procs[i], nrecv++));
    }
  }

  allswap.nsend = nsend;
  allswap.nrecv = nrecv;
  allswap.sproc = new int[nsend];
  allswap.rproc = new int[nrecv];
  allswap.spart = new Particle *[nsend];
  allswap.rpart = new Particle *[nrecv];
  allswap.scount = new int[nsend];
  allswap.rcount = new int[nrecv];
  allswap.smax = new int[nsend];
  allswap.request = new MPI_Request[nsend];
  allswap.status = new MPI_Status[nrecv];

  std::map<int, int>::iterator it_map;
  for (int i = 0; i < nsend; i++)
  {
    for (it_map = allswap.shash.begin(); it_map != allswap.shash.end(); it_map++)
    {
      allswap.sproc[it_map->second] = it_map->first;
    }
    allswap.scount[i] = 0;
    allswap.smax[i] = PARTNUM;
    allswap.spart[i] = (Particle *)memory->smalloc(PARTNUM * sizeof(Particle), "app:particle");
  }

  for (int i = 0; i < nrecv; i++)
  {
    for (it_map = allswap.rhash.begin(); it_map != allswap.rhash.end(); it_map++)
    {
      allswap.rproc[it_map->second] = it_map->first;
    }
    allswap.rcount[i] = 0;
    allswap.rpart[i] = new Particle[PARTNUM];
  }

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;
  double boxxhi = domain->boxxhi;
  double boxyhi = domain->boxyhi;
  double boxzhi = domain->boxzhi;
  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;
  double xlattice = domain->lattice->xlattice;
  double ylattice = domain->lattice->ylattice;
  double zlattice = domain->lattice->zlattice;
  double unitdis = 2 / domain->lattice->getLatconst();

  double xpb, xnb, ypb, ynb, zpb, znb;
  xnb = subxlo + 1.99 * xlattice;
  xpb = subxhi - 2.01 * xlattice;
  ynb = subylo + 1.99 * ylattice;
  ypb = subyhi - 2.01 * ylattice;
  znb = subzlo + 1.99 * zlattice;
  zpb = subzhi - 2.01 * zlattice;

  double x, y, z;
  int xdir = -1, ydir = -1, zdir = -1, xydir = -1, xzdir = -1, yzdir = -1, xyzdir = -1;

  std::set<int> tmp; // to remove identical proc
  std::set<int>::iterator it_set;
  for (int i = 0; i < nlocal + nghost; i++)
  {
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];

    if (x >= xpb)
      xdir = XP;
    else if (x < xnb)
      xdir = XN;
    else
      xdir = -1;
    if (y >= ypb)
      ydir = YP;
    else if (y < ynb)
      ydir = YN;
    else
      ydir = -1;
    if (z >= zpb)
      zdir = ZP;
    else if (z < znb)
      zdir = ZN;
    else
      zdir = -1;

    if (subxlo == boxxlo && x >= boxxhi - 2.01 * xlattice)
      xdir = XN;
    if (subxhi == boxxhi && x < boxxlo + 1.99 * xlattice)
      xdir = XP;
    if (subylo == boxylo && y >= boxyhi - 2.01 * ylattice)
      ydir = YN;
    if (subyhi == boxyhi && y < boxylo + 1.99 * ylattice)
      ydir = YP;
    if (subzlo == boxzlo && z >= boxzhi - 2.01 * zlattice)
      zdir = ZN;
    if (subzhi == boxzhi && z < boxzlo + 1.99 * zlattice)
      zdir = ZP;

    if (xdir == XP && ydir == YP)
      xydir = XPYP;
    else if (xdir == XP && ydir == YN)
      xydir = XPYN;
    else if (xdir == XN && ydir == YP)
      xydir = XNYP;
    else if (xdir == XN && ydir == YN)
      xydir = XNYN;
    else
      xydir = -1;
    if (xdir == XP && zdir == ZP)
      xzdir = XPZP;
    else if (xdir == XP && zdir == ZN)
      xzdir = XPZN;
    else if (xdir == XN && zdir == ZP)
      xzdir = XNZP;
    else if (xdir == XN && zdir == ZN)
      xzdir = XNZN;
    else
      xzdir = -1;
    if (ydir == YP && zdir == ZP)
      yzdir = YPZP;
    else if (ydir == YP && zdir == ZN)
      yzdir = YPZN;
    else if (ydir == YN && zdir == ZP)
      yzdir = YNZP;
    else if (ydir == YN && zdir == ZN)
      yzdir = YNZN; //?
    else
      yzdir = -1;

    if (zdir == ZP)
    {
      if (xydir == XPYP)
        xyzdir = XPYPZP;
      else if (xydir == XPYN)
        xyzdir = XPYNZP;
      else if (xydir == XNYP)
        xyzdir = XNYPZP;
      else if (xydir == XNYN)
        xyzdir = XNYNZP;
      else
        xyzdir = -1;
    }
    else if (zdir == ZN)
    {
      if (xydir == XPYP)
        xyzdir = XPYPZN;
      else if (xydir == XPYN)
        xyzdir = XPYNZN;
      else if (xydir == XNYP)
        xyzdir = XNYPZN;
      else if (xydir == XNYN)
        xyzdir = XNYNZN;
      else
        xyzdir = -1;
    }
    else
      xyzdir = -1;

    // insert into needer
    int dirs[7] = {xdir, ydir, zdir, xydir, xzdir, yzdir, xyzdir};
    int *proc = new int[7];
    for (int j = 0; j < 7; j++)
    {
      proc[j] = -1;
      if (dirs[j] != -1 && dirs[j] != me)
      {
        tmp.insert(dirs[j]);
      }
    }
    int k = 0;
    for (it_set = tmp.begin(); it_set != tmp.end(); it_set++)
    {
      proc[k++] = *it_set;
    }
    tmp.clear();
    if (k > 0)
      needer.insert(std::pair<int, int *>(i, proc));

    // insert into send queue
    if (i < nlocal && k > 0)
    {
      int xx = round(x * unitdis);
      int yy = round(y * unitdis);
      int zz = round(z * unitdis);
      Particle particle = {xx, yy, zz, iarray[0][i], 0, darray[0][i], darray[1][i], darray[2][i]};
      for (int j = 0; j < 7 && proc[j] != -1; j++)
      {
        int p = proc[j];
        int index = allswap.shash[p];
        if (allswap.scount[index] == allswap.smax[index])
        {
          allswap.smax[index] += PARTNUM;
          allswap.spart[index] = (Particle *)memory->srealloc(allswap.spart[index], allswap.smax[index] * sizeof(Particle), "app:particle");
        }
        allswap.spart[index][allswap.scount[index]] = particle;
        allswap.scount[index]++;
      }
    }
  }
}

void CommLattice::perform_all()
{

  AppVacancy *appvacancy = (AppVacancy *)app;

  int i, j;

  for (i = 0; i < allswap.nsend; i++)
  {
    for (j = 0; j < allswap.scount[i]; j++)
    {
      Particle *particle = &(allswap.spart[i][j]);
      int x = particle->x;
      int y = particle->y;
      int z = particle->z;
      int it;
      double ev, es, er;
      appvacancy->get_it_e(x, y, z, &it, &ev, &er, &es);
      particle->it = it;
      particle->ev = ev;
      particle->er = er;
      particle->es = es;
    }
  }
  MPI_Status status[allswap.nsend];
  for (i = 0; i < allswap.nsend; i++)
  {
    MPI_Isend(allswap.spart[i], allswap.scount[i], ctype, allswap.sproc[i], 0, world, &allswap.request[i]);
  }

  for (i = 0; i < allswap.nrecv; i++)
  {
    MPI_Probe(allswap.rproc[i], 0, world, &allswap.status[i]);
    MPI_Get_count(&allswap.status[i], ctype, &allswap.rcount[i]);
    delete[] allswap.rpart[i];
    allswap.rpart[i] = new Particle[allswap.rcount[i]];
    MPI_Recv(allswap.rpart[i], allswap.rcount[i], ctype, allswap.rproc[i], 0, world, &allswap.status[i]);
    for (int j = 0; j < allswap.rcount[i]; j++)
    {
      int x = allswap.rpart[i][j].x;
      int y = allswap.rpart[i][j].y;
      int z = allswap.rpart[i][j].z;
      int it = allswap.rpart[i][j].it;
      double ev = allswap.rpart[i][j].ev;
      double er = allswap.rpart[i][j].er;
      double es = allswap.rpart[i][j].es;

      appvacancy->set_it_e(x, y, z, it, ev, er, es);
    }
  }
  MPI_Waitall(allswap.nsend, allswap.request, status);
}

void CommLattice::perform_priority(int iset)
{
  AppVacancy *appvacancy = (AppVacancy *)app;

  int i, j, k, flag_Iprobe; // wubd 2016-12-29
  int flag;
  int pos;

  int n;
  int *ptr;
  get_Priority_target(iset, &n, &ptr);

  MPI_Status status[n];
  MPI_Request request[n];

  for (i = 0; i < n; i++)
  {

    j = 0;
    while (cur_setswap->sproc[j] != ptr[i])
    {
      j++;
    }
    MPI_Isend(cur_setswap->spart[j], cur_setswap->scount[j], ctype, cur_setswap->sproc[j], 0, world, &request[i]);
    cur_setswap->scount[j] = 0;
  }

  get_Priority_source(iset, &n, &ptr);
  for (k = 0; k < n; k++)
  {
    i = 0;
    while (cur_setswap->rproc[i] != ptr[k])
    {
      i++;
    }

    MPI_Probe(cur_setswap->rproc[i], 0, world, &status[k]);
    MPI_Get_count(&status[k], ctype, &cur_setswap->rcount[i]);
    delete[] cur_setswap->rpart[i];
    cur_setswap->rpart[i] = new Particle[cur_setswap->rcount[i]];
    MPI_Recv(cur_setswap->rpart[i], cur_setswap->rcount[i], ctype, cur_setswap->rproc[i], 0, world, &status[k]);
    for (int j = 0; j < cur_setswap->rcount[i]; j++)
    {
      int x = cur_setswap->rpart[i][j].x;
      int y = cur_setswap->rpart[i][j].y;
      int z = cur_setswap->rpart[i][j].z;
      int it = cur_setswap->rpart[i][j].it;
      double ev = cur_setswap->rpart[i][j].ev;
      double er = cur_setswap->rpart[i][j].er;
      double es = cur_setswap->rpart[i][j].es;
      pos = appvacancy->set_it_e_set(x, y, z, it, ev, er, es, flag);
      if ((cur_setswap->rpart[i][j].change_type == 1) && flag && (it == 0))
      {
        int mysector = appvacancy->in_which_set(pos);
        cvector_push_back(listVacanciesUpdated[mysector], pos);
      }
    }
  }
  MPI_Waitall(n, request, status);
}

void CommLattice::create_set()
{
  setswap = new Swap[nset];
  int xleft, xright, yleft, yright, zleft, zright;

  xleft = domain->procneigh[0][0];
  xright = domain->procneigh[0][1];
  yleft = domain->procneigh[1][0];
  yright = domain->procneigh[1][1];
  zleft = domain->procneigh[2][0];
  zright = domain->procneigh[2][1];

  int xdir, ydir, zdir, xydir, xzdir, yzdir, xyzdir, xdelta, zdelta;

  // get sectors' neighbor processors, they are processors to send.
  for (int iset = 0; iset < nset; iset++)
  {
    if (iset % 2 == 0)
      xdir = xleft;
    else
      xdir = xright;
    if (iset < 4)
      zdir = zleft;
    else
      zdir = zright;
    if (iset == 0 || iset == 1 || iset == 4 || iset == 5)
      ydir = yleft;
    else
      ydir = yright;
    xdelta = xdir - me;
    xydir = ydir + xdelta;
    zdelta = zdir - me;
    xzdir = xdir + zdelta;
    yzdir = ydir + zdelta;
    xyzdir = xydir + zdelta;

    int tmp[7] = {xdir, ydir, zdir, xydir, xzdir, yzdir, xyzdir};
    std::map<int, bool> sproc;
    sproc[xdir] = true;
    sproc[ydir] = true;
    sproc[zdir] = true;
    sproc[xydir] = true;
    sproc[xzdir] = true;
    sproc[yzdir] = true;
    sproc[xyzdir] = true;

    setswap[iset].nsend = 0;
    setswap[iset].sproc = new int[sproc.size()];
    std::map<int, bool>::iterator it;
    for (it = sproc.begin(); it != sproc.end(); it++)
    {
      if (it->first != me)
      {
        setswap[iset].sproc[setswap[iset].nsend++] = it->first;
      }
    }
  }
  // the processors to send is the neighbors of the diagonal sector.
  for (int iset = 0; iset < nset; iset++)
  {
    int nsend = setswap[iset].nsend;
    int nrecv = setswap[7 - iset].nsend;
    setswap[iset].nrecv = nrecv;
    setswap[iset].rproc = new int[nrecv];
    for (int i = 0; i < nrecv; i++)
    {
      setswap[iset].rproc[i] = setswap[7 - iset].sproc[i];
    }

    setswap[iset].spart = new Particle *[nsend];
    setswap[iset].rpart = new Particle *[nrecv];
    setswap[iset].scount = new int[nsend];
    setswap[iset].smax = new int[nsend];
    setswap[iset].rcount = new int[nrecv];
    setswap[iset].request = new MPI_Request[nsend];
    setswap[iset].status = new MPI_Status[nrecv];

    for (int j = 0; j < nsend; j++)
    {
      setswap[iset].shash.insert(std::pair<int, int>(setswap[iset].sproc[j], j));
      setswap[iset].scount[j] = 0;
      setswap[iset].smax[j] = PARTNUM;
      setswap[iset].spart[j] = (Particle *)memory->smalloc(PARTNUM * sizeof(Particle), "app:particle");
    }
    for (int j = 0; j < nrecv; j++)
    {
      setswap[iset].rhash.insert(std::pair<int, int>(setswap[iset].rproc[j], j));
      setswap[iset].rcount[j] = 0;
      setswap[iset].rpart[j] = new Particle[PARTNUM];
    }
  }
}

void CommLattice::set_iset(int i)
{
  cur_setswap = &setswap[i];
}
void CommLattice::add_to_send(int site, int x, int y, int z, int change_type)
{
  int i, p, index;
  int *proc;
  std::map<int, int>::iterator loc;
  Particle particle = {x, y, z, iarray[0][site], change_type, darray[0][site], darray[1][site], darray[2][site]};
  if (needer.find(site) != needer.end())
  {
    proc = needer[site];
    for (i = 0; i < nset - 1 && proc[i] != -1; i++)
    {
      p = proc[i];
      loc = cur_setswap->shash.find(p);
      if (loc == cur_setswap->shash.end())
      {
        error->all(FLERR, "cur_setswap error!");
      }
      index = loc->second;
      if (cur_setswap->scount[index] == cur_setswap->smax[index])
      {
        cur_setswap->smax[index] += PARTNUM;
        cur_setswap->spart[index] = (Particle *)memory->srealloc(cur_setswap->spart[index], cur_setswap->smax[index] * sizeof(Particle), "app:particle");
      }
      cur_setswap->spart[index][cur_setswap->scount[index]] = particle;
      cur_setswap->scount[index]++;
    }
  }
}

void CommLattice::perform_set()
{
  int flag, pos;
  AppVacancy *appvacancy = (AppVacancy *)app;

  int i, flag_Iprobe; // wubd 2016-12-29

  MPI_Status status[cur_setswap->nsend];
  for (i = 0; i < cur_setswap->nsend; i++)
  {
    MPI_Isend(cur_setswap->spart[i], cur_setswap->scount[i], ctype, cur_setswap->sproc[i], 0, world, &cur_setswap->request[i]);
    cur_setswap->scount[i] = 0;
  }
  for (i = 0; i < cur_setswap->nrecv; i++)
  {
    MPI_Probe(cur_setswap->rproc[i], 0, world, &cur_setswap->status[i]);
    MPI_Get_count(&cur_setswap->status[i], ctype, &cur_setswap->rcount[i]);
    delete[] cur_setswap->rpart[i];
    cur_setswap->rpart[i] = new Particle[cur_setswap->rcount[i]];
    MPI_Recv(cur_setswap->rpart[i], cur_setswap->rcount[i], ctype, cur_setswap->rproc[i], 0, world, &cur_setswap->status[i]);
    for (int j = 0; j < cur_setswap->rcount[i]; j++)
    {
      int x = cur_setswap->rpart[i][j].x;
      int y = cur_setswap->rpart[i][j].y;
      int z = cur_setswap->rpart[i][j].z;
      int it = cur_setswap->rpart[i][j].it;
      double ev = cur_setswap->rpart[i][j].ev;
      double er = cur_setswap->rpart[i][j].er;
      double es = cur_setswap->rpart[i][j].es;
      pos = appvacancy->set_it_e_set(x, y, z, it, ev, er, es, flag);

      if ((cur_setswap->rpart[i][j].change_type == 1) && flag && (it == 0))
      {
        int mysector = appvacancy->in_which_set(pos);
        cvector_push_back(listVacanciesUpdated[mysector], pos);
      }
    }
  }
  MPI_Waitall(cur_setswap->nsend, cur_setswap->request, status);
}

void CommLattice::new_type_particle(MPI_Datatype *ctype)
{
  int blockcounts[2];
  MPI_Datatype oldtypes[2];
  MPI_Aint offsets[2];

  blockcounts[0] = 5;
  blockcounts[1] = 3;

  offsets[0] = 0;
  offsets[1] = sizeof(int) * 5;

  oldtypes[0] = MPI_INT;
  oldtypes[1] = MPI_DOUBLE;

  MPI_Type_struct(2, blockcounts, offsets, oldtypes, ctype);
  MPI_Type_commit(ctype);
}

/* ----------------------------------------------------------------------
   free all stl map data from memory wubd 2016-12-29
------------------------------------------------------------------------- */

void CommLattice::destroy_map(std::map<int, int *> *doclean)
{
  std::map<int, int *>::iterator iter;
  for (iter = (*doclean).begin(); iter != (*doclean).end(); iter++)
  {
    delete iter->second;
  }
  (*doclean).clear();
}

void CommLattice::free_swap(Swap *swap)
{
  delete[] swap->sproc;
  delete[] swap->scount;
  delete[] swap->smax;
  swap->shash.clear();
  for (int j = 0; j < swap->nsend; j++)
    memory->destroy(swap->spart[j]);
  delete[] swap->spart;

  delete[] swap->rproc;
  delete[] swap->rcount;
  swap->rhash.clear();
  for (int j = 0; j < swap->nrecv; j++)
    memory->destroy(swap->rpart[j]);
  delete[] swap->rpart;

  delete[] swap->request;
  delete[] swap->status;
}

void CommLattice::set_PriorityTargetList(int *a)
{
  for (int i = 0; i < 24; i++)
  {
    PriorityTargetList[i] = a[i];
  }
}

void CommLattice::set_PrioritySourceList(int *a)
{
  for (int i = 0; i < 24; i++)
  {
    PrioritySourceList[i] = a[i];
  }
}

void CommLattice::get_Priority_target(int mysectorID, int *n, int **ptr)
{
  int flag = mysectorID % 2;
  (*ptr) = PriorityTargetList + mysectorID % 2 + (mysectorID / 2) * 4 + (mysectorID / 4) * 4;
  *n = 1;
  if (flag)
  {
    *n = 3;
  }
  flag = mysectorID % 4;
  if (flag == 3)
    *n = 7;
}
void CommLattice::get_Priority_source(int mysectorID, int *n, int **ptr)
{
  int flag = mysectorID % 2;
  (*ptr) = PrioritySourceList + mysectorID % 2 + (mysectorID / 2) * 4 + (mysectorID / 4) * 4;
  *n = 1;
  if (flag)
  {
    *n = 3;
  }
  flag = mysectorID % 4;
  if (flag == 3)
    *n = 7;
}

/* ----------------------------------------------------------------------
   create send portion of a Swap communication pattern
   create from list of sites I send and procs who own them
------------------------------------------------------------------------- */
// deleted by bdwu 170301
// void CommLattice::create_send_from_list(int nsite, Site *buf, Swap *swap)

/* ----------------------------------------------------------------------
   create recv portion of a Swap communication pattern
   start with list of sites I will send
   circulate list to all procs
------------------------------------------------------------------------- */
// deleted by bdwu 170301
//  void CommLattice::create_recv_from_send(int nsite, int maxsite,
//           Site *buf, Swap *swap)

/* ----------------------------------------------------------------------
   communicate site values via Swap instructions
   use site array = iarray[0] as source/destination
------------------------------------------------------------------------- */
// deleted by bdwu 170228
// void CommLattice::perform_swap_site(Swap *swap)

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use iarray as source/destination, all integer data
------------------------------------------------------------------------- */
// deleted by bdwu 170228
//  void CommLattice::perform_swap_int(Swap *swap)

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use darray as source/destination, all double data
------------------------------------------------------------------------- */
// deleted by bdwu 170228
//  void CommLattice::perform_swap_double(Swap *swap)