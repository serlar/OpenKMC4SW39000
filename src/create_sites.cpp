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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_sites.h"
#include "app.h"
#include "app_lattice.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "potential.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "timer.h"
#include "error.h"

#include <map>
#include <time.h>

using namespace SPPARKS_NS;

// same as in lattice.cpp

enum
{
  NONE,
  LINE_2N,
  SQ_4N,
  SQ_8N,
  TRI,
  SC_6N,
  SC_26N,
  FCC,
  BCC,
  DIAMOND,
  FCC_OCTA_TETRA,
  RANDOM_1D,
  RANDOM_2D,
  RANDOM_3D
};

enum
{
  BOX,
  REGION
};
enum
{
  DUMMY,
  IARRAY,
  DARRAY
};

#define DELTALOCAL 10000
#define DELTABUF 10000
#define EPSILON 0.0001

/* ---------------------------------------------------------------------- */

CreateSites::CreateSites(SPPARKS *spk) : Pointers(spk) {}

/* ---------------------------------------------------------------------- */

void CreateSites::command(int narg, char **arg)
{
  if (app == NULL)
    error->all(FLERR, "Create_sites command before app_style set");
  if (domain->box_exist == 0)
    error->all(FLERR, "Create_sites command before simulation box is defined");
  if (app->sites_exist == 1)
    error->all(FLERR, "Cannot create sites after sites already exist");
  if (domain->lattice == NULL)
    error->all(FLERR, "Cannot create sites with undefined lattice");

  if (narg < 1)
    error->all(FLERR, "Illegal create_sites command");

  int iarg;
  if (strcmp(arg[0], "box") == 0)
  {
    style = BOX;
    iarg = 1;
  }
  else if (strcmp(arg[0], "region") == 0)
  {
    style = REGION;
    if (narg < 2)
      error->all(FLERR, "Illegal create_sites command");
    nregion = domain->find_region(arg[1]);
    if (nregion == -1)
      error->all(FLERR, "Create_sites region ID does not exist");
    iarg = 2;
  }
  else
    error->all(FLERR, "Illegal create_sites command");

  // parse optional args

  valueflag = DUMMY;
  nbasis = domain->lattice->nbasis;
  basisflag = new int[nbasis + 1];
  basis_ivalue = new int[nbasis + 1];
  basis_dvalue = new double[nbasis + 1];
  for (int i = 1; i <= nbasis; i++)
    basisflag[i] = 0;

  while (iarg < narg)
  {
    if (strcmp(arg[iarg], "value") == 0)
    {
      if (iarg + 3 > narg)
        error->all(FLERR, "Illegal create_sites command");
      valueflag = 1;
      if (strcmp(arg[iarg + 1], "site") == 0)
      {
        valueflag = IARRAY;
        valueindex = 0;
        if (app->iarray == NULL)
          error->all(FLERR, "Creating a quantity application does not support");
      }
      else if (arg[iarg + 1][0] == 'i')
      {
        valueflag = IARRAY;
        valueindex = atoi(&arg[iarg + 1][1]);
        if (valueindex < 1 || valueindex > app->ninteger)
          error->all(FLERR, "Creating a quantity application does not support");
        valueindex--;
      }
      else if (arg[iarg + 1][0] == 'd')
      {
        valueflag = DARRAY;
        valueindex = atoi(&arg[iarg + 1][1]);
        if (valueindex < 1 || valueindex > app->ndouble)
          error->all(FLERR, "Creating a quantity application does not support");
        valueindex--;
      }
      if (valueflag == IARRAY)
        ivalue = atoi(arg[iarg + 2]);
      else
        dvalue = atof(arg[iarg + 2]);
      iarg += 3;
    }
    else if (strcmp(arg[iarg], "basis") == 0)
    {
      if (iarg + 3 > narg)
        error->all(FLERR, "Illegal create_sites command");
      if (valueflag == DUMMY)
        error->all(FLERR, "Must use value option before basis option "
                          "in create_sites command");
      int ilo, ihi;
      if (nbasis == 0)
        error->all(FLERR, "Cannot use create_sites basis with random lattice");
      potential->bounds(arg[iarg + 1], nbasis, ilo, ihi);
      int count = 0;
      for (int i = ilo; i <= ihi; i++)
      {
        basisflag[i] = 1;
        if (valueflag == IARRAY)
          basis_ivalue[i] = atoi(arg[iarg + 2]);
        else if (valueflag == DARRAY)
          basis_dvalue[i] = atof(arg[iarg + 2]);
        count++;
      }
      if (count == 0)
        error->all(FLERR, "Illegal create_sites command");
      iarg += 3;
    }
    else
      error->all(FLERR, "Illegal create_sites command");
  }

  // create sites, either on-lattice or off-lattice

  if (domain->me == 0)
  {
    if (screen)
      fprintf(screen, "Creating sites ...\n");
    if (logfile)
      fprintf(logfile, "Creating sites ...\n");
  }

  app->sites_exist = 1;

  if (app->appclass == App::LATTICE)
  {
    applattice = (AppLattice *)app;
    latticeflag = 1;
  } /*else if (app->appclass == App::OFF_LATTICE) { //delete by bdwu 0703
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  }*/

  int dimension = domain->dimension;
  latstyle = domain->lattice->style;

  if (latstyle == BCC)
  {

    xlattice = domain->lattice->xlattice;
    ylattice = domain->lattice->ylattice;
    zlattice = domain->lattice->zlattice;

    time_t timer;
    struct tm *tblock;
    timer = time(NULL);
    tblock = localtime(&timer);

    if (domain->me == 0)
    {
      if (screen)
        fprintf(screen, "Begining structured_lattice ...%s", asctime(tblock));
      if (logfile)
        fprintf(logfile, "Begining structured_lattice ...%s", asctime(tblock));
    }
    structured_lattice();
  }

  time_t timer;
  struct tm *tblock;
  timer = time(NULL);
  tblock = localtime(&timer);
  if (domain->me == 0)
  {
    if (screen)
      fprintf(screen, "Begining transport ...%s", asctime(tblock));
    if (logfile)
      fprintf(logfile, "Begining transport ...%s", asctime(tblock));
  }
  applattice->transport();
  timer = time(NULL);
  tblock = localtime(&timer);
  if (domain->me == 0)
  {
    if (screen)
      fprintf(screen, "End transport ...%s", asctime(tblock));
    if (logfile)
      fprintf(logfile, "End transport ...%s", asctime(tblock));
  }

  delete[] basisflag;
  delete[] basis_ivalue;
  delete[] basis_dvalue;
}

/* ----------------------------------------------------------------------
   generate sites on structured lattice that fits in simulation box
   loop over entire lattice
   if style = REGION, require a site be in region as well
   each proc keeps sites in its sub-domain
 ------------------------------------------------------------------------- */

void CreateSites::structured_lattice()
{
  int dimension = domain->dimension;
  int nonperiodic = domain->nonperiodic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

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

  double **basis = domain->lattice->basis;
  int **iarray = app->iarray;
  double **darray = app->darray;

  // in periodic dims:
  // check that simulation box is integer multiple of lattice spacing

  nx = static_cast<int>(domain->xprd / xlattice);
  if (dimension >= 2)
    ny = static_cast<int>(domain->yprd / ylattice);
  else
    ny = 1;
  if (dimension == 3)
    nz = static_cast<int>(domain->zprd / zlattice);
  else
    nz = 1;

  if (xperiodic &&
      fabs(nx * xlattice - domain->xprd) > EPSILON)
    error->all(FLERR, "Periodic box is not a multiple of lattice spacing");
  if (dimension > 1 && yperiodic &&
      fabs(ny * ylattice - domain->yprd) > EPSILON)
    error->all(FLERR, "Periodic box is not a multiple of lattice spacing");
  if (dimension > 2 && zperiodic &&
      fabs(nz * zlattice - domain->zprd) > EPSILON)
    error->all(FLERR, "Periodic box is not a multiple of lattice spacing");

  // set domain->nx,ny,nz iff style = BOX and system is fully periodic
  // else site IDs may be non-contiguous and/or ordered irregularly

  if (style == BOX && nonperiodic == 0)
  {
    domain->nx = nx;
    domain->ny = ny;
    domain->nz = nz;
  }

  // if dim is periodic:
  //   lattice origin = lower box boundary
  //   loop bounds = 0 to N-1
  // if dim is non-periodic:
  //   lattice origin = 0.0
  //   loop bounds = enough to tile box completely, with all basis atoms

  if (xperiodic)
  {
    xorig = boxxlo;
    xlo = 0;
    xhi = nx - 1;
  }
  else
  {
    xorig = 0.0;
    xlo = static_cast<int>(boxxlo / xlattice);
    while ((xlo + 1) * xlattice > boxxlo)
      xlo--;
    xlo++;
    xhi = static_cast<int>(boxxhi / xlattice);
    while (xhi * xlattice <= boxxhi)
      xhi++;
    xhi--;
  }

  if (yperiodic)
  {
    yorig = boxylo;
    ylo = 0;
    yhi = ny - 1;
  }
  else
  {
    yorig = 0.0;
    ylo = static_cast<int>(boxylo / ylattice);
    while ((ylo + 1) * ylattice > boxylo)
      ylo--;
    ylo++;
    yhi = static_cast<int>(boxyhi / ylattice);
    while (yhi * ylattice <= boxyhi)
      yhi++;
    yhi--;
  }

  if (zperiodic)
  {
    zorig = boxzlo;
    zlo = 0;
    zhi = nz - 1;
  }
  else
  {
    zorig = 0.0;
    zlo = static_cast<int>(boxzlo / zlattice);
    while ((zlo + 1) * zlattice > boxzlo)
      zlo--;
    zlo++;
    zhi = static_cast<int>(boxzhi / zlattice);
    while (zhi * zlattice <= boxzhi)
      zhi++;
    zhi--;
  }

  // generate xyz coords and store them with site ID
  // tile the simulation box from origin, respecting PBC
  // site IDs should be contiguous if style = BOX and fully periodic
  // for non-periodic dims, check if site is within global box
  // for style = REGION, check if site is within region
  // if non-periodic or style = REGION, IDs may not be contiguous

  int i, j, k, m, nlocal;
  double x, y, z;

  int maxlocal = 0;
  siteijk = NULL;
  double **ghosts = NULL;
  tagint *ghostids = NULL;

  tagint n = 0;

  tagint localcount = 0;
  tagint ghostcount = 0;

  if (domain->me == 0)
  {
    fprintf(screen, " structured_lattice: before 4-level-loop(mijk)\n");
    fprintf(logfile, " structured_lattice: before 4-level-loop(mijk)\n");
    timer->current();
  }
  // replace line 405-459 with line 323-402
  int xxlo = floor(subxlo / xlattice - 2.0);
  int xxhi = ceil(subxhi / xlattice + 2.0);
  int yylo = floor(subylo / ylattice - 2.0);
  int yyhi = ceil(subyhi / ylattice + 2.0);
  int zzlo = floor(subzlo / zlattice - 2.0);
  int zzhi = ceil(subzhi / zlattice + 2.0);
  if (subxlo == boxxlo && subxhi == boxxhi)
  {
    xxlo = 0;
    xxhi = xhi;
  }
  if (subylo == boxylo && subyhi == boxyhi)
  {
    yylo = 0;
    yyhi = yhi;
  }
  if (subzlo == boxzlo && subzhi == boxzhi)
  {
    zzlo = 0;
    zzhi = zhi;
  }

  int ii, jj, kk;

  for (k = zzlo; k <= zzhi; k++)
    for (j = yylo; j <= yyhi; j++)
      for (i = xxlo; i <= xxhi; i++)
        for (m = 0; m < nbasis; m++)
        {
          ii = i;
          jj = j;
          kk = k;
          if (kk < 0)
            kk += nz;
          if (kk > zhi)
            kk -= nz;
          if (jj < 0)
            jj += ny;
          if (jj > yhi)
            jj -= ny;
          if (ii < 0)
            ii += nx;
          if (ii > xhi)
            ii -= nx;

          tagint gid = (kk - zlo) * (yhi - ylo + 1) * (xhi - xlo + 1) * nbasis + (jj - ylo) * (xhi - xlo + 1) * nbasis + (ii - xlo) * nbasis + m + 1;

          x = (ii + basis[m][0]) * xlattice + xorig;
          y = (jj + basis[m][1]) * ylattice + yorig;
          z = (kk + basis[m][2]) * zlattice + zorig;

          /** add by smz 170106 ==BEGIN==**/
          if ((x >= subxlo - 2.01 * xlattice && x < subxhi + 1.99 * xlattice || subxlo == boxxlo && x >= boxxhi - 2.01 * xlattice || subxhi == boxxhi && x < boxxlo + 1.99 * xlattice) &&
              (y >= subylo - 2.01 * ylattice && y < subyhi + 1.99 * ylattice || subylo == boxylo && y >= boxyhi - 2.01 * ylattice || subyhi == boxyhi && y < boxylo + 1.99 * ylattice) &&
              (z >= subzlo - 2.01 * zlattice && z < subzhi + 1.99 * zlattice || subzlo == boxzlo && z >= boxzhi - 2.01 * zlattice || subzhi == boxzhi && z < boxzlo + 1.99 * zlattice))
          {
            if (x >= subxlo && x < subxhi &&
                y >= subylo && y < subyhi &&
                z >= subzlo && z < subzhi)
            {
              applattice->add_site(gid, x, y, z);
              nlocal = app->nlocal;
              if (valueflag == IARRAY)
              {
                if (basisflag[m + 1])
                  iarray[valueindex][nlocal - 1] = basis_ivalue[m + 1];
                else
                  iarray[valueindex][nlocal - 1] = ivalue; // run
              }
              else if (valueflag == DARRAY)
              {
                if (basisflag[m + 1])
                  darray[valueindex][nlocal - 1] = basis_dvalue[m + 1];
                else
                  darray[valueindex][nlocal - 1] = dvalue;
              }
              localcount++;
            }
            else
            {
              ghostcount++;
              if (ghostcount > maxlocal)
              {
                maxlocal += DELTALOCAL;
                memory->grow(ghosts, maxlocal, 3, "create:ghosts");
                memory->grow(ghostids, maxlocal, "create:ghostids");
              }
              ghostids[ghostcount - 1] = gid;
              ghosts[ghostcount - 1][0] = x;
              ghosts[ghostcount - 1][1] = y;
              ghosts[ghostcount - 1][2] = z;
            }
          }
        }

  if (domain->me == 0)
  {
    fprintf(screen, " structured_lattice: before add_ghost\n");
    fprintf(logfile, " structured_lattice: before add_ghost\n");
    timer->current();
  }
  for (n = 0; n < ghostcount; n++)
  {
    applattice->add_ghost(ghostids[n], ghosts[n][0], ghosts[n][1], ghosts[n][2]);
  }
  memory->destroy(ghosts);
  memory->destroy(ghostids);
  /** add by smz 170106 ==END==**/
  if (domain->me == 0)
  {
    fprintf(screen, " structured_lattice: before MPI_Allreduce\n");
    fprintf(logfile, " structured_lattice: before MPI_Allreduce\n");
    timer->current();
  }
  tagint nbig = app->nlocal;
  app->nglobal = nbasis * nx * ny * nz;
  // MPI_Allreduce(&nbig,&app->nglobal,1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (domain->me == 0)
  {
    if (screen)
      fprintf(screen, "  " TAGINT_FORMAT " sites\n", app->nglobal);
    if (logfile)
      fprintf(logfile, "  " TAGINT_FORMAT " sites\n", app->nglobal);
  }

  // for style = BOX and periodic system, check if nglobal is correct

  if (style == BOX && domain->nonperiodic == 0)
  {
    nbig = nbasis;
    nbig = nbig * nx * ny * nz;
    if (style == BOX && app->nglobal != nbig)
      error->all(FLERR, "Did not create correct number of sites");
  }
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets(double **basis)
{
  // deleted by bdwu 170301
  // if (latstyle == BCC)
  for (int m = 0; m < nbasis; m++)
    offsets_3d(m, basis, sqrt(3.0) / 2.0 * xlattice, sqrt(3.0) / 2.0 * xlattice,
               maxneigh, cmap[m]);
}

/* ---------------------------------------------------------------------- */
// offsets_2dd elete by bdwu 160703

/* ---------------------------------------------------------------------- */

void CreateSites::offsets_3d(int ibasis, double **basis,
                             double cutlo, double cuthi,
                             int ntarget, int **cmapone)
{
  int i, j, k, m, n;
  double x0, y0, z0, delx, dely, delz, r;

  n = 0;
  x0 = basis[ibasis][0] * xlattice;
  y0 = basis[ibasis][1] * ylattice;
  z0 = basis[ibasis][2] * zlattice;
  for (i = -1; i <= 1; i++)
  {
    for (j = -1; j <= 1; j++)
    {
      for (k = -1; k <= 1; k++)
      {
        for (m = 0; m < nbasis; m++)
        {
          delx = (i + basis[m][0]) * xlattice - x0;
          dely = (j + basis[m][1]) * ylattice - y0;
          delz = (k + basis[m][2]) * zlattice - z0;
          r = sqrt(delx * delx + dely * dely + delz * delz);
          if (r > cutlo - EPSILON && r < cuthi + EPSILON)
          {
            if (n == ntarget)
              error->all(FLERR, "Incorrect lattice neighbor count");
            cmapone[n][0] = i;
            cmapone[n][1] = j;
            cmapone[n][2] = k;
            cmapone[n][3] = m;
            n++;
          }
        }
      }
    }
  }
  if (n != ntarget)
    error->all(FLERR, "Incorrect lattice neighbor count");
}

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
   respect non-periodic boundaries
   only called for on-lattice models
 ------------------------------------------------------------------------- */
// deleted by bdwu 170301
// void CreateSites::structured_connectivity()

/* ----------------------------------------------------------------------
   generate random sites
   each proc keeps those in its sub-domain
 ------------------------------------------------------------------------- */
// delete by bdwu 160703

/* ----------------------------------------------------------------------
   infer connectivity from neighbors within cutoff distance
   respect non-periodic boundaries
   only called for on-lattice models
 ------------------------------------------------------------------------- */
// delete by bdwu 160703

/* ----------------------------------------------------------------------
   set maxneigh and initialize idneigh when lattice created via read_sites
   called from read_sites when it reads in sites and neighbors
 ------------------------------------------------------------------------- */

void CreateSites::read_sites(AppLattice *apl)
{ /* delete by smz 161017
   int i,j;

   maxneigh = apl->maxneigh;
   memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

   int *numneigh = apl->numneigh;
   int **neighbor = apl->neighbor;
   int nlocal = app->nlocal;

   for (i = 0; i < nlocal; i++)
     for (j = 0; j < numneigh[i]; j++)
       idneigh[i][j] = neighbor[i][j];
       */
}

/* ----------------------------------------------------------------------
   create ghosts sites around local sub-domain
   only called for on-lattice models
   pass applattice as pointer so can call from ReadSites
   numneigh and global neighbor IDs of each owned site are known as input
   add ghost sites for delpropensity layers
   form neigh list for each layer of ghost sites one layer at a time
   when done, delpropensity-1 layers have a full numneigh and neigh list
     last delpropensity layers has a partial numneigh and neigh list
   convert neighbor IDs from global indices to local indices
 ------------------------------------------------------------------------- */
// deleted by bdwu 170301
// void CreateSites::ghosts_from_connectivity(AppLattice *apl, int delpropensity)
