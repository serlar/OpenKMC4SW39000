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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dump.h"
#include "app.h"
#include "app_lattice.h"
 
#include "domain.h"
 
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

Dump *Dump::dumpptr;

#define BIG 1.0e20
#define EPSILON 1.0e-6

enum{ASCEND,DESCEND};

/* ---------------------------------------------------------------------- */

Dump::Dump(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // parse dump args

  if (narg < 4) error->all(FLERR,"Illegal dump command");

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  delta = atof(arg[2]);
  if (delta <= 0.0) error->all(FLERR,"Illegal dump command");

  n = strlen(arg[3]) + 1;
  filename = new char[n];
  strcpy(filename,arg[3]);

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;
  multiname = NULL;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    multiname = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",filename,me,ptr+1);
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;

  // on-lattice or off-lattice app

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  } else
    error->all(FLERR,"Dump command can only be used for spatial applications");

  // dump params

  boxxlo = domain->boxxlo;
  boxxhi = domain->boxxhi;
  boxylo = domain->boxylo;
  boxyhi = domain->boxyhi;
  boxzlo = domain->boxzlo;
  boxzhi = domain->boxzhi;

  logfreq = 0;
  delay = 0.0;
  flush_flag = 1;
  padflag = 0;
  sort_flag = 0;

  maxbuf = maxids = maxsort = maxproc = 0;
  buf = bufsort = NULL;
  ids = idsort = NULL;
  index = proclist = NULL;
  // irregular = NULL;

  firstflag = 1;
  idump = 0;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;
  delete [] multiname;

  memory->destroy(buf);
  memory->destroy(bufsort);
  memory->destroy(ids);
  memory->destroy(idsort);
  memory->destroy(index);
  memory->destroy(proclist);
  // delete irregular;

  if (multiproc) MPI_Comm_free(&clustercomm);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();

  if (!sort_flag) {
    memory->destroy(bufsort);
    memory->destroy(ids);
    memory->destroy(idsort);
    memory->destroy(index);
    memory->destroy(proclist);
    // delete irregular;

    maxids = maxsort = maxproc = 0;
    bufsort = NULL;
    ids = idsort = NULL;
    index = proclist = NULL;
    // irregular = NULL;
  }

  if (sort_flag) { 
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of site values
------------------------------------------------------------------------- */

void Dump::write(double time)
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // nmine = # of dump lines this proc will contribute to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_SPK_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm


  // write timestep header

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_SPK_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader,time);

  firstflag = 0;
  idump++;

  // grow communication buffer if necessary

  if (nmax > maxbuf) {
    maxbuf = nmax;
    memory->sfree(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // insure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(NULL);
  if (sort_flag) sort();

  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines,nchars;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= size_one;
      } else nlines = nme;
      
      write_data(nlines,buf);
    }
    if (flush_flag) fflush(fp);
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
  }

  // write timestep footer, if needed

  write_footer();

  // if file per timestep, close file if I am filewriter

  if (multifile) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   derived classes may override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,idump,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,idump,ptr+1);
    }
    *ptr = '*';
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (compressed) {
#ifdef SPPARKS_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
#ifdef _WIN32
      fp = _popen(gzip,"wb");
#else
      fp = popen(gzip,"w");
#endif
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    //} else if (append_flag) {
    // fp = fopen(filecurrent,"a");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ----------------------------------------------------------------------
   parallel sort of buf across all procs
   changes nme, reorders datums in buf, grows buf if necessary
------------------------------------------------------------------------- */

void Dump::sort()
{
  
}

/* ----------------------------------------------------------------------
   compare two site IDs
   called via qsort() in sort() method
   is a static method so access data via dumpptr
------------------------------------------------------------------------- */

int Dump::idcompare(const void *pi, const void *pj)
{
  tagint *idsort = dumpptr->idsort;

  int i = *((int *) pi);
  int j = *((int *) pj);

  if (idsort[i] < idsort[j]) return -1;
  if (idsort[i] > idsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via dumpptr
   sort in ASCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare(const void *pi, const void *pj)
{
  double *bufsort = dumpptr->bufsort;
  int size_one = dumpptr->size_one;
  int sortcolm1 = dumpptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] < bufsort[j]) return -1;
  if (bufsort[i] > bufsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via dumpptr
   sort in DESCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare_reverse(const void *pi, const void *pj)
{
  double *bufsort = dumpptr->bufsort;
  int size_one = dumpptr->size_one;
  int sortcolm1 = dumpptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] > bufsort[j]) return -1;
  if (bufsort[i] < bufsort[j]) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      delay = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"delta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      delta = atof(arg[iarg+1]);
      if (delta <= 0.0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
	error->all(FLERR,"Cannot use dump_modify fileper "
                   "without % in dump file name");
      int nper = atoi(arg[iarg+1]);
      if (nper <= 0) error->all(FLERR,"Illegal dump_modify command");
      
      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      int icluster = fileproc/nper;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      idump = atoi(arg[iarg+1]);
      if (idump < 0.0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (nrepeat < 0) error->all(FLERR,"Illegal dump_modify command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"loglinfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (nrepeat < 0) error->all(FLERR,"Illegal dump_modify command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 2;
      iarg += 3;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
	error->all(FLERR,"Cannot use dump_modify nfile "
                   "without % in dump file name");
      int nfile = atoi(arg[iarg+1]);
      if (nfile <= 0) error->all(FLERR,"Illegal dump_modify command");
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext = 
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      padflag = atoi(arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"off") == 0) sort_flag = 0;
      else if (strcmp(arg[iarg+1],"id") == 0) {
        sort_flag = 1;
        sortcol = 0;
        sortorder = ASCEND;
      } else {
        sort_flag = 1;
        sortcol = atoi(arg[iarg+1]);
        sortorder = ASCEND;
        if (sortcol == 0) error->all(FLERR,"Illegal dump_modify command");
        if (sortcol < 0) {
          sortorder = DESCEND;
          sortcol = -sortcol;
        }
        sortcolm1 = sortcol - 1;
      }
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += n;
    }
  }
}
