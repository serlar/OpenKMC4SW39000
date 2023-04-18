/* ----------------------------------------------------------------------
   Vacancy tries to build a model for vacancy diffusion. The calculation are based on the
   diffusion model and LAKIMOCA.
   By SMZ, 2015-12.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_vacancy.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "timer.h"

#define CLOCKRATE 2.0E12

using namespace SPPARKS_NS;

enum
{
	VACANT,
	Fe,
	Cu,
	Si,
	Mn,
	Ni,
	TOP
};
enum
{
	DEPOSITION,
	NNHOP,
	SCHWOEBEL
};
enum
{
	NOSWEEP,
	RANDOM,
	RASTER,
	COLOR,
	COLOR_STRICT
}; // from app_lattice.cpp
enum
{
	ZERO,
	CM_EAM,
	CM_FS,
	CM_EAM3,
	CM_PPAIR2
};

#define DELTAEVENT 100000
#define kb (1.38066E-23 / 1.6E-19)
#define nu 6.E12

static inline unsigned long rpcc()
{
	unsigned long time;
	asm("rtc %0"
		: "=r"(time)
		:);
	return time;
}
/* ---------------------------------------------------------------------- */

AppVacancy::AppVacancy(SPPARKS *spk, int narg, char **arg)
	: AppLattice(spk, narg, arg)
{
	// these can be changed by model choice, see below

	ninteger = 1;
	ndouble = 3;
	delpropensity = 3;
	delevent = 3;
	allow_kmc = 1;
	allow_rejection = 1;
	allow_masking = 0;
	numrandom = 1;
	a = 0;
	tag_COHESIVE_MODEL = ZERO;

	create_arrays();

	maxevent = 0;
	events = NULL;
	firstevent = NULL;

	barrier = NULL;

	allocated = 0;
}

/* ---------------------------------------------------------------------- */

AppVacancy::~AppVacancy()
{
	int i, j;

	memory->sfree(events);
	memory->destroy(firstevent);
	memory->destroy(barrier);
	wait_send.clear();
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppVacancy::input_app(char *command, int narg, char **arg)
{
	if (sites_exist == 0)
	{
		char str[128];
		sprintf(str, "Cannot use %s command until sites exist", command);
		error->all(FLERR, str);
	}

	if (!allocated)
		allocate_data();
	allocated = 1;
	if (strcmp(command, "nbody") == 0)
	{
		a = domain->lattice->getLatconst();
		if (strcmp(arg[0], "eam") == 0)
		{
			tag_COHESIVE_MODEL = CM_EAM;
			read_pot_file(arg[1]);
		}
		else if (strcmp(arg[0], "fs") == 0)
		{
			tag_COHESIVE_MODEL = CM_FS;
			read_pot_file(arg[1]);
		}
		else if (strcmp(arg[0], "eam3") == 0)
		{
			tag_COHESIVE_MODEL = CM_EAM3;
			read_pot_file(arg[1]);
		}
		else if (strcmp(arg[0], "ppair2") == 0)
		{
			tag_COHESIVE_MODEL = CM_PPAIR2;
			init_potential_pair_file(arg[1]);
		}
		else
			error->all(FLERR, "Wrong cohesive model specified\nsupported model are: eam, fs, eam3\n");
	}
	else if (strcmp(command, "barrier") == 0)
	{
		if (narg != 2)
			error->all(FLERR, "Illegal barrier command");
		int type = -1;
		if (strcmp(arg[0], "Fe") == 0)
			type = Fe;
		else if (strcmp(arg[0], "Ni") == 0)
			type = Ni;
		else if (strcmp(arg[0], "Si") == 0)
			type = Si;
		else if (strcmp(arg[0], "Mn") == 0)
			type = Mn;
		else if (strcmp(arg[0], "Cu") == 0)
			type = Cu;
		else
			error->all(FLERR, "Illegal barrier command");
		barrier[type] = atof(arg[1]);
		if (barrier[type] < 0.1)
		{
			printf("Warning: very low value of barrier\n");
			printf("Warning: I hope you know what you are doing\n");
		}
	}
	else
		error->all(FLERR, "Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppVacancy::grow_app()
{
	lattice = iarray[0];
	E_V = darray[0];
	E_R = darray[1];
	E_S = darray[2];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppVacancy::init_app()
{
	if (!allocated)
		allocate_data();
	allocated = 1;

	if (maxneigh == 0)
		dt_sweep = 2655216.00; // deleted by bdwu 170301
	else
		dt_sweep = 1.0 / maxneigh; // deleted by bdwu 170301

	// site validity

	int flag = 0;
	for (int i = 0; i < nlocal; i++)
		if (lattice[i] < VACANT || lattice[i] >= TOP)
			flag = 1;
	int flagall;
	MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
	if (flagall)
		error->all(FLERR, "One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppVacancy::setup_app()
{
	// clear event list
	nevents = 0;
	for (int i = 0; i < nlocal; i++)
		firstevent[i] = -1;
	for (int i = 0; i < maxevent; i++)
		events[i].next = i + 1;
	freeevent = 0;

	// add:
	find_border();
	comm->perform_all();
	compute_energy();
}

void AppVacancy::transport()
{

	unitdis = 2 / domain->lattice->getLatconst();

	boxxlo = round(domain->boxxlo * unitdis);
	boxylo = round(domain->boxylo * unitdis);
	boxzlo = round(domain->boxzlo * unitdis);
	boxxhi = round(domain->boxxhi * unitdis);
	boxyhi = round(domain->boxyhi * unitdis);
	boxzhi = round(domain->boxzhi * unitdis);
	boxx = boxxhi - boxxlo;
	boxy = boxyhi - boxylo;
	boxz = boxzhi - boxzlo;

	subxlo = round(domain->subxlo * unitdis);
	subylo = round(domain->subylo * unitdis);
	subzlo = round(domain->subzlo * unitdis);
	subxhi = round(domain->subxhi * unitdis);
	subyhi = round(domain->subyhi * unitdis);
	subzhi = round(domain->subzhi * unitdis);
	subx = subxhi - subxlo;
	suby = subyhi - subylo;
	subz = subzhi - subzlo;
	subxlen = subx + 9;
	subylen = suby + 9;
	subzlen = subz + 9;

	if (subxlo == boxxlo && subxhi == boxxhi)
	{
		xghost = 0;
	}
	else
	{
		xghost = 4;
	}

	if (subylo == boxylo && subyhi == boxyhi)
	{
		yghost = 0;
	}
	else
	{
		yghost = 4;
	}

	if (subzlo == boxzlo && subzhi == boxzhi)
	{
		zghost = 0;
	}
	else
	{
		zghost = 4;
	}
}
void AppVacancy::get_neighbor(int initsite, int *neighbors)
{
	int i, j, k;
	int ip1, im1, jp1, jm1, kp1, km1;

	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);

	inc_px(i, 1, &ip1);
	inc_py(j, 1, &jp1);
	inc_pz(k, 1, &kp1);
	inc_px(i, -1, &im1);
	inc_py(j, -1, &jm1);
	inc_pz(k, -1, &km1);

	box2sub(&ip1, &jp1, &kp1);
	box2sub(&im1, &jm1, &km1);

	neighbors[0] = get_POS_ID(im1, jm1, km1);
	neighbors[1] = get_POS_ID(im1, jp1, km1);
	neighbors[2] = get_POS_ID(ip1, jm1, km1);
	neighbors[3] = get_POS_ID(ip1, jp1, km1);
	neighbors[4] = get_POS_ID(im1, jm1, kp1);
	neighbors[5] = get_POS_ID(im1, jp1, kp1);
	neighbors[6] = get_POS_ID(ip1, jm1, kp1);
	neighbors[7] = get_POS_ID(ip1, jp1, kp1);
}
void AppVacancy::find_border()
{

	int i, j, m;
	int nlayer = 1;
	int ntotal = nlocal + nghost;
	int nsites;
	int *site2i;
	// flag sites with -1 that are not in sector
	// flag sites with 0 that are in sector
	int *flag;
	int *neighbor = new int[8];
	memory->create(flag, ntotal, "app:flag");
	for (int isector = 0; isector < nset; isector++)
	{

		nsites = set[isector].nlocal;
		site2i = set[isector].site2i;

		for (i = 0; i < ntotal; i++)
			flag[i] = -1;
		for (m = 0; m < nsites; m++)
			flag[site2i[m]] = 0;

		// flag sector sites with -1 that have non-sector neighbor up to nlayer away
		for (int ilayer = 0; ilayer < nlayer; ilayer++)
		{
			for (m = 0; m < nsites; m++)
			{
				i = site2i[m];
				if (flag[i])
					continue;
				get_neighbor(i, neighbor);
				for (j = 0; j < 8; j++)
				{
					if (flag[neighbor[j]] < 0)
						break;
				}
				if (j < 8)
					flag[i] = 1;
			}
			for (m = 0; m < nsites; m++)
			{
				i = site2i[m];
				if (flag[i] > 0)
					flag[i] = -1;
			}
		}

		// nborder = # of border sites
		// allocate border and fill with site indices

		int nborder = 0;
		for (m = 0; m < nsites; m++)
		{
			i = site2i[m];
			if (flag[i] < 0)
				nborder++;
		}

		int *border;
		memory->create(border, nborder, "app:border");

		nborder = 0;
		for (m = 0; m < nsites; m++)
		{
			i = site2i[m];
			if (flag[i] < 0)
				border[nborder++] = i;
		}

		set[isector].border = border;
		set[isector].nborder = nborder;
	}
	delete[] neighbor;
	memory->destroy(flag);
}

void AppVacancy::box2sub(int *x, int *y, int *z)
{
	// int a=*x, b=*y, c=*z;
	if (*x < subxlo - 4)
		*x = *x + boxxhi - (subxlo - 4);
	else if (*x > subxhi + 4)
		*x = *x - boxxhi - (subxlo - 4);
	else
		*x = *x - (subxlo - 4);

	if (*y < subylo - 4)
		*y = *y + boxyhi - (subylo - 4);
	else if (*y > subyhi + 4)
		*y = *y - boxyhi - (subylo - 4);
	else
		*y = *y - (subylo - 4);

	if (*z < subzlo - 4)
		*z = *z + boxzhi - (subzlo - 4);
	else if (*z > subzhi + 4)
		*z = *z - boxzhi - (subzlo - 4);
	else
		*z = *z - (subzlo - 4);
}

void AppVacancy::count_vac(int nlocal, int *site2i,
						   cvector_vector_type(int) * vacant)
{
	int i, site;
	for (site = 0; site < nlocal; site++)
	{
		i = site2i[site];
		if (lattice[i] == 0)
		{
			cvector_push_back((*vacant), i);
		}
	}
}
/* ----------------------------------------------------------------------
   to compute EV, ER, ES using eatom
------------------------------------------------------------------------- */
void AppVacancy::compute_energy()
{
	int i;
	for (i = 0; i < nlocal; i++)
	{
		compute_e(i);
	}
}

/* ----------------------------------------------------------------------
  only local site call this function
 ---------------------------------------------------------------------- */
void AppVacancy::compute_e(int initsite)
{
	double e_v, e_r, e_s;
	int ip, im, jp, jm, kp, km;
	int ip1, im1, jp1, jm1, kp1, km1;
	int ip2, im2, jp2, jm2, kp2, km2;
	int it_nn1[8], it_nn2[6], it_nn3[12];
	int i, j, k;

	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);

	inc_px(i, 1, &ip1);
	inc_py(j, 1, &jp1);
	inc_pz(k, 1, &kp1);
	inc_px(i, -1, &im1);
	inc_py(j, -1, &jm1);
	inc_pz(k, -1, &km1);

	inc_px(i, 2, &ip2);
	inc_py(j, 2, &jp2);
	inc_pz(k, 2, &kp2);
	inc_px(i, -2, &im2);
	inc_py(j, -2, &jm2);
	inc_pz(k, -2, &km2);

	box2sub(&i, &j, &k);
	box2sub(&ip1, &jp1, &kp1);
	box2sub(&im1, &jm1, &km1);
	box2sub(&ip2, &jp2, &kp2);
	box2sub(&im2, &jm2, &km2);

	it_nn1[0] = lattice[get_POS_ID(im1, jm1, km1)];
	it_nn1[1] = lattice[get_POS_ID(im1, jp1, km1)];
	it_nn1[2] = lattice[get_POS_ID(ip1, jm1, km1)];
	it_nn1[3] = lattice[get_POS_ID(ip1, jp1, km1)];
	it_nn1[4] = lattice[get_POS_ID(im1, jm1, kp1)];
	it_nn1[5] = lattice[get_POS_ID(im1, jp1, kp1)];
	it_nn1[6] = lattice[get_POS_ID(ip1, jm1, kp1)];
	it_nn1[7] = lattice[get_POS_ID(ip1, jp1, kp1)];

	it_nn2[0] = lattice[get_POS_ID(im2, j, k)];
	it_nn2[1] = lattice[get_POS_ID(ip2, j, k)];
	it_nn2[2] = lattice[get_POS_ID(i, jm2, k)];
	it_nn2[3] = lattice[get_POS_ID(i, jp2, k)];
	it_nn2[4] = lattice[get_POS_ID(i, j, km2)];
	it_nn2[5] = lattice[get_POS_ID(i, j, kp2)];

	it_nn3[0] = lattice[get_POS_ID(i, jm2, km2)];
	it_nn3[1] = lattice[get_POS_ID(i, jp2, km2)];
	it_nn3[2] = lattice[get_POS_ID(i, jm2, kp2)];
	it_nn3[3] = lattice[get_POS_ID(i, jp2, kp2)];
	it_nn3[4] = lattice[get_POS_ID(im2, j, kp2)];
	it_nn3[5] = lattice[get_POS_ID(im2, j, km2)];
	it_nn3[6] = lattice[get_POS_ID(ip2, j, kp2)];
	it_nn3[7] = lattice[get_POS_ID(ip2, j, km2)];
	it_nn3[8] = lattice[get_POS_ID(im2, jm2, k)];
	it_nn3[9] = lattice[get_POS_ID(im2, jp2, k)];
	it_nn3[10] = lattice[get_POS_ID(ip2, jm2, k)];
	it_nn3[11] = lattice[get_POS_ID(ip2, jp2, k)];

	eatom(lattice[initsite], it_nn1, it_nn2, it_nn3, &e_v, &e_r, &e_s);
	set_e(initsite, e_v, e_r, e_s);
}

int AppVacancy::type(int i)
{
	return lattice[i];
}

int AppVacancy::get_cu_alone(int initsite)
{

	int ip, im, jp, jm, kp, km;
	int ip1, im1, jp1, jm1, kp1, km1;
	int it_nn1[8];
	int i, j, k;

	int n_cu_alone = 0;

	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);

	inc_px(i, 1, &ip1);
	inc_py(j, 1, &jp1);
	inc_pz(k, 1, &kp1);
	inc_px(i, -1, &im1);
	inc_py(j, -1, &jm1);
	inc_pz(k, -1, &km1);

	box2sub(&i, &j, &k);
	box2sub(&ip1, &jp1, &kp1);
	box2sub(&im1, &jm1, &km1);

	it_nn1[0] = lattice[get_POS_ID(im1, jm1, km1)];
	it_nn1[1] = lattice[get_POS_ID(im1, jp1, km1)];
	it_nn1[2] = lattice[get_POS_ID(ip1, jm1, km1)];
	it_nn1[3] = lattice[get_POS_ID(ip1, jp1, km1)];
	it_nn1[4] = lattice[get_POS_ID(im1, jm1, kp1)];
	it_nn1[5] = lattice[get_POS_ID(im1, jp1, kp1)];
	it_nn1[6] = lattice[get_POS_ID(ip1, jm1, kp1)];
	it_nn1[7] = lattice[get_POS_ID(ip1, jp1, kp1)];

	if (lattice[initsite] == 2 && it_nn1[0] == 1 && it_nn1[1] == 1 &&
		it_nn1[2] == 1 && it_nn1[3] == 1 && it_nn1[4] == 1 && it_nn1[5] == 1 &&
		it_nn1[6] == 1 && it_nn1[7] == 1)
	{
		n_cu_alone++;
	}
	return n_cu_alone;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppVacancy::site_energy(int i)
{
	int it;

	it = lattice[i];
	return single_energy(it, E_V[i], E_R[i], E_S[i]);
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
   i--- lattice index
   site_propensity---Ñ°ÕÒ²¢½¨Á¢site i¿ÉÄÜ·¢ÉúµÄÊÂ¼þ£¬¼ÆËã¸÷ÊÂ¼þ¸ÅÂÊ£¬·µ»Ø×Ü¸ÅÂÊ
------------------------------------------------------------------------- */

double AppVacancy::site_propensity(int i, double *vac_propensity_dir,
								   int *vac_jumpsite_dir, int *num_dir)
{
	/** if VACANT, proball = 0.
	 * else, proall = sum up of probs of jumping to 8 neighbors
	 * */
	double probone, proball, dene = 0;
	int jumpsite;
	int num;

	if (lattice[i] != VACANT)
	{
		return 0.0;
	}
	proball = 0.0;

	num = 0;
	for (int di = -1; di <= 1; di += 2)
	{
		for (int dj = -1; dj <= 1; dj += 2)
		{
			for (int dk = -1; dk <= 1; dk += 2)
			{
				dene = calcul_de(i, di, dj, dk, &jumpsite); // change of energy

				if (lattice[jumpsite] == VACANT)
				{
					continue;
				}

				probone = nu * exp(-(0.5 * dene + barrier[lattice[jumpsite]]) *
								   t_inverse / kb);
				proball += probone;

				vac_propensity_dir[num] = probone;
				vac_jumpsite_dir[num] = jumpsite;
				num = num + 1;
			}
		}
	}
	*num_dir = num - 1;
	return proball;
}

// added by xulei
int AppVacancy::get_POS_ID(int x, int y, int z)
{

	int x1, x2, y1, y2, z1, z2;
	int pos;

	if (z < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 + nlocal;
		return pos;
	}

	// The top of the ghost zone
	if (z >= 4 + subzhi - subzlo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2;
		return pos;
	}

	// The middle front end of the ghost zone
	if (y < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 + nlocal;
		return pos;
	}

	// The middle back end of the ghost zone
	if (y >= 4 + subyhi - subylo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (z1 - zghost / 2 + 1) * (subxhi - subxlo) * (subyhi - subylo) / 2 + nlocal;
		return pos;
	}

	// The middle left side of the ghost zone
	if (x < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (y1 - yghost / 2) * (subxhi - subxlo) - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 +
			  nlocal;
		return pos;
	}

	// The middle right side of the ghost zone
	if (x >= 4 + subxhi - subxlo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (y1 - yghost / 2 + 1) * (subxhi - subxlo) - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 +
			  nlocal;
		return pos;
	}

	// for local area
	x1 = (x - 4) / 2;
	x2 = (x - 4) % 2;
	y1 = (y - 4) / 2;
	y2 = (y - 4) % 2;
	z1 = (z - 4) / 2;
	z2 = (z - 4) % 2;

	pos = x1 * 2 + y1 * (subxhi - subxlo) + z1 * (subxhi - subxlo) / 2 * (subyhi - subylo) + x2;

	return pos;
}

/*
	add a flag in get_POS_ID_v2,
	That flag is 1 means this pos is in the local area;
				0 means this pos is in the ghost area.
*/
int AppVacancy::get_POS_ID_v2(int x, int y, int z, int &flag)
{

	int x1, x2, y1, y2, z1, z2;
	int pos;

	flag = 0;

	if (z < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 + nlocal;
		return pos;
	}

	// The top of the ghost zone
	if (z >= 4 + subzhi - subzlo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2;
		return pos;
	}

	// The middle front end of the ghost zone
	if (y < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 + nlocal;
		return pos;
	}

	// The middle back end of the ghost zone
	if (y >= 4 + subyhi - subylo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (z1 - zghost / 2 + 1) * (subxhi - subxlo) * (subyhi - subylo) / 2 +
			  nlocal;
		return pos;
	}

	// The middle left side of the ghost zone
	if (x < 4)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (y1 - yghost / 2) * (subxhi - subxlo) - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 +
			  nlocal;
		return pos;
	}

	// The middle right side of the ghost zone
	if (x >= 4 + subxhi - subxlo)
	{
		x1 = (x - 4 + xghost) / 2;
		x2 = x % 2;
		y1 = (y - 4 + yghost) / 2;
		y2 = y % 2;
		z1 = (z - 4 + zghost) / 2;
		z2 = z % 2;

		pos = x1 * 2 + y1 * (subxhi - subxlo + xghost * 2) +
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) +
			  x2 - (y1 - yghost / 2 + 1) * (subxhi - subxlo) - (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 +
			  nlocal;
		return pos;
	}

	// for local area
	if (x == 4 || y == 4 || z == 4 || x == (3 + subxhi - subxlo) || y == (3 + subyhi - subylo) || z == (3 + subzhi - subzlo))
		flag = 1;

	x1 = (x - 4) / 2;
	x2 = (x - 4) % 2;
	y1 = (y - 4) / 2;
	y2 = (y - 4) % 2;
	z1 = (z - 4) / 2;
	z2 = (z - 4) % 2;

	pos = x1 * 2 + y1 * (subxhi - subxlo) + z1 * (subxhi - subxlo) / 2 * (subyhi - subylo) + x2;

	return pos;
}

double AppVacancy::calcul_de(int initsite, int di, int dj, int dk, int *jumpsite)
{
	int i, j, k;
	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);

	int di_p1, dj_p1, dk_p1, di_m1, dj_m1, dk_m1;
	int di_m2, di_p2, dj_m2, dj_p2, dk_m2, dk_p2;
	int di_p3, dj_p3, dk_p3;
	int new_di_p1, new_dj_p1, new_dk_p1, new_di_m1, new_dj_m1, new_dk_m1;
	int new_di_m2, new_di_p2, new_dj_m2, new_dj_p2, new_dk_m2, new_dk_p2;
	int it;

	inc_px(i, di, &di_p1);
	inc_py(j, dj, &dj_p1);
	inc_pz(k, dk, &dk_p1);
	box2sub(&di_p1, &dj_p1, &dk_p1);

	it = lattice[get_POS_ID(di_p1, dj_p1, dk_p1)];
	*jumpsite = get_POS_ID(di_p1, dj_p1, dk_p1);

	if (it == 0)
	{
		return 0.0;
	}

	inc_px(i, -di, &di_m1);
	inc_py(j, -dj, &dj_m1);
	inc_pz(k, -dk, &dk_m1);

	inc_px(i, 2 * di, &di_p2);
	inc_py(j, 2 * dj, &dj_p2);
	inc_pz(k, 2 * dk, &dk_p2);
	inc_px(i, -2 * di, &di_m2);
	inc_py(j, -2 * dj, &dj_m2);
	inc_pz(k, -2 * dk, &dk_m2);

	inc_px(i, 3 * di, &di_p3);
	inc_py(j, 3 * dj, &dj_p3);
	inc_pz(k, 3 * dk, &dk_p3);

	box2sub(&i, &j, &k);
	box2sub(&di_m1, &dj_m1, &dk_m1);
	box2sub(&di_p2, &dj_p2, &dk_p2);
	box2sub(&di_m2, &dj_m2, &dk_m2);
	box2sub(&di_p3, &dj_p3, &dk_p3);

	double dene = 0.;

	// 1.energy of jumpsite
	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		dene += -site_energy(*jumpsite);
	}
	else
	{
		dene += -site_energy(get_POS_ID(i, j, k)) - site_energy(*jumpsite);
	}

	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		// 2.energy changes of neighbors
		///* 1nn vacancy */
		dene = dene + calcul_dene_v1(it, di_m1, dj_m1, dk_m1, 1, +1.0); // 1

		dene = dene + calcul_dene_v2(it, di_m1, dj_p1, dk_p1, 1, 2); // 2
		dene = dene + calcul_dene_v2(it, di_p1, dj_m1, dk_m1, 1, 3); // 3

		dene = dene + calcul_dene_v2(it, di_p1, dj_p1, dk_m1, 1, 2); // 2
		dene = dene + calcul_dene_v2(it, di_m1, dj_m1, dk_p1, 1, 3); // 3

		dene = dene + calcul_dene_v2(it, di_p1, dj_m1, dk_p1, 1, 2); // 2
		dene = dene + calcul_dene_v2(it, di_m1, dj_p1, dk_m1, 1, 3); // 3

		///* 2nn vacancy and 1nn atom it */
		dene = dene + calcul_dene_v2(it, di_p2, j, k, 2, 1); // 4
		dene = dene + calcul_dene_v2(it, i, dj_p2, k, 2, 1); // 4
		dene = dene + calcul_dene_v2(it, i, j, dk_p2, 2, 1); // 4

		///* 1nn atom it */
		dene = dene + calcul_dene_v1(it, di_p2, dj_p2, dk_p2, 1, -1.0); // 5

		dene = dene + calcul_dene_v2(it, i, dj_p2, dk_p2, 3, 1); // 6
		dene = dene + calcul_dene_v2(it, di_p2, j, dk_p2, 3, 1); // 6
		dene = dene + calcul_dene_v2(it, di_p2, dj_p2, k, 3, 1); // 6

		///* 2nn vacancy */
		dene = dene + calcul_dene_v1(it, di_m2, j, k, 2, +1.0); // 7
		dene = dene + calcul_dene_v1(it, i, dj_m2, k, 2, +1.0); // 7
		dene = dene + calcul_dene_v1(it, i, j, dk_m2, 2, +1.0); // 7

		///* 2nn atom it */
		dene = dene + calcul_dene_v1(it, di_p3, dj_p1, dk_p1, 2, -1.0); // 8
		dene = dene + calcul_dene_v1(it, di_p1, dj_p3, dk_p1, 2, -1.0); // 8
		dene = dene + calcul_dene_v1(it, di_p1, dj_p1, dk_p3, 2, -1.0); // 8

		///* 3nn vacancy */
		dene = dene + calcul_dene_v1(it, di_m2, dj_m2, k, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, di_p2, dj_m2, k, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, di_m2, dj_p2, k, 3, +1.0); // 9

		dene = dene + calcul_dene_v1(it, di_m2, j, dk_p2, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, di_p2, j, dk_m2, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, di_m2, j, dk_m2, 3, +1.0); // 9

		dene = dene + calcul_dene_v1(it, i, dj_m2, dk_m2, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, i, dj_m2, dk_p2, 3, +1.0); // 9
		dene = dene + calcul_dene_v1(it, i, dj_p2, dk_m2, 3, +1.0); // 9

		///* 3nn atom it */
		dene = dene + calcul_dene_v1(it, di_p3, dj_p3, dk_p1, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_p3, dj_m1, dk_p1, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_m1, dj_p3, dk_p1, 3, -1.0); // 10

		dene = dene + calcul_dene_v1(it, di_p1, dj_p3, dk_p3, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_p1, dj_p3, dk_m1, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_p1, dj_m1, dk_p3, 3, -1.0); // 10

		dene = dene + calcul_dene_v1(it, di_p3, dj_p1, dk_p3, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_p3, dj_p1, dk_m1, 3, -1.0); // 10
		dene = dene + calcul_dene_v1(it, di_m1, dj_p1, dk_p3, 3, -1.0); // 10
	}

	double e_v, e_r, e_s;
	int it_nn1[8], it_nn2[6], it_nn3[12];

	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		it_nn1[0] = lattice[get_POS_ID(di_m1, dj_m1, dk_m1)];
		it_nn1[1] = lattice[get_POS_ID(di_m1, dj_p1, dk_m1)];
		it_nn1[2] = lattice[get_POS_ID(di_p1, dj_m1, dk_m1)];
		it_nn1[3] = lattice[get_POS_ID(di_p1, dj_p1, dk_m1)];
		it_nn1[4] = lattice[get_POS_ID(di_m1, dj_m1, dk_p1)];
		it_nn1[5] = lattice[get_POS_ID(di_m1, dj_p1, dk_p1)];
		it_nn1[6] = lattice[get_POS_ID(di_p1, dj_m1, dk_p1)];
		it_nn1[7] = 0;

		it_nn2[0] = lattice[get_POS_ID(di_m2, j, k)];
		it_nn2[1] = lattice[get_POS_ID(di_p2, j, k)];
		it_nn2[2] = lattice[get_POS_ID(i, dj_m2, k)];
		it_nn2[3] = lattice[get_POS_ID(i, dj_p2, k)];
		it_nn2[4] = lattice[get_POS_ID(i, j, dk_m2)];
		it_nn2[5] = lattice[get_POS_ID(i, j, dk_p2)];

		it_nn3[0] = lattice[get_POS_ID(i, dj_m2, dk_m2)];
		it_nn3[1] = lattice[get_POS_ID(i, dj_p2, dk_m2)];
		it_nn3[2] = lattice[get_POS_ID(i, dj_m2, dk_p2)];
		it_nn3[3] = lattice[get_POS_ID(i, dj_p2, dk_p2)];
		it_nn3[4] = lattice[get_POS_ID(di_m2, j, dk_p2)];
		it_nn3[5] = lattice[get_POS_ID(di_m2, j, dk_m2)];
		it_nn3[6] = lattice[get_POS_ID(di_p2, j, dk_p2)];
		it_nn3[7] = lattice[get_POS_ID(di_p2, j, dk_m2)];
		it_nn3[8] = lattice[get_POS_ID(di_m2, dj_m2, k)];
		it_nn3[9] = lattice[get_POS_ID(di_m2, dj_p2, k)];
		it_nn3[10] = lattice[get_POS_ID(di_p2, dj_m2, k)];
		it_nn3[11] = lattice[get_POS_ID(di_p2, dj_p2, k)];

		eatom(it, it_nn1, it_nn2, it_nn3, &e_v, &e_r, &e_s);
		delete1(it, &e_v, &e_r, &e_s);
		dene += single_energy(it, e_v, e_r, e_s);
	}
	else
	{
		it_nn1[0] = lattice[get_POS_ID(di_m1, dj_m1, dk_m1)];
		it_nn1[1] = lattice[get_POS_ID(di_m1, dj_p1, dk_m1)];
		it_nn1[2] = lattice[get_POS_ID(di_p1, dj_m1, dk_m1)];
		it_nn1[3] = lattice[get_POS_ID(di_p1, dj_p1, dk_m1)];
		it_nn1[4] = lattice[get_POS_ID(di_m1, dj_m1, dk_p1)];
		it_nn1[5] = lattice[get_POS_ID(di_m1, dj_p1, dk_p1)];
		it_nn1[6] = lattice[get_POS_ID(di_p1, dj_m1, dk_p1)];
		it_nn1[7] = 0;

		it_nn2[0] = lattice[get_POS_ID(di_m2, j, k)];
		it_nn2[1] = lattice[get_POS_ID(di_p2, j, k)];
		it_nn2[2] = lattice[get_POS_ID(i, dj_m2, k)];
		it_nn2[3] = lattice[get_POS_ID(i, dj_p2, k)];
		it_nn2[4] = lattice[get_POS_ID(i, j, dk_m2)];
		it_nn2[5] = lattice[get_POS_ID(i, j, dk_p2)];

		eatom(it, it_nn1, it_nn2, it_nn3, &e_v, &e_r, &e_s);
		dene = dene + 0.5 * e_v;

		i = round(xyz[initsite][0] * unitdis);
		j = round(xyz[initsite][1] * unitdis);
		k = round(xyz[initsite][2] * unitdis);
		inc_px(i, di, &di_p1);
		inc_py(j, dj, &dj_p1);
		inc_pz(k, dk, &dk_p1);

		inc_px(di_p1, di, &new_di_p1);
		inc_py(dj_p1, dj, &new_dj_p1);
		inc_pz(dk_p1, dk, &new_dk_p1);
		inc_px(di_p1, -di, &new_di_m1);
		inc_py(dj_p1, -dj, &new_dj_m1);
		inc_pz(dk_p1, -dk, &new_dk_m1);

		inc_px(di_p1, 2 * di, &new_di_p2);
		inc_py(dj_p1, 2 * dj, &new_dj_p2);
		inc_pz(dk_p1, 2 * dk, &new_dk_p2);
		inc_px(di_p1, -2 * di, &new_di_m2);
		inc_py(dj_p1, -2 * dj, &new_dj_m2);
		inc_pz(dk_p1, -2 * dk, &new_dk_m2);

		box2sub(&di_p1, &dj_p1, &dk_p1);
		box2sub(&new_di_p1, &new_dj_p1, &new_dk_p1);
		box2sub(&new_di_m1, &new_dj_m1, &new_dk_m1);
		box2sub(&new_di_p2, &new_dj_p2, &new_dk_p2);
		box2sub(&new_di_m2, &new_dj_m2, &new_dk_m2);

		it_nn1[0] = it;

		it_nn1[1] = lattice[get_POS_ID(new_di_m1, new_dj_p1, new_dk_m1)];
		it_nn1[2] = lattice[get_POS_ID(new_di_p1, new_dj_m1, new_dk_m1)];
		it_nn1[3] = lattice[get_POS_ID(new_di_p1, new_dj_p1, new_dk_m1)];
		it_nn1[4] = lattice[get_POS_ID(new_di_m1, new_dj_m1, new_dk_p1)];
		it_nn1[5] = lattice[get_POS_ID(new_di_m1, new_dj_p1, new_dk_p1)];
		it_nn1[6] = lattice[get_POS_ID(new_di_p1, new_dj_m1, new_dk_p1)];
		it_nn1[7] = lattice[get_POS_ID(new_di_p1, new_dj_p1, new_dk_p1)];

		it_nn2[0] = lattice[get_POS_ID(new_di_m2, dj_p1, dk_p1)];
		it_nn2[1] = lattice[get_POS_ID(new_di_p2, dj_p1, dk_p1)];
		it_nn2[2] = lattice[get_POS_ID(di_p1, new_dj_m2, dk_p1)];
		it_nn2[3] = lattice[get_POS_ID(di_p1, new_dj_p2, dk_p1)];
		it_nn2[4] = lattice[get_POS_ID(di_p1, dj_p1, new_dk_m2)];
		it_nn2[5] = lattice[get_POS_ID(di_p1, dj_p1, new_dk_p2)];

		eatom(0, it_nn1, it_nn2, it_nn3, &e_v, &e_r, &e_s);
		dene = dene + 0.5 * e_v;
	}

	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		return dene;
	}
	else
	{
		return dene * 2.0;
	}
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
   * i -- lattice index
------------------------------------------------------------------------- */

void AppVacancy::site_event(int ii, int *notinset, class RandomPark *random)
{
	int i_last_time, j_last_time, vac, jsite;
	int nsite = 0, *index;
	double random_number, threshhold;
	double tmp_sum = 0.0;
	int num;
	index = new int[cvector_size(*vacant) + 1];

	static unsigned int counter = 0;
	counter++;
	i_last_time = (*vacant)[ii];
	random_number = random->uniform();
	threshhold = propensity[ii] * random_number;
	num = num_dir[ii];
	while (tmp_sum < threshhold)
	{
		tmp_sum = tmp_sum + propensity_dir[ii][num];
		num = num - 1;
	}
	j_last_time = jumpsite_dir[ii][num + 1];

	do_jump(i_last_time, j_last_time);

	int flag;
	if (sectorflag == 0)
	{
		flag = 1;
	}
	else
	{
		flag = is_in_set(j_last_time);
	}

	unsigned long tst_kernel, ted_kernel;
	if (flag == 1)
	{
		(*vacant)[ii] = j_last_time;
	}
	else
	{ // if jsite is not in the set, update all propensity of other vacancies.

		(*notinset)++;
		int last_idx = cvector_size((*vacant)) - 1;
		vac = (*vacant)[last_idx];
		(*vacant)[ii] = vac;
		index[nsite++] = ii;
		propensity[ii] = propensity[last_idx];
		num_dir[ii] = num_dir[last_idx];
		for (int k_dir = 0; k_dir < 8; k_dir++)
		{
			propensity_dir[ii][k_dir] = propensity_dir[last_idx][k_dir];
			jumpsite_dir[ii][k_dir] = jumpsite_dir[last_idx][k_dir];
		}
		index[nsite++] = last_idx;
		propensity[last_idx] = 0;
		num_dir[last_idx] = 0;
		for (int k_dir = 0; k_dir < 8; k_dir++)
		{
			propensity_dir[last_idx][k_dir] = 0;
			jumpsite_dir[last_idx][k_dir] = 0;
		}

		cvector_erase(*vacant, last_idx);
	}
	int vac_size = cvector_size((*vacant));

	set_arg_of_kernel_func(index, &nsite, i_last_time, j_last_time, *vacant, vac_size);

	kernel_func(&argOfKernel);

	if (nset > 1)
		update_wait_send_map(i_last_time, j_last_time);
	if (flag == 0 && j_last_time < nlocal)
	{
		int mysector = in_which_set(j_last_time);
		cvector_push_back(comm->listVacanciesUpdated[mysector], j_last_time);
	}

	call_athread_join();

	solve->update(nsite, index, propensity);

	delete[] index;
}

/* ---------------------------------------------------------------------- */

void AppVacancy::do_jump(int initsite, int jumpsite)
{
	int i, j, k, i2, j2, k2, di, dj, dk;

	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);
	i2 = round(xyz[jumpsite][0] * unitdis);
	j2 = round(xyz[jumpsite][1] * unitdis);
	k2 = round(xyz[jumpsite][2] * unitdis);
	di = i2 - i;
	if (di == -boxx + 1)
		di = 1;
	if (di == boxx - 1)
		di = -1;
	dj = j2 - j;
	if (dj == -boxy + 1)
		dj = 1;
	if (dj == boxy - 1)
		dj = -1;
	dk = k2 - k;
	if (dk == -boxz + 1)
		dk = 1;
	if (dk == boxz - 1)
		dk = -1;
	if (abs(di) != 1 || abs(dj) != 1 || abs(dk) != 1)
		error->all(FLERR, "The jumpsite is not the neighbor of initsite.\n");

	int di_p1, dj_p1, dk_p1, di_m1, dj_m1, dk_m1;
	int di_m2, di_p2, dj_m2, dj_p2, dk_m2, dk_p2;
	int di_p3, dj_p3, dk_p3;
	int it;

	inc_px(i, di, &di_p1);
	inc_py(j, dj, &dj_p1);
	inc_pz(k, dk, &dk_p1);
	inc_px(i, -di, &di_m1);
	inc_py(j, -dj, &dj_m1);
	inc_pz(k, -dk, &dk_m1);

	inc_px(i, 2 * di, &di_p2);
	inc_py(j, 2 * dj, &dj_p2);
	inc_pz(k, 2 * dk, &dk_p2);
	inc_px(i, -2 * di, &di_m2);
	inc_py(j, -2 * dj, &dj_m2);
	inc_pz(k, -2 * dk, &dk_m2);

	inc_px(i, 3 * di, &di_p3);
	inc_py(j, 3 * dj, &dj_p3);
	inc_pz(k, 3 * dk, &dk_p3);

	box2sub(&i, &j, &k);
	box2sub(&di_p1, &dj_p1, &dk_p1);
	box2sub(&di_m1, &dj_m1, &dk_m1);
	box2sub(&di_p2, &dj_p2, &dk_p2);
	box2sub(&di_m2, &dj_m2, &dk_m2);
	box2sub(&di_p3, &dj_p3, &dk_p3);

	it = lattice[get_POS_ID(di_p1, dj_p1, dk_p1)];

	/* 1nn vacancy */
	do_jump_v1(it, di_m1, dj_m1, dk_m1, 1, +1.0); // 1 //problem

	do_jump_v2(it, di_m1, dj_p1, dk_p1, 1, 2); // 2
	do_jump_v2(it, di_p1, dj_m1, dk_m1, 1, 3); // 3

	do_jump_v2(it, di_p1, dj_p1, dk_m1, 1, 2); // 2
	do_jump_v2(it, di_m1, dj_m1, dk_p1, 1, 3); // 3

	do_jump_v2(it, di_p1, dj_m1, dk_p1, 1, 2); // 2
	do_jump_v2(it, di_m1, dj_p1, dk_m1, 1, 3); // 3

	/* 2nn vacancy and 1nn atom it */
	do_jump_v2(it, di_p2, j, k, 2, 1); // 4
	do_jump_v2(it, i, dj_p2, k, 2, 1); // 4
	do_jump_v2(it, i, j, dk_p2, 2, 1); // 4

	/* 1nn atom it */
	do_jump_v1(it, di_p2, dj_p2, dk_p2, 1, -1.0); // 5

	do_jump_v2(it, i, dj_p2, dk_p2, 3, 1); // 6
	do_jump_v2(it, di_p2, j, dk_p2, 3, 1); // 6
	do_jump_v2(it, di_p2, dj_p2, k, 3, 1); // 6

	/* 2nn vacancy */
	do_jump_v1(it, di_m2, j, k, 2, +1.0); // 7
	do_jump_v1(it, i, dj_m2, k, 2, +1.0); // 7
	do_jump_v1(it, i, j, dk_m2, 2, +1.0); // 7

	/* 2nn atom it */
	do_jump_v1(it, di_p3, dj_p1, dk_p1, 2, -1.0); // 8
	do_jump_v1(it, di_p1, dj_p3, dk_p1, 2, -1.0); // 8
	do_jump_v1(it, di_p1, dj_p1, dk_p3, 2, -1.0); // 8

	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		/* 3nn vacancy */
		do_jump_v1(it, di_m2, dj_m2, k, 3, +1.0); // 9
		do_jump_v1(it, di_p2, dj_m2, k, 3, +1.0); // 9
		do_jump_v1(it, di_m2, dj_p2, k, 3, +1.0); // 9

		do_jump_v1(it, di_m2, j, dk_p2, 3, +1.0); // 9
		do_jump_v1(it, di_p2, j, dk_m2, 3, +1.0); // 9
		do_jump_v1(it, di_m2, j, dk_m2, 3, +1.0); // 9

		do_jump_v1(it, i, dj_m2, dk_m2, 3, +1.0); // 9
		do_jump_v1(it, i, dj_m2, dk_p2, 3, +1.0); // 9
		do_jump_v1(it, i, dj_p2, dk_m2, 3, +1.0); // 9

		/* 3nn atom it */
		do_jump_v1(it, di_p3, dj_p3, dk_p1, 3, -1.0); // 10
		do_jump_v1(it, di_p3, dj_m1, dk_p1, 3, -1.0); // 10
		do_jump_v1(it, di_m1, dj_p3, dk_p1, 3, -1.0); // 10

		do_jump_v1(it, di_p1, dj_p3, dk_p3, 3, -1.0); // 10
		do_jump_v1(it, di_p1, dj_p3, dk_m1, 3, -1.0); // 10
		do_jump_v1(it, di_p1, dj_m1, dk_p3, 3, -1.0); // 10

		do_jump_v1(it, di_p3, dj_p1, dk_p3, 3, -1.0); // 10
		do_jump_v1(it, di_p3, dj_p1, dk_m1, 3, -1.0); // 10
		do_jump_v1(it, di_m1, dj_p1, dk_p3, 3, -1.0); // 10
	}

	lattice[initsite] = lattice[jumpsite];
	lattice[jumpsite] = VACANT;

	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		compute_e(initsite);
		set_e(jumpsite, 0, 0, 0);
	}
	else
	{
		compute_e(initsite);
		compute_e(jumpsite);
	}

	// if (nset > 1)
	// {
	// 	add_to_wait(i, j, k);
	// 	add_to_wait4TypeChange(di_p1, dj_p1, dk_p1);

	// 	add_to_wait(di_m1, dj_m1, dk_m1);
	// 	add_to_wait(di_m1, dj_p1, dk_p1);
	// 	add_to_wait(di_p1, dj_m1, dk_m1);
	// 	add_to_wait(di_p1, dj_p1, dk_m1);
	// 	add_to_wait(di_m1, dj_m1, dk_p1);
	// 	add_to_wait(di_p1, dj_m1, dk_p1);
	// 	add_to_wait(di_m1, dj_p1, dk_m1);

	// 	add_to_wait(di_p2, j, k);
	// 	add_to_wait(i, dj_p2, k);
	// 	add_to_wait(i, j, dk_p2);
	// 	add_to_wait(di_m2, j, k);
	// 	add_to_wait(i, dj_m2, k);
	// 	add_to_wait(i, j, dk_m2);

	// 	add_to_wait(di_p2, dj_p2, dk_p2);
	// 	add_to_wait(i, dj_p2, dk_p2);
	// 	add_to_wait(di_p2, j, dk_p2);
	// 	add_to_wait(di_p2, dj_p2, k);
	// 	add_to_wait(di_m2, dj_m2, k);
	// 	add_to_wait(di_m2, j, dk_m2);
	// 	add_to_wait(i, dj_m2, dk_m2);

	// 	add_to_wait(di_p3, dj_p1, dk_p1);
	// 	add_to_wait(di_p1, dj_p3, dk_p1);
	// 	add_to_wait(di_p1, dj_p1, dk_p3);
	// 	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	// 	{
	// 		add_to_wait(di_p2, dj_m2, k);
	// 		add_to_wait(di_m2, dj_p2, k);
	// 		add_to_wait(di_m2, j, dk_p2);
	// 		add_to_wait(di_p2, j, dk_m2);
	// 		add_to_wait(i, dj_m2, dk_p2);
	// 		add_to_wait(i, dj_p2, dk_m2);
	// 		add_to_wait(di_p3, dj_p3, dk_p1);
	// 		add_to_wait(di_p3, dj_m1, dk_p1);
	// 		add_to_wait(di_m1, dj_p3, dk_p1);
	// 		add_to_wait(di_p1, dj_p3, dk_p3);
	// 		add_to_wait(di_p1, dj_p3, dk_m1);
	// 		add_to_wait(di_p1, dj_m1, dk_p3);
	// 		add_to_wait(di_p3, dj_p1, dk_p3);
	// 		add_to_wait(di_p3, dj_p1, dk_m1);
	// 		add_to_wait(di_m1, dj_p1, dk_p3);
	// 	}
	// }
}

void AppVacancy::update_wait_send_map(int initsite, int jumpsite)
{
	int i, j, k, i2, j2, k2, di, dj, dk;

	i = round(xyz[initsite][0] * unitdis);
	j = round(xyz[initsite][1] * unitdis);
	k = round(xyz[initsite][2] * unitdis);
	i2 = round(xyz[jumpsite][0] * unitdis);
	j2 = round(xyz[jumpsite][1] * unitdis);
	k2 = round(xyz[jumpsite][2] * unitdis);

	di = i2 - i;
	if (di == -boxx + 1)
		di = 1;
	if (di == boxx - 1)
		di = -1;
	dj = j2 - j;
	if (dj == -boxy + 1)
		dj = 1;
	if (dj == boxy - 1)
		dj = -1;
	dk = k2 - k;
	if (dk == -boxz + 1)
		dk = 1;
	if (dk == boxz - 1)
		dk = -1;

	int di_p1, dj_p1, dk_p1, di_m1, dj_m1, dk_m1;
	int di_m2, di_p2, dj_m2, dj_p2, dk_m2, dk_p2;
	int di_p3, dj_p3, dk_p3;

	inc_px(i, di, &di_p1);
	inc_py(j, dj, &dj_p1);
	inc_pz(k, dk, &dk_p1);
	inc_px(i, -di, &di_m1);
	inc_py(j, -dj, &dj_m1);
	inc_pz(k, -dk, &dk_m1);

	inc_px(i, 2 * di, &di_p2);
	inc_py(j, 2 * dj, &dj_p2);
	inc_pz(k, 2 * dk, &dk_p2);
	inc_px(i, -2 * di, &di_m2);
	inc_py(j, -2 * dj, &dj_m2);
	inc_pz(k, -2 * dk, &dk_m2);

	inc_px(i, 3 * di, &di_p3);
	inc_py(j, 3 * dj, &dj_p3);
	inc_pz(k, 3 * dk, &dk_p3);

	box2sub(&i, &j, &k);
	box2sub(&di_p1, &dj_p1, &dk_p1);
	box2sub(&di_m1, &dj_m1, &dk_m1);
	box2sub(&di_p2, &dj_p2, &dk_p2);
	box2sub(&di_m2, &dj_m2, &dk_m2);
	box2sub(&di_p3, &dj_p3, &dk_p3);

	add_to_wait(i, j, k);
	add_to_wait4TypeChange(di_p1, dj_p1, dk_p1);

	add_to_wait(di_m1, dj_m1, dk_m1);
	add_to_wait(di_m1, dj_p1, dk_p1);
	add_to_wait(di_p1, dj_m1, dk_m1);
	add_to_wait(di_p1, dj_p1, dk_m1);
	add_to_wait(di_m1, dj_m1, dk_p1);
	add_to_wait(di_p1, dj_m1, dk_p1);
	add_to_wait(di_m1, dj_p1, dk_m1);

	add_to_wait(di_p2, j, k);
	add_to_wait(i, dj_p2, k);
	add_to_wait(i, j, dk_p2);
	add_to_wait(di_m2, j, k);
	add_to_wait(i, dj_m2, k);
	add_to_wait(i, j, dk_m2);

	add_to_wait(di_p2, dj_p2, dk_p2);
	add_to_wait(i, dj_p2, dk_p2);
	add_to_wait(di_p2, j, dk_p2);
	add_to_wait(di_p2, dj_p2, k);
	add_to_wait(di_m2, dj_m2, k);
	add_to_wait(di_m2, j, dk_m2);
	add_to_wait(i, dj_m2, dk_m2);

	add_to_wait(di_p3, dj_p1, dk_p1);
	add_to_wait(di_p1, dj_p3, dk_p1);
	add_to_wait(di_p1, dj_p1, dk_p3);
	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		add_to_wait(di_p2, dj_m2, k);
		add_to_wait(di_m2, dj_p2, k);
		add_to_wait(di_m2, j, dk_p2);
		add_to_wait(di_p2, j, dk_m2);
		add_to_wait(i, dj_m2, dk_p2);
		add_to_wait(i, dj_p2, dk_m2);
		add_to_wait(di_p3, dj_p3, dk_p1);
		add_to_wait(di_p3, dj_m1, dk_p1);
		add_to_wait(di_m1, dj_p3, dk_p1);
		add_to_wait(di_p1, dj_p3, dk_p3);
		add_to_wait(di_p1, dj_p3, dk_m1);
		add_to_wait(di_p1, dj_m1, dk_p3);
		add_to_wait(di_p3, dj_p1, dk_p3);
		add_to_wait(di_p3, dj_p1, dk_m1);
		add_to_wait(di_m1, dj_p1, dk_p3);
	}
}

void AppVacancy::add_to_wait(int x, int y, int z)
{
	int site = get_POS_ID(x, y, z);
	wait_send.insert(std::pair<int, int>(site, 0));
}

void AppVacancy::add_to_wait4TypeChange(int x, int y, int z)
{
	int site = get_POS_ID(x, y, z);
	auto ret = wait_send.insert(std::pair<int, int>(site, 1));
	if (!ret.second)
		ret.first->second = 1;
}

void AppVacancy::add_to_send()
{
	int i, p, index;
	int *proc;
	int x, y, z;
	int site;
	std::map<int, int>::iterator wait;
	int a = 1;
	for (wait = wait_send.begin(); wait != wait_send.end(); wait++)
	{
		site = wait->first;
		x = round(xyz[site][0] * unitdis);
		y = round(xyz[site][1] * unitdis);
		z = round(xyz[site][2] * unitdis);
		comm->add_to_send(site, x, y, z, wait->second);
	}
	wait_send.clear();
}

void AppVacancy::communicate()
{
	add_to_send();
	comm->perform_set();
}

void AppVacancy::communicate_priority(int iset)
{
	add_to_send();
	comm->perform_priority(iset);
}

void AppVacancy::communicate_perform_set(int iset)
{
	// printf("xulei communicate_perform_set %d\n", iset);
	comm->set_iset(iset);
	comm->perform_set();
	if (iset < 7)
		comm->set_iset(iset + 1);
}

void AppVacancy::set_it(int x, int y, int z, int it)
{
	box2sub(&x, &y, &z);
	lattice[get_POS_ID(x, y, z)] = it;
}

void AppVacancy::set_it_e(int x, int y, int z, int it, double ev, double er, double es)
{
	box2sub(&x, &y, &z);
	int pos = get_POS_ID(x, y, z);
	lattice[pos] = it;
	E_V[pos] = ev;
	E_R[pos] = er;
	E_S[pos] = es;
}

int AppVacancy::set_it_e_set(int x, int y, int z, int it, double ev, double er, double es, int &flag)
{
	box2sub(&x, &y, &z);
	int pos = get_POS_ID_v2(x, y, z, flag);
	lattice[pos] = it;
	E_V[pos] = ev;
	E_R[pos] = er;
	E_S[pos] = es;
	return pos;
}

void AppVacancy::get_it_e(int x, int y, int z, int *it, double *ev, double *er, double *es)
{
	box2sub(&x, &y, &z);
	int pos = get_POS_ID(x, y, z);
	*it = lattice[pos];
	*ev = E_V[pos];
	*er = E_R[pos];
	*es = E_S[pos];
}
void AppVacancy::set_e(int i, double ev, double er, double es)
{
	E_V[i] = ev;
	E_R[i] = er;
	E_S[i] = es;
}

void AppVacancy::get_e(int i, double *ev, double *er, double *es)
{
	*ev = E_V[i];
	*er = E_R[i];
	*es = E_S[i];
}

void AppVacancy::inc_px(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxxhi)
		*ii = *ii - boxxhi + boxxlo;
	else if (*ii < boxxlo)
		*ii = boxxhi + *ii - boxxlo;
}

void AppVacancy::inc_py(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxyhi)
		*ii = *ii - boxyhi + boxylo;
	else if (*ii < boxylo)
		*ii = boxyhi + *ii - boxylo;
}

void AppVacancy::inc_pz(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxzhi)
		*ii = *ii - boxzhi + boxzlo;
	else if (*ii < boxzlo)
		*ii = boxzhi + *ii - boxzlo;
}

double AppVacancy::compute_linedis(int site1, int site2)
{
	int i, j, k, i2, j2, k2, x, y, z;
	i = round(xyz[site1][0] * unitdis);
	j = round(xyz[site1][1] * unitdis);
	k = round(xyz[site1][2] * unitdis);
	i2 = round(xyz[site2][0] * unitdis);
	j2 = round(xyz[site2][1] * unitdis);
	k2 = round(xyz[site2][2] * unitdis);
	x = abs(i - i2);
	y = abs(j - j2);
	z = abs(k - k2);
	if (x * 2 > boxx)
		x = boxx - x;
	if (y * 2 > boxy)
		y = boxy - y;
	if (z * 2 > boxz)
		z = boxz - z;
	int dis = x * x + y * y + z * z;
	return dis;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppVacancy::clear_events(int i)
{
	int next;
	int index = firstevent[i];
	while (index >= 0)
	{
		next = events[index].next;
		events[index].next = freeevent;
		freeevent = index;
		nevents--;
		index = next;
	}
	firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
   i ---ÆðÊ¼OCCUPIEDµã
   destination ---ÌøÔ¾µ½µÄVACENTÄ¿µÄµã
   propensity ---ÌøÔ¾µÄ¸ÅÂÊ
   eventflag ---±êÊ¶ÄÄÖÖÊÂ¼þDEPOSITION/NNHOP/SCHWOEBEL
------------------------------------------------------------------------- */

void AppVacancy::add_event(int i, int destination,
						   double propensity)
{
	// grow event list and setup free list
	if (nevents == maxevent)
	{
		maxevent += DELTAEVENT;
		events =
			(Event *)memory->srealloc(events, maxevent * sizeof(Event), "app:events");
		for (int m = nevents; m < maxevent; m++)
			events[m].next = m + 1;
		freeevent = nevents;
	}

	int next = events[freeevent].next;
	events[freeevent].propensity = propensity;
	events[freeevent].destination = destination;
	events[freeevent].next = firstevent[i];
	firstevent[i] = freeevent;
	freeevent = next;
	nevents++;
}

void AppVacancy::init_potential_pair_file(char *arg)
{
	FILE *fin;
	int i, j, k;
	float tmp;
	char SS[200];

	fin = fopen(arg, "r");

	for (i = 0; i < NSPECIES; i++)
		for (j = 0; j < NSPECIES; j++)
		{
			vij1[i][j] = vij2[i][j] = vij3[i][j] = vij4[i][j] = 0.;
			rho1[i][j] = rho2[i][j] = rho3[i][j] = rho4[i][j] = 0.;
			for (k = 1; k < 5; k++)
				vij[k][i][j] = rho[k][i][j] = 0.;
		}

	while (fgets(SS, 200, fin) != NULL)
	{
		if (SS[0] == '#')
			printf("%s", SS);
		else
		{

			sscanf(SS, "%d %d %d %f", &i, &j, &k, &tmp);
			vij[i][j][k] = vij[i][k][j] = tmp;
			switch (i)
			{
			case 1:
				vij1[j][k] = vij1[k][j] = tmp;
				break;
			case 2:
				vij2[j][k] = vij2[k][j] = tmp;
				break;
			case 3:
				vij3[j][k] = vij3[k][j] = tmp;
				break;
			}
		}
	}

	fclose(fin);
}

void AppVacancy::read_pot_file(char *arg)
{
	char SS[100];
	char S1[100];
	char S2[100];
	char S3[100];
	FILE *fpot;
	float r_min, r_max, dr, value;
	int it1, it2, npt;
	int i, j, k;
	int scale_embed, scale_pair, scale_sembed;
	int ok;

	double PAIR[NPT_EAM];
	double RRRR[NPT_EAM];

	scale_pair = scale_embed = scale_sembed = 1;

	if ((fpot = fopen(arg, "r")) == NULL)
	{
		error->all(FLERR, "Open potentiel file failed.");
	}

	/*--Initialisation de pair et densite a 0. */
	for (i = 0; i < NSPECIES; i++)
		for (j = 0; j < NSPECIES; j++)
		{
			vij1[i][j] = vij2[i][j] = vij3[i][j] = vij4[i][j] = 0.;
			rho1[i][j] = rho2[i][j] = rho3[i][j] = rho4[i][j] = 0.;
			srho1[i][j] = srho2[i][j] = srho3[i][j] = srho4[i][j] = 0.;
			for (k = 0; k < 5; k++)
				vij[k][i][j] = rho[k][i][j] = srho[k][i][j] = 0.;
		}

	while (fgets(SS, 100, fpot) != NULL)
	{

		if (SS[0] != '*')
		{

			sscanf(SS, "%s %s", S1, S2);
			npt = it1 = 0;

			/*--         Potential pair it1 it2 npt r_min  r_max */
			if ((S2[0] == 'P') || (S2[0] == 'p'))
			{
				sscanf(SS, "%s %s %d %d %d %g %g", S1, S2, &it1, &it2, &npt, &r_min, &r_max);
				if (npt >= NPT_EAM)
				{
					printf("ERROR: parameter NPT_EAM in vlkmc.h file too small\n");
					printf("       NPT_EAM must be >= %d\n", npt + 1);
					exit(1);
				}
				for (i = 1; i <= npt; i++)
				{
					ok = fscanf(fpot, "%g", &value);
					PAIR[i] = value;
					RRRR[i] = r_min + (double)(i - 1) * (r_max - r_min) / (double)(npt - 1);
				}
				ok = fscanf(fpot, "\n");
				value = interpol(npt, RRRR, PAIR, a * 0.5 * sqrt(3.0));
				vij1[it1][it2] = vij1[it2][it1] = value;
				vij[1][it1][it2] = vij[1][it2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a);
				vij2[it1][it2] = vij2[it2][it1] = value;
				vij[2][it1][it2] = vij[2][it2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a * sqrt(2.0));
				vij3[it1][it2] = vij3[it2][it1] = value;
				vij[3][it1][it2] = vij[3][it2][it1] = value;
				if (r_max > a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5))
				{
					value = interpol(npt, RRRR, PAIR, a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5));
					vij4[it1][it2] = vij4[it2][it1] = value;
					vij[4][it1][it2] = vij[4][it2][it1] = value;
				}
			}

			/*--         Potential dens it1 npt r_min  r_max */
			if ((S2[0] == 'D') || (S2[0] == 'd'))
			{
				sscanf(SS, "%s %s %d %d %g %g", S1, S2, &it1, &npt, &r_min, &r_max);
				if (npt >= NPT_EAM)
				{
					printf("ERROR: parameter NPT_EAM in vlkmc.h file too small\n");
					printf("       NPT_EAM must be >= %d\n", npt + 1);
					exit(1);
				}
				for (i = 1; i <= npt; i++)
				{
					ok = fscanf(fpot, "%g", &value);
					PAIR[i] = value;
					RRRR[i] = r_min + (double)(i - 1) * (r_max - r_min) / (double)(npt - 1);
				}
				ok = fscanf(fpot, "\n");
				value = interpol(npt, RRRR, PAIR, a * 0.5 * sqrt(3.0));
				rho1[1][it1] = rho1[2][it1] = value;
				rho[1][1][it1] = rho[1][2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a);
				rho2[1][it1] = rho2[2][it1] = value;
				rho[2][1][it1] = rho[2][2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a * sqrt(2.0));
				rho3[1][it1] = rho3[2][it1] = value;
				rho[3][1][it1] = rho[3][2][it1] = value;
				if (r_max > a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5))
				{
					value = interpol(npt, RRRR, PAIR, a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5));
					rho4[1][it1] = rho4[2][it1] = value;
					rho[4][1][it1] = rho[4][2][it1] = value;
				}
			}

			/*--         Potential sdens it1 it2 npt r_min  r_max */
			if (((S2[0] == 'S') || (S2[0] == 's')) && ((S2[1] == 'D') || (S2[1] == 'd')))
			{
				sscanf(SS, "%s %s %d %d %d %g %g", S1, S2, &it1, &it2, &npt, &r_min, &r_max);
				if (npt >= NPT_EAM)
				{
					printf("ERROR: parameter NPT_EAM in vlkmc.h file too small\n");
					printf("       NPT_EAM must be >= %d\n", npt + 1);
					exit(1);
				}
				for (i = 1; i <= npt; i++)
				{
					ok = fscanf(fpot, "%g", &value);
					PAIR[i] = value;
					RRRR[i] = r_min + (double)(i - 1) * (r_max - r_min) / (double)(npt - 1);
				}
				ok = fscanf(fpot, "\n");
				value = interpol(npt, RRRR, PAIR, a * 0.5 * sqrt(3.0));
				srho1[it1][it2] = srho1[it2][it1] = value;
				srho[1][it1][it2] = srho[1][it2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a);
				srho2[it1][it2] = srho2[it2][it1] = value;
				srho[2][it1][it2] = srho[2][it2][it1] = value;
				value = interpol(npt, RRRR, PAIR, a * sqrt(2.0));
				srho3[it1][it2] = srho3[it2][it1] = value;
				srho[3][it1][it2] = srho[3][it2][it1] = value;
				if (r_max > a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5))
				{
					value = interpol(npt, RRRR, PAIR, a * sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5));
					srho4[it1][it2] = srho4[it2][it1] = value;
					srho[4][it1][it2] = srho[4][it2][it1] = value;
				}
			}

			/*--         Potential embed it1 npt r_min  r_max */
			if ((S2[0] == 'E') || (S2[0] == 'e'))
			{
				sscanf(SS, "%s %s %d %d %g %g", S1, S2, &it1, &npt, &r_min, &r_max);
				if (npt >= NPT_EMBED)
				{
					printf("ERROR: parameter NPT_EMBED in vlkmc.h file too small\n");
					printf("       NPT_EMBED must be >= %d\n", npt + 1);
					exit(1);
				}
				FEMBED[it1][0] = npt;
				FEMBED_dr[it1] = (r_max - r_min) / (double)(npt - 1);
				FEMBED_rmin[it1] = r_min;
				FEMBED_rmax[it1] = r_max;
				for (i = 1; i <= npt; i++)
				{
					ok = fscanf(fpot, "%g", &value);
					FEMBED[it1][i] = value;
				}
				ok = fscanf(fpot, "\n");
			}

			/*--         Potential embed it1 npt r_min  r_max */
			if (((S2[0] == 'S') || (S2[0] == 's')) && ((S2[1] == 'E') || (S2[1] == 'e')))
			{
				sscanf(SS, "%s %s %d %d %g %g", S1, S2, &it1, &npt, &r_min, &r_max);
				if (npt >= NPT_EMBED)
				{
					printf("ERROR: parameter NPT_EMBED in vlkmc.h file too small\n");
					printf("       NPT_EMBED must be >= %d\n", npt + 1);
					exit(1);
				}
				SEMBED[it1][0] = npt;
				SEMBED_dr[it1] = (r_max - r_min) / (double)(npt - 1);
				SEMBED_rmin[it1] = r_min;
				SEMBED_rmax[it1] = r_max;
				for (i = 1; i <= npt; i++)
				{
					ok = fscanf(fpot, "%g", &value);
					SEMBED[it1][i] = value;
				}
				ok = fscanf(fpot, "\n");
			}

			/*--         Potential scale value {pair|embed} */
			if ((S2[0] == 'S') || (S2[0] == 's'))
			{
				sscanf(SS, "%s %s %g %s", S1, S2, &value, S3);
				if ((S3[0] == 'E') || (S3[0] == 'e'))
				{
					scale_embed = 0;
				}
				if ((S3[0] == 'S') || (S3[0] == 's'))
				{
					scale_sembed = 0;
				}
				if ((S3[0] == 'P') || (S3[0] == 'p'))
				{
					scale_pair = 0;
				}
			}
		}

	} /* while */
	fclose(fpot);

	/*--Conversion (eventuelle) pair K -> eV */
	if (scale_pair == 1)
	{
		for (i = 1; i <= 2; i++)
			for (j = 1; j <= 2; j++)
			{
				for (k = 1; k <= 4; k++)
					vij[k][i][j] = vij[k][i][j] / 11604.49987;
				vij1[i][j] = vij1[i][j] / 11604.49987;
				vij2[i][j] = vij2[i][j] / 11604.49987;
				vij3[i][j] = vij3[i][j] / 11604.49987;
				vij4[i][j] = vij4[i][j] / 11604.49987;
			}
	}

	/*--Conversion (eventuelle) embed K -> eV */
	if (scale_embed == 1)
		for (i = 1; i <= 2; i++)
			for (j = 1; j <= FEMBED[i][0]; j++)
				FEMBED[i][j] = FEMBED[i][j] / 11604.49987;

	/*--Conversion (eventuelle) sembed K -> eV */
	if (scale_sembed == 1 && tag_COHESIVE_MODEL == CM_EAM3)
		for (i = 1; i <= 2; i++)
			for (j = 1; j <= SEMBED[i][0]; j++)
				SEMBED[i][j] = SEMBED[i][j] / 11604.49987;

	/*-- EAM potential with  RHO_AB = SQRT( RHO_A * RHO_B ) --*/

	if (tag_COHESIVE_MODEL == CM_FS)
	{
		it1 = 1;
		it2 = 2;
		value = sqrt(rho1[it1][it1] * rho1[it2][it2]);
		rho1[it1][it2] = rho1[it2][it1] = value;
		rho[1][it1][it2] = rho[1][it2][it1] = value;
		value = sqrt(rho2[it1][it1] * rho2[it2][it2]);
		rho2[it1][it2] = rho2[it2][it1] = value;
		rho[2][it1][it2] = rho[2][it2][it1] = value;
		value = sqrt(rho3[it1][it1] * rho3[it2][it2]);
		rho3[it1][it2] = rho3[it2][it1] = value;
		rho[3][it1][it2] = rho[3][it2][it1] = value;
		value = sqrt(rho4[it1][it1] * rho4[it2][it2]);
		rho4[it1][it2] = rho4[it2][it1] = value;
		rho[4][it1][it2] = rho[4][it2][it1] = value;
	}
}

double AppVacancy::interpol(int npt, double RRRR[], double PAIR[], double r)
{
	int i;
	double value;

	i = 1;

	while ((r > RRRR[i]) && (i <= npt))
	{
		i++;
	}
	i--;

	if (i == npt)
	{
		value = 0.0;
	}
	else
	{
		value = PAIR[i] + (PAIR[i + 1] - PAIR[i]) * (r - RRRR[i]) / (RRRR[i + 1] - RRRR[i]);
	}
	return (value);
}

/*--------------------------------------------------------------------*/
/*-- Evaluation fonction entourage pour atome de type it et densite dens */

double AppVacancy::eval_embed(int it, double dens)
{

	double value;
	int i;

	if (dens < FEMBED_rmin[it])
	{
		char s[120];
		if (dens < 0)
		{
			sprintf(s, "Error in eval_embed: dens %g less than 0\n", dens);
			error->all(FLERR, s);
		}
		else
		{
			i = 1;					 // calculate slope at FEMBED_rmin[it],slope k = (FEMBED[it][2]-FEMBED[it][1])/FEMBED_dr[it],so i = 1
			value = FEMBED_rmin[it]; // r(1)
			value = FEMBED[it][i] + (FEMBED[it][i + 1] - FEMBED[it][i]) * (dens - value) / FEMBED_dr[it];
			return value;
		}
	}
	else if (dens > FEMBED_rmax[it])
	{
		char s[120];
		sprintf(s, "Error in eval_embed: dens too large %g (limit %g)\n", dens, FEMBED_rmax[it]);
		error->all(FLERR, s);
	}

	value = (dens - FEMBED_rmin[it]) / FEMBED_dr[it]; // bdwu 1129
													  // when r_min < r < r_max, calculate how many interval strides between dens and FEMBED_rmin[it]
	i = 1 + (int)value;								  // get the index of start point of the inerval where dens is

	value = FEMBED_rmin[it] + (double)(i - 1.) * FEMBED_dr[it]; // get the value of start point of the inerval where dens is  -->  r(i)

	value = FEMBED[it][i] + (FEMBED[it][i + 1] - FEMBED[it][i]) * (dens - value) / FEMBED_dr[it];

	return (value);
}

/*----------------------------------------------------------------*/
/*-- edapted from eval_embed() subroutine ---*/

double AppVacancy::eval_sembed(int it, double sdens)
{
	double value;
	int i;

	if (sdens < SEMBED_rmin[it])
	{
		char s[120];
		sprintf(s, "me:%d Error in eval_sembed: sdens too small %g (limit %g)\n", me, sdens, SEMBED_rmin[it]);
		error->all(FLERR, s);
	}
	else if (sdens > SEMBED_rmax[it])
	{
		char s[120];
		sprintf(s, "Error in eval_sembed: sdens too large %g (limit %g)\n", sdens, SEMBED_rmax[it]);
		error->all(FLERR, s);
	}

	value = (sdens - SEMBED_rmin[it]) / SEMBED_dr[it];
	i = 1 + (int)value;

	value = SEMBED_rmin[it] + (double)(i - 1.) * SEMBED_dr[it];

	value = SEMBED[it][i] +
			(SEMBED[it][i + 1] - SEMBED[it][i]) * (sdens - value) / SEMBED_dr[it];
	if (i < 1)
		printf("sembed %d  %g %g\n", i, sdens, value);

	return (value);
}

void AppVacancy::eatom(int it, int *it_nn1, int *it_nn2, int *it_nn3,
					   double *ev, double *er, double *es)
{
	double e_v, e_r, e_s;
	int i, j, k, l;
	e_v = e_r = e_s = 0.0;
	for (i = 0; i < 8; i++)
	{
		e_v = e_v + vij1[it][it_nn1[i]];
		e_r = e_r + rho1[it][it_nn1[i]];
	}
	for (j = 0; j < 6; j++)
	{
		e_v = e_v + vij2[it][it_nn2[j]];
		e_r = e_r + rho2[it][it_nn2[j]];
	}
	if (tag_COHESIVE_MODEL != CM_PPAIR2)
	{
		for (k = 0; k < 12; k++)
		{
			e_v = e_v + vij3[it][it_nn3[k]];
			e_r = e_r + rho3[it][it_nn3[k]];
		}
	}
	if (tag_COHESIVE_MODEL == CM_EAM3)
	{
		for (i = 0; i < 8; i++)
		{
			e_s = e_s + srho1[it][it_nn1[i]];
		}
		for (j = 0; j < 6; j++)
		{
			e_s = e_s + srho2[it][it_nn2[j]];
		}
		for (k = 0; k < 12; k++)
		{
			e_s = e_s + srho3[it][it_nn3[k]];
		}
	}

	*ev = e_v;
	*er = e_r;
	*es = e_s;
}

double AppVacancy::single_energy(int it, double e_v, double e_r, double e_s)
{
	double ene = 0.0;

	switch (tag_COHESIVE_MODEL)
	{
	case CM_EAM:
		if (it > 0)
			ene = 0.5 * e_v + eval_embed(it, e_r);
		else
			ene = 0.;
		break;
	case CM_FS:
		if (it > 0)
			ene = 0.5 * e_v - sqrt(e_r);
		else
			ene = 0.;
		break;
	case CM_EAM3:
		if (it > 0)
			ene = 0.5 * e_v + eval_embed(it, e_r) + eval_sembed(it, e_s);
		else
			ene = 0.;
		break;
	case CM_PPAIR2:
		ene = 0.5 * e_v;
		break;
	}
	return ene;
}

/*  it---type of jumpsite */
double AppVacancy::delete1(int it, double *e_v, double *e_r, double *e_s)
{
	*e_v = *e_v - vij1[it][it];
	*e_r = *e_r - rho1[it][it];
	*e_s = *e_s - srho1[it][it];
}

double AppVacancy::calcul_dene_v1(int it, int ii, int jj, int kk, int nn, double sign)
{
	int pos = get_POS_ID(ii, jj, kk);
	int it2 = lattice[pos];
	double er2 = E_R[pos];
	double dene = 0.;

	switch (tag_COHESIVE_MODEL)
	{
	case CM_EAM:
		if (it2 > 0)
			dene = 0.5 * sign * vij[nn][it2][it] - eval_embed(it2, er2) + eval_embed(it2, er2 + sign * rho[nn][it2][it]);
		break;
	case CM_FS:
		if (it2 > 0)
			dene = 0.5 * sign * vij[nn][it][it2] + sqrt(er2) - sqrt(er2 + sign * rho[nn][it][it2]);
		break;
	case CM_EAM3:
		if (it2 > 0)
		{
			double es2 = E_S[pos];
			dene = 0.5 * sign * vij[nn][it2][it] - eval_embed(it2, er2) +
				   eval_embed(it2, er2 + sign * rho[nn][it2][it]) -
				   eval_sembed(it2, es2) +
				   eval_sembed(it2, es2 + sign * srho[nn][it2][it]);
		}
		break;
	case CM_PPAIR2:
		dene = 0.5 * sign * (vij[nn][it][it2] - vij[nn][0][it2]);
		break;
	default:
		printf("not supported dene_v1\n");
		break;
	}

	return (dene);
}

double AppVacancy::calcul_dene_v2(int it, int ii, int jj, int kk, int nn1, int nn2)
{

	double dene = 0.;
	int pos = get_POS_ID(ii, jj, kk);
	int it2 = lattice[pos];
	double er2 = E_R[pos];

	switch (tag_COHESIVE_MODEL)
	{
	case CM_EAM:
		if (it2 > 0)
			dene = 0.5 * (vij[nn1][it2][it] - vij[nn2][it2][it]) -
				   eval_embed(it2, er2) +
				   eval_embed(it2, er2 + rho[nn1][it2][it] - rho[nn2][it2][it]);
		break;
	case CM_FS:
		if (it2 > 0)
			dene = 0.5 * (vij[nn1][it][it2] - vij[nn2][it][it2]) +
				   sqrt(er2) -
				   sqrt(er2 + rho[nn1][it][it2] - rho[nn2][it][it2]);
		break;
	case CM_EAM3:
		if (it2 > 0)
		{
			double es2 = E_S[pos];
			dene = 0.5 * (vij[nn1][it2][it] - vij[nn2][it2][it]) -
				   eval_embed(it2, er2) +
				   eval_embed(it2, er2 + rho[nn1][it2][it] - rho[nn2][it2][it]) -
				   eval_sembed(it2, es2) +
				   eval_sembed(it2, es2 + srho[nn1][it2][it] - srho[nn2][it2][it]);
		}
		break;

	case CM_PPAIR2:
		dene = 0.5 * (vij[nn1][it][it2] - vij[nn2][it][it2] - vij[nn1][0][it2] + vij[nn2][0][it2]);
		break;
	default:
		printf("not supported dene_v2\n");
		break;
	}
	return (dene);
}
/*-----------------------------------------------------------*/

void AppVacancy::do_jump_v1(int it, int ii, int jj, int kk, int nn, double sign)
{
	int pos = get_POS_ID(ii, jj, kk);
	int it2 = lattice[pos];
	switch (tag_COHESIVE_MODEL)
	{
	case CM_EAM:
		E_V[pos] += sign * vij[nn][it][it2];
		E_R[pos] += sign * rho[nn][it2][it];
		break;
	case CM_FS:
		E_V[pos] += sign * vij[nn][it][it2];
		E_R[pos] += sign * rho[nn][it2][it];
		break;
	case CM_EAM3:
		E_V[pos] += sign * vij[nn][it][it2];
		E_R[pos] += sign * rho[nn][it2][it];
		E_S[pos] += sign * srho[nn][it2][it];
		break;
	case CM_PPAIR2:
		E_V[pos] += sign * (vij[nn][it][it2] - vij[nn][0][it2]);
		break;
	default:
		printf("not supported do_jump_v1\n");
		break;
	}
}

void AppVacancy::do_jump_v2(int it, int ii, int jj, int kk, int nn1, int nn2)
{
	int pos = get_POS_ID(ii, jj, kk);
	int it2 = lattice[pos];
	switch (tag_COHESIVE_MODEL)
	{
	case CM_EAM:
		E_V[pos] += vij[nn1][it][it2] - vij[nn2][it][it2];
		E_R[pos] += rho[nn1][it2][it] - rho[nn2][it2][it];
		break;
	case CM_FS:
		E_V[pos] += vij[nn1][it][it2] - vij[nn2][it][it2];
		E_R[pos] += rho[nn1][it2][it] - rho[nn2][it2][it];
		break;
	case CM_EAM3:
		E_V[pos] += vij[nn1][it][it2] - vij[nn2][it][it2];
		E_R[pos] += rho[nn1][it2][it] - rho[nn2][it2][it];
		E_S[pos] += srho[nn1][it2][it] - srho[nn2][it2][it];
		break;
	case CM_PPAIR2:
		E_V[pos] += vij[nn1][it][it2] - vij[nn2][it][it2] - vij[nn1][0][it2] + vij[nn2][0][it2];
		break;
	}
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppVacancy::allocate_data()
{
	memory->create(barrier, TOP - VACANT, "app:barrier");
	for (int i = 0; i < TOP; i++)
		barrier[i] = 0.68;

	memory->create(firstevent, nlocal, "app:firstevent");
}

int AppVacancy::is_in_set(int i)
{
	int iwhich, jwhich, kwhich, msector;

	if (xyz[i][0] < xmid)
		iwhich = 0;
	else
		iwhich = 1;
	if (xyz[i][1] < ymid)
		jwhich = 0;
	else
		jwhich = 1;
	if (xyz[i][2] < zmid)
		kwhich = 0;
	else
		kwhich = 1;

	if (nset == 2)
		msector = iwhich + 1;
	else if (nset == 4)
		msector = 2 * jwhich + iwhich + 1;
	else
		msector = 4 * kwhich + 2 * jwhich + iwhich + 1;

	if (i > nlocal || current_set + 1 != msector)
		return 0;
	else
		return 1;
}

int AppVacancy::in_which_set(int i)
{
	int iwhich, jwhich, kwhich, msector;

	if (xyz[i][0] < xmid)
		iwhich = 0;
	else
		iwhich = 1;
	if (xyz[i][1] < ymid)
		jwhich = 0;
	else
		jwhich = 1;
	if (xyz[i][2] < zmid)
		kwhich = 0;
	else
		kwhich = 1;

	if (nset == 2)
		msector = iwhich;
	else if (nset == 4)
		msector = 2 * jwhich + iwhich;
	else
		msector = 4 * kwhich + 2 * jwhich + iwhich;

	return msector;
}

void AppVacancy::update_site_propensity(cvector_vector_type(int) vacant,
										double *propensity,
										double **propensity_dir,
										int **jumpsite_dir, int *num_dir)
{
	int vac_size = cvector_size((vacant));

	argOfSite_propensity.c_vacant = &vacant[0];
	argOfSite_propensity.vac_size = vac_size;
	int i, j;
	for (i = 0; i < 6; i++)
	{
		argOfSite_propensity.c_barrier[i] = barrier[i];
		for (j = 0; j < 6; j++)
		{
			argOfSite_propensity.vij1[i][j] = vij1[i][j];
			argOfSite_propensity.vij2[i][j] = vij2[i][j];
		}
	}
	argOfSite_propensity.xyz = xyz;
	argOfSite_propensity.unitdis = unitdis;
	argOfSite_propensity.lattice = lattice;
	argOfSite_propensity.E_V = E_V;
	argOfSite_propensity.boxxlo = boxxlo;
	argOfSite_propensity.boxxhi = boxxhi;
	argOfSite_propensity.subxlo = subxlo;
	argOfSite_propensity.subxhi = subxhi;
	argOfSite_propensity.boxylo = boxylo;
	argOfSite_propensity.boxyhi = boxyhi;
	argOfSite_propensity.subylo = subylo;
	argOfSite_propensity.subyhi = subyhi;
	argOfSite_propensity.boxzlo = boxzlo;
	argOfSite_propensity.boxzhi = boxzhi;
	argOfSite_propensity.subzlo = subzlo;
	argOfSite_propensity.subzhi = subzhi;
	argOfSite_propensity.xghost = xghost;
	argOfSite_propensity.yghost = yghost;
	argOfSite_propensity.zghost = zghost;
	argOfSite_propensity.nlocal = nlocal;
	argOfSite_propensity.t_inverse = t_inverse;
	argOfSite_propensity.propensity = propensity;
	argOfSite_propensity.propensity_dir = propensity_dir;
	argOfSite_propensity.jumpsite_dir = jumpsite_dir;
	argOfSite_propensity.num_dir = num_dir;

	kernel_site_propensity(&argOfSite_propensity);
}

void AppVacancy::call_athread_join()
{
	c_call_athread_join();
}

void AppVacancy::set_arg_of_kernel_func(int *index, int *nsite,
										int isite, int jsite,
										int *vacant, int vac_size)
{

	argOfKernel.c_vacant = vacant;
	argOfKernel.t_inverse = t_inverse;
	argOfKernel.vac_size = vac_size;
	argOfKernel.site1[0] = round(xyz[isite][0] * unitdis);
	argOfKernel.site1[1] = round(xyz[isite][1] * unitdis);
	argOfKernel.site1[2] = round(xyz[isite][2] * unitdis);
	argOfKernel.site2[0] = round(xyz[jsite][0] * unitdis);
	argOfKernel.site2[1] = round(xyz[jsite][1] * unitdis);
	argOfKernel.site2[2] = round(xyz[jsite][2] * unitdis);
	int i, j;
	for (i = 0; i < 6; i++)
	{
		argOfKernel.c_barrier[i] = barrier[i];
		for (j = 0; j < 6; j++)
		{
			argOfKernel.vij1[i][j] = vij1[i][j];
			argOfKernel.vij2[i][j] = vij2[i][j];
		}
	}
	argOfKernel.nsite = nsite;
	argOfKernel.index = index;
	argOfKernel.boxx = boxx;
	argOfKernel.boxy = boxy;
	argOfKernel.boxz = boxz;
	argOfKernel.boxxlo = boxxlo;
	argOfKernel.boxxhi = boxxhi;
	argOfKernel.subxlo = subxlo;
	argOfKernel.subxhi = subxhi;
	argOfKernel.boxylo = boxylo;
	argOfKernel.boxyhi = boxyhi;
	argOfKernel.subylo = subylo;
	argOfKernel.subyhi = subyhi;
	argOfKernel.boxzlo = boxzlo;
	argOfKernel.boxzhi = boxzhi;
	argOfKernel.subzlo = subzlo;
	argOfKernel.subzhi = subzhi;
	argOfKernel.unitdis = unitdis;
	argOfKernel.xyz = xyz;
	argOfKernel.lattice = lattice;
	argOfKernel.E_V = E_V;
	argOfKernel.propensity = propensity;
	argOfKernel.propensity_dir = propensity_dir;
	argOfKernel.jumpsite_dir = jumpsite_dir;
	argOfKernel.num_dir = num_dir;
	argOfKernel.xghost = xghost;
	argOfKernel.yghost = yghost;
	argOfKernel.zghost = zghost;
	argOfKernel.nlocal = nlocal;
}
