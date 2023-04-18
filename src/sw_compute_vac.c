#include <stdio.h>
#include <math.h>
#include <string.h>
#include "crts.h"
#include "swarg.h"
#include <simd.h>

__thread_local_share int nsite;
__thread_local_share int share_k_site[100][3];
__thread_local_share int share_k[100];
__thread_local volatile unsigned long get_reply, put_reply;
__thread_local volatile int my_id;
__thread_local int subxlo, subxhi, subylo, subyhi, subzlo, subzhi;
__thread_local int boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi;
__thread_local int xghost, yghost, zghost;
__thread_local int nlocal;

static inline int get_POS_ID(int x, int y, int z)
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
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) + x2 + nlocal;
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
			  z1 * (subxhi - subxlo + xghost * 2) / 2 * (subyhi - subylo + yghost * 2) + x2;
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
			  x2 - (y1 - yghost / 2) * (subxhi - subxlo) -
			  (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 + nlocal;
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
			  x2 - (y1 - yghost / 2 + 1) * (subxhi - subxlo) -
			  (z1 - zghost / 2) * (subxhi - subxlo) * (subyhi - subylo) / 2 + nlocal;
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

	// printf("pos : %d \n", pos);
	return pos;
}

void c_inc_px(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxxhi)
		*ii = *ii - boxxhi + boxxlo;
	else if (*ii < boxxlo)
		*ii = boxxhi + *ii - boxxlo;
}

void c_inc_py(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxyhi)
		*ii = *ii - boxyhi + boxylo;
	else if (*ii < boxylo)
		*ii = boxyhi + *ii - boxylo;
}

void c_inc_pz(int i, int di, int *ii)
{
	*ii = i + di;
	if (*ii >= boxzhi)
		*ii = *ii - boxzhi + boxzlo;
	else if (*ii < boxzlo)
		*ii = boxzhi + *ii - boxzlo;
}

void c_box2sub(int *x, int *y, int *z)
{
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

double c_calcul_de(double **xyz, double unitdis, int *initsite,
				   int di, int dj, int dk, int *jumpsite, int *lattice,
				   double *E_V, double vij1[6][6], double vij2[6][6])
{
	int i, j, k;
	i = initsite[0];
	j = initsite[1];
	k = initsite[2];

	int di_p1, dj_p1, dk_p1, di_m1, dj_m1, dk_m1;
	int di_m2, di_p2, dj_m2, dj_p2, dk_m2, dk_p2;
	int new_di_p1, new_dj_p1, new_dk_p1, new_di_m1, new_dj_m1, new_dk_m1;
	int new_di_m2, new_di_p2, new_dj_m2, new_dj_p2, new_dk_m2, new_dk_p2;
	int it;

	c_inc_px(i, di, &di_p1);
	c_inc_py(j, dj, &dj_p1);
	c_inc_pz(k, dk, &dk_p1);
	c_box2sub(&di_p1, &dj_p1, &dk_p1);

	int pos = get_POS_ID(di_p1, dj_p1, dk_p1);
	it = lattice[pos];
	*jumpsite = pos;

	if (it == 0)
	{
		return 0.0;
	}

	c_inc_px(i, -di, &di_m1);
	c_inc_py(j, -dj, &dj_m1);
	c_inc_pz(k, -dk, &dk_m1);

	c_inc_px(i, 2 * di, &di_p2);
	c_inc_py(j, 2 * dj, &dj_p2);
	c_inc_pz(k, 2 * dk, &dk_p2);
	c_inc_px(i, -2 * di, &di_m2);
	c_inc_py(j, -2 * dj, &dj_m2);
	c_inc_pz(k, -2 * dk, &dk_m2);

	c_box2sub(&i, &j, &k);
	c_box2sub(&di_m1, &dj_m1, &dk_m1);
	c_box2sub(&di_p2, &dj_p2, &dk_p2);
	c_box2sub(&di_m2, &dj_m2, &dk_m2);

	double dene = 0.;

	// 1. energy of inital : Ei
	dene += -0.5 * E_V[get_POS_ID(i, j, k)] -
			0.5 * E_V[get_POS_ID(di_p1, dj_p1, dk_p1)];

	// 2. energy of the final : Ef
	double e_v;
	int it_nn1[8], it_nn2[6];

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

	e_v = 0.0;
	for (i = 0; i < 8; i++)
	{
		e_v = e_v + vij1[it][it_nn1[i]];
	}
	for (j = 0; j < 6; j++)
	{
		e_v = e_v + vij2[it][it_nn2[j]];
	}
	dene = dene + 0.5 * e_v;

	i = initsite[0];
	j = initsite[1];
	k = initsite[2];

	c_inc_px(i, di, &di_p1);
	c_inc_py(j, dj, &dj_p1);
	c_inc_pz(k, dk, &dk_p1);

	c_inc_px(di_p1, di, &new_di_p1);
	c_inc_py(dj_p1, dj, &new_dj_p1);
	c_inc_pz(dk_p1, dk, &new_dk_p1);
	c_inc_px(di_p1, -di, &new_di_m1);
	c_inc_py(dj_p1, -dj, &new_dj_m1);
	c_inc_pz(dk_p1, -dk, &new_dk_m1);

	c_inc_px(di_p1, 2 * di, &new_di_p2);
	c_inc_py(dj_p1, 2 * dj, &new_dj_p2);
	c_inc_pz(dk_p1, 2 * dk, &new_dk_p2);
	c_inc_px(di_p1, -2 * di, &new_di_m2);
	c_inc_py(dj_p1, -2 * dj, &new_dj_m2);
	c_inc_pz(dk_p1, -2 * dk, &new_dk_m2);

	c_box2sub(&di_p1, &dj_p1, &dk_p1);
	c_box2sub(&new_di_p1, &new_dj_p1, &new_dk_p1);
	c_box2sub(&new_di_m1, &new_dj_m1, &new_dk_m1);
	c_box2sub(&new_di_p2, &new_dj_p2, &new_dk_p2);
	c_box2sub(&new_di_m2, &new_dj_m2, &new_dk_m2);

	it_nn1[0] = it;
	it_nn1[1] = it_nn2[3];
	it_nn1[2] = it_nn2[1];
	it_nn1[4] = it_nn2[5];

	it_nn2[0] = it_nn1[5];
	it_nn2[2] = it_nn1[6];
	it_nn2[4] = it_nn1[3];

	it_nn1[3] = lattice[get_POS_ID(new_di_p1, new_dj_p1, new_dk_m1)];
	it_nn1[5] = lattice[get_POS_ID(new_di_m1, new_dj_p1, new_dk_p1)];
	it_nn1[6] = lattice[get_POS_ID(new_di_p1, new_dj_m1, new_dk_p1)];
	it_nn1[7] = lattice[get_POS_ID(new_di_p1, new_dj_p1, new_dk_p1)];
	it_nn2[1] = lattice[get_POS_ID(new_di_p2, dj_p1, dk_p1)];
	it_nn2[3] = lattice[get_POS_ID(di_p1, new_dj_p2, dk_p1)];
	it_nn2[5] = lattice[get_POS_ID(di_p1, dj_p1, new_dk_p2)];

	it = 0;
	e_v = 0.0;
	for (i = 0; i < 8; i++)
	{
		e_v = e_v + vij1[it][it_nn1[i]];
	}
	for (j = 0; j < 6; j++)
	{
		e_v = e_v + vij2[it][it_nn2[j]];
	}
	dene = dene + 0.5 * e_v;

	return dene * 2.0;
}

double c_site_propensity(double *c_barrier, double **xyz, double unitdis,
						 int *i, int *lattice, double *E_V, double vij1[6][6], double vij2[6][6],
						 double t_inverse,
						 double *vac_propensity_dir, int *vac_jumpsite_dir, int *num_dir)
{

#define DELTAEVENT 100000
#define kb (1.38066E-23 / 1.6E-19)
#define nu 6.E12

	double probone, proball, dene = 0;
	int jumpsite;
	int num;
	proball = 0.0;
	num = 0;
	int di, dj, dk;
	for (di = -1; di <= 1; di += 2)
	{
		for (dj = -1; dj <= 1; dj += 2)
		{
			for (dk = -1; dk <= 1; dk += 2)
			{
				dene = c_calcul_de(xyz, unitdis, i, di, dj, dk, &jumpsite,
								   lattice, E_V, vij1, vij2); // change of energy

				if (lattice[jumpsite] == 0)
				{
					continue;
				}
				probone = nu * exp(-(0.5 * dene + c_barrier[lattice[jumpsite]]) * t_inverse / kb);
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

void sw_site_propensity(struct _sw_site_propensity_arg *marg)
{
	my_id = CRTS_smng_get_tid();
	struct _sw_site_propensity_arg sarg;

	int sw_k_site[3];
	int sw_nsite;
	int sw_nsite_st;
	int sw_task;

	get_reply = 0;
	CRTS_dma_iget(&sarg, marg, sizeof(struct _sw_site_propensity_arg), &get_reply);
	while (get_reply != 1)
		;
	int i, j, k;

	subxlo = sarg.subxlo;
	subxhi = sarg.subxhi;
	subylo = sarg.subylo;
	subyhi = sarg.subyhi;
	subzlo = sarg.subzlo;
	subzhi = sarg.subzhi;
	boxxlo = sarg.boxxlo;
	boxxhi = sarg.boxxhi;
	boxylo = sarg.boxylo;
	boxyhi = sarg.boxyhi;
	boxzlo = sarg.boxzlo;
	boxzhi = sarg.boxzhi;
	xghost = sarg.xghost;
	yghost = sarg.yghost;
	zghost = sarg.zghost;
	nlocal = sarg.nlocal;

	sw_nsite_st = nsite;

	long load = sarg.vac_size / NTHREAD;
	long rest = sarg.vac_size % NTHREAD;
	long stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest)
		load++;
	long edid = stid + load;

	long batch = load;
	while (batch > 0)
	{
		if (batch > DMAMAX)
			batch = DMAMAX;

		for (i = 0; i < batch; i++)
		{
			k = i + stid;
			int vacant = sarg.c_vacant[k];
			sw_k_site[0] = (int)round(sarg.xyz[vacant][0] * sarg.unitdis);
			sw_k_site[1] = (int)round(sarg.xyz[vacant][1] * sarg.unitdis);
			sw_k_site[2] = (int)round(sarg.xyz[vacant][2] * sarg.unitdis);

			sarg.propensity[k] =
				c_site_propensity(sarg.c_barrier, sarg.xyz, sarg.unitdis,
								  sw_k_site, sarg.lattice, sarg.E_V, sarg.vij1, sarg.vij2,
								  sarg.t_inverse, sarg.propensity_dir[k],
								  sarg.jumpsite_dir[k], sarg.num_dir + k);
		}
		stid += batch;
		batch = edid - stid;
	}
}

void sw_compute_vac(struct _swarg *marg)
{
	my_id = CRTS_smng_get_tid();
	struct _swarg sarg;
	int xyz1[3], xyz2[3];
	int sw_k_site[3];
	int dis1, dis2;
	int sw_nsite;
	int sw_nsite_st;
	int sw_task;

	get_reply = 0;
	CRTS_dma_iget(&sarg, marg, sizeof(struct _swarg), &get_reply);
	while (get_reply != 1)
		;
	int i, j, k;

	subxlo = sarg.subxlo;
	subxhi = sarg.subxhi;
	subylo = sarg.subylo;
	subyhi = sarg.subyhi;
	subzlo = sarg.subzlo;
	subzhi = sarg.subzhi;
	boxxlo = sarg.boxxlo;
	boxxhi = sarg.boxxhi;
	boxylo = sarg.boxylo;
	boxyhi = sarg.boxyhi;
	boxzlo = sarg.boxzlo;
	boxzhi = sarg.boxzhi;
	xghost = sarg.xghost;
	yghost = sarg.yghost;
	zghost = sarg.zghost;
	nlocal = sarg.nlocal;

	if (my_id == 0)
	{
		nsite = *(sarg.nsite);
	}

	CRTS_ssync_array();

	sw_nsite_st = nsite;

	long load = sarg.vac_size / NTHREAD;
	long rest = sarg.vac_size % NTHREAD;
	long stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest)
		load++;
	long edid = stid + load;

	long debug_stid = stid;

	long batch = load;

	intv16 Vx, Vy, Vz;
	intv16 x01, y01, z01, x02, y02, z02;
	x01 = sarg.site1[0];
	x02 = sarg.site2[0];
	y01 = sarg.site1[1];
	y02 = sarg.site2[1];
	z01 = sarg.site1[2];
	z02 = sarg.site2[2];
	intv16 sign;
	intv16 Vtemp1;
	intv16 Vtemp2;
	intv16 dx1, dx2, dy1, dy2, dz1, dz2;
	intv16 Vboxx = sarg.boxx;
	intv16 Vboxy = sarg.boxy;
	intv16 Vboxz = sarg.boxz;
	intv16 d1, d2;
	intv16 result;
	int sw_compute_vac[16];
	int sw_vac_site_int[3][32] __attribute__((aligned(512)));

	double prefetch_buffer[16][3];
	double prefetch_buffer2[16][3];

	for (int i = 0; i < 16 && i < batch; ++i)
	{
		int vac = sarg.c_vacant[stid + i];
		for (int j = 0; j < 3; j++)
		{
			prefetch_buffer[i][j] = sarg.xyz[vac][j];
		}
	}

	while (batch > 0)
	{
		if (batch > 16)
		{
			for (int i = 0; i < 16; i++)
			{

				for (int j = 0; j < 3; j++)
				{
					prefetch_buffer2[i][j] = prefetch_buffer[i][j];
				}
				if ((stid + i + 16) < edid)
				{
					int vac = sarg.c_vacant[stid + i + 16];
					for (int j = 0; j < 3; j++)
					{
						prefetch_buffer[i][j] = sarg.xyz[vac][j];
					}
				}
			}
			batch = 16;
		}
		else
		{
			for (int i = 0; i < batch; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					prefetch_buffer2[i][j] = prefetch_buffer[i][j];
				}
			}
		}

		if (batch == 16)
		{

			for (int j = 0; j < 16; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					sw_vac_site_int[i][j] = round(prefetch_buffer2[j][i] * sarg.unitdis);
				}
			}

			simd_load(Vx, &sw_vac_site_int[0][0]);
			simd_load(Vy, &sw_vac_site_int[1][0]);
			simd_load(Vz, &sw_vac_site_int[2][0]);

			dx1 = Vx - x01;
			Vtemp1 = simd_vsrlwi(dx1, 31);
			dx1 = simd_vxorw(Vtemp1, dx1 + Vtemp1);
			sign = simd_vcmpltw(Vboxx, 2 * dx1);
			dx1 = dx1 + sign * (Vboxx - 2 * dx1);

			dy1 = Vy - y01;
			Vtemp1 = simd_vsrlwi(dy1, 31);
			dy1 = simd_vxorw(Vtemp1, dy1 + Vtemp1);
			sign = simd_vcmpltw(Vboxy, 2 * dy1);
			dy1 = dy1 + sign * (Vboxy - 2 * dy1);

			dz1 = Vz - z01;
			Vtemp1 = simd_vsrlwi(dz1, 31);
			dz1 = simd_vxorw(Vtemp1, dz1 + Vtemp1);
			sign = simd_vcmpltw(Vboxz, 2 * dz1);
			dz1 = dz1 + sign * (Vboxz - 2 * dz1);

			dx2 = Vx - x02;
			Vtemp1 = simd_vsrlwi(dx2, 31);
			dx2 = simd_vxorw(Vtemp1, dx2 + Vtemp1);
			sign = simd_vcmpltw(Vboxx, 2 * dx2);
			dx2 = dx2 + sign * (Vboxx - 2 * dx2);

			dy2 = Vy - y02;
			Vtemp1 = simd_vsrlwi(dy2, 31);
			dy2 = simd_vxorw(Vtemp1, dy2 + Vtemp1);
			sign = simd_vcmpltw(Vboxy, 2 * dy2);
			dy2 = dy2 + sign * (Vboxy - 2 * dy2);

			dz2 = Vz - z02;
			Vtemp1 = simd_vsrlwi(dz2, 31);
			dz2 = simd_vxorw(Vtemp1, dz2 + Vtemp1);
			sign = simd_vcmpltw(Vboxz, 2 * dz2);
			dz2 = dz2 + sign * (Vboxz - 2 * dz2);

			d1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
			d2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
			result = simd_vcmpltwi(d1, 27) + simd_vcmpltwi(d2, 27);

			simd_store(result, &sw_compute_vac[0]);
			for (int i = 0; i < 16; ++i)
			{
				if (sw_compute_vac[i])
				{
					CRTS_smutex_lock_node();
					sw_nsite = nsite;
					nsite++;
					CRTS_smutex_unlock_node();
					share_k_site[sw_nsite - sw_nsite_st][0] = sw_vac_site_int[0][i];
					share_k_site[sw_nsite - sw_nsite_st][1] = sw_vac_site_int[1][i];
					share_k_site[sw_nsite - sw_nsite_st][2] = sw_vac_site_int[2][i];
					sarg.index[sw_nsite] = stid + i;
					share_k[sw_nsite - sw_nsite_st] = stid + i;
				}
			}
		}
		else
		{
			for (i = 0; i < batch; i++)
			{
				for (j = 0; j < 3; j++)
				{
					sw_k_site[j] = round(prefetch_buffer2[i][j] * sarg.unitdis);
					xyz1[j] = abs(sarg.site1[j] - sw_k_site[j]);
					xyz2[j] = abs(sarg.site2[j] - sw_k_site[j]);
				}

				if (xyz1[0] * 2 > sarg.boxx)
					xyz1[0] = sarg.boxx - xyz1[0];
				if (xyz1[1] * 2 > sarg.boxy)
					xyz1[1] = sarg.boxy - xyz1[1];
				if (xyz1[2] * 2 > sarg.boxz)
					xyz1[2] = sarg.boxz - xyz1[2];
				if (xyz2[0] * 2 > sarg.boxx)
					xyz2[0] = sarg.boxx - xyz2[0];
				if (xyz2[1] * 2 > sarg.boxy)
					xyz2[1] = sarg.boxy - xyz2[1];
				if (xyz2[2] * 2 > sarg.boxz)
					xyz2[2] = sarg.boxz - xyz2[2];

				dis1 = 0;
				dis2 = 0;
				for (j = 0; j < 3; j++)
				{
					dis1 = dis1 + xyz1[j] * xyz1[j];
					dis2 = dis2 + xyz2[j] * xyz2[j];
				}

				if (dis1 < 27 || dis2 < 27)
				{
					CRTS_smutex_lock_node();
					sw_nsite = nsite;
					nsite++;
					CRTS_smutex_unlock_node();
					share_k_site[sw_nsite - sw_nsite_st][0] = sw_k_site[0];
					share_k_site[sw_nsite - sw_nsite_st][1] = sw_k_site[1];
					share_k_site[sw_nsite - sw_nsite_st][2] = sw_k_site[2];
					sarg.index[sw_nsite] = stid + i;
					share_k[sw_nsite - sw_nsite_st] = stid + i;
				}
			}
		}
		stid += batch;
		batch = edid - stid;
	}
	CRTS_ssync_array();

	load = (nsite - sw_nsite_st) / NTHREAD;
	rest = (nsite - sw_nsite_st) % NTHREAD;
	stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest)
		load++;
	edid = stid + load;

	batch = load;
	while (batch > 0)
	{
		if (batch > DMAMAX2)
			batch = DMAMAX2;
		for (i = 0; i < batch; i++)
		{
			k = share_k[i + stid];
			sw_k_site[0] = share_k_site[i + stid][0];
			sw_k_site[1] = share_k_site[i + stid][1];
			sw_k_site[2] = share_k_site[i + stid][2];

			sarg.propensity[k] = c_site_propensity(sarg.c_barrier, sarg.xyz, sarg.unitdis,
												   sw_k_site, sarg.lattice, sarg.E_V, sarg.vij1, sarg.vij2,
												   sarg.t_inverse, sarg.propensity_dir[k],
												   sarg.jumpsite_dir[k], sarg.num_dir + k);
		}
		stid += batch;
		batch = edid - stid;
	}

	if (my_id == 0)
		*(sarg.nsite) = nsite;
}
