#define NTHREAD 64
#define DMAMAX 512
#define DMAMAX2 64
#define DMAMAX3 32
#define USE_SIMD
//#define REG_COM

struct _sw_site_propensity_arg
{
	int *c_vacant;
	int vac_size;

	double c_barrier[6];
	double vij1[6][6];
	double vij2[6][6];

	double **xyz;
	double unitdis;
	int *lattice;
	double *E_V;

	int boxxlo;
	int boxxhi;
	int subxlo;
	int subxhi;
	int boxylo;
	int boxyhi;
	int subylo;
	int subyhi;
	int boxzlo;
	int boxzhi;
	int subzlo;
	int subzhi;
	int xghost, yghost, zghost;
	int nlocal;

	double t_inverse;

	double *propensity;
	double **propensity_dir;
	int **jumpsite_dir;
	int *num_dir;
};

struct _swarg
{
	int *c_vacant;
	double t_inverse;
	int vac_size;

	int site1[3];
	int site2[3];

	double c_barrier[6];
	double vij1[6][6];
	double vij2[6][6];

	int *nsite;
	int *index;

	int boxx;
	int boxy;
	int boxz;
	int boxxlo;
	int boxxhi;
	int subxlo;
	int subxhi;
	int boxylo;
	int boxyhi;
	int subylo;
	int subyhi;
	int boxzlo;
	int boxzhi;
	int subzlo;
	int subzhi;
	double unitdis;
	double **xyz;
	int *lattice;
	double *E_V;
	double *propensity;
	double **propensity_dir;
	int **jumpsite_dir;
	int *num_dir;
	int xghost, yghost, zghost;
	int nlocal;
};