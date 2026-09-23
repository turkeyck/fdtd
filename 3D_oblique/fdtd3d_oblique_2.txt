/*
 * fdtd3d_oblique.c -- 3D Yee FDTD, oblique plane wave, PBC in x/z, CPML in y, TF/SF at y = y0.
 *
 * Normalized units: c = eps0 = mu0 = 1, lambda0 = 1 (omega0 = 2*pi).
 * Equations, index conventions and sign derivations: derivation.md (section numbers cited below).
 *
 * Build: make            (gcc, C99, optional OpenMP)
 * Usage: ./fdtd3d_oblique key=value ...   (see parse_args for the list and defaults)
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#define MKDIR(p) mkdir(p, 0755)
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.14159265358979323846
#define MAXPLANES 64

/* ------------------------------------------------------------------ parameters */
typedef struct {
    int nl;                 /* cells per lambda0 */
    double S, Lx, Lz;
    int m, n;
    char pol;               /* 's' or 'p' */
    int npml, sf, tf;       /* near PML, SF length, TF length (cells); Ny = 2*npml + sf + tf */
    double pml_m, sig_fac, kappa_max, alpha_max;
    double eps2;            /* relative permittivity of the half space y >= y1 */
    int y1;                 /* interface offset from j0 in cells (<0: no medium) */
    char ifmode;            /* 'a': arithmetic mean for tangential E on the interface plane */
    char inc;               /* 'a': aux modal line, 'n': analytic*g(t), '0': none */
    int ky_cont;            /* analytic mode only: 1 -> continuous ky (control) */
    char ramp;              /* 'r': raised cosine, 'e': erf */
    double ramp_T;          /* raised-cosine duration (periods) */
    double erf_t0, erf_tau; /* erf ramp centre and width (periods) */
    double off_t;           /* start of the turn-off ramp (periods); <0: never */
    int na;                 /* aux injection point j_a = j0 - na */
    long nsteps, dft0, dft1;
    int nyplanes, yplanes[MAXPLANES];
    int zk, xi;             /* x-y slice at k = zk, y-z slice at i = xi (<0: none) */
    char init;              /* '0': zero, 'a': analytic fill + Dirichlet y (debug), 'b': gaussian blob */
    int energy_every, snap_every;
    char snapcomp[4];
    int dump;               /* dump full final fields */
    int auxspan;            /* aux-line DFT recorded up to j0 + auxspan (0: up to Ny) */
    int div_every;          /* log max |div E|, |div H| over the TF interior every N steps (0: off) */
    char out[512];
} Params;

/* ------------------------------------------------------------------ globals */
static Params P;
static int Nx, Ny, Nz, SY, SZ, j0, j1;
static double D, dt, W0, kx, ky, kz, kyc;
static double Kt[3], wt, E0v[3], H0v[3], Kt_src[3];
static size_t NTOT;
static double *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
static double *ceEx, *ceEy, *ceEz, ch;
static int NP, *pmlE, *pmlH;
static double *bE, *cE, *ikE, *bH, *cH, *ikH, *sigE, *sigH, *kapE, *kapH, *alpE, *alpH;
static double *psiExy, *psiEzy, *psiHxy, *psiHzy;
static double *cosP1, *sinP1, *cosP2, *sinP2;  /* P1: x=(i+1/2)D, z=kD (Ex,Hz); P2: x=iD, z=(k+1/2)D (Ez,Hx) */
static double *ExI, *EzI, *HxI, *HzI;           /* incident samples on the TF/SF plane */

#define ID(i, j, k) ((((size_t)(i) + 1) * (size_t)SY + (size_t)(j)) * (size_t)SZ + (size_t)(k) + 1)

/* Yee offsets (derivation.md §1): x, y, z in cells; time in steps. Order Ex,Ey,Ez,Hx,Hy,Hz. */
static const double OFF[6][4] = {
    {0.5, 0.0, 0.0, 0.0}, {0.0, 0.5, 0.0, 0.0}, {0.0, 0.0, 0.5, 0.0},
    {0.0, 0.5, 0.5, 0.5}, {0.5, 0.0, 0.5, 0.5}, {0.5, 0.5, 0.0, 0.5}};
static const char *CNAME[6] = {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz"};

static void die(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(2);
}

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (!p) die("out of memory");
    return p;
}

/* ------------------------------------------------------------------ args */
static void parse_list(const char *s, int *arr, int *cnt) {
    *cnt = 0;
    char buf[2048];
    strncpy(buf, s, sizeof buf - 1);
    buf[sizeof buf - 1] = 0;
    for (char *t = strtok(buf, ","); t && *cnt < MAXPLANES; t = strtok(NULL, ",")) arr[(*cnt)++] = atoi(t);
}

static void parse_args(int argc, char **argv) {
    Params p = {0};
    p.nl = 20; p.S = 0.5; p.Lx = 2.0; p.Lz = 3.0; p.m = 1; p.n = 1; p.pol = 's';
    p.npml = 20; p.sf = 20; p.tf = 200;
    p.pml_m = 3.0; p.sig_fac = 0.8; p.kappa_max = 5.0; p.alpha_max = 0.05 * 2.0 * PI;
    p.eps2 = 1.0; p.y1 = -1; p.ifmode = 'a';
    p.inc = 'a'; p.ky_cont = 0; p.ramp = 'e'; p.ramp_T = 10.0; p.erf_t0 = 20.0; p.erf_tau = 4.0; p.off_t = -1.0;
    p.na = 40; p.nsteps = 1000; p.dft0 = -1; p.dft1 = -1; p.nyplanes = 0; p.zk = -1; p.xi = -1;
    p.init = '0'; p.energy_every = 10; p.snap_every = 0; strcpy(p.snapcomp, "Ez"); p.dump = 0;
    strcpy(p.out, "out");
    for (int a = 1; a < argc; a++) {
        char *eq = strchr(argv[a], '=');
        if (!eq) { fprintf(stderr, "bad arg %s\n", argv[a]); exit(2); }
        *eq = 0;
        const char *k = argv[a], *v = eq + 1;
        if (!strcmp(k, "nl")) p.nl = atoi(v);
        else if (!strcmp(k, "S")) p.S = atof(v);
        else if (!strcmp(k, "Lx")) p.Lx = atof(v);
        else if (!strcmp(k, "Lz")) p.Lz = atof(v);
        else if (!strcmp(k, "m")) p.m = atoi(v);
        else if (!strcmp(k, "n")) p.n = atoi(v);
        else if (!strcmp(k, "pol")) p.pol = v[0];
        else if (!strcmp(k, "npml")) p.npml = atoi(v);
        else if (!strcmp(k, "sf")) p.sf = atoi(v);
        else if (!strcmp(k, "tf")) p.tf = atoi(v);
        else if (!strcmp(k, "pml_m")) p.pml_m = atof(v);
        else if (!strcmp(k, "sig_fac")) p.sig_fac = atof(v);
        else if (!strcmp(k, "kappa_max")) p.kappa_max = atof(v);
        else if (!strcmp(k, "alpha_max")) p.alpha_max = atof(v);
        else if (!strcmp(k, "eps2")) p.eps2 = atof(v);
        else if (!strcmp(k, "y1")) p.y1 = atoi(v);
        else if (!strcmp(k, "ifmode")) p.ifmode = v[0];
        else if (!strcmp(k, "inc")) p.inc = v[0];
        else if (!strcmp(k, "ky")) p.ky_cont = (v[0] == 'c');
        else if (!strcmp(k, "ramp")) p.ramp = v[0];
        else if (!strcmp(k, "ramp_T")) p.ramp_T = atof(v);
        else if (!strcmp(k, "erf_t0")) p.erf_t0 = atof(v);
        else if (!strcmp(k, "erf_tau")) p.erf_tau = atof(v);
        else if (!strcmp(k, "off_t")) p.off_t = atof(v);
        else if (!strcmp(k, "na")) p.na = atoi(v);
        else if (!strcmp(k, "nsteps")) p.nsteps = atol(v);
        else if (!strcmp(k, "dft0")) p.dft0 = atol(v);
        else if (!strcmp(k, "dft1")) p.dft1 = atol(v);
        else if (!strcmp(k, "yplanes")) parse_list(v, p.yplanes, &p.nyplanes);
        else if (!strcmp(k, "zk")) p.zk = atoi(v);
        else if (!strcmp(k, "xi")) p.xi = atoi(v);
        else if (!strcmp(k, "init")) p.init = v[0];
        else if (!strcmp(k, "energy_every")) p.energy_every = atoi(v);
        else if (!strcmp(k, "snap")) p.snap_every = atoi(v);
        else if (!strcmp(k, "snapcomp")) { strncpy(p.snapcomp, v, 3); p.snapcomp[3] = 0; }
        else if (!strcmp(k, "dump")) p.dump = atoi(v);
        else if (!strcmp(k, "auxspan")) p.auxspan = atoi(v);
        else if (!strcmp(k, "div_every")) p.div_every = atoi(v);
        else if (!strcmp(k, "out")) { strncpy(p.out, v, sizeof p.out - 1); }
        else { fprintf(stderr, "unknown key %s\n", k); exit(2); }
    }
    P = p;
}

/* ------------------------------------------------------------------ theory (derivation.md §2) */
static double ky_disc(double kxx, double kzz, double nidx) {
    double q = pow(nidx / P.S * sin(W0 * dt / 2.0), 2) - pow(sin(kxx * D / 2.0), 2) - pow(sin(kzz * D / 2.0), 2);
    if (!(q > 0.0 && q <= 1.0)) {
        fprintf(stderr, "sin^2(ky D/2) = %.6g outside (0,1]\n", q);
        die("no propagating discrete ky for this (m, n, Delta, dt)");
    }
    return 2.0 / D * asin(sqrt(q));
}

static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(const double *a) { return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

/* E0, H0 for the wave vector whose y component is kyy (s_hat ∝ K~ x y_hat, p_hat ∝ s_hat x K~). */
static void amplitudes(double kyy, double *K, double *E, double *H) {
    double yh[3] = {0, 1, 0}, s[3], p[3];
    K[0] = 2.0 / D * sin(kx * D / 2.0);
    K[1] = 2.0 / D * sin(kyy * D / 2.0);
    K[2] = 2.0 / D * sin(kz * D / 2.0);
    cross(K, yh, s);
    double ns = norm3(s);
    if (ns == 0.0) die("normal incidence (m = n = 0) is not supported: s/p basis degenerate");
    for (int c = 0; c < 3; c++) s[c] /= ns;
    cross(s, K, p);
    double np = norm3(p);
    for (int c = 0; c < 3; c++) p[c] /= np;
    for (int c = 0; c < 3; c++) E[c] = (P.pol == 's') ? s[c] : p[c];
    cross(K, E, H);
    for (int c = 0; c < 3; c++) H[c] /= wt; /* mu0 = 1 */
}

/* ------------------------------------------------------------------ time envelope */
static double ramp_on(double t) {
    double T0 = 1.0; /* lambda0 / c */
    if (P.ramp == 'e') return 0.5 * erfc(-(t - P.erf_t0 * T0) / (P.erf_tau * T0));
    double Tr = P.ramp_T * T0;
    if (t <= 0.0) return 0.0;
    if (t >= Tr) return 1.0;
    return 0.5 * (1.0 - cos(PI * t / Tr));
}

static double envelope(double t) {
    double g = ramp_on(t);
    if (P.off_t >= 0.0) g *= 1.0 - ramp_on(t - P.off_t);
    return g;
}

/* complex amplitude of component c at global y-index j (its own y offset) and time step tn (own offset) */
static double complex inc_amp(int c, long j, double tn, const double *E, const double *H, double kyy) {
    double y = ((double)j + OFF[c][1]) * D;
    double t = (tn + OFF[c][3]) * dt;
    double a = (c < 3) ? E[c] : H[c - 3];
    return envelope(t) * a * cexp(I * (kyy * y - W0 * t));
}

/* ------------------------------------------------------------------ aux modal line (derivation.md §6) */
static long AL, Ajlo, Aua;
static double complex *aEx, *aEy, *aEz, *aHx, *aHy, *aHz;
static double complex ikxD, ikzD;

static void aux_init(void) {
    long ja = j0 - P.na;
    long half = P.nsteps / 2 + 20;
    Ajlo = ja - half;
    long jhi = j0 + half + (P.auxspan > 0 ? P.auxspan : Ny - j0); /* no end reflection reaches a recorded node */
    AL = jhi - Ajlo;
    Aua = ja - Ajlo;
    aEx = xcalloc(AL + 1, sizeof *aEx); aEz = xcalloc(AL + 1, sizeof *aEz); aHy = xcalloc(AL + 1, sizeof *aHy);
    aEy = xcalloc(AL, sizeof *aEy); aHx = xcalloc(AL, sizeof *aHx); aHz = xcalloc(AL, sizeof *aHz);
    ikxD = I * Kt[0] * D;
    ikzD = I * Kt[2] * D;
}

static void aux_update_H(long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 0) lo = 0;
    if (hi > AL - 1) hi = AL - 1;
    double c = dt / D; /* mu0 = 1 */
    for (long u = lo; u <= hi; u++) {
        aHx[u] -= c * ((aEz[u + 1] - aEz[u]) - ikzD * aEy[u]);
        aHz[u] -= c * (ikxD * aEy[u] - (aEx[u + 1] - aEx[u]));
    }
    for (long u = lo; u <= hi + 1 && u <= AL; u++) aHy[u] -= c * (ikzD * aEx[u] - ikxD * aEz[u]);
    /* 1D TF/SF at j_a: same signs as derivation.md §3 (a), (b) */
    long ja = Ajlo + Aua;
    aHx[Aua - 1] += c * inc_amp(2, ja, (double)n, E0v, H0v, ky);
    aHz[Aua - 1] -= c * inc_amp(0, ja, (double)n, E0v, H0v, ky);
}

static void aux_update_E(long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 1) lo = 1;
    if (hi > AL - 1) hi = AL - 1;
    double c = dt / D; /* eps0 = 1; the aux line is vacuum */
    for (long u = lo; u <= hi; u++) {
        aEx[u] += c * ((aHz[u] - aHz[u - 1]) - ikzD * aHy[u]);
        aEz[u] += c * (ikxD * aHy[u] - (aHx[u] - aHx[u - 1]));
    }
    for (long u = lo - 1; u <= hi && u < AL; u++) aEy[u] += c * (ikzD * aHx[u] - ikxD * aHz[u]);
    long ja = Ajlo + Aua;
    aEx[Aua] -= c * inc_amp(5, ja - 1, (double)n, E0v, H0v, ky);
    aEz[Aua] += c * inc_amp(3, ja - 1, (double)n, E0v, H0v, ky);
}

/* ------------------------------------------------------------------ setup */
static void setup(void) {
    D = 1.0 / P.nl;
    dt = P.S * D;
    W0 = 2.0 * PI;
    if (P.S >= 1.0 / sqrt(3.0)) die("Courant number S must be < 1/sqrt(3)");
    Nx = (int)lround(P.Lx * P.nl);
    Nz = (int)lround(P.Lz * P.nl);
    if (fabs(Nx * D - P.Lx) > 1e-12 * P.Lx || fabs(Nz * D - P.Lz) > 1e-12 * P.Lz)
        die("Lx and Lz must be integer multiples of Delta");
    kx = 2.0 * PI * P.m / P.Lx;
    kz = 2.0 * PI * P.n / P.Lz;
    double kt = hypot(kx, kz);
    if (kt >= 0.95 * W0) die("|k_t| >= 0.95 k0: insufficient propagation margin");
    wt = 2.0 / dt * sin(W0 * dt / 2.0);
    ky = ky_disc(kx, kz, 1.0);
    kyc = sqrt(W0 * W0 - kt * kt);
    amplitudes(ky, Kt, E0v, H0v);
    for (int c = 0; c < 3; c++) Kt_src[c] = Kt[c];
    if (P.inc == 'n' && P.ky_cont) { /* control C1: continuous ky everywhere in the source */
        amplitudes(kyc, Kt_src, E0v, H0v);
        ky = kyc;
    }
    if (P.pol != 's' && P.pol != 'p') die("pol must be s or p");
    if (P.tf < 10 * P.nl) fprintf(stderr, "WARNING: TF length %d cells < 10 lambda0\n", P.tf);

    Ny = 2 * P.npml + P.sf + P.tf;
    j0 = P.npml + P.sf;
    j1 = (P.y1 >= 0) ? j0 + P.y1 : -1;
    if (P.inc != '0' && P.init == 'a') die("init=analytic is a debug mode without a source");
    if (P.inc != '0' && (P.sf < 2 || P.tf < 2)) die("sf and tf must be >= 2 cells");
    if (j1 >= Ny - P.npml) die("interface must lie in the physical TF region");
    SY = Ny + 1;
    SZ = Nz + 2;
    NTOT = (size_t)(Nx + 2) * SY * SZ;
    Ex = xcalloc(NTOT, sizeof(double)); Ey = xcalloc(NTOT, sizeof(double)); Ez = xcalloc(NTOT, sizeof(double));
    Hx = xcalloc(NTOT, sizeof(double)); Hy = xcalloc(NTOT, sizeof(double)); Hz = xcalloc(NTOT, sizeof(double));

    /* material (derivation.md §8): Ex,Ez at integer j; Ey at j+1/2 never lies on an integer-y interface */
    ceEx = xcalloc(SY, sizeof(double)); ceEy = xcalloc(SY, sizeof(double)); ceEz = xcalloc(SY, sizeof(double));
    for (int j = 0; j <= Ny; j++) {
        double et = 1.0, en = 1.0;
        if (j1 >= 0) {
            if (j > j1) et = P.eps2;
            else if (j == j1) et = (P.ifmode == 'a') ? 0.5 * (1.0 + P.eps2) : P.eps2;
            if (j >= j1) en = P.eps2; /* Ey at (j+1/2) >= j1 */
        }
        ceEx[j] = ceEz[j] = dt / (et * D);
        ceEy[j] = dt / (en * D);
    }
    ch = dt / D;

    /* CPML (derivation.md §7) */
    pmlE = xcalloc(SY, sizeof(int)); pmlH = xcalloc(SY, sizeof(int));
    bE = xcalloc(SY, sizeof(double)); cE = xcalloc(SY, sizeof(double)); ikE = xcalloc(SY, sizeof(double));
    bH = xcalloc(SY, sizeof(double)); cH = xcalloc(SY, sizeof(double)); ikH = xcalloc(SY, sizeof(double));
    sigE = xcalloc(SY, sizeof(double)); sigH = xcalloc(SY, sizeof(double));
    kapE = xcalloc(SY, sizeof(double)); kapH = xcalloc(SY, sizeof(double));
    alpE = xcalloc(SY, sizeof(double)); alpH = xcalloc(SY, sizeof(double));
    NP = 2 * (P.npml + 1);
    double d = P.npml * D;
    double n_far = (j1 >= 0) ? sqrt(P.eps2) : 1.0;
    for (int j = 0; j <= Ny; j++) {
        for (int h = 0; h < 2; h++) {
            double y = (j + 0.5 * h) * D;
            double rho = 0.0, nloc = 1.0;
            int slot = -1;
            if (P.npml > 0 && y < P.npml * D) { rho = P.npml * D - y; slot = j; }
            else if (P.npml > 0 && y > (Ny - P.npml) * D) { rho = y - (Ny - P.npml) * D; slot = j - (Ny - P.npml) + P.npml + 1; nloc = n_far; }
            if (h == 1 && j == Ny) slot = -1;
            double r = (slot >= 0) ? rho / d : 0.0;
            double smax = P.sig_fac * (P.pml_m + 1.0) / (D * nloc); /* eta0 = 1, eta = 1/n */
            double sg = smax * pow(r, P.pml_m);
            double kp = 1.0 + (P.kappa_max - 1.0) * pow(r, P.pml_m);
            double al = P.alpha_max * (1.0 - r);
            double b = exp(-(sg / kp + al) * dt);
            double cc = (sg > 0.0) ? sg * (b - 1.0) / (kp * (sg + kp * al)) : 0.0;
            if (h == 0) { pmlE[j] = slot; bE[j] = b; cE[j] = cc; ikE[j] = 1.0 / kp; sigE[j] = sg; kapE[j] = kp; alpE[j] = al; }
            else { pmlH[j] = slot; bH[j] = b; cH[j] = cc; ikH[j] = 1.0 / kp; sigH[j] = sg; kapH[j] = kp; alpH[j] = al; }
        }
    }
    size_t np = (size_t)Nx * NP * Nz;
    psiExy = xcalloc(np, sizeof(double)); psiEzy = xcalloc(np, sizeof(double));
    psiHxy = xcalloc(np, sizeof(double)); psiHzy = xcalloc(np, sizeof(double));

    /* incident-plane phase tables */
    cosP1 = xcalloc((size_t)Nx * Nz, sizeof(double)); sinP1 = xcalloc((size_t)Nx * Nz, sizeof(double));
    cosP2 = xcalloc((size_t)Nx * Nz, sizeof(double)); sinP2 = xcalloc((size_t)Nx * Nz, sizeof(double));
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            double p1 = kx * (i + 0.5) * D + kz * k * D, p2 = kx * i * D + kz * (k + 0.5) * D;
            cosP1[i * Nz + k] = cos(p1); sinP1[i * Nz + k] = sin(p1);
            cosP2[i * Nz + k] = cos(p2); sinP2[i * Nz + k] = sin(p2);
        }
    ExI = xcalloc((size_t)Nx * Nz, sizeof(double)); EzI = xcalloc((size_t)Nx * Nz, sizeof(double));
    HxI = xcalloc((size_t)Nx * Nz, sizeof(double)); HzI = xcalloc((size_t)Nx * Nz, sizeof(double));
    if (P.inc == 'a') aux_init();
}

/* ------------------------------------------------------------------ analytic fill (debug, Stage 1) */
static double analytic_val(int c, int i, int j, int k, double tn) {
    double x = (i + OFF[c][0]) * D, y = (j + OFF[c][1]) * D, z = (k + OFF[c][2]) * D, t = (tn + OFF[c][3]) * dt;
    double a = (c < 3) ? E0v[c] : H0v[c - 3];
    return a * cos(kx * x + ky * y + kz * z - W0 * t);
}

static double *comp_ptr(int c) {
    double *arr[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    return arr[c];
}

static void fill_analytic(double tnE, double tnH) {
    for (int c = 0; c < 6; c++) {
        double *F = comp_ptr(c);
        int jmax = (OFF[c][1] == 0.0) ? Ny : Ny - 1;
        double tn = (c < 3) ? tnE : tnH;
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j <= jmax; j++)
                for (int k = 0; k < Nz; k++) F[ID(i, j, k)] = analytic_val(c, i, j, k, tn);
    }
}

static void dirichlet_y(double tn) {
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            Ex[ID(i, 0, k)] = analytic_val(0, i, 0, k, tn); Ex[ID(i, Ny, k)] = analytic_val(0, i, Ny, k, tn);
            Ez[ID(i, 0, k)] = analytic_val(2, i, 0, k, tn); Ez[ID(i, Ny, k)] = analytic_val(2, i, Ny, k, tn);
        }
}

/* ------------------------------------------------------------------ PBC ghost layers (derivation.md §4) */
static void ghost_fwd(double *F) { /* F[Nx] <- F[0], F[.][.][Nz] <- F[.][.][0] */
    for (int j = 0; j <= Ny; j++)
        for (int k = 0; k < Nz; k++) F[ID(Nx, j, k)] = F[ID(0, j, k)];
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++) F[ID(i, j, Nz)] = F[ID(i, j, 0)];
}

static void ghost_bwd(double *F) { /* F[-1] <- F[Nx-1], F[.][.][-1] <- F[.][.][Nz-1] */
    for (int j = 0; j <= Ny; j++)
        for (int k = 0; k < Nz; k++) F[ID(-1, j, k)] = F[ID(Nx - 1, j, k)];
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++) F[ID(i, j, -1)] = F[ID(i, j, Nz - 1)];
}

/* ------------------------------------------------------------------ updates */
static void update_H(void) {
    const size_t sI = (size_t)SY * SZ, sJ = SZ;
    ghost_fwd(Ex); ghost_fwd(Ey); ghost_fwd(Ez);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            size_t b = ID(i, j, 0);
            double *hy = Hy + b;
            const double *ex = Ex + b, *ez = Ez + b, *ey = Ey + b;
            for (int k = 0; k < Nz; k++) hy[k] -= ch * ((ex[k + 1] - ex[k]) - (ez[k + sI] - ez[k]));
            if (j == Ny) continue;
            double *hx = Hx + b, *hz = Hz + b;
            int s = pmlH[j];
            if (s < 0) {
                for (int k = 0; k < Nz; k++) {
                    hx[k] -= ch * ((ez[k + sJ] - ez[k]) - (ey[k + 1] - ey[k]));
                    hz[k] -= ch * ((ey[k + sI] - ey[k]) - (ex[k + sJ] - ex[k]));
                }
            } else {
                double bb = bH[j], cc = cH[j], ik = ikH[j];
                double *px = psiHxy + ((size_t)i * NP + s) * Nz, *pz = psiHzy + ((size_t)i * NP + s) * Nz;
                for (int k = 0; k < Nz; k++) {
                    double dEz = ez[k + sJ] - ez[k], dEx = ex[k + sJ] - ex[k];
                    px[k] = bb * px[k] + cc * dEz;
                    pz[k] = bb * pz[k] + cc * dEx;
                    hx[k] -= ch * ((ik * dEz + px[k]) - (ey[k + 1] - ey[k]));
                    hz[k] -= ch * ((ey[k + sI] - ey[k]) - (ik * dEx + pz[k]));
                }
            }
        }
    }
}

/* E update. If doW, also accumulates the Yee-conserved energy
   W^{n+1/2} = 1/2 sum( eps E^n . E^{n+1} + mu H^{n+1/2} . H^{n+1/2} ) Delta^3
   (exactly constant in a lossless PEC/PBC cavity), per i into partT/partP (total / outside PML). */
static double *partT, *partP;
static void update_E(int doW) {
    const size_t sI = (size_t)SY * SZ, sJ = SZ;
    ghost_bwd(Hx); ghost_bwd(Hy); ghost_bwd(Hz);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < Nx; i++) {
        double wT = 0.0, wP = 0.0;
        for (int j = 0; j <= Ny; j++) {
            size_t b = ID(i, j, 0);
            const double *hx = Hx + b, *hy = Hy + b, *hz = Hz + b;
            double wI = 0.0, wH = 0.0; /* integer-y nodes, half-y nodes */
            if (j < Ny) {
                double *ey = Ey + b, c = ceEy[j], epsy = dt / (c * D);
                for (int k = 0; k < Nz; k++) {
                    double old = ey[k];
                    ey[k] += c * ((hx[k] - hx[k - 1]) - (hz[k] - hz[k - sI]));
                    if (doW) wH += epsy * old * ey[k] + hx[k] * hx[k] + hz[k] * hz[k];
                }
            }
            if (doW)
                for (int k = 0; k < Nz; k++) wI += hy[k] * hy[k];
            if (j > 0 && j < Ny) { /* j = 0, Ny: PEC, tangential E fixed (Dirichlet in debug mode) */
                double *ex = Ex + b, *ez = Ez + b, cx = ceEx[j], cz = ceEz[j];
                double epsx = dt / (cx * D), epsz = dt / (cz * D);
                int s = pmlE[j];
                if (s < 0) {
                    for (int k = 0; k < Nz; k++) {
                        double ox = ex[k], oz = ez[k];
                        ex[k] += cx * ((hz[k] - hz[k - sJ]) - (hy[k] - hy[k - 1]));
                        ez[k] += cz * ((hy[k] - hy[k - sI]) - (hx[k] - hx[k - sJ]));
                        if (doW) wI += epsx * ox * ex[k] + epsz * oz * ez[k];
                    }
                } else {
                    double bb = bE[j], cc = cE[j], ik = ikE[j];
                    double *px = psiExy + ((size_t)i * NP + s) * Nz, *pz = psiEzy + ((size_t)i * NP + s) * Nz;
                    for (int k = 0; k < Nz; k++) {
                        double ox = ex[k], oz = ez[k];
                        double dHz = hz[k] - hz[k - sJ], dHx = hx[k] - hx[k - sJ];
                        px[k] = bb * px[k] + cc * dHz;
                        pz[k] = bb * pz[k] + cc * dHx;
                        ex[k] += cx * ((ik * dHz + px[k]) - (hy[k] - hy[k - 1]));
                        ez[k] += cz * ((hy[k] - hy[k - sI]) - (ik * dHx + pz[k]));
                        if (doW) wI += epsx * ox * ex[k] + epsz * oz * ez[k];
                    }
                }
            } else if (doW && P.init == 'a') {
                double *ex = Ex + b, *ez = Ez + b;
                for (int k = 0; k < Nz; k++) wI += ex[k] * ex[k] + ez[k] * ez[k];
            }
            wT += wI + wH;
            if (j >= P.npml && j <= Ny - P.npml) wP += wI;
            if (j >= P.npml && j <= Ny - P.npml - 1) wP += wH;
        }
        if (doW) { partT[i] = wT; partP[i] = wP; }
    }
}

static void energy_sum(double *Wt, double *Wp) {
    double a = 0.0, b = 0.0;
    for (int i = 0; i < Nx; i++) { a += partT[i]; b += partP[i]; }
    *Wt = 0.5 * a * D * D * D;
    *Wp = 0.5 * b * D * D * D;
}

/* incident field on the TF/SF plane from complex amplitudes (derivation.md §6.2) */
static void incident_E(long n) { /* Ex, Ez at j0, time n */
    double complex ax, az;
    if (P.inc == 'a') { ax = aEx[j0 - Ajlo]; az = aEz[j0 - Ajlo]; }
    else { ax = inc_amp(0, j0, (double)n, E0v, H0v, ky); az = inc_amp(2, j0, (double)n, E0v, H0v, ky); }
    for (int q = 0; q < Nx * Nz; q++) {
        ExI[q] = creal(ax) * cosP1[q] - cimag(ax) * sinP1[q];
        EzI[q] = creal(az) * cosP2[q] - cimag(az) * sinP2[q];
    }
}

static void incident_H(long n) { /* Hx, Hz at j0 - 1/2, time n + 1/2 */
    double complex ax, az;
    if (P.inc == 'a') { ax = aHx[j0 - 1 - Ajlo]; az = aHz[j0 - 1 - Ajlo]; }
    else { ax = inc_amp(3, j0 - 1, (double)n, E0v, H0v, ky); az = inc_amp(5, j0 - 1, (double)n, E0v, H0v, ky); }
    for (int q = 0; q < Nx * Nz; q++) {
        HxI[q] = creal(ax) * cosP2[q] - cimag(ax) * sinP2[q];
        HzI[q] = creal(az) * cosP1[q] - cimag(az) * sinP1[q];
    }
}

/* TF/SF corrections, derivation.md §3 (a)-(d) */
static void tfsf_H(void) {
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            Hx[ID(i, j0 - 1, k)] += ch * EzI[i * Nz + k];
            Hz[ID(i, j0 - 1, k)] -= ch * ExI[i * Nz + k];
        }
}

static void tfsf_E(void) {
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            Ex[ID(i, j0, k)] -= ceEx[j0] * HzI[i * Nz + k];
            Ez[ID(i, j0, k)] += ceEz[j0] * HxI[i * Nz + k];
        }
}

/* ------------------------------------------------------------------ DFT */
typedef struct {
    char name[16];
    int type;         /* 0: y-plane (Nx x Nz), 1: x-y slice at k (Nx x SY), 2: y-z slice at i (SY x Nz) */
    int idx;
    size_t cnt;
    double *re[6], *im[6];
} Slice;
static Slice *SL;
static int NSL;

static void add_slice(const char *name, int type, int idx) {
    Slice *s = &SL[NSL++];
    snprintf(s->name, sizeof s->name, "%s", name);
    s->type = type;
    s->idx = idx;
    s->cnt = (type == 0) ? (size_t)Nx * Nz : (type == 1) ? (size_t)Nx * SY : (size_t)SY * Nz;
    for (int c = 0; c < 6; c++) { s->re[c] = xcalloc(s->cnt, sizeof(double)); s->im[c] = xcalloc(s->cnt, sizeof(double)); }
}

static double complex *auxdft[6];
static long aux_dft_lo, aux_dft_hi;

static void dft_setup(void) {
    SL = xcalloc(MAXPLANES + 2, sizeof(Slice));
    NSL = 0;
    char nm[16];
    for (int p = 0; p < P.nyplanes; p++) {
        if (P.yplanes[p] < 0 || P.yplanes[p] >= Ny) die("y-plane index out of range");
        snprintf(nm, sizeof nm, "y%d", P.yplanes[p]);
        add_slice(nm, 0, P.yplanes[p]);
    }
    if (P.zk >= 0) add_slice("xy", 1, P.zk);
    if (P.xi >= 0) add_slice("yz", 2, P.xi);
    if (P.inc == 'a') {
        aux_dft_lo = 0 - Ajlo;
        aux_dft_hi = (P.auxspan > 0 ? j0 + P.auxspan : Ny) - Ajlo;
        if (aux_dft_lo < 0) aux_dft_lo = 0;
        if (aux_dft_hi > AL - 1) aux_dft_hi = AL - 1;
        for (int c = 0; c < 6; c++) auxdft[c] = xcalloc(aux_dft_hi - aux_dft_lo + 1, sizeof(double complex));
    }
}

static void dft_accumulate(int first, int last, double tsec) {
    double cr = cos(W0 * tsec), ci = sin(W0 * tsec);
    for (int s = 0; s < NSL; s++) {
        Slice *sl = &SL[s];
        for (int c = first; c <= last; c++) {
            double *F = comp_ptr(c);
            size_t q = 0;
            if (sl->type == 0) {
                for (int i = 0; i < Nx; i++)
                    for (int k = 0; k < Nz; k++, q++) { double f = F[ID(i, sl->idx, k)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci; }
            } else if (sl->type == 1) {
                for (int i = 0; i < Nx; i++)
                    for (int j = 0; j <= Ny; j++, q++) { double f = F[ID(i, j, sl->idx)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci; }
            } else {
                for (int j = 0; j <= Ny; j++)
                    for (int k = 0; k < Nz; k++, q++) { double f = F[ID(sl->idx, j, k)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci; }
            }
        }
    }
    if (P.inc == 'a') {
        double complex e = cexp(I * W0 * tsec);
        double complex *arr[6] = {aEx, aEy, aEz, aHx, aHy, aHz};
        for (int c = first; c <= last; c++)
            for (long u = aux_dft_lo; u <= aux_dft_hi; u++) auxdft[c][u - aux_dft_lo] += arr[c][u] * e;
    }
}

/* ------------------------------------------------------------------ diagnostics */
static void sf_max(double *mE, double *mH) {
    double me = 0.0, mh = 0.0;
    int jlo = P.npml + 1;           /* integer-y nodes strictly outside the near PML */
    for (int i = 0; i < Nx; i++)
        for (int j = P.npml; j < j0; j++)
            for (int k = 0; k < Nz; k++) {
                size_t q = ID(i, j, k);
                double a;
                a = fabs(Ey[q]); if (a > me) me = a;           /* half-integer nodes y = j+1/2 in (npml, j0) */
                a = fabs(Hx[q]); if (a > mh) mh = a;
                a = fabs(Hz[q]); if (a > mh) mh = a;
                if (j >= jlo) {
                    a = fabs(Ex[q]); if (a > me) me = a;
                    a = fabs(Ez[q]); if (a > me) me = a;
                    a = fabs(Hy[q]); if (a > mh) mh = a;
                }
            }
    *mE = me;
    *mH = mh;
}

/* max |div E| at nodes (i,j,k) and max |div H| at cell centres over the TF interior
   j in [j0+2, Ny-npml-2] (excludes the TF/SF plane +-1 cell and the PML); explicit PBC wrap. */
static void div_max(double *dE, double *dH) {
    double me = 0.0, mh = 0.0;
    int jlo = j0 + 2, jhi = Ny - P.npml - 2;
    if (j1 >= 0 && j1 - 1 < jhi) jhi = j1 - 1; /* stay on the vacuum side of an interface */
#pragma omp parallel for schedule(static) reduction(max : me, mh)
    for (int i = 0; i < Nx; i++) {
        int im = (i + Nx - 1) % Nx, ip = (i + 1) % Nx;
        for (int j = jlo; j <= jhi; j++)
            for (int k = 0; k < Nz; k++) {
                int km = (k + Nz - 1) % Nz, kp = (k + 1) % Nz;
                double de = (Ex[ID(i, j, k)] - Ex[ID(im, j, k)]) + (Ey[ID(i, j, k)] - Ey[ID(i, j - 1, k)]) +
                            (Ez[ID(i, j, k)] - Ez[ID(i, j, km)]);
                double dh = (Hx[ID(ip, j, k)] - Hx[ID(i, j, k)]) + (Hy[ID(i, j + 1, k)] - Hy[ID(i, j, k)]) +
                            (Hz[ID(i, j, kp)] - Hz[ID(i, j, k)]);
                if (fabs(de) > me) me = fabs(de);
                if (fabs(dh) > mh) mh = fabs(dh);
            }
    }
    *dE = me / D;
    *dH = mh / D;
}

/* ------------------------------------------------------------------ output */
static FILE *xfopen(const char *name, const char *mode) {
    char path[1024];
    snprintf(path, sizeof path, "%s/%s", P.out, name);
    FILE *f = fopen(path, mode);
    if (!f) die("cannot open output file");
    return f;
}

static void write_dft(void) {
    long N = P.dft1 - P.dft0;
    if (N <= 0) return;
    double sc = 2.0 / (double)N; /* phasor A with f = Re[A exp(-i w t)] */
    for (int s = 0; s < NSL; s++) {
        char fn[64];
        snprintf(fn, sizeof fn, "dft_%s.bin", SL[s].name);
        FILE *f = xfopen(fn, "wb");
        for (int c = 0; c < 6; c++)
            for (size_t q = 0; q < SL[s].cnt; q++) {
                double v[2] = {SL[s].re[c][q] * sc, SL[s].im[c][q] * sc};
                fwrite(v, sizeof(double), 2, f);
            }
        fclose(f);
    }
    if (P.inc == 'a') {
        FILE *f = xfopen("dft_aux.bin", "wb");
        for (int c = 0; c < 6; c++)
            for (long u = aux_dft_lo; u <= aux_dft_hi; u++) {
                double complex z = auxdft[c][u - aux_dft_lo] / (double)N;
                double v[2] = {creal(z), cimag(z)};
                fwrite(v, sizeof(double), 2, f);
            }
        fclose(f);
    }
}

static void dump_fields(void) {
    FILE *f = xfopen("fields_final.bin", "wb");
    for (int c = 0; c < 6; c++) {
        double *F = comp_ptr(c);
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j <= Ny; j++) fwrite(&F[ID(i, j, 0)], sizeof(double), Nz, f);
    }
    fclose(f);
}

static void write_meta(double runtime, long nsnap) {
    FILE *f = xfopen("meta.json", "w");
    fprintf(f, "{\n");
    fprintf(f, "  \"units\": \"c = eps0 = mu0 = 1, lambda0 = 1\",\n");
    fprintf(f, "  \"nl\": %d, \"Delta\": %.17g, \"dt\": %.17g, \"S\": %.17g, \"omega0\": %.17g,\n", P.nl, D, dt, P.S, W0);
    fprintf(f, "  \"Lx\": %.17g, \"Lz\": %.17g, \"m\": %d, \"n\": %d, \"pol\": \"%c\",\n", P.Lx, P.Lz, P.m, P.n, P.pol);
    fprintf(f, "  \"Nx\": %d, \"Ny\": %d, \"Nz\": %d, \"j0\": %d, \"j1\": %d, \"npml\": %d, \"sf\": %d, \"tf\": %d,\n",
            Nx, Ny, Nz, j0, j1, P.npml, P.sf, P.tf);
    fprintf(f, "  \"kx\": %.17g, \"ky\": %.17g, \"kz\": %.17g, \"ky_disc\": %.17g, \"ky_cont\": %.17g,\n",
            kx, ky, kz, ky_disc(kx, kz, 1.0), kyc);
    fprintf(f, "  \"Ktilde\": [%.17g, %.17g, %.17g], \"Ktilde_src\": [%.17g, %.17g, %.17g], \"omega_tilde\": %.17g,\n",
            Kt[0], Kt[1], Kt[2], Kt_src[0], Kt_src[1], Kt_src[2], wt);
    fprintf(f, "  \"E0\": [%.17g, %.17g, %.17g], \"H0\": [%.17g, %.17g, %.17g],\n",
            E0v[0], E0v[1], E0v[2], H0v[0], H0v[1], H0v[2]);
    fprintf(f, "  \"offsets\": {");
    for (int c = 0; c < 6; c++)
        fprintf(f, "\"%s\": [%g, %g, %g, %g]%s", CNAME[c], OFF[c][0], OFF[c][1], OFF[c][2], OFF[c][3], c < 5 ? ", " : "},\n");
    fprintf(f, "  \"pml\": {\"m\": %g, \"sig_fac\": %g, \"kappa_max\": %g, \"alpha_max\": %.17g},\n",
            P.pml_m, P.sig_fac, P.kappa_max, P.alpha_max);
    fprintf(f, "  \"eps2\": %.17g, \"y1\": %d, \"ifmode\": \"%c\",\n", P.eps2, P.y1, P.ifmode);
    fprintf(f, "  \"inc\": \"%c\", \"ky_mode\": \"%s\", \"ramp\": \"%c\", \"ramp_T\": %g, \"erf_t0\": %g, \"erf_tau\": %g, \"off_t\": %g, \"na\": %d,\n",
            P.inc, P.ky_cont ? "continuous" : "discrete", P.ramp, P.ramp_T, P.erf_t0, P.erf_tau, P.off_t, P.na);
    fprintf(f, "  \"init\": \"%c\", \"nsteps\": %ld, \"dft0\": %ld, \"dft1\": %ld,\n", P.init, P.nsteps, P.dft0, P.dft1);
    fprintf(f, "  \"aux\": {\"jlo\": %ld, \"L\": %ld, \"ja\": %ld, \"dft_jlo\": %ld, \"dft_jhi\": %ld},\n",
            P.inc == 'a' ? Ajlo : 0, P.inc == 'a' ? AL : 0, P.inc == 'a' ? Ajlo + Aua : 0,
            P.inc == 'a' ? aux_dft_lo + Ajlo : 0, P.inc == 'a' ? aux_dft_hi + Ajlo : -1);
    fprintf(f, "  \"slices\": [");
    for (int s = 0; s < NSL; s++) {
        const char *shape = SL[s].type == 0 ? "Nx,Nz" : SL[s].type == 1 ? "Nx,Ny+1" : "Ny+1,Nz";
        fprintf(f, "{\"name\": \"%s\", \"type\": %d, \"idx\": %d, \"shape\": \"%s\"}%s", SL[s].name, SL[s].type, SL[s].idx, shape,
                s < NSL - 1 ? ", " : "");
    }
    fprintf(f, "],\n");
    fprintf(f, "  \"snap_every\": %d, \"snapcomp\": \"%s\", \"nsnap\": %ld, \"zk\": %d, \"xi\": %d, \"dump\": %d,\n",
            P.snap_every, P.snapcomp, nsnap, P.zk, P.xi, P.dump);
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif
    fprintf(f, "  \"threads\": %d, \"runtime_s\": %.3f\n}\n", nth, runtime);
    fclose(f);
}

/* ------------------------------------------------------------------ main */
int main(int argc, char **argv) {
    parse_args(argc, argv);
    setup();
    MKDIR(P.out);
    if (P.dft0 < 0 || P.dft1 < 0) { P.dft0 = 0; P.dft1 = 0; }
    dft_setup();

    if (P.init == 'a') {
        if (P.npml != 0) die("init=analytic requires npml=0");
        fill_analytic(0.0, -1.0); /* E^0, H^{-1/2} */
    } else if (P.init == 'b') {
        /* Divergence-free packet: Ez depends only on (x, y), Ex only on (y, z), so div E = 0 exactly
           and no static charge is left behind. Carrier ky0 = 5 under a Gaussian of width 1.2, so the
           near-grazing (ky ~ 0) content is ~exp(-(ky0 w)^2/4) ~ 1e-4 in amplitude. */
        double yc = 0.5 * Ny * D, w = 1.2, ky0 = 5.0;
        for (int i = 0; i < Nx; i++)
            for (int j = 1; j < Ny; j++)
                for (int k = 0; k < Nz; k++) {
                    double y = j * D, g = cos(ky0 * (y - yc)) * exp(-(y - yc) * (y - yc) / (w * w));
                    Ez[ID(i, j, k)] = cos(2.0 * PI * i * D / P.Lx) * g;
                    Ex[ID(i, j, k)] = cos(2.0 * PI * k * D / P.Lz) * g;
                }
    }

    partT = xcalloc(Nx, sizeof(double));
    partP = xcalloc(Nx, sizeof(double));
    FILE *lg = xfopen("log.csv", "w");
    fprintf(lg, "n,t,maxE_SF,maxH_SF,W_total,W_phys,divE_TF,divH_TF\n");
    FILE *sn = NULL;
    int snapc = 2;
    for (int c = 0; c < 6; c++) if (!strcmp(P.snapcomp, CNAME[c])) snapc = c;
    if (P.snap_every > 0) {
        if (P.zk < 0) die("snapshots need zk >= 0");
        sn = xfopen("snapshots.bin", "wb");
    }
    long nsnap = 0;

    clock_t c0 = clock();
    time_t w0 = time(NULL);
    for (long n = 0; n < P.nsteps; n++) {
        /* H: n-1/2 -> n+1/2 */
        if (P.inc != '0') incident_E(n);
        update_H();
        if (P.inc != '0') tfsf_H();
        if (P.inc == 'a') aux_update_H(n);
        if (n >= P.dft0 && n < P.dft1) dft_accumulate(3, 5, (n + 0.5) * dt);
        /* E: n -> n+1 */
        if (P.inc != '0') incident_H(n);
        int doW = (P.energy_every > 0 && (n + 1) % P.energy_every == 0);
        update_E(doW);
        if (P.inc != '0') tfsf_E();
        if (P.inc == 'a') aux_update_E(n);
        if (P.init == 'a') dirichlet_y((double)(n + 1));
        if (n + 1 >= P.dft0 && n + 1 < P.dft1) dft_accumulate(0, 2, (n + 1) * dt);

        double mE, mH, Wt = NAN, Wp = NAN, dvE = NAN, dvH = NAN;
        sf_max(&mE, &mH);
        if (doW) energy_sum(&Wt, &Wp); /* W at t = (n + 1/2) dt */
        if (P.div_every > 0 && (n + 1) % P.div_every == 0) div_max(&dvE, &dvH);
        fprintf(lg, "%ld,%.10g,%.9e,%.9e,%.17e,%.17e,%.6e,%.6e\n", n + 1, (n + 1) * dt, mE, mH, Wt, Wp, dvE, dvH);
        if (!isfinite(mE) || !isfinite(mH)) {
            fprintf(stderr, "non-finite field at step %ld\n", n + 1);
            die("instability detected");
        }
        if (sn && (n + 1) % P.snap_every == 0) {
            double *F = comp_ptr(snapc);
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j <= Ny; j++) fwrite(&F[ID(i, j, P.zk)], sizeof(double), 1, sn);
            nsnap++;
        }
        if (P.nsteps >= 10 && (n + 1) % (P.nsteps / 10) == 0) {
            fprintf(stderr, "  step %ld/%ld (%.0f s)\n", n + 1, P.nsteps, difftime(time(NULL), w0));
        }
    }
    fclose(lg);
    if (sn) fclose(sn);
    double runtime = difftime(time(NULL), w0);
    (void)c0;
    write_dft();
    if (P.dump) dump_fields();
    write_meta(runtime, nsnap);
    fprintf(stderr, "done: %ld steps, grid %d x %d x %d, %.0f s\n", P.nsteps, Nx, Ny + 1, Nz, runtime);
    return 0;
}
