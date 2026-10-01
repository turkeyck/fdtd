/*
 * fdtd3d_oblique.c -- 3D Yee FDTD, oblique plane wave, PBC in x/z, CPML in y, TF/SF at y = y0,
 * on a tensor-product (possibly nonuniform) grid.
 *
 * Normalized units: c = eps0 = mu0 = 1, lambda0 = 1 (omega0 = 2*pi).
 * Uniform-grid equations: derivation.md; nonuniform generalization: docs/derivation_nonuniform.md (§ cited below).
 * mesh=uniform (default) builds the legacy equally spaced grid internally; mesh=file reads a grid_gen.py JSON.
 * Both paths use the same update code: the uniform grid is the special case h = d = Delta.
 *
 * Build: make            (gcc, C99, optional OpenMP)
 * Usage: ./fdtd3d_oblique key=value ...   (see parse_args for the list and defaults)
 */
#include <complex.h>
#include <ctype.h>
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
#define MAXDUMP 64

/* ------------------------------------------------------------------ parameters */
typedef struct {
    int nl;                 /* cells per lambda0 (mesh=uniform) */
    double S, Lx, Lz;
    int m, n;
    char pol;               /* 's' or 'p' */
    int npml, sf, tf;       /* near PML, SF length, TF length (cells); Ny = 2*npml + sf + tf (mesh=uniform) */
    double pml_m, sig_fac, kappa_max, alpha_max;
    double eps2;            /* relative permittivity of the half space y >= y1 (mesh=uniform) */
    int y1;                 /* interface offset from j0 in cells (<0: no medium) */
    char ifmode;            /* 'a': arithmetic mean for tangential E on the interface plane, 's': staircase */
    char inc;               /* 'a': aux modal line, 'n': analytic*g(t), '0': none */
    int ky_cont;            /* analytic mode only: 1 -> continuous ky (control) */
    char ramp;              /* 'r': raised cosine, 'e': erf */
    double ramp_T;          /* raised-cosine duration (periods) */
    double erf_t0, erf_tau; /* erf ramp centre and width (periods) */
    double off_t;           /* start of the turn-off ramp (periods); <0: never */
    int na;                 /* aux injection point j_a = j0 - na (mesh=uniform) */
    long nsteps, dft0, dft1;
    int nyplanes, yplanes[MAXPLANES];
    int zk, xi;             /* x-y slice at k = zk, y-z slice at i = xi (<0: none) */
    char init;              /* '0': zero, 'a': analytic + Dirichlet y (debug), 'b': gaussian blob, 'r': random */
    int energy_every, snap_every;
    char snapcomp[4];
    int dump;               /* dump full final fields */
    int auxspan;            /* aux-line DFT recorded up to j0 + auxspan (0: up to Ny) */
    int div_every;          /* log max |div E|, |div H| over the TF interior every N steps (0: off) */
    char out[512];
    /* nonuniform-grid additions (SPEC_nonuniform §17.1) */
    char mesh;              /* 'u': uniform (legacy layout), 'f': grid file */
    char grid[512];
    double dtfac;           /* dt = dtfac * Courant-formula dt */
    double dt_set;          /* >0: dt given explicitly */
    int auxref;             /* 1: exact 1D reduction of the whole main grid (aux_ref line) */
    int ref_every;          /* full-domain main vs aux_ref deviation every N steps (0: off) */
    unsigned long seed;     /* init=r */
    char divop;             /* 'n': nonuniform divergence operator, 'u': (wrong) uniform operator, control only */
    int ndump;
    long dump_at[MAXDUMP];
    char mode;              /* 't': time stepping (default), 'e': power iteration for lambda_max (gate 0-4) */
    int jsrc;               /* inc=j: y index of the current sheet (default j0) */
    int proj;               /* 1: DFT of every y row projected on the transverse Floquet phase (dft_proj.bin) */
    int planar;             /* 1: full-volume DFT kept in memory; per (component, y plane) RMS phase residual
                               against kx x + kz z + const written to dft_planar.bin (gate 1-2) */
    long eig_maxit;
    double eig_tol;
} Params;

/* ------------------------------------------------------------------ globals */
static Params P;
static int Nx, Ny, Nz, SY, SZ, j0, j1, NPlo, NPhi, JA;
static double D, dt, dt_cour, W0, kx, ky, kz, kyc, Da;
static double Kt[3], wt, E0v[3], H0v[3], Kt_src[3];
static size_t NTOT;
static double *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
/* geometry (docs/derivation_nonuniform.md §0): primal nodes, dual (midpoint) coordinates, primal spacings h,
   dual spacings d (x, z periodic; y: d at the two walls = adjacent h) */
static double *xn, *xd, *hxv, *dxv, *yn, *yd, *hyv, *dyv, *zn, *zd, *hzv, *dzv;
static double *epsT, *epsN;
static int xnonuni, znonuni, grating;
/* per-(i, j) materials (stage B, x-y dependent eps; exact copies of epsT/epsN for layered media) and coefficients */
static double *eTx, *eTz, *eNy, *ayx, *ayz, *cEy2;
static int gr_jlo, gr_jhi, gr_ilo, gr_ihi;
static double gr_ridge, gr_groove;
static char grid_hash[80] = "";
/* update coefficients (derivation_nonuniform §1) */
static double *cHy, *ay, *cEyn, *cHz, *zfac;
static double *rExz, *rHxz, *rEzx, *rHzx, *rEyk, *rHyk;
static int NP, *pmlE, *pmlH;
static double *bE, *cE, *ikE, *bH, *cH, *ikH, *sigE, *sigH, *kapE, *kapH, *alpE, *alpH;
static double *psiExy, *psiEzy, *psiHxy, *psiHzy;
static double *cosP1, *sinP1, *cosP2, *sinP2, *cosP3, *sinP3, *cosP4, *sinP4;
/* P1: (x_{i+1/2}, z_k) Ex,Hz; P2: (x_i, z_{k+1/2}) Ez,Hx; P3: (x_i, z_k) Ey; P4: (x_{i+1/2}, z_{k+1/2}) Hy */
static double *ExI, *EzI, *HxI, *HzI;
static double complex incEx, incEz, incHx, incHz;

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
    void *p = calloc(n ? n : 1, s);
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

static void parse_llist(const char *s, long *arr, int *cnt) {
    *cnt = 0;
    char buf[2048];
    strncpy(buf, s, sizeof buf - 1);
    buf[sizeof buf - 1] = 0;
    for (char *t = strtok(buf, ","); t && *cnt < MAXDUMP; t = strtok(NULL, ",")) arr[(*cnt)++] = atol(t);
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
    p.mode = 't'; p.jsrc = -1; p.eig_maxit = 200000; p.eig_tol = 1e-10;
    p.mesh = 'u'; p.dtfac = 1.0; p.dt_set = 0.0; p.auxref = 0; p.ref_every = 0; p.seed = 12345; p.divop = 'n';
    int got_eps = 0;
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
        else if (!strcmp(k, "eps2")) { p.eps2 = atof(v); got_eps = 1; }
        else if (!strcmp(k, "y1")) { p.y1 = atoi(v); got_eps = 1; }
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
        else if (!strcmp(k, "mesh")) p.mesh = (v[0] == 'f') ? 'f' : 'u';
        else if (!strcmp(k, "grid")) { strncpy(p.grid, v, sizeof p.grid - 1); p.mesh = 'f'; }
        else if (!strcmp(k, "dtfac")) p.dtfac = atof(v);
        else if (!strcmp(k, "dt")) p.dt_set = atof(v);
        else if (!strcmp(k, "auxref")) p.auxref = atoi(v);
        else if (!strcmp(k, "ref_every")) p.ref_every = atoi(v);
        else if (!strcmp(k, "seed")) p.seed = strtoul(v, NULL, 10);
        else if (!strcmp(k, "divop")) p.divop = (v[0] == 'u') ? 'u' : 'n';
        else if (!strcmp(k, "dump_at")) parse_llist(v, p.dump_at, &p.ndump);
        else if (!strcmp(k, "mode")) p.mode = (v[0] == 'e') ? 'e' : 't';
        else if (!strcmp(k, "proj")) p.proj = atoi(v);
        else if (!strcmp(k, "jsrc")) p.jsrc = atoi(v);
        else if (!strcmp(k, "planar")) p.planar = atoi(v);
        else if (!strcmp(k, "eig_maxit")) p.eig_maxit = atol(v);
        else if (!strcmp(k, "eig_tol")) p.eig_tol = atof(v);
        else { fprintf(stderr, "unknown key %s\n", k); exit(2); }
    }
    if (p.planar) p.proj = 1;
    if (p.mesh == 'f' && !p.grid[0]) die("mesh=file needs grid=<path.json>");
    if (p.mesh == 'f' && got_eps) die("eps2/y1 cannot be combined with mesh=file (materials come from the grid file)");
    P = p;
}

/* ------------------------------------------------------------------ minimal JSON reader (grid files only) */
typedef struct {
    char key[160];
    double *v;
    int n;
} JEnt;
static JEnt *JE;
static int NJE, CJE;
static const char *jpos;

static void jws(void) { while (*jpos && isspace((unsigned char)*jpos)) jpos++; }

static void jadd(const char *key, double *v, int n) {
    if (NJE == CJE) {
        CJE = CJE ? 2 * CJE : 64;
        JE = realloc(JE, (size_t)CJE * sizeof *JE);
        if (!JE) die("out of memory");
    }
    snprintf(JE[NJE].key, sizeof JE[NJE].key, "%s", key);
    JE[NJE].v = v;
    JE[NJE].n = n;
    NJE++;
}

static void jstring(char *out, size_t cap) {
    if (*jpos != '"') die("grid file: expected string");
    jpos++;
    size_t q = 0;
    while (*jpos && *jpos != '"') {
        if (*jpos == '\\' && jpos[1]) jpos++;
        if (out && q + 1 < cap) out[q++] = *jpos;
        jpos++;
    }
    if (*jpos != '"') die("grid file: unterminated string");
    jpos++;
    if (out) out[q] = 0;
}

static void jvalue(const char *key);

static void jobject(const char *prefix) {
    jpos++;
    jws();
    if (*jpos == '}') { jpos++; return; }
    for (;;) {
        char k[128], child[160];
        jws();
        jstring(k, sizeof k);
        jws();
        if (*jpos != ':') die("grid file: expected ':'");
        jpos++;
        if (prefix[0]) snprintf(child, sizeof child, "%s.%s", prefix, k);
        else snprintf(child, sizeof child, "%s", k);
        jvalue(child);
        jws();
        if (*jpos == ',') { jpos++; continue; }
        if (*jpos == '}') { jpos++; return; }
        die("grid file: expected ',' or '}'");
    }
}

static void jarray(const char *key) {
    jpos++;
    int cap = 64, n = 0, numeric = 1;
    double *buf = xcalloc((size_t)cap, sizeof(double));
    jws();
    if (*jpos == ']') { jpos++; jadd(key, buf, 0); return; }
    for (;;) {
        jws();
        if (*jpos == '-' || *jpos == '+' || isdigit((unsigned char)*jpos) || *jpos == '.') {
            char *end;
            double x = strtod(jpos, &end);
            if (end == jpos) die("grid file: bad number");
            jpos = end;
            if (n == cap) {
                cap *= 2;
                if (cap > 20000000) die("grid file: array too long");
                buf = realloc(buf, (size_t)cap * sizeof(double));
                if (!buf) die("out of memory");
            }
            buf[n++] = x;
        } else {
            numeric = 0;
            jvalue("");
        }
        jws();
        if (*jpos == ',') { jpos++; continue; }
        if (*jpos == ']') { jpos++; break; }
        die("grid file: expected ',' or ']'");
    }
    if (numeric && key[0]) jadd(key, buf, n);
    else free(buf);
}

static void jvalue(const char *key) {
    jws();
    if (*jpos == '{') jobject(key);
    else if (*jpos == '[') jarray(key);
    else if (*jpos == '"') {
        if (!strcmp(key, "hash")) jstring(grid_hash, sizeof grid_hash);
        else jstring(NULL, 0);
    } else if (!strncmp(jpos, "true", 4)) { jpos += 4; double *v = xcalloc(1, sizeof(double)); *v = 1; if (key[0]) jadd(key, v, 1); }
    else if (!strncmp(jpos, "false", 5)) { jpos += 5; double *v = xcalloc(1, sizeof(double)); if (key[0]) jadd(key, v, 1); }
    else if (!strncmp(jpos, "null", 4)) { jpos += 4; }
    else {
        char *end;
        double x = strtod(jpos, &end);
        if (end == jpos) die("grid file: unexpected character");
        jpos = end;
        double *v = xcalloc(1, sizeof(double));
        *v = x;
        if (key[0]) jadd(key, v, 1);
    }
}

static double *jget(const char *key, int *n, int required) {
    for (int q = 0; q < NJE; q++)
        if (!strcmp(JE[q].key, key)) { if (n) *n = JE[q].n; return JE[q].v; }
    if (required) {
        fprintf(stderr, "ERROR: grid file: missing key \"%s\"\n", key);
        exit(2);
    }
    if (n) *n = 0;
    return NULL;
}

static void read_grid_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) die("grid file: cannot open");
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz <= 0 || sz > 400000000L) die("grid file: bad size");
    char *txt = xcalloc((size_t)sz + 1, 1);
    if (fread(txt, 1, (size_t)sz, f) != (size_t)sz) die("grid file: read error");
    fclose(f);
    jpos = txt;
    jws();
    if (*jpos != '{') die("grid file: not a JSON object");
    jobject("");
    free(txt);
}

/* ------------------------------------------------------------------ theory (derivation_nonuniform §5) */
/* discrete ky in a uniform sub-region of spacing h and index nidx: sin^2(ky h/2) = (h/2)^2 (n^2 w~^2 - Kx~^2 - Kz~^2) */
static double ky_local(double h, double nidx) {
    double dx = hxv[0], dz = hzv[0];
    double a = (P.mesh == 'u' && P.dt_set <= 0.0 && P.dtfac == 1.0) ? nidx / P.S * sin(W0 * dt / 2.0)
                                                                     : nidx * (h / dt) * sin(W0 * dt / 2.0);
    double q = pow(a, 2) - pow(sin(kx * dx / 2.0) * (h / dx), 2) - pow(sin(kz * dz / 2.0) * (h / dz), 2);
    if (!(q > 0.0 && q <= 1.0)) {
        fprintf(stderr, "sin^2(ky h/2) = %.6g outside (0,1] for h = %.6g\n", q, h);
        die("no propagating discrete ky for this (m, n, h, dt)");
    }
    return 2.0 / h * asin(sqrt(q));
}

static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(const double *a) { return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

/* E0, H0 for the wave vector whose y component is kyy on the aux-source spacing Da */
static void amplitudes(double kyy, double *K, double *E, double *H) {
    double yh[3] = {0, 1, 0}, s[3], p[3];
    K[0] = 2.0 / hxv[0] * sin(kx * hxv[0] / 2.0);
    K[1] = 2.0 / Da * sin(kyy * Da / 2.0);
    K[2] = 2.0 / hzv[0] * sin(kz * hzv[0] / 2.0);
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

/* global y geometry, extended uniformly beyond [0, Ny] with the end spacings (aux line, §7.2) */
static double gh(long j) { return j < 0 ? hyv[0] : (j >= Ny ? hyv[Ny - 1] : hyv[j]); }
static double gy(long j) {
    if (j < 0) return yn[0] + (double)j * hyv[0];
    if (j > Ny) return yn[Ny] + (double)(j - Ny) * hyv[Ny - 1];
    return yn[j];
}
static double gyd(long j) {
    if (j < 0) return yn[0] + ((double)j + 0.5) * hyv[0];
    if (j >= Ny) return yn[Ny] + ((double)(j - Ny) + 0.5) * hyv[Ny - 1];
    return yd[j];
}
static double gd(long j) {
    double a = gh(j - 1), b = gh(j);
    return a == b ? a : 0.5 * (a + b);
}

/* complex amplitude of component c at global y-index j (its own y offset) and time step tn (own offset) */
static double complex inc_amp(int c, long j, double tn, const double *E, const double *H, double kyy) {
    double y = (OFF[c][1] == 0.0) ? gy(j) : gyd(j);
    double t = (tn + OFF[c][3]) * dt;
    double a = (c < 3) ? E[c] : H[c - 3];
    return envelope(t) * a * cexp(I * (kyy * y - W0 * t));
}

/* ------------------------------------------------------------------ aux modal line (derivation_nonuniform §3) */
static long AL, Ajlo, Aua;
static double complex *aEx, *aEy, *aEz, *aHx, *aHy, *aHz;
static double *acH, *acE;                       /* dt/h (H side, per dual u), dt/d (E side, per primal u) */
static double complex *ikxH, *ikzH, *ikxE, *ikzE; /* i K~ h and i K~ d: acH*ikH = dt i K~ (y and phasor terms split) */

static void aux_init(void) {
    long ja = JA;
    long half = P.nsteps / 2 + 20;
    Ajlo = ja - half;
    long jhi = j0 + half + (P.auxspan > 0 ? P.auxspan : Ny - j0); /* no end reflection reaches a recorded node */
    AL = jhi - Ajlo;
    Aua = ja - Ajlo;
    aEx = xcalloc(AL + 1, sizeof *aEx); aEz = xcalloc(AL + 1, sizeof *aEz); aHy = xcalloc(AL + 1, sizeof *aHy);
    aEy = xcalloc(AL, sizeof *aEy); aHx = xcalloc(AL, sizeof *aHx); aHz = xcalloc(AL, sizeof *aHz);
    acH = xcalloc(AL + 1, sizeof(double)); acE = xcalloc(AL + 1, sizeof(double));
    ikxH = xcalloc(AL + 1, sizeof *ikxH); ikzH = xcalloc(AL + 1, sizeof *ikzH);
    ikxE = xcalloc(AL + 1, sizeof *ikxE); ikzE = xcalloc(AL + 1, sizeof *ikzE);
    for (long u = 0; u <= AL; u++) {
        long j = Ajlo + u;
        double h = gh(j), d = gd(j);
        acH[u] = dt / h;
        acE[u] = dt / d;
        ikxH[u] = I * Kt[0] * h; ikzH[u] = I * Kt[2] * h;
        ikxE[u] = I * Kt[0] * d; ikzE[u] = I * Kt[2] * d;
    }
}

static void aux_update_H(long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 0) lo = 0;
    if (hi > AL - 1) hi = AL - 1;
    for (long u = lo; u <= hi; u++) {
        aHx[u] -= acH[u] * ((aEz[u + 1] - aEz[u]) - ikzH[u] * aEy[u]);
        aHz[u] -= acH[u] * (ikxH[u] * aEy[u] - (aEx[u + 1] - aEx[u]));
    }
    for (long u = lo; u <= hi + 1 && u <= AL; u++) aHy[u] -= acE[u] * (ikzE[u] * aEx[u] - ikxE[u] * aEz[u]);
    /* 1D TF/SF at j_a (derivation.md §3 (a), (b)); the source region is uniform, spacing Da */
    long ja = Ajlo + Aua;
    aHx[Aua - 1] += acH[Aua - 1] * inc_amp(2, ja, (double)n, E0v, H0v, ky);
    aHz[Aua - 1] -= acH[Aua - 1] * inc_amp(0, ja, (double)n, E0v, H0v, ky);
}

static void aux_update_E(long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 1) lo = 1;
    if (hi > AL - 1) hi = AL - 1;
    for (long u = lo; u <= hi; u++) { /* eps0 = 1; the aux line is vacuum */
        aEx[u] += acE[u] * ((aHz[u] - aHz[u - 1]) - ikzE[u] * aHy[u]);
        aEz[u] += acE[u] * (ikxE[u] * aHy[u] - (aHx[u] - aHx[u - 1]));
    }
    for (long u = lo - 1; u <= hi && u < AL; u++) aEy[u] += acH[u] * (ikzH[u] * aHx[u] - ikxH[u] * aHz[u]);
    long ja = Ajlo + Aua;
    aEx[Aua] -= acE[Aua] * inc_amp(5, ja - 1, (double)n, E0v, H0v, ky);
    aEz[Aua] += acE[Aua] * inc_amp(3, ja - 1, (double)n, E0v, H0v, ky);
}

/* ------------------------------------------------------------------ inc=m: transverse discrete-mode injection (D22)
   On a nonuniform periodic axis the sampled exp(i k x) is not a discrete eigenfunction: the operator D_d D_p
   (D_p u = (u_{i+1} - u_i)/h_i at dual nodes, D_d g = (g_i - g_{i-1})/d_i at primal nodes, the weights of the main
   update) splits the +-m pair into two real modes D_d D_p phi = -kappa^2 phi. A product mode phi_a(x) phi_b(z)
   separates exactly: every transverse difference acts as i kappa when primal-node components carry phi and
   dual-node components carry psi = D_p phi / (i kappa). Each product mode is fed by its own 1D aux line with
   K~ = (kappa_x, K~y, kappa_z) and weighted by the projection of exp(i(kx x + kz z)) on it. A uniform axis keeps the
   exact exponential (one mode, kappa = K~). The injected field is then an exact discrete solution. */
typedef struct {
    int nm;                          /* number of kept modes: 1 (uniform axis or m = 0) or 2 */
    double kap[2], res[2], rest;     /* kappa, relative eigen-residual, power of exp(ikx) outside the kept modes */
    double complex c[2];             /* weights <phi_q, exp(i k x)>_d (phi orthonormal under the d weights) */
    double complex *pp[2], *pd[2];   /* profiles at primal / dual nodes */
} AxisModes;
static AxisModes MX, MZ;

typedef struct {
    double complex *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
    double complex *ikxH, *ikzH, *ikxE, *ikzE;
    double K[3], E0[3], H0[3], ky;
    double complex w;                /* product weight c_a c_b */
    double complex *T1, *T2;         /* w * profile at (x dual, z primal) [Ex, Hz] and (x primal, z dual) [Ez, Hx] */
    int a, b;
} MAux;
static MAux MA[4];
static int NMA;
/* y geometry of the modal aux lines: the main-grid nodes up to j0, then uniform with the j0 spacing (j0 +- 5 is
   uniform by the grid rules). The incident is then a pure forward discrete wave: reflections of the main grid's own
   y grading beyond j0 belong to the scattered field (an aux copy of the whole main grid would carry them in the
   "incident" and remove them from the SF region, biasing R measured there by ~2|r||echo|). */
static double *acHm, *acEm;
static double ghm(long j) { return j <= j0 ? gh(j) : hyv[j0]; }
static double gdm(long j) {
    double a = ghm(j - 1), b = ghm(j);
    return a == b ? a : 0.5 * (a + b);
}

static void lu_factor(int N, double *A, int *piv) { /* in place, partial pivoting, row-major */
    for (int c = 0; c < N; c++) {
        int p = c;
        for (int r = c + 1; r < N; r++) if (fabs(A[(size_t)r * N + c]) > fabs(A[(size_t)p * N + c])) p = r;
        piv[c] = p;
        if (p != c)
            for (int q = 0; q < N; q++) {
                double t = A[(size_t)c * N + q]; A[(size_t)c * N + q] = A[(size_t)p * N + q]; A[(size_t)p * N + q] = t;
            }
        double d = A[(size_t)c * N + c];
        if (d == 0.0) die("inc=m: singular shifted mode matrix");
        for (int r = c + 1; r < N; r++) {
            double f = A[(size_t)r * N + c] / d;
            A[(size_t)r * N + c] = f;
            for (int q = c + 1; q < N; q++) A[(size_t)r * N + q] -= f * A[(size_t)c * N + q];
        }
    }
}

static void lu_solve(int N, const double *A, const int *piv, double *b) {
    /* whole rows (with their multipliers) were swapped in lu_factor: permute first, then substitute */
    for (int c = 0; c < N; c++)
        if (piv[c] != c) { double t = b[c]; b[c] = b[piv[c]]; b[piv[c]] = t; }
    for (int c = 0; c < N; c++)
        for (int r = c + 1; r < N; r++) b[r] -= A[(size_t)r * N + c] * b[c];
    for (int r = N - 1; r >= 0; r--) {
        double s = b[r];
        for (int q = r + 1; q < N; q++) s -= A[(size_t)r * N + q] * b[q];
        b[r] = s / A[(size_t)r * N + r];
    }
}

/* (A u)_i = -[(u_{i+1} - u_i)/h_i - (u_i - u_{i-1})/h_{i-1}] = -d_i (D_d D_p u)_i  (periodic) */
static void apply_A(int N, const double *h, const double *u, double *out) {
    for (int i = 0; i < N; i++) {
        int ip = (i + 1) % N, im = (i + N - 1) % N;
        out[i] = -((u[ip] - u[i]) / h[i] - (u[i] - u[im]) / h[im]);
    }
}

static void axis_modes(int N, const double *h, const double *d, const double *xn, const double *xd, double L, int m,
                       int nonuni, AxisModes *M) {
    double k = 2.0 * PI * m / L;
    for (int q = 0; q < 2; q++) { M->pp[q] = xcalloc(N, sizeof(double complex)); M->pd[q] = xcalloc(N, sizeof(double complex)); }
    M->rest = 0.0;
    if (!nonuni || m == 0) { /* exact: the discrete exponential (or the constant) is an eigenfunction */
        M->nm = 1;
        M->kap[0] = (m == 0) ? 0.0 : 2.0 / h[0] * sin(k * h[0] / 2.0);
        M->c[0] = 1.0;
        M->res[0] = 0.0;
        for (int i = 0; i < N; i++) { M->pp[0][i] = cexp(I * k * xn[i]); M->pd[0][i] = cexp(I * k * xd[i]); }
        return;
    }
    /* shift-invert subspace iteration on A v = lambda W v (W = diag d) for the pair nearest sigma = k^2 */
    double sigma = k * k;
    double *B = xcalloc((size_t)N * N, sizeof(double));
    int *piv = xcalloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        int ip = (i + 1) % N, im = (i + N - 1) % N;
        B[(size_t)i * N + i] += 1.0 / h[i] + 1.0 / h[im] - sigma * d[i];
        B[(size_t)i * N + ip] -= 1.0 / h[i];
        B[(size_t)i * N + im] -= 1.0 / h[im];
    }
    lu_factor(N, B, piv);
    double *V[2], *T = xcalloc(N, sizeof(double)), *AV[2], lam[2] = {0.0, 0.0};
    for (int q = 0; q < 2; q++) { V[q] = xcalloc(N, sizeof(double)); AV[q] = xcalloc(N, sizeof(double)); }
    for (int i = 0; i < N; i++) { V[0][i] = cos(k * xn[i]); V[1][i] = sin(k * xn[i]); }
    for (int it = 0; it < 60; it++) {
        for (int q = 0; q < 2; q++) {
            for (int i = 0; i < N; i++) V[q][i] *= d[i];
            lu_solve(N, B, piv, V[q]);
        }
        for (int q = 0; q < 2; q++) { /* W-orthonormalize (modified Gram-Schmidt, twice) */
            for (int rep = 0; rep < 2; rep++)
                for (int p = 0; p < q; p++) {
                    double s = 0.0;
                    for (int i = 0; i < N; i++) s += d[i] * V[p][i] * V[q][i];
                    for (int i = 0; i < N; i++) V[q][i] -= s * V[p][i];
                }
            double nn = 0.0;
            for (int i = 0; i < N; i++) nn += d[i] * V[q][i] * V[q][i];
            nn = sqrt(nn);
            for (int i = 0; i < N; i++) V[q][i] /= nn;
        }
        for (int q = 0; q < 2; q++) apply_A(N, h, V[q], AV[q]);
        double a00 = 0.0, a01 = 0.0, a11 = 0.0; /* Rayleigh-Ritz in the 2D subspace */
        for (int i = 0; i < N; i++) { a00 += V[0][i] * AV[0][i]; a01 += V[0][i] * AV[1][i]; a11 += V[1][i] * AV[1][i]; }
        double th = 0.5 * atan2(2.0 * a01, a00 - a11), cs = cos(th), sn = sin(th);
        for (int i = 0; i < N; i++) {
            double v0 = V[0][i], v1 = V[1][i];
            V[0][i] = cs * v0 + sn * v1; V[1][i] = -sn * v0 + cs * v1;
        }
        double rmax = 0.0;
        for (int q = 0; q < 2; q++) {
            apply_A(N, h, V[q], AV[q]);
            double num = 0.0, den = 0.0;
            for (int i = 0; i < N; i++) num += V[q][i] * AV[q][i];
            lam[q] = num;                       /* V W-orthonormal */
            double r = 0.0;
            for (int i = 0; i < N; i++) { r = fmax(r, fabs(AV[q][i] - lam[q] * d[i] * V[q][i])); den = fmax(den, fabs(AV[q][i])); }
            M->res[q] = r / den;
            rmax = fmax(rmax, M->res[q]);
        }
        if (it >= 2 && rmax < 1e-13) break;
    }
    M->nm = 2;
    double pw = 0.0;
    for (int q = 0; q < 2; q++) {
        if (!(lam[q] > 0.0)) die("inc=m: non-positive mode eigenvalue");
        M->kap[q] = sqrt(lam[q]);
        double complex cq = 0.0;
        for (int i = 0; i < N; i++) cq += d[i] * V[q][i] * cexp(I * k * xn[i]);
        M->c[q] = cq;
        pw += creal(cq * conj(cq));
        for (int i = 0; i < N; i++) {
            M->pp[q][i] = V[q][i];
            M->pd[q][i] = (V[q][(i + 1) % N] - V[q][i]) / h[i] / (I * M->kap[q]);
        }
    }
    M->rest = 1.0 - pw / L;                     /* sum_i d_i |exp(ikx)|^2 = L */
    for (int q = 0; q < 2; q++) { free(V[q]); free(AV[q]); }
    free(T); free(B); free(piv);
}

static void modal_init(void) {
    axis_modes(Nx, hxv, dxv, xn, xd, P.Lx, P.m, xnonuni, &MX);
    axis_modes(Nz, hzv, dzv, zn, zd, P.Lz, P.n, znonuni, &MZ);
    for (int q = 0; q < MX.nm; q++) if (MX.res[q] > 1e-10) die("inc=m: x mode solve did not converge");
    for (int q = 0; q < MZ.nm; q++) if (MZ.res[q] > 1e-10) die("inc=m: z mode solve did not converge");
    acHm = xcalloc(AL + 1, sizeof(double)); acEm = xcalloc(AL + 1, sizeof(double));
    for (long u = 0; u <= AL; u++) { acHm[u] = dt / ghm(Ajlo + u); acEm[u] = dt / gdm(Ajlo + u); }
    NMA = 0;
    for (int a = 0; a < MX.nm; a++)
        for (int b = 0; b < MZ.nm; b++) {
            MAux *A = &MA[NMA++];
            A->a = a; A->b = b;
            A->K[0] = MX.kap[a]; A->K[2] = MZ.kap[b];
            double q = pow(Da * wt / 2.0, 2) - pow(A->K[0] * Da / 2.0, 2) - pow(A->K[2] * Da / 2.0, 2);
            if (!(q > 0.0 && q <= 1.0)) die("inc=m: no propagating discrete ky for a mode");
            A->ky = 2.0 / Da * asin(sqrt(q));
            A->K[1] = 2.0 / Da * sqrt(q);
            double yh[3] = {0, 1, 0}, s[3], p[3];
            cross(A->K, yh, s);
            double ns = norm3(s);
            if (ns == 0.0) die("normal incidence (m = n = 0) is not supported: s/p basis degenerate");
            for (int c = 0; c < 3; c++) s[c] /= ns;
            cross(s, A->K, p);
            double np = norm3(p);
            for (int c = 0; c < 3; c++) p[c] /= np;
            for (int c = 0; c < 3; c++) A->E0[c] = (P.pol == 's') ? s[c] : p[c];
            cross(A->K, A->E0, A->H0);
            for (int c = 0; c < 3; c++) A->H0[c] /= wt;
            A->w = MX.c[a] * MZ.c[b];
            A->Ex = xcalloc(AL + 1, sizeof(double complex)); A->Ez = xcalloc(AL + 1, sizeof(double complex));
            A->Hy = xcalloc(AL + 1, sizeof(double complex)); A->Ey = xcalloc(AL, sizeof(double complex));
            A->Hx = xcalloc(AL, sizeof(double complex)); A->Hz = xcalloc(AL, sizeof(double complex));
            A->ikxH = xcalloc(AL + 1, sizeof(double complex)); A->ikzH = xcalloc(AL + 1, sizeof(double complex));
            A->ikxE = xcalloc(AL + 1, sizeof(double complex)); A->ikzE = xcalloc(AL + 1, sizeof(double complex));
            for (long u = 0; u <= AL; u++) {
                long j = Ajlo + u;
                double h = ghm(j), dd = gdm(j);
                A->ikxH[u] = I * A->K[0] * h; A->ikzH[u] = I * A->K[2] * h;
                A->ikxE[u] = I * A->K[0] * dd; A->ikzE[u] = I * A->K[2] * dd;
            }
            A->T1 = xcalloc((size_t)Nx * Nz, sizeof(double complex));
            A->T2 = xcalloc((size_t)Nx * Nz, sizeof(double complex));
            for (int i = 0; i < Nx; i++)
                for (int k = 0; k < Nz; k++) {
                    A->T1[(size_t)i * Nz + k] = A->w * MX.pd[a][i] * MZ.pp[b][k];
                    A->T2[(size_t)i * Nz + k] = A->w * MX.pp[a][i] * MZ.pd[b][k];
                }
        }
}

static void maux_update_H(MAux *A, long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 0) lo = 0;
    if (hi > AL - 1) hi = AL - 1;
    for (long u = lo; u <= hi; u++) {
        A->Hx[u] -= acHm[u] * ((A->Ez[u + 1] - A->Ez[u]) - A->ikzH[u] * A->Ey[u]);
        A->Hz[u] -= acHm[u] * (A->ikxH[u] * A->Ey[u] - (A->Ex[u + 1] - A->Ex[u]));
    }
    for (long u = lo; u <= hi + 1 && u <= AL; u++) A->Hy[u] -= acEm[u] * (A->ikzE[u] * A->Ex[u] - A->ikxE[u] * A->Ez[u]);
    long ja = Ajlo + Aua;
    A->Hx[Aua - 1] += acHm[Aua - 1] * inc_amp(2, ja, (double)n, A->E0, A->H0, A->ky);
    A->Hz[Aua - 1] -= acHm[Aua - 1] * inc_amp(0, ja, (double)n, A->E0, A->H0, A->ky);
}

static void maux_update_E(MAux *A, long n) {
    long lo = Aua - n - 3, hi = Aua + n + 3;
    if (lo < 1) lo = 1;
    if (hi > AL - 1) hi = AL - 1;
    for (long u = lo; u <= hi; u++) {
        A->Ex[u] += acEm[u] * ((A->Hz[u] - A->Hz[u - 1]) - A->ikzE[u] * A->Hy[u]);
        A->Ez[u] += acEm[u] * (A->ikxE[u] * A->Hy[u] - (A->Hx[u] - A->Hx[u - 1]));
    }
    for (long u = lo - 1; u <= hi && u < AL; u++) A->Ey[u] += acHm[u] * (A->ikzH[u] * A->Hx[u] - A->ikxH[u] * A->Hz[u]);
    long ja = Ajlo + Aua;
    A->Ex[Aua] -= acEm[Aua] * inc_amp(5, ja - 1, (double)n, A->E0, A->H0, A->ky);
    A->Ez[Aua] += acEm[Aua] * inc_amp(3, ja - 1, (double)n, A->E0, A->H0, A->ky);
}

static FILE *xfopen(const char *name, const char *mode);

static void modal_write_axis(FILE *f, const char *name, const AxisModes *M, int N) {
    fprintf(f, "  \"%s\": {\"nm\": %d, \"rest_power\": %.17g, \"modes\": [", name, M->nm, M->rest);
    for (int q = 0; q < M->nm; q++) {
        fprintf(f, "{\"kappa\": %.17g, \"residual\": %.3e, \"c\": [%.17g, %.17g], \"primal\": [", M->kap[q], M->res[q],
                creal(M->c[q]), cimag(M->c[q]));
        for (int i = 0; i < N; i++) fprintf(f, "[%.17g, %.17g]%s", creal(M->pp[q][i]), cimag(M->pp[q][i]), i < N - 1 ? ", " : "");
        fprintf(f, "], \"dual\": [");
        for (int i = 0; i < N; i++) fprintf(f, "[%.17g, %.17g]%s", creal(M->pd[q][i]), cimag(M->pd[q][i]), i < N - 1 ? ", " : "");
        fprintf(f, "]}%s", q < M->nm - 1 ? ", " : "");
    }
    fprintf(f, "]}");
}

static void modal_write(void) { /* modes.json: profiles and per-product-mode aux data (B1-3 analysis) */
    FILE *f = xfopen("modes.json", "w");
    fprintf(f, "{\n");
    modal_write_axis(f, "x", &MX, Nx);
    fprintf(f, ",\n");
    modal_write_axis(f, "z", &MZ, Nz);
    fprintf(f, ",\n  \"product\": [");
    for (int t = 0; t < NMA; t++) {
        const MAux *A = &MA[t];
        fprintf(f, "{\"a\": %d, \"b\": %d, \"K\": [%.17g, %.17g, %.17g], \"ky\": %.17g, \"w\": [%.17g, %.17g], "
                   "\"E0\": [%.17g, %.17g, %.17g], \"H0\": [%.17g, %.17g, %.17g]}%s",
                A->a, A->b, A->K[0], A->K[1], A->K[2], A->ky, creal(A->w), cimag(A->w), A->E0[0], A->E0[1], A->E0[2],
                A->H0[0], A->H0[1], A->H0[2], t < NMA - 1 ? ", " : "");
    }
    fprintf(f, "]\n}\n");
    fclose(f);
}

/* ------------------------------------------------------------------ aux_ref: exact 1D reduction of the main grid
   (same y nodes, eps, CPML; TF/SF at j0 fed by the aux line). Main field = Re[aux_ref * exp(i(kx x + kz z))]. */
static double complex *rEx, *rEy, *rEz, *rHx, *rHy, *rHz, *rpEx, *rpEz, *rpHx, *rpHz;

static void auxref_init(void) {
    rEx = xcalloc(SY, sizeof *rEx); rEy = xcalloc(SY, sizeof *rEy); rEz = xcalloc(SY, sizeof *rEz);
    rHx = xcalloc(SY, sizeof *rHx); rHy = xcalloc(SY, sizeof *rHy); rHz = xcalloc(SY, sizeof *rHz);
    rpEx = xcalloc(NP, sizeof *rpEx); rpEz = xcalloc(NP, sizeof *rpEz);
    rpHx = xcalloc(NP, sizeof *rpHx); rpHz = xcalloc(NP, sizeof *rpHz);
}

static void auxref_update_H(void) {
    const double complex iKx = I * Kt[0], iKz = I * Kt[2];
    for (int j = 0; j <= Ny; j++) rHy[j] -= dt * (iKz * rEx[j] - iKx * rEz[j]);
    for (int j = 0; j < Ny; j++) {
        double complex dEz = rEz[j + 1] - rEz[j], dEx = rEx[j + 1] - rEx[j];
        int s = pmlH[j];
        if (s >= 0) {
            rpHx[s] = bH[j] * rpHx[s] + cH[j] * dEz;
            rpHz[s] = bH[j] * rpHz[s] + cH[j] * dEx;
            dEz = ikH[j] * dEz + rpHx[s];
            dEx = ikH[j] * dEx + rpHz[s];
        }
        rHx[j] -= cHy[j] * (dEz - iKz * hyv[j] * rEy[j]);
        rHz[j] -= cHy[j] * (iKx * hyv[j] * rEy[j] - dEx);
    }
    rHx[j0 - 1] += cHy[j0 - 1] * incEz;
    rHz[j0 - 1] -= cHy[j0 - 1] * incEx;
}

static void auxref_update_E(void) {
    const double complex iKx = I * Kt[0], iKz = I * Kt[2];
    for (int j = 0; j < Ny; j++) rEy[j] += (dt / epsN[j]) * (iKz * rHx[j] - iKx * rHz[j]);
    for (int j = 1; j < Ny; j++) {
        double complex dHz = rHz[j] - rHz[j - 1], dHx = rHx[j] - rHx[j - 1];
        int s = pmlE[j];
        if (s >= 0) {
            rpEx[s] = bE[j] * rpEx[s] + cE[j] * dHz;
            rpEz[s] = bE[j] * rpEz[s] + cE[j] * dHx;
            dHz = ikE[j] * dHz + rpEx[s];
            dHx = ikE[j] * dHx + rpEz[s];
        }
        rEx[j] += ay[j] * (dHz - iKz * dyv[j] * rHy[j]);
        rEz[j] += ay[j] * (iKx * dyv[j] * rHy[j] - dHx);
    }
    rEx[j0] -= ay[j0] * incHz;
    rEz[j0] += ay[j0] * incHx;
}

/* ------------------------------------------------------------------ geometry */
static double maxreldiff(const double *h, int lo, int hi, double ref) {
    double m = 0.0;
    for (int q = lo; q < hi; q++) m = fmax(m, fabs(h[q] / ref - 1.0));
    return m;
}

static void periodic_axis(int N, const double *h, double **node, double **dual, double **dd) {
    *node = xcalloc(N + 1, sizeof(double)); *dual = xcalloc(N, sizeof(double)); *dd = xcalloc(N, sizeof(double));
    double s = 0.0, c = 0.0; /* Kahan summation */
    (*node)[0] = 0.0;
    for (int q = 0; q < N; q++) {
        double yv = h[q] - c, t = s + yv;
        c = (t - s) - yv;
        s = t;
        (*node)[q + 1] = s;
    }
    for (int q = 0; q < N; q++) {
        (*dual)[q] = 0.5 * ((*node)[q] + (*node)[q + 1]);
        double a = h[(q + N - 1) % N], b = h[q];
        (*dd)[q] = a == b ? a : 0.5 * (a + b);
    }
}

static void build_uniform_mesh(void) {
    D = 1.0 / P.nl;
    Nx = (int)lround(P.Lx * P.nl);
    Nz = (int)lround(P.Lz * P.nl);
    if (fabs(Nx * D - P.Lx) > 1e-12 * P.Lx || fabs(Nz * D - P.Lz) > 1e-12 * P.Lz)
        die("Lx and Lz must be integer multiples of Delta");
    Ny = 2 * P.npml + P.sf + P.tf;
    j0 = P.npml + P.sf;
    NPlo = NPhi = P.npml;
    JA = j0 - P.na;
    j1 = (P.y1 >= 0) ? j0 + P.y1 : -1;
    SY = Ny + 1;
    hxv = xcalloc(Nx, sizeof(double)); dxv = xcalloc(Nx, sizeof(double));
    xn = xcalloc(Nx + 1, sizeof(double)); xd = xcalloc(Nx, sizeof(double));
    hzv = xcalloc(Nz, sizeof(double)); dzv = xcalloc(Nz, sizeof(double));
    zn = xcalloc(Nz + 1, sizeof(double)); zd = xcalloc(Nz, sizeof(double));
    hyv = xcalloc(SY, sizeof(double)); dyv = xcalloc(SY, sizeof(double));
    yn = xcalloc(SY, sizeof(double)); yd = xcalloc(SY, sizeof(double));
    for (int i = 0; i <= Nx; i++) xn[i] = i * D;
    for (int i = 0; i < Nx; i++) { xd[i] = (i + 0.5) * D; hxv[i] = D; dxv[i] = D; }
    for (int k = 0; k <= Nz; k++) zn[k] = k * D;
    for (int k = 0; k < Nz; k++) { zd[k] = (k + 0.5) * D; hzv[k] = D; dzv[k] = D; }
    for (int j = 0; j <= Ny; j++) { yn[j] = j * D; dyv[j] = D; }
    for (int j = 0; j < Ny; j++) { yd[j] = (j + 0.5) * D; hyv[j] = D; }
    /* material (derivation.md §8): Ex,Ez at integer j; Ey at j+1/2 never lies on an integer-y interface */
    epsT = xcalloc(SY, sizeof(double)); epsN = xcalloc(SY, sizeof(double));
    for (int j = 0; j <= Ny; j++) {
        double et = 1.0, en = 1.0;
        if (j1 >= 0) {
            if (j > j1) et = P.eps2;
            else if (j == j1) et = (P.ifmode == 'a') ? 0.5 * (1.0 + P.eps2) : P.eps2;
            if (j >= j1) en = P.eps2; /* Ey at (j+1/2) >= j1 */
        }
        epsT[j] = et;
        epsN[j] = en;
    }
}

static void build_file_mesh(void) {
    read_grid_file(P.grid);
    int nx, ny, nz, ne, nt, nn;
    double *hx = jget("x.h", &nx, 1), *hy = jget("y.h", &ny, 1), *hz = jget("z.h", &nz, 1);
    double *ey = jget("eps_y", &ne, 1);
    double *et = jget("eps_t", &nt, 0), *en = jget("eps_n", &nn, 0);
    if (ne != ny) die("grid file: len(eps_y) != len(y.h)");
    Nx = nx; Ny = ny; Nz = nz; SY = Ny + 1;
    NPlo = (int)lround(*jget("zones.npml_lo", NULL, 1));
    NPhi = (int)lround(*jget("zones.npml_hi", NULL, 1));
    j0 = (int)lround(*jget("zones.j0", NULL, 1));
    JA = (int)lround(*jget("zones.ja", NULL, 1));
    j1 = -1;
    hxv = xcalloc(Nx, sizeof(double)); hzv = xcalloc(Nz, sizeof(double)); hyv = xcalloc(SY, sizeof(double));
    memcpy(hxv, hx, Nx * sizeof(double)); memcpy(hzv, hz, Nz * sizeof(double)); memcpy(hyv, hy, Ny * sizeof(double));
    periodic_axis(Nx, hxv, &xn, &xd, &dxv);
    periodic_axis(Nz, hzv, &zn, &zd, &dzv);
    yn = xcalloc(SY, sizeof(double)); yd = xcalloc(SY, sizeof(double)); dyv = xcalloc(SY, sizeof(double));
    double s = 0.0, c = 0.0;
    for (int j = 0; j < Ny; j++) {
        double v = hyv[j] - c, t = s + v;
        c = (t - s) - v;
        s = t;
        yn[j + 1] = s;
    }
    for (int j = 0; j < Ny; j++) yd[j] = 0.5 * (yn[j] + yn[j + 1]);
    dyv[0] = hyv[0];
    dyv[Ny] = hyv[Ny - 1];
    for (int j = 1; j < Ny; j++) dyv[j] = hyv[j - 1] == hyv[j] ? hyv[j] : 0.5 * (hyv[j - 1] + hyv[j]);
    P.Lx = xn[Nx];
    P.Lz = zn[Nz];
    D = hxv[0];
    for (int i = 1; i < Nx; i++) if (fabs(hxv[i] / hxv[0] - 1.0) > 1e-14) xnonuni = 1;
    for (int k = 1; k < Nz; k++) if (fabs(hzv[k] / hzv[0] - 1.0) > 1e-14) znonuni = 1;
    /* materials: eps_n per dual node (cell), eps_t per primal node (h-weighted average, derivation_nonuniform §1) */
    epsT = xcalloc(SY, sizeof(double)); epsN = xcalloc(SY, sizeof(double));
    for (int j = 0; j < Ny; j++) epsN[j] = (en && nn == Ny) ? en[j] : ey[j];
    epsN[Ny] = epsN[Ny - 1];
    if (et) {
        if (nt != Ny + 1) die("grid file: len(eps_t) != Ny+1");
        for (int j = 0; j <= Ny; j++) epsT[j] = et[j];
    } else {
        epsT[0] = epsN[0];
        epsT[Ny] = epsN[Ny - 1];
        for (int j = 1; j < Ny; j++) {
            double e1 = epsN[j - 1], e2 = epsN[j], h1 = hyv[j - 1], h2 = hyv[j];
            epsT[j] = (e1 == e2) ? e1 : (h1 == h2 ? 0.5 * (e1 + e2) : (h1 * e1 + h2 * e2) / (h1 + h2));
        }
    }
    /* structural checks (the sha256 is verified by grid_gen.check / fdtd_io before a run) */
    int u = 5;
    if (NPlo < 0 || NPhi < 0 || j0 - u < NPlo || j0 + u > Ny - NPhi) die("grid file: zones inconsistent");
    if (maxreldiff(hyv, 0, NPlo + u, hyv[0]) > 1e-12) die("grid file: near PML (+5 cells) not uniform");
    if (maxreldiff(hyv, Ny - NPhi - u, Ny, hyv[Ny - 1]) > 1e-12) die("grid file: far PML (+5 cells) not uniform");
    if (maxreldiff(hyv, j0 - u, j0 + u, hyv[j0]) > 1e-12) {
        fprintf(stderr, "ERROR: j0+-5 not uniform (max rel diff %.1e)\n", maxreldiff(hyv, j0 - u, j0 + u, hyv[j0]));
        exit(2);
    }
    {
        int lo = JA - u < 0 ? 0 : JA - u, hi = JA + u > Ny ? Ny : JA + u;
        double ref = JA - u < 0 ? hyv[0] : hyv[lo];
        if (hi > lo && maxreldiff(hyv, lo, hi, ref) > 1e-12) die("grid file: aux source j_a+-5 not uniform");
    }
    for (int j = 1; j < Ny; j++) if (hyv[j] <= 0.0) die("grid file: non-positive spacing");
    double *gj = jget("grating.j_lo", NULL, 0);
    if (gj) { /* lamellar grating: ridge cells i in [i_lo, i_hi) of the layer cells j in [j_lo, j_hi) */
        grating = 1;
        gr_jlo = (int)lround(*gj);
        gr_jhi = (int)lround(*jget("grating.j_hi", NULL, 1));
        gr_ilo = (int)lround(*jget("grating.i_lo", NULL, 1));
        gr_ihi = (int)lround(*jget("grating.i_hi", NULL, 1));
        gr_ridge = *jget("grating.eps_ridge", NULL, 1);
        gr_groove = *jget("grating.eps_groove", NULL, 1);
        if (gr_jlo < 1 || gr_jhi > Ny - 1 || gr_ilo < 0 || gr_ihi > Nx || gr_ilo >= gr_ihi) die("grid file: bad grating");
    }
}

/* eps of cell (i, j) = [x_i, x_i+1] x [y_j, y_j+1]; periodic in i */
static double eps_cell2(int i, int j) {
    i = (i % Nx + Nx) % Nx;
    if (grating && j >= gr_jlo && j < gr_jhi) return (i >= gr_ilo && i < gr_ihi) ? gr_ridge : gr_groove;
    return epsN[j];
}

/* stage B materials: tangential components average the adjacent cells with length weights
   (Ex over y, Ey over x, Ez over the four cells); layered media reproduce epsT/epsN exactly */
static void build_materials(void) {
    size_t n = (size_t)Nx * SY;
    eTx = xcalloc(n, sizeof(double)); eTz = xcalloc(n, sizeof(double)); eNy = xcalloc(n, sizeof(double));
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++) {
            size_t q = (size_t)i * SY + j;
            if (!grating || j < gr_jlo - 1 || j > gr_jhi) {
                eTx[q] = epsT[j];
                eTz[q] = epsT[j];
                eNy[q] = epsN[j];
                continue;
            }
            int jm = j > 0 ? j - 1 : 0, jp = j < Ny ? j : Ny - 1;
            double wm = (j > 0) ? hyv[jm] : 0.0, wp = (j < Ny) ? hyv[jp] : 0.0;
            eTx[q] = (wm * eps_cell2(i, jm) + wp * eps_cell2(i, jp)) / (wm + wp);
            double a = hxv[(i + Nx - 1) % Nx], b = hxv[i];
            eTz[q] = (wm * (a * eps_cell2(i - 1, jm) + b * eps_cell2(i, jm)) + wp * (a * eps_cell2(i - 1, jp) + b * eps_cell2(i, jp)))
                     / ((wm + wp) * (a + b));
            eNy[q] = (j < Ny) ? (a * eps_cell2(i - 1, j) + b * eps_cell2(i, j)) / (a + b) : epsN[j];
        }
}

/* ------------------------------------------------------------------ setup */
static void setup(void) {
    W0 = 2.0 * PI;
    if (P.mesh == 'u') build_uniform_mesh();
    else build_file_mesh();
    if ((xnonuni || znonuni) && (P.inc == 'a' || P.auxref)) {
        fprintf(stderr, "ERROR: phasor aux line requires uniform x and z (x nonuniform=%d, z nonuniform=%d); "
                        "use inc=p|j|m\n", xnonuni, znonuni);
        exit(2);
    }
    if (P.S >= 1.0 / sqrt(3.0)) die("Courant number S must be < 1/sqrt(3)");
    /* time step (SPEC_nonuniform §7.5): exactly S*Delta when the three minimum spacings coincide */
    {
        double mx = hxv[0], my = hyv[0], mz = hzv[0];
        for (int i = 0; i < Nx; i++) mx = fmin(mx, hxv[i]);
        for (int j = 0; j < Ny; j++) my = fmin(my, hyv[j]);
        for (int k = 0; k < Nz; k++) mz = fmin(mz, hzv[k]);
        if (P.mesh == 'u') dt_cour = P.S * D;
        else if (mx == my && my == mz) dt_cour = P.S * mx;
        else dt_cour = P.S * sqrt(3.0) / sqrt(1.0 / (mx * mx) + 1.0 / (my * my) + 1.0 / (mz * mz));
        dt = (P.dt_set > 0.0) ? P.dt_set : dt_cour * P.dtfac;
        if (dt > dt_cour / P.S / sqrt(3.0) * (1.0 + 1e-12))
            fprintf(stderr, "WARNING: dt exceeds the minimum-spacing Courant bound by %.2f%%\n",
                    100.0 * (dt / (dt_cour / P.S / sqrt(3.0)) - 1.0));
    }
    /* source-free runs (energy / stability / eigen tests) need no plane wave: kx = kz = 0, no amplitudes */
    int wave = !(P.mesh == 'f' && P.inc == '0' && P.init != 'a');   /* legacy mesh=uniform runs unchanged */
    kx = wave ? 2.0 * PI * P.m / P.Lx : 0.0;
    kz = wave ? 2.0 * PI * P.n / P.Lz : 0.0;
    double kt = hypot(kx, kz);
    if (wave && kt >= 0.95 * W0) die("|k_t| >= 0.95 k0: insufficient propagation margin");
    wt = 2.0 / dt * sin(W0 * dt / 2.0);
    Da = gh(JA);
    ky = ky_local(Da, 1.0);
    kyc = sqrt(W0 * W0 - kt * kt);
    if (wave) amplitudes(ky, Kt, E0v, H0v);
    for (int c = 0; c < 3; c++) Kt_src[c] = Kt[c];
    if (P.inc == 'p' || P.inc == 'j') { /* stage B: textbook continuous wave (continuous k, omega) */
        double kc[3] = {kx, kyc, kz}, yh[3] = {0, 1, 0}, s[3], pp[3];
        cross(kc, yh, s);
        double ns = norm3(s);
        for (int c = 0; c < 3; c++) s[c] /= ns;
        cross(s, kc, pp);
        double np = norm3(pp);
        for (int c = 0; c < 3; c++) pp[c] /= np;
        for (int c = 0; c < 3; c++) E0v[c] = (P.pol == 's') ? s[c] : pp[c];
        cross(kc, E0v, H0v);
        for (int c = 0; c < 3; c++) { H0v[c] /= W0; Kt_src[c] = kc[c]; }
        ky = kyc;
    }
    if (P.inc == 'n' && P.ky_cont) { /* control C1: continuous ky everywhere in the source */
        amplitudes(kyc, Kt_src, E0v, H0v);
        ky = kyc;
    }
    if (P.pol != 's' && P.pol != 'p') die("pol must be s or p");
    if (P.mesh == 'u' && P.tf < 10 * P.nl) fprintf(stderr, "WARNING: TF length %d cells < 10 lambda0\n", P.tf);
    if (P.inc != '0' && P.init == 'a') die("init=analytic is a debug mode without a source");
    if (P.inc != '0' && (j0 - NPlo < 2 || Ny - NPhi - j0 < 2)) die("sf and tf must be >= 2 cells");
    if (j1 >= Ny - NPhi) die("interface must lie in the physical TF region");
    if (P.auxref && P.inc != 'a') die("auxref=1 needs inc=a");
    if (P.inc == 'j' && P.jsrc < 0) P.jsrc = j0;
    if (P.inc == 'j' && (P.jsrc <= NPlo || P.jsrc >= Ny - NPhi)) die("jsrc must lie outside the PML");
    SZ = Nz + 2;
    NTOT = (size_t)(Nx + 2) * SY * SZ;
    Ex = xcalloc(NTOT, sizeof(double)); Ey = xcalloc(NTOT, sizeof(double)); Ez = xcalloc(NTOT, sizeof(double));
    Hx = xcalloc(NTOT, sizeof(double)); Hy = xcalloc(NTOT, sizeof(double)); Hz = xcalloc(NTOT, sizeof(double));

    /* coefficients (derivation_nonuniform §1): y coefficient factored out, the other term scaled by an exact
       spacing ratio (== 1.0 on a uniform grid, so the legacy arithmetic is reproduced) */
    cHy = xcalloc(SY, sizeof(double)); ay = xcalloc(SY, sizeof(double)); cEyn = xcalloc(SY, sizeof(double));
    cHz = xcalloc(Nz, sizeof(double)); zfac = xcalloc(Nz, sizeof(double));
    for (int j = 0; j < Ny; j++) cHy[j] = dt / hyv[j];
    for (int j = 0; j <= Ny; j++) {
        ay[j] = dt / (epsT[j] * dyv[j]);
        cEyn[j] = dt / (epsN[j] * dzv[0]);
    }
    for (int k = 0; k < Nz; k++) { cHz[k] = dt / hzv[k]; zfac[k] = dzv[0] / dzv[k]; }
    build_materials();
    ayx = xcalloc((size_t)Nx * SY, sizeof(double)); ayz = xcalloc((size_t)Nx * SY, sizeof(double));
    cEy2 = xcalloc((size_t)Nx * SY, sizeof(double));
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++) {
            size_t q = (size_t)i * SY + j;
            ayx[q] = (eTx[q] == epsT[j]) ? ay[j] : dt / (eTx[q] * dyv[j]);
            ayz[q] = (eTz[q] == epsT[j]) ? ay[j] : dt / (eTz[q] * dyv[j]);
            cEy2[q] = (eNy[q] == epsN[j]) ? cEyn[j] : dt / (eNy[q] * dzv[0]);
        }
    if (grating && P.auxref) die("auxref=1 needs a layered (x-independent) structure");
    rExz = xcalloc((size_t)SY * Nz, sizeof(double)); rHxz = xcalloc((size_t)SY * Nz, sizeof(double));
    rEzx = xcalloc((size_t)Nx * SY, sizeof(double)); rHzx = xcalloc((size_t)Nx * SY, sizeof(double));
    rEyk = xcalloc((size_t)Nx * Nz, sizeof(double)); rHyk = xcalloc((size_t)Nx * Nz, sizeof(double));
    for (int j = 0; j <= Ny; j++)
        for (int k = 0; k < Nz; k++) {
            rExz[(size_t)j * Nz + k] = dyv[j] / dzv[k];
            rHxz[(size_t)j * Nz + k] = (j < Ny) ? hyv[j] / hzv[k] : 0.0;
        }
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            rEzx[(size_t)i * SY + j] = dyv[j] / dxv[i];
            rHzx[(size_t)i * SY + j] = (j < Ny) ? hyv[j] / hxv[i] : 0.0;
        }
        for (int k = 0; k < Nz; k++) {
            rEyk[(size_t)i * Nz + k] = dzv[k] / dxv[i];
            rHyk[(size_t)i * Nz + k] = hzv[k] / hxv[i];
        }
    }

    /* CPML (derivation.md §7): profiles by physical distance from the PML inner face (SPEC_nonuniform §7.4) */
    pmlE = xcalloc(SY, sizeof(int)); pmlH = xcalloc(SY, sizeof(int));
    bE = xcalloc(SY, sizeof(double)); cE = xcalloc(SY, sizeof(double)); ikE = xcalloc(SY, sizeof(double));
    bH = xcalloc(SY, sizeof(double)); cH = xcalloc(SY, sizeof(double)); ikH = xcalloc(SY, sizeof(double));
    sigE = xcalloc(SY, sizeof(double)); sigH = xcalloc(SY, sizeof(double));
    kapE = xcalloc(SY, sizeof(double)); kapH = xcalloc(SY, sizeof(double));
    alpE = xcalloc(SY, sizeof(double)); alpH = xcalloc(SY, sizeof(double));
    NP = (NPlo + 1) + (NPhi + 1);
    double ylo = yn[NPlo], yhi = yn[Ny - NPhi];
    double dlo = (P.mesh == 'u') ? P.npml * D : ylo - yn[0];
    double dhi = (P.mesh == 'u') ? P.npml * D : yn[Ny] - yhi;
    double nlo = sqrt(epsN[0]), nhi = sqrt(epsN[Ny - 1]);
    double Dlo = hyv[0], Dhi = hyv[Ny - 1];
    if (P.mesh == 'u') { nlo = 1.0; nhi = (j1 >= 0) ? sqrt(P.eps2) : 1.0; }
    for (int j = 0; j <= Ny; j++) {
        for (int h = 0; h < 2; h++) {
            double y = (h == 0) ? yn[j] : (j < Ny ? yd[j] : yn[Ny]);
            double rho = 0.0, nloc = 1.0, dd = 1.0, Dl = D;
            int slot = -1;
            if (NPlo > 0 && y < ylo) { rho = ylo - y; slot = j; nloc = nlo; dd = dlo; Dl = Dlo; }
            else if (NPhi > 0 && y > yhi) { rho = y - yhi; slot = j - (Ny - NPhi) + NPlo + 1; nloc = nhi; dd = dhi; Dl = Dhi; }
            if (h == 1 && j == Ny) slot = -1;
            double r = (slot >= 0) ? rho / dd : 0.0;
            double smax = P.sig_fac * (P.pml_m + 1.0) / (Dl * nloc); /* eta0 = 1, eta = 1/n */
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
    size_t nq = (size_t)Nx * Nz;
    cosP1 = xcalloc(nq, sizeof(double)); sinP1 = xcalloc(nq, sizeof(double));
    cosP2 = xcalloc(nq, sizeof(double)); sinP2 = xcalloc(nq, sizeof(double));
    cosP3 = xcalloc(nq, sizeof(double)); sinP3 = xcalloc(nq, sizeof(double));
    cosP4 = xcalloc(nq, sizeof(double)); sinP4 = xcalloc(nq, sizeof(double));
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            double p1, p2;
            if (P.mesh == 'u') { p1 = kx * (i + 0.5) * D + kz * k * D; p2 = kx * i * D + kz * (k + 0.5) * D; }
            else { p1 = kx * xd[i] + kz * zn[k]; p2 = kx * xn[i] + kz * zd[k]; }
            double p3 = kx * xn[i] + kz * zn[k], p4 = kx * xd[i] + kz * zd[k];
            size_t q = (size_t)i * Nz + k;
            cosP1[q] = cos(p1); sinP1[q] = sin(p1);
            cosP2[q] = cos(p2); sinP2[q] = sin(p2);
            cosP3[q] = cos(p3); sinP3[q] = sin(p3);
            cosP4[q] = cos(p4); sinP4[q] = sin(p4);
        }
    ExI = xcalloc(nq, sizeof(double)); EzI = xcalloc(nq, sizeof(double));
    HxI = xcalloc(nq, sizeof(double)); HzI = xcalloc(nq, sizeof(double));
    if (P.inc == 'a') aux_init();
    if (P.inc == 'm') { aux_init(); modal_init(); } /* shared y geometry of the aux line, one line per product mode */
    if (P.auxref) auxref_init();
}

/* ------------------------------------------------------------------ analytic fill (debug, Stage 1) */
static double comp_x(int c, int i) { return OFF[c][0] == 0.0 ? xn[i] : xd[i]; }
static double comp_y(int c, int j) { return OFF[c][1] == 0.0 ? yn[j] : yd[j]; }
static double comp_z(int c, int k) { return OFF[c][2] == 0.0 ? zn[k] : zd[k]; }

static double analytic_val(int c, int i, int j, int k, double tn) {
    double x, y, z;
    if (P.mesh == 'u') { x = (i + OFF[c][0]) * D; y = (j + OFF[c][1]) * D; z = (k + OFF[c][2]) * D; }
    else { x = comp_x(c, i); y = comp_y(c, j); z = comp_z(c, k); }
    double t = (tn + OFF[c][3]) * dt;
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

/* ------------------------------------------------------------------ updates (derivation_nonuniform §1)
   If doW, update_H also accumulates the conserved energy of derivation_nonuniform §2 at t = n dt:
   W_mod = 1/2 sum( eps W_E |E^n|^2 + mu W_H H^{n+1/2}.H^{n-1/2} ). */
static double Wmod_acc;
static void update_H(int doW) {
    const size_t sI = (size_t)SY * SZ, sJ = SZ;
    ghost_fwd(Ex); ghost_fwd(Ey); ghost_fwd(Ez);
    double wsum = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : wsum)
    for (int i = 0; i < Nx; i++) {
        const double *rHy_i = rHyk + (size_t)i * Nz, *rHz_i = rHzx + (size_t)i * SY;
        for (int j = 0; j <= Ny; j++) {
            size_t b = ID(i, j, 0);
            double *hy = Hy + b;
            const double *ex = Ex + b, *ez = Ez + b, *ey = Ey + b;
            double wHy = 0.0, wEx = 0.0, wEz = 0.0, wEy = 0.0, wHx = 0.0, wHz = 0.0;
            for (int k = 0; k < Nz; k++) {
                double o = hy[k];
                hy[k] -= cHz[k] * ((ex[k + 1] - ex[k]) - rHy_i[k] * (ez[k + sI] - ez[k]));
                if (doW) {
                    wHy += hzv[k] * o * hy[k];
                    wEx += dzv[k] * ex[k] * ex[k];
                    wEz += hzv[k] * ez[k] * ez[k];
                }
            }
            if (doW) wsum += hxv[i] * dyv[j] * wHy + dyv[j] * (eTx[(size_t)i * SY + j] * hxv[i] * wEx + eTz[(size_t)i * SY + j] * dxv[i] * wEz);
            if (j == Ny) continue;
            double *hx = Hx + b, *hz = Hz + b;
            const double *rHx_j = rHxz + (size_t)j * Nz;
            double chy = cHy[j], rz = rHz_i[j];
            int s = pmlH[j];
            if (s < 0) {
                for (int k = 0; k < Nz; k++) {
                    double ox = hx[k], oz = hz[k];
                    hx[k] -= chy * ((ez[k + sJ] - ez[k]) - rHx_j[k] * (ey[k + 1] - ey[k]));
                    hz[k] -= chy * (rz * (ey[k + sI] - ey[k]) - (ex[k + sJ] - ex[k]));
                    if (doW) {
                        wHx += hzv[k] * ox * hx[k];
                        wHz += dzv[k] * oz * hz[k];
                        wEy += dzv[k] * ey[k] * ey[k];
                    }
                }
            } else {
                double bb = bH[j], cc = cH[j], ik = ikH[j];
                double *px = psiHxy + ((size_t)i * NP + s) * Nz, *pz = psiHzy + ((size_t)i * NP + s) * Nz;
                for (int k = 0; k < Nz; k++) {
                    double ox = hx[k], oz = hz[k];
                    double dEz = ez[k + sJ] - ez[k], dEx = ex[k + sJ] - ex[k];
                    px[k] = bb * px[k] + cc * dEz;
                    pz[k] = bb * pz[k] + cc * dEx;
                    hx[k] -= chy * ((ik * dEz + px[k]) - rHx_j[k] * (ey[k + 1] - ey[k]));
                    hz[k] -= chy * (rz * (ey[k + sI] - ey[k]) - (ik * dEx + pz[k]));
                    if (doW) {
                        wHx += hzv[k] * ox * hx[k];
                        wHz += dzv[k] * oz * hz[k];
                        wEy += dzv[k] * ey[k] * ey[k];
                    }
                }
            }
            if (doW) wsum += hyv[j] * (dxv[i] * wHx + hxv[i] * wHz + eNy[(size_t)i * SY + j] * dxv[i] * wEy);
        }
    }
    if (doW) Wmod_acc = 0.5 * wsum;
}

/* E update. If doW, also accumulates the Yee-conserved energy
   W^{n+1/2} = 1/2 sum( eps W_E E^n . E^{n+1} + mu W_H H^{n+1/2} . H^{n+1/2} )
   (exactly constant in a lossless PEC/PBC cavity), per i into partT/partP (total / outside PML). */
static double *partT, *partP;
static void update_E(int doW) {
    const size_t sI = (size_t)SY * SZ, sJ = SZ;
    ghost_bwd(Hx); ghost_bwd(Hy); ghost_bwd(Hz);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < Nx; i++) {
        double wT = 0.0, wP = 0.0;
        const double *rEy_i = rEyk + (size_t)i * Nz, *rEz_i = rEzx + (size_t)i * SY;
        for (int j = 0; j <= Ny; j++) {
            size_t b = ID(i, j, 0);
            const double *hx = Hx + b, *hy = Hy + b, *hz = Hz + b;
            double wI = 0.0, wH = 0.0; /* integer-y nodes, half-y nodes */
            if (j < Ny) {
                double *ey = Ey + b, c = cEy2[(size_t)i * SY + j];
                double sx = 0.0, sz = 0.0, sE = 0.0;
                for (int k = 0; k < Nz; k++) {
                    double old = ey[k];
                    ey[k] += c * zfac[k] * ((hx[k] - hx[k - 1]) - rEy_i[k] * (hz[k] - hz[k - sI]));
                    if (doW) {
                        sE += dzv[k] * old * ey[k];
                        sx += hzv[k] * hx[k] * hx[k];
                        sz += dzv[k] * hz[k] * hz[k];
                    }
                }
                if (doW) wH = hyv[j] * (eNy[(size_t)i * SY + j] * dxv[i] * sE + dxv[i] * sx + hxv[i] * sz);
            }
            if (doW) {
                double sy = 0.0;
                for (int k = 0; k < Nz; k++) sy += hzv[k] * hy[k] * hy[k];
                wI += hxv[i] * dyv[j] * sy;
            }
            if (j > 0 && j < Ny) { /* j = 0, Ny: PEC, tangential E fixed (Dirichlet in debug mode) */
                double *ex = Ex + b, *ez = Ez + b, a = ayx[(size_t)i * SY + j], a2 = ayz[(size_t)i * SY + j], rzx = rEz_i[j];
                const double *rxz = rExz + (size_t)j * Nz;
                double sx = 0.0, sz = 0.0;
                int s = pmlE[j];
                if (s < 0) {
                    for (int k = 0; k < Nz; k++) {
                        double ox = ex[k], oz = ez[k];
                        ex[k] += a * ((hz[k] - hz[k - sJ]) - rxz[k] * (hy[k] - hy[k - 1]));
                        ez[k] += a2 * (rzx * (hy[k] - hy[k - sI]) - (hx[k] - hx[k - sJ]));
                        if (doW) { sx += dzv[k] * ox * ex[k]; sz += hzv[k] * oz * ez[k]; }
                    }
                } else {
                    double bb = bE[j], cc = cE[j], ik = ikE[j];
                    double *px = psiExy + ((size_t)i * NP + s) * Nz, *pz = psiEzy + ((size_t)i * NP + s) * Nz;
                    for (int k = 0; k < Nz; k++) {
                        double ox = ex[k], oz = ez[k];
                        double dHz = hz[k] - hz[k - sJ], dHx = hx[k] - hx[k - sJ];
                        px[k] = bb * px[k] + cc * dHz;
                        pz[k] = bb * pz[k] + cc * dHx;
                        ex[k] += a * ((ik * dHz + px[k]) - rxz[k] * (hy[k] - hy[k - 1]));
                        ez[k] += a2 * (rzx * (hy[k] - hy[k - sI]) - (ik * dHx + pz[k]));
                        if (doW) { sx += dzv[k] * ox * ex[k]; sz += hzv[k] * oz * ez[k]; }
                    }
                }
                if (doW) wI += dyv[j] * (eTx[(size_t)i * SY + j] * hxv[i] * sx + eTz[(size_t)i * SY + j] * dxv[i] * sz);
            } else if (doW && P.init == 'a') {
                double *ex = Ex + b, *ez = Ez + b;
                double sx = 0.0, sz = 0.0;
                for (int k = 0; k < Nz; k++) { sx += dzv[k] * ex[k] * ex[k]; sz += hzv[k] * ez[k] * ez[k]; }
                wI += dyv[j] * (hxv[i] * sx + dxv[i] * sz);
            }
            wT += wI + wH;
            if (j >= NPlo && j <= Ny - NPhi) wP += wI;
            if (j >= NPlo && j <= Ny - NPhi - 1) wP += wH;
        }
        if (doW) { partT[i] = wT; partP[i] = wP; }
    }
}

static void energy_sum(double *Wt, double *Wp) {
    double a = 0.0, b = 0.0;
    for (int i = 0; i < Nx; i++) { a += partT[i]; b += partP[i]; }
    *Wt = 0.5 * a;
    *Wp = 0.5 * b;
}

/* incident field on the TF/SF plane from complex amplitudes (derivation.md §6.2) */
static void incident_E(long n) { /* Ex, Ez at j0, time n */
    double complex ax, az;
    if (P.inc == 'm') { /* sum of the product modes: Re[aux_t(j0) * w_t * profile_t] */
        for (int q = 0; q < Nx * Nz; q++) { ExI[q] = 0.0; EzI[q] = 0.0; }
        for (int t = 0; t < NMA; t++) {
            const MAux *A = &MA[t];
            ax = A->Ex[j0 - Ajlo]; az = A->Ez[j0 - Ajlo];
            for (int q = 0; q < Nx * Nz; q++) { ExI[q] += creal(ax * A->T1[q]); EzI[q] += creal(az * A->T2[q]); }
        }
        return;
    }
    if (P.inc == 'a') { ax = aEx[j0 - Ajlo]; az = aEz[j0 - Ajlo]; }
    else { ax = inc_amp(0, j0, (double)n, E0v, H0v, ky); az = inc_amp(2, j0, (double)n, E0v, H0v, ky); }
    incEx = ax; incEz = az;
    for (int q = 0; q < Nx * Nz; q++) {
        ExI[q] = creal(ax) * cosP1[q] - cimag(ax) * sinP1[q];
        EzI[q] = creal(az) * cosP2[q] - cimag(az) * sinP2[q];
    }
}

static void incident_H(long n) { /* Hx, Hz at j0 - 1/2, time n + 1/2 */
    double complex ax, az;
    if (P.inc == 'm') {
        for (int q = 0; q < Nx * Nz; q++) { HxI[q] = 0.0; HzI[q] = 0.0; }
        for (int t = 0; t < NMA; t++) {
            const MAux *A = &MA[t];
            ax = A->Hx[j0 - 1 - Ajlo]; az = A->Hz[j0 - 1 - Ajlo];
            for (int q = 0; q < Nx * Nz; q++) { HxI[q] += creal(ax * A->T2[q]); HzI[q] += creal(az * A->T1[q]); }
        }
        return;
    }
    if (P.inc == 'a') { ax = aHx[j0 - 1 - Ajlo]; az = aHz[j0 - 1 - Ajlo]; }
    else { ax = inc_amp(3, j0 - 1, (double)n, E0v, H0v, ky); az = inc_amp(5, j0 - 1, (double)n, E0v, H0v, ky); }
    incHx = ax; incHz = az;
    for (int q = 0; q < Nx * Nz; q++) {
        HxI[q] = creal(ax) * cosP2[q] - cimag(ax) * sinP2[q];
        HzI[q] = creal(az) * cosP1[q] - cimag(az) * sinP1[q];
    }
}

/* inc=j: electric surface current J_s = Re[J0 e^{i(kx x + kz z)} g(t) e^{-i w t}] (J0 = tangential E0) on the
   primal plane jsrc, applied as a volume current J_s / d_j in the E update at t = (n + 1/2) dt */
static void source_J(long n) {
    double t = (n + 0.5) * dt;
    double complex g = envelope(t) * cexp(-I * W0 * t);
    double complex jx = g * E0v[0], jz = g * E0v[2];
    int j = P.jsrc;
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            size_t q = (size_t)i * Nz + k, m = (size_t)i * SY + j;
            Ex[ID(i, j, k)] -= ayx[m] * (creal(jx) * cosP1[q] - cimag(jx) * sinP1[q]);
            Ez[ID(i, j, k)] -= ayz[m] * (creal(jz) * cosP2[q] - cimag(jz) * sinP2[q]);
        }
}

/* TF/SF corrections, derivation.md §3 (a)-(d); coefficients at the local spacing (derivation_nonuniform §4) */
static void tfsf_H(void) {
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            Hx[ID(i, j0 - 1, k)] += cHy[j0 - 1] * EzI[i * Nz + k];
            Hz[ID(i, j0 - 1, k)] -= cHy[j0 - 1] * ExI[i * Nz + k];
        }
}

static void tfsf_E(void) {
    for (int i = 0; i < Nx; i++)
        for (int k = 0; k < Nz; k++) {
            Ex[ID(i, j0, k)] -= ayx[(size_t)i * SY + j0] * HzI[i * Nz + k];
            Ez[ID(i, j0, k)] += ayz[(size_t)i * SY + j0] * HxI[i * Nz + k];
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

static double complex *auxdft[6], *refdft[6], *projdft[6];
static double *vre[6], *vim[6];
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
    if (P.auxref)
        for (int c = 0; c < 6; c++) refdft[c] = xcalloc(SY, sizeof(double complex));
    if (P.proj)
        for (int c = 0; c < 6; c++) projdft[c] = xcalloc(SY, sizeof(double complex));
    if (P.planar)
        for (int c = 0; c < 6; c++) {
            vre[c] = xcalloc((size_t)Nx * SY * Nz, sizeof(double));
            vim[c] = xcalloc((size_t)Nx * SY * Nz, sizeof(double));
        }
}

static void dft_accumulate(int first, int last, double tsec) {
    double cr = cos(W0 * tsec), ci = sin(W0 * tsec);
    for (int s = 0; s < NSL; s++) {
        Slice *sl = &SL[s];
        for (int c = first; c <= last; c++) {
            double *F = comp_ptr(c);
            int jmax = (OFF[c][1] == 0.0) ? Ny : Ny - 1; /* half-y components have no sample at j = Ny */
            size_t q = 0;
            if (sl->type == 0) {
                for (int i = 0; i < Nx; i++)
                    for (int k = 0; k < Nz; k++, q++) { double f = F[ID(i, sl->idx, k)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci; }
            } else if (sl->type == 1) {
                for (int i = 0; i < Nx; i++)
                    for (int j = 0; j <= Ny; j++, q++) {
                        if (j > jmax) continue;
                        double f = F[ID(i, j, sl->idx)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci;
                    }
            } else {
                for (int j = 0; j <= Ny; j++)
                    for (int k = 0; k < Nz; k++, q++) {
                        if (j > jmax) continue;
                        double f = F[ID(sl->idx, j, k)]; sl->re[c][q] += f * cr; sl->im[c][q] += f * ci;
                    }
            }
        }
    }
    double complex e = cexp(I * W0 * tsec);
    if (P.inc == 'a') {
        double complex *arr[6] = {aEx, aEy, aEz, aHx, aHy, aHz};
        for (int c = first; c <= last; c++)
            for (long u = aux_dft_lo; u <= aux_dft_hi; u++) auxdft[c][u - aux_dft_lo] += arr[c][u] * e;
    }
    if (P.auxref) {
        double complex *arr[6] = {rEx, rEy, rEz, rHx, rHy, rHz};
        for (int c = first; c <= last; c++)
            for (int j = 0; j <= Ny; j++) refdft[c][j] += arr[c][j] * e;
    }
    if (P.proj) { /* reduced phasor a_c(j): F = Re[a e^{i(phi_c - w t)}]; x-z Floquet projection with the quadrature
                     weights of the component's sample positions (primal: dual spacing, dual: primal spacing):
                     a = 2/(Lx Lz Nt) sum_t sum_xz w F e^{-i phi_c} e^{iwt} (uniform grid: w = Delta^2) */
        const double *cp[6] = {cosP1, cosP3, cosP2, cosP2, cosP4, cosP1}, *sp[6] = {sinP1, sinP3, sinP2, sinP2, sinP4, sinP1};
        for (int c = first; c <= last; c++) {
            double *F = comp_ptr(c);
            const double *wx = (OFF[c][0] == 0.0) ? dxv : hxv, *wz = (OFF[c][2] == 0.0) ? dzv : hzv;
#pragma omp parallel for schedule(static)
            for (int j = 0; j <= Ny; j++) {
                double sr = 0.0, si = 0.0;
                for (int i = 0; i < Nx; i++)
                    for (int k = 0; k < Nz; k++) {
                        size_t q = (size_t)i * Nz + k;
                        double f = F[ID(i, j, k)] * wx[i] * wz[k];
                        sr += f * cp[c][q];
                        si -= f * sp[c][q];
                    }
                projdft[c][j] += (sr + I * si) * e;
            }
        }
    }
    if (P.planar)
        for (int c = first; c <= last; c++) {
            double *F = comp_ptr(c), *R = vre[c], *M = vim[c];
#pragma omp parallel for schedule(static)
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j <= Ny; j++)
                    for (int k = 0; k < Nz; k++) {
                        size_t q = ((size_t)i * SY + j) * Nz + k;
                        double f = F[ID(i, j, k)];
                        R[q] += f * cr;
                        M[q] += f * ci;
                    }
        }
}

/* ------------------------------------------------------------------ diagnostics */
static void sf_max(double *mE, double *mH) {
    double me = 0.0, mh = 0.0;
    int jlo = NPlo + 1;           /* integer-y nodes strictly outside the near PML */
    for (int i = 0; i < Nx; i++)
        for (int j = NPlo; j < j0; j++)
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

/* max |F_main - Re[aux_ref e^{i phi}]| over y-index range [jlo, jhi] (E and H separately) */
static void ref_dev(int jlo, int jhi, double *dE, double *dH) {
    double me = 0.0, mh = 0.0;
    const double *cp[6] = {cosP1, cosP3, cosP2, cosP2, cosP4, cosP1}, *sp[6] = {sinP1, sinP3, sinP2, sinP2, sinP4, sinP1};
    double complex *ar[6] = {rEx, rEy, rEz, rHx, rHy, rHz};
    for (int c = 0; c < 6; c++) {
        double *F = comp_ptr(c);
        int top = (OFF[c][1] == 0.0) ? jhi : (jhi < Ny ? jhi : Ny - 1);
        double m = 0.0; /* serial: small region every step; a parallel region per component costs more */
        for (int i = 0; i < Nx; i++)
            for (int j = jlo; j <= top; j++) {
                double ar_ = creal(ar[c][j]), ai_ = cimag(ar[c][j]);
                for (int k = 0; k < Nz; k++) {
                    size_t q = (size_t)i * Nz + k;
                    double v = fabs(F[ID(i, j, k)] - (ar_ * cp[c][q] - ai_ * sp[c][q]));
                    if (v > m) m = v;
                }
            }
        if (c < 3) me = fmax(me, m); else mh = fmax(mh, m);
    }
    *dE = me;
    *dH = mh;
}

/* max |div E| at nodes (i,j,k) and max |div H| at cell centres over the TF interior
   j in [j0+2, Ny-npml_hi-2] (excludes the TF/SF plane +-1 cell and the PML); explicit PBC wrap.
   divop=n: nonuniform operator (derivation_nonuniform §2); divop=u: uniform operator with Delta = hx[0] (control). */
static void div_max(double *dE, double *dH) {
    double me = 0.0, mh = 0.0;
    int jlo = j0 + 2, jhi = Ny - NPhi - 2;
    if (j1 >= 0 && j1 - 1 < jhi) jhi = j1 - 1; /* stay on the vacuum side of an interface */
    for (int j = j0 + 1; j < Ny; j++) if (epsN[j] != epsN[j - 1]) { if (j - 1 < jhi) jhi = j - 1; break; }
    int nu = (P.divop == 'n');
#pragma omp parallel for schedule(static) reduction(max : me, mh)
    for (int i = 0; i < Nx; i++) {
        int im = (i + Nx - 1) % Nx, ip = (i + 1) % Nx;
        for (int j = jlo; j <= jhi; j++)
            for (int k = 0; k < Nz; k++) {
                int km = (k + Nz - 1) % Nz, kp = (k + 1) % Nz;
                double a1 = Ex[ID(i, j, k)] - Ex[ID(im, j, k)], a2 = Ey[ID(i, j, k)] - Ey[ID(i, j - 1, k)],
                       a3 = Ez[ID(i, j, k)] - Ez[ID(i, j, km)];
                double b1 = Hx[ID(ip, j, k)] - Hx[ID(i, j, k)], b2 = Hy[ID(i, j + 1, k)] - Hy[ID(i, j, k)],
                       b3 = Hz[ID(i, j, kp)] - Hz[ID(i, j, k)];
                double de, dh;
                if (nu) {
                    de = a1 / dxv[i] + a2 / dyv[j] + a3 / dzv[k];
                    dh = b1 / hxv[i] + b2 / hyv[j] + b3 / hzv[k];
                } else {
                    de = (a1 + a2 + a3) / D;
                    dh = (b1 + b2 + b3) / D;
                }
                if (fabs(de) > me) me = fabs(de);
                if (fabs(dh) > mh) mh = fabs(dh);
            }
    }
    *dE = me;
    *dH = mh;
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
        for (int c = 0; c < 6; c++) {
            int half = (OFF[c][1] != 0.0);
            for (size_t q = 0; q < SL[s].cnt; q++) {
                double v[2] = {SL[s].re[c][q] * sc, SL[s].im[c][q] * sc};
                if (half && SL[s].type == 1 && (q % SY) == (size_t)Ny) v[0] = v[1] = NAN;
                if (half && SL[s].type == 2 && (q / Nz) == (size_t)Ny) v[0] = v[1] = NAN;
                fwrite(v, sizeof(double), 2, f);
            }
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
    if (P.proj) {
        FILE *f = xfopen("dft_proj.bin", "wb");
        double pc = 2.0 / (P.Lx * P.Lz * (double)N);
        for (int c = 0; c < 6; c++)
            for (int j = 0; j <= Ny; j++) {
                double complex z = projdft[c][j] * pc;
                double v[2] = {creal(z), cimag(z)};
                if (OFF[c][1] != 0.0 && j == Ny) v[0] = v[1] = NAN;
                fwrite(v, sizeof(double), 2, f);
            }
        fclose(f);
    }
    if (P.planar) { /* RMS over (x, z) of arg(F_dft e^{-i phi_c} / a_c(j)) for every component and y plane */
        FILE *f = xfopen("dft_planar.bin", "wb");
        const double *cp[6] = {cosP1, cosP3, cosP2, cosP2, cosP4, cosP1}, *sp[6] = {sinP1, sinP3, sinP2, sinP2, sinP4, sinP1};
        for (int c = 0; c < 6; c++)
            for (int j = 0; j <= Ny; j++) {
                double complex a = projdft[c][j];
                double s2 = 0.0;
                for (int i = 0; i < Nx; i++)
                    for (int k = 0; k < Nz; k++) {
                        size_t q = ((size_t)i * SY + j) * Nz + k, t = (size_t)i * Nz + k;
                        double complex z = (vre[c][q] + I * vim[c][q]) * (cp[c][t] - I * sp[c][t]);
                        double ph = carg(z * conj(a));
                        s2 += ph * ph;
                    }
                double v = (cabs(a) > 0.0 && !(OFF[c][1] != 0.0 && j == Ny)) ? sqrt(s2 / ((double)Nx * Nz)) : NAN;
                fwrite(&v, sizeof(double), 1, f);
            }
        fclose(f);
    }
    if (P.auxref) {
        FILE *f = xfopen("dft_auxref.bin", "wb");
        for (int c = 0; c < 6; c++)
            for (int j = 0; j <= Ny; j++) {
                double complex z = refdft[c][j] / (double)N;
                double v[2] = {creal(z), cimag(z)};
                if (OFF[c][1] != 0.0 && j == Ny) v[0] = v[1] = NAN;
                fwrite(v, sizeof(double), 2, f);
            }
        fclose(f);
    }
}

static void dump_fields(const char *name) {
    FILE *f = xfopen(name, "wb");
    for (int c = 0; c < 6; c++) {
        double *F = comp_ptr(c);
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j <= Ny; j++) fwrite(&F[ID(i, j, 0)], sizeof(double), Nz, f);
    }
    fclose(f);
}

static void dump_auxref(long step) {
    char fn[64];
    snprintf(fn, sizeof fn, "auxref_n%ld.bin", step);
    FILE *f = xfopen(fn, "wb");
    double complex *arr[6] = {rEx, rEy, rEz, rHx, rHy, rHz};
    for (int c = 0; c < 6; c++)
        for (int j = 0; j <= Ny; j++) {
            double v[2] = {creal(arr[c][j]), cimag(arr[c][j])};
            fwrite(v, sizeof(double), 2, f);
        }
    fclose(f);
}

static void jarr(FILE *f, const char *name, const double *a, int n, int last) {
    fprintf(f, "  \"%s\": [", name);
    for (int q = 0; q < n; q++) fprintf(f, "%.17g%s", a[q], q < n - 1 ? ", " : "");
    fprintf(f, "]%s\n", last ? "" : ",");
}

static void write_grid_used(void) {
    FILE *f = xfopen("grid_used.json", "w");
    fprintf(f, "{\n  \"mesh\": \"%s\", \"grid_file\": \"%s\", \"grid_hash\": \"%s\",\n", P.mesh == 'u' ? "uniform" : "file",
            P.grid, grid_hash);
    fprintf(f, "  \"Nx\": %d, \"Ny\": %d, \"Nz\": %d, \"j0\": %d, \"ja\": %d, \"npml_lo\": %d, \"npml_hi\": %d,\n",
            Nx, Ny, Nz, j0, JA, NPlo, NPhi);
    fprintf(f, "  \"valid_j\": {\"Ex\": [0, %d], \"Ey\": [0, %d], \"Ez\": [0, %d], \"Hx\": [0, %d], \"Hy\": [0, %d], \"Hz\": [0, %d]},\n",
            Ny, Ny - 1, Ny, Ny - 1, Ny, Ny - 1);
    jarr(f, "x", xn, Nx + 1, 0); jarr(f, "x_dual", xd, Nx, 0); jarr(f, "hx", hxv, Nx, 0); jarr(f, "dx", dxv, Nx, 0);
    jarr(f, "y", yn, Ny + 1, 0); jarr(f, "y_dual", yd, Ny, 0); jarr(f, "hy", hyv, Ny, 0); jarr(f, "dy", dyv, Ny + 1, 0);
    jarr(f, "z", zn, Nz + 1, 0); jarr(f, "z_dual", zd, Nz, 0); jarr(f, "hz", hzv, Nz, 0); jarr(f, "dz", dzv, Nz, 0);
    jarr(f, "eps_t", epsT, Ny + 1, 0); jarr(f, "eps_n", epsN, Ny, 1);
    fprintf(f, "}\n");
    fclose(f);
}

static void write_meta(double runtime, long nsnap) {
    FILE *f = xfopen("meta.json", "w");
    double hymin = hyv[0], hymax = hyv[0], rmax = 1.0;
    for (int j = 0; j < Ny; j++) { hymin = fmin(hymin, hyv[j]); hymax = fmax(hymax, hyv[j]); }
    for (int j = 1; j < Ny; j++) rmax = fmax(rmax, fmax(hyv[j] / hyv[j - 1], hyv[j - 1] / hyv[j]));
    fprintf(f, "{\n");
    fprintf(f, "  \"units\": \"c = eps0 = mu0 = 1, lambda0 = 1\",\n");
    fprintf(f, "  \"nl\": %d, \"Delta\": %.17g, \"dt\": %.17g, \"S\": %.17g, \"omega0\": %.17g,\n", P.nl, D, dt, P.S, W0);
    fprintf(f, "  \"mesh\": \"%s\", \"grid_file\": \"%s\", \"grid_hash\": \"%s\", \"Delta_is_reference\": %s,\n",
            P.mesh == 'u' ? "uniform" : "file", P.grid, grid_hash, P.mesh == 'u' ? "false" : "true");
    fprintf(f, "  \"dt_courant\": %.17g, \"dt_penalty\": %.17g, \"Delta_y_min\": %.17g, \"Delta_y_max\": %.17g, "
               "\"r_max_y\": %.17g, \"x_nonuniform\": %d, \"z_nonuniform\": %d, \"auxref\": %d, \"divop\": \"%c\",\n",
            dt_cour, (P.S * hymax) / dt, hymin, hymax, rmax, xnonuni, znonuni, P.auxref, P.divop);
    fprintf(f, "  \"Lx\": %.17g, \"Lz\": %.17g, \"m\": %d, \"n\": %d, \"pol\": \"%c\",\n", P.Lx, P.Lz, P.m, P.n, P.pol);
    fprintf(f, "  \"Nx\": %d, \"Ny\": %d, \"Nz\": %d, \"j0\": %d, \"j1\": %d, \"npml\": %d, \"npml_lo\": %d, \"npml_hi\": %d, "
               "\"sf\": %d, \"tf\": %d, \"ja\": %d,\n",
            Nx, Ny, Nz, j0, j1, NPlo, NPlo, NPhi, j0 - NPlo, Ny - NPhi - j0, JA);
    fprintf(f, "  \"kx\": %.17g, \"ky\": %.17g, \"kz\": %.17g, \"ky_disc\": %.17g, \"ky_cont\": %.17g, \"Delta_a\": %.17g,\n",
            kx, ky, kz, ky_local(Da, 1.0), kyc, Da);
    fprintf(f, "  \"Ktilde\": [%.17g, %.17g, %.17g], \"Ktilde_src\": [%.17g, %.17g, %.17g], \"omega_tilde\": %.17g,\n",
            Kt[0], Kt[1], Kt[2], Kt_src[0], Kt_src[1], Kt_src[2], wt);
    fprintf(f, "  \"E0\": [%.17g, %.17g, %.17g], \"H0\": [%.17g, %.17g, %.17g],\n",
            E0v[0], E0v[1], E0v[2], H0v[0], H0v[1], H0v[2]);
    fprintf(f, "  \"offsets\": {");
    for (int c = 0; c < 6; c++)
        fprintf(f, "\"%s\": [%g, %g, %g, %g]%s", CNAME[c], OFF[c][0], OFF[c][1], OFF[c][2], OFF[c][3], c < 5 ? ", " : "},\n");
    fprintf(f, "  \"pml\": {\"m\": %g, \"sig_fac\": %g, \"kappa_max\": %g, \"alpha_max\": %.17g},\n",
            P.pml_m, P.sig_fac, P.kappa_max, P.alpha_max);
    fprintf(f, "  \"eps2\": %.17g, \"y1\": %d, \"ifmode\": \"%c\", \"grating\": %d, \"jsrc\": %d,\n", P.eps2, P.y1, P.ifmode,
            grating, P.jsrc);
    fprintf(f, "  \"inc\": \"%c\", \"ky_mode\": \"%s\", \"ramp\": \"%c\", \"ramp_T\": %g, \"erf_t0\": %g, \"erf_tau\": %g, \"off_t\": %g, \"na\": %d,\n",
            P.inc, P.ky_cont ? "continuous" : "discrete", P.ramp, P.ramp_T, P.erf_t0, P.erf_tau, P.off_t, j0 - JA);
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
    fprintf(f, "  \"dump_at\": [");
    for (int q = 0; q < P.ndump; q++) fprintf(f, "%ld%s", P.dump_at[q], q < P.ndump - 1 ? ", " : "");
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

static double urand(unsigned long long *s) { /* xorshift64*, uniform in [-1, 1) */
    *s ^= *s >> 12; *s ^= *s << 25; *s ^= *s >> 27;
    return (double)((*s * 2685821657736338717ULL) >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}

/* ------------------------------------------------------------------ power iteration (gate 0-4)
   M = eps^-1 C_H mu^-1 C_E is self-adjoint and positive semi-definite in <.,.>_{eps W_E} (derivation_nonuniform §2).
   One application uses the solver's own update kernels: H = 0, H -= dt C_E v; E = 0, E += dt eps^-1 C_H H,
   so M v = -E / dt^2. CPML is switched off (PEC walls at j = 0, Ny; x, z periodic); the leapfrog is stable iff
   dt < dt_max = 2 / sqrt(lambda_max). Writes eig.json. */
static double enorm2(const double *X, const double *Y, const double *Z, const double *X2, const double *Y2,
                     const double *Z2) { /* <E, E2>_{eps W_E} */
    double s = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : s)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++)
            for (int k = 0; k < Nz; k++) {
                size_t q = ID(i, j, k);
                size_t m = (size_t)i * SY + j;
                s += dyv[j] * (eTx[m] * hxv[i] * dzv[k] * X[q] * X2[q] + eTz[m] * dxv[i] * hzv[k] * Z[q] * Z2[q]);
                if (j < Ny) s += eNy[m] * dxv[i] * hyv[j] * dzv[k] * Y[q] * Y2[q];
            }
    return s;
}

static void power_iteration(void) {
    for (int j = 0; j <= Ny; j++) pmlE[j] = pmlH[j] = -1;
    double *vx = xcalloc(NTOT, sizeof(double)), *vy = xcalloc(NTOT, sizeof(double)), *vz = xcalloc(NTOT, sizeof(double));
    unsigned long long s = 0x9E3779B97F4A7C15ULL ^ (unsigned long long)P.seed;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j <= Ny; j++)
            for (int k = 0; k < Nz; k++) {
                size_t q = ID(i, j, k);
                double a = urand(&s), b = urand(&s), c = urand(&s);
                vx[q] = (j == 0 || j == Ny) ? 0.0 : a;
                vz[q] = (j == 0 || j == Ny) ? 0.0 : c;
                vy[q] = (j == Ny) ? 0.0 : b;
            }
    double nv = sqrt(enorm2(vx, vy, vz, vx, vy, vz));
    for (size_t q = 0; q < NTOT; q++) { vx[q] /= nv; vy[q] /= nv; vz[q] /= nv; }
    double lam = 0.0, lam_old = 0.0;
    long it;
    FILE *lg = xfopen("eig_log.csv", "w");
    fprintf(lg, "it,rayleigh\n");
    for (it = 1; it <= P.eig_maxit; it++) {
        memset(Hx, 0, NTOT * sizeof(double)); memset(Hy, 0, NTOT * sizeof(double)); memset(Hz, 0, NTOT * sizeof(double));
        memcpy(Ex, vx, NTOT * sizeof(double)); memcpy(Ey, vy, NTOT * sizeof(double)); memcpy(Ez, vz, NTOT * sizeof(double));
        update_H(0);
        memset(Ex, 0, NTOT * sizeof(double)); memset(Ey, 0, NTOT * sizeof(double)); memset(Ez, 0, NTOT * sizeof(double));
        update_E(0);
        double sc = -1.0 / (dt * dt);
        for (size_t q = 0; q < NTOT; q++) { Ex[q] *= sc; Ey[q] *= sc; Ez[q] *= sc; } /* M v */
        lam = enorm2(vx, vy, vz, Ex, Ey, Ez);                                         /* <v, Mv>, |v| = 1 */
        double nm = sqrt(enorm2(Ex, Ey, Ez, Ex, Ey, Ez));
        for (size_t q = 0; q < NTOT; q++) { vx[q] = Ex[q] / nm; vy[q] = Ey[q] / nm; vz[q] = Ez[q] / nm; }
        if (it % 100 == 0) fprintf(lg, "%ld,%.17e\n", it, lam);
        if (it > 10 && fabs(lam - lam_old) < P.eig_tol * lam) break;
        lam_old = lam;
    }
    fclose(lg);
    double dtmax = 2.0 / sqrt(lam);
    FILE *f = xfopen("eig.json", "w");
    fprintf(f, "{\n  \"lambda_max\": %.17g, \"dt_max\": %.17g, \"dt_courant\": %.17g, \"S\": %.17g,\n"
               "  \"iterations\": %ld, \"last_rel_change\": %.3e, \"grid_hash\": \"%s\", \"Nx\": %d, \"Ny\": %d, \"Nz\": %d\n}\n",
            lam, dtmax, dt_cour / P.S / sqrt(3.0), P.S, it, fabs(lam - lam_old) / lam, grid_hash, Nx, Ny, Nz);
    fclose(f);
    fprintf(stderr, "power iteration: lambda_max = %.15g, dt_max = %.15g (Courant bound %.15g), %ld iterations\n",
            lam, dtmax, dt_cour / P.S / sqrt(3.0), it);
}

/* ------------------------------------------------------------------ main */

int main(int argc, char **argv) {
    parse_args(argc, argv);
    setup();
    MKDIR(P.out);
    if (P.mode == 'e') {
        power_iteration();
        write_grid_used();
        write_meta(0.0, 0);
        return 0;
    }
    if (P.dft0 < 0 || P.dft1 < 0) { P.dft0 = 0; P.dft1 = 0; }
    dft_setup();

    if (P.init == 'a') {
        if (NPlo != 0 || NPhi != 0) die("init=analytic requires npml=0");
        fill_analytic(0.0, -1.0); /* E^0, H^{-1/2} */
    } else if (P.init == 'b') {
        /* Divergence-free packet: Ez depends only on (x, y), Ex only on (y, z), so div E = 0 exactly
           and no static charge is left behind. Carrier ky0 = 5 under a Gaussian of width 1.2, so the
           near-grazing (ky ~ 0) content is ~exp(-(ky0 w)^2/4) ~ 1e-4 in amplitude. */
        double yc = 0.5 * (yn[0] + yn[Ny]), w = 1.2, ky0 = 5.0;
        if (P.mesh == 'u') yc = 0.5 * Ny * D;
        for (int i = 0; i < Nx; i++)
            for (int j = 1; j < Ny; j++)
                for (int k = 0; k < Nz; k++) {
                    double y = yn[j], g = cos(ky0 * (y - yc)) * exp(-(y - yc) * (y - yc) / (w * w));
                    Ez[ID(i, j, k)] = cos(2.0 * PI * xn[i] / P.Lx) * g;
                    Ex[ID(i, j, k)] = cos(2.0 * PI * zn[k] / P.Lz) * g;
                }
    } else if (P.init == 'r') {
        /* random E^0, H^{-1/2} (PEC walls: tangential E = 0 at j = 0, Ny) -- energy / stability tests */
        unsigned long long s = 0x9E3779B97F4A7C15ULL ^ (unsigned long long)P.seed;
        for (int c = 0; c < 6; c++) {
            double *F = comp_ptr(c);
            int half = (OFF[c][1] != 0.0);
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j <= Ny; j++)
                    for (int k = 0; k < Nz; k++) {
                        double v = urand(&s);
                        if (half && j == Ny) v = 0.0;
                        if ((c == 0 || c == 2) && (j == 0 || j == Ny)) v = 0.0;
                        F[ID(i, j, k)] = v;
                    }
        }
    }

    partT = xcalloc(Nx, sizeof(double));
    partP = xcalloc(Nx, sizeof(double));
    FILE *lg = xfopen("log.csv", "w");
    fprintf(lg, "n,t,maxE_SF,maxH_SF,W_total,W_phys,divE_TF,divH_TF,W_mod,leakE_ref,leakH_ref,devE_ref,devH_ref\n");
    FILE *sn = NULL;
    int snapc = 2;
    for (int c = 0; c < 6; c++) if (!strcmp(P.snapcomp, CNAME[c])) snapc = c;
    if (P.snap_every > 0) {
        if (P.zk < 0) die("snapshots need zk >= 0");
        sn = xfopen("snapshots.bin", "wb");
    }
    long nsnap = 0;
    write_grid_used();
    if (P.inc == 'm') modal_write();

    time_t w0 = time(NULL);
    for (long n = 0; n < P.nsteps; n++) {
        int doW = (P.energy_every > 0 && (n + 1) % P.energy_every == 0);
        /* H: n-1/2 -> n+1/2 */
        int tfsf = (P.inc == 'a' || P.inc == 'n' || P.inc == 'p' || P.inc == 'm');
        if (tfsf) incident_E(n);
        update_H(doW);
        double Wm = doW ? Wmod_acc : NAN;     /* energy at t = n dt */
        if (tfsf) tfsf_H();
        if (P.auxref) auxref_update_H();
        if (P.inc == 'a') aux_update_H(n);
        if (P.inc == 'm') for (int t = 0; t < NMA; t++) maux_update_H(&MA[t], n);
        if (n >= P.dft0 && n < P.dft1) dft_accumulate(3, 5, (n + 0.5) * dt);
        /* E: n -> n+1 */
        if (tfsf) incident_H(n);
        update_E(doW);
        if (tfsf) tfsf_E();
        if (P.inc == 'j') source_J(n);
        if (P.auxref) auxref_update_E();
        if (P.inc == 'a') aux_update_E(n);
        if (P.inc == 'm') for (int t = 0; t < NMA; t++) maux_update_E(&MA[t], n);
        if (P.init == 'a') dirichlet_y((double)(n + 1));
        if (n + 1 >= P.dft0 && n + 1 < P.dft1) dft_accumulate(0, 2, (n + 1) * dt);

        double mE, mH, Wt = NAN, Wp = NAN, dvE = NAN, dvH = NAN, lE = NAN, lH = NAN, gE = NAN, gH = NAN;
        sf_max(&mE, &mH);
        if (doW) energy_sum(&Wt, &Wp); /* W at t = (n + 1/2) dt */
        if (P.div_every > 0 && (n + 1) % P.div_every == 0) div_max(&dvE, &dvH);
        if (P.auxref) {
            ref_dev(NPlo, j0 - 1, &lE, &lH);   /* scattered-field region (half-y rows up to j0 - 1/2) */
            if (P.ref_every > 0 && (n + 1) % P.ref_every == 0) ref_dev(0, Ny, &gE, &gH);
        }
        fprintf(lg, "%ld,%.10g,%.9e,%.9e,%.17e,%.17e,%.6e,%.6e,%.17e,%.6e,%.6e,%.6e,%.6e\n", n + 1, (n + 1) * dt, mE, mH,
                Wt, Wp, dvE, dvH, Wm, lE, lH, gE, gH);
        if (!isfinite(mE) || !isfinite(mH) || (doW && !isfinite(Wm))) {
            fflush(lg);
            fprintf(stderr, "non-finite field at step %ld\n", n + 1);
            die("instability detected");
        }
        if (sn && (n + 1) % P.snap_every == 0) {
            double *F = comp_ptr(snapc);
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j <= Ny; j++) fwrite(&F[ID(i, j, P.zk)], sizeof(double), 1, sn);
            nsnap++;
        }
        for (int q = 0; q < P.ndump; q++)
            if (P.dump_at[q] == n + 1) {
                char fn[64];
                snprintf(fn, sizeof fn, "fields_n%ld.bin", n + 1);
                dump_fields(fn);
                if (P.auxref) dump_auxref(n + 1);
            }
        if (P.nsteps >= 10 && (n + 1) % (P.nsteps / 10) == 0) {
            fprintf(stderr, "  step %ld/%ld (%.0f s)\n", n + 1, P.nsteps, difftime(time(NULL), w0));
        }
    }
    fclose(lg);
    if (sn) fclose(sn);
    double runtime = difftime(time(NULL), w0);
    write_dft();
    if (P.dump) dump_fields("fields_final.bin");
    write_meta(runtime, nsnap);
    fprintf(stderr, "done: %ld steps, grid %d x %d x %d, %.0f s\n", P.nsteps, Nx, Ny + 1, Nz, runtime);
    return 0;
}
