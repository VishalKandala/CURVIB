#include "petscvec.h"
#include "petscdmda.h"
#include "petscksp.h"
#include "petscsnes.h"
#include <unistd.h>
#include<string.h>
#include<stdio.h> 

#include <stdlib.h>

// Seokkoo
#if defined (__CPLUSPLUS__) || defined (__cplusplus)
#include <vector>
#include <algorithm>
#include <assert.h>
#endif

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  PetscReal t, f;
} FlowWave;

typedef struct {
	PetscScalar	x, y;
} Cmpnts2;

typedef struct {
  PetscReal	x, y;
} Cpt2D;

typedef struct {
  Vec	Ubcs; // An ugly hack, waste of memory
} BCS;

typedef struct {
  PetscInt	i1, j1, k1;
  PetscInt	i2, j2, k2;
  PetscInt	i3, j3, k3;
  PetscReal	cr1, cr2, cr3; // coefficients
  PetscReal	d_i; // distance to interception point on grid cells
  PetscInt	imode; // interception mode

  PetscInt	ni, nj, nk;
  PetscReal	d_s; // shortest distance to solid surfaces
  Cmpnts	pmin;
  PetscInt	cell; // shortest distance surface cell
  PetscReal	cs1, cs2, cs3;

  PetscInt	i11, j11, k11;
  PetscInt	i22, j22, k22;
  PetscInt	i33, j33, k33;
  PetscReal	cr11, cr22, cr33; // coefficients
  PetscReal	d_ii; // distance to interception point on grid cells
  PetscInt	iimode; // interception mode
  PetscReal	cs11, cs22, cs33;

  PetscInt	ii1, jj1, kk1;
  PetscInt	ii2, jj2, kk2;
  PetscInt	ii3, jj3, kk3;
  PetscReal	ct1, ct2, ct3; // coefficients
  //PetscReal	d_s; // distance to interception point on grid cells
  PetscInt	smode; // interception mode

  PetscInt	ii11, jj11, kk11;
  PetscInt	ii22, jj22, kk22;
  PetscInt	ii33, jj33, kk33;
  PetscReal	ct11, ct22, ct33; // coefficients
  PetscReal	d_ss; // distance to interception point on grid cells
  PetscInt	ssmode; // interception mode

  
/*   PetscInt      bi1, bj1, bk1; */
/*   PetscInt      bi2, bj2, bk2; */
/*   PetscInt      bi3, bj3, bk3; */
/*   PetscInt      bi4, bj4, bk4; */

/*   PetscReal     bcr1,bcr2,bcr3,bcr4; */
} IBMInfo;

typedef struct {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;	// number of vertices and number of elements
  PetscInt	*nv1, *nv2, *nv3;	// Node index of each triangle
  PetscReal	*nf_x, *nf_y, *nf_z;	// Normal direction
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM surface nodes
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
  PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270];
  Cmpnts	*u, *uold, *urm1;
  
  // Added 4/1/06 iman
  PetscReal     *dA ;         // area of an element
  // Added 4/2/06 iman
  PetscReal     *nt_x, *nt_y, *nt_z; //tangent dir
  PetscReal     *ns_x, *ns_y, *ns_z; //azimuthal dir
  // Added 6/4/06 iman
  //Cmpnts        *cent; // center of the element 
  PetscReal     *cent_x,*cent_y,*cent_z;

  // for forces on each element
  PetscReal      *pres, *tau0, *tauN;
  PetscReal      *Bvel_u, *Bvel_v, *Bvel_w;
  PetscReal       x_min,x_max,y_min,y_max,z_min,z_max;
  // for radius check
  Cmpnts *qvec;
  PetscReal *radvec; 
} IBMNodes;
typedef struct {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;	// number of vertices and number of elements
  PetscInt	*nv1, *nv2, *nv3, *nv4;	// Node index of each tetra
  
  PetscReal	*x_bp, *y_bp, *z_bp;	// Coordinates of IBM volume nodes
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
  
  Cmpnts	*u, *uold, *urm1;
  PetscReal      V;
  
  PetscReal     *dV0 ;         // volume of an element

  PetscReal     *cent_x,*cent_y,*cent_z;
  PetscReal     x_c,y_c,z_c;
  PetscReal     J[3][3];
  PetscReal     I_inv[3][3];
} IBMVNodes;
typedef struct node{
  PetscInt Node;
  struct node *next;
} node;

typedef struct list{
  node *head;
} List;


typedef struct list_node {
  PetscInt	index;
  struct list_node *next;
} Node_List;

typedef struct IBMListNode {
  IBMInfo ibm_intp;
  struct IBMListNode* next;
} IBMListNode;

typedef struct IBMList {
  IBMListNode *head;
} IBMList;


typedef struct UserCtx {
  DM da;	/* Data structure for scalars (include the grid geometry
		   informaion, to obtain the grid information, 
		   use DAGetCoordinates) */
  DM fda;	// Data Structure for vectors
  DMDALocalInfo info;

  Vec	Cent;	// Coordinates of cell centers
  Vec   Centx,Centy,Centz;
  Vec 	Csi, Eta, Zet, Aj;
  Vec 	ICsi, IEta, IZet, IAj;
  Vec 	JCsi, JEta, JZet, JAj;
  Vec 	KCsi, KEta, KZet, KAj;
  //  Vec   Area,lArea, Volume, lVolume;
  Vec 	Ucont,Vcont;	// Contravariant velocity components
  Vec 	Ucat, Wcat;	// Cartesian velocity components
  Vec	Ucat_o;
  Vec	Ucont_o, Ucont_rm1, Rhs, dUcont, pUcont;
  Vec 	P, P_o;
  ///////////////////////////
  Vec   Qnew;
  Vec   Ql;
  Vec   shrr;
  Vec   ItfcQ;//, ItfcP;
  Vec	nhostQ;//, nhostP;
  Vec   lItfcQ;//, lItfcP;
  PetscReal LVOUT;


  PetscReal cdisx,cdisy,cdisz;
  ///////////////////////////
  Vec	Phi;
  Vec	GridSpace;
  Vec	Nvert;
  Vec	Nvert_o;
  Vec   Itfc;//, ItfcP;
  BCS	Bcs;

  Vec	lUcont, lUcat, lP, lPhi, lVcont,lWcat;
  Vec   lUcont_o, lUcont_rm1;//, ldUcont;
  Vec	lCsi, lEta, lZet, lAj;
  Vec	lICsi, lIEta, lIZet, lIAj;
  Vec	lJCsi, lJEta, lJZet, lJAj;
  Vec	lKCsi, lKEta, lKZet, lKAj;
  Vec	lGridSpace;
  Vec	lNvert, lNvert_o, lNFace;
  Vec	lCent;
  Vec   lItfc;//, lItfcP;

  Vec   lMAreaCsi, lMAreaEta, lMAreaZet; //Moments of Area

  Vec	inletU;

  Vec	nhostU;//, nhostP;


  Vec	DUold;

  Vec	Forcing;
  Vec	Ucont_MG;

  Vec	Dt, Nu_t, CS;

  AO	ao;

  PetscReal	ren;	// Reynolds number
  PetscReal	dt; 	// time step
  PetscReal	st;	// Strouhal number

  // added for new poisson solver
  PetscReal     FluxIntpSum;

  // Added Iman
  Vec           psuedot;
  PetscReal     cfl, vnn;

  PetscReal	r[101], tin[101], uinr[101][1001];

  PetscInt      ip[10],jp[10],kp[10],dispx[10],dispy[10],dispz[10],dispnn[10];
  PetscReal     *itfchostx, *itfchosty, *itfchostz;
  PetscReal     FluxInSum, FluxOutSum, FluxIntfcSum;
  PetscReal                AreaOutSum, AreaIntfcSum;
/*   PetscReal	lFluxOutSum[6];  */
/*   PetscReal	lAreaSum[6]; */

  PetscErrorCode aotopetsc;
  PetscBool assignedA;

  PetscInt _this;
  PetscInt *idx_from;

  PetscInt bctype[6],inttype[7];
  PetscInt itfcptsnumber;
  PetscInt *itfcI, *itfcJ, *itfcK;
  PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
  PetscInt	IM, JM, KM; // dimensions of grid
  PetscReal Max_X,Max_Y,Max_Z,Min_X,Min_Y,Min_Z;

  PetscInt ibmnumber;
  IBMInfo  *ibm_intp;
  Mat	A, C;
  KSP   ksp;

  IBMNodes *ibm;

  DM	*da_f, *da_c;
  struct UserCtx *user_f, *user_c;
  Vec	*lNvert_c;

  Vec	B;
  Vec	Rhsp, X, R;
  
  Mat	MR, MP;
  MatNullSpace nullsp;

  
  /* Variables for multi-nullspace case */
  PetscInt *KSKE;
  PetscBool multinullspace;

  IBMList *ibmlist;

  PetscInt thislevel, mglevels;

  PetscInt isc, jsc, ksc;
  PetscInt cgrid;  

  FlowWave *inflow;
  PetscInt number_flowwave;

  // New
  Vec Ucat_sum;			// sum of u, v, w over time
  Vec Ucat_cross_sum;	// sum of uv, vw, wu
  Vec Ucat_square_sum;	// sum of uu, vv, ww
  Vec P_sum;				// sum of p
  Vec lSx, lSy, lSz, lS;  
  Vec lLM, lMM, lNM;
  //  std::vector<double> k_area;
  Vec	lNu_t, lF1;		// eddy viscosity, f1 constant for SST model
  Vec	lCs;
  Vec	K_Omega, lK_Omega, K_Omega_o, lK_Omega_o;//, K_Omega_rm1;
  Vec	Distance;
  Vec   lItfcKO;
  // for hypre
  int	rhs_count;
  Vec	Gid, Gidm;				// seokkoo, Global ID for local array containg ghost points
  Vec	Phi2, B2, Ucont2;
  int	local_Phi2_size, p_global_begin;
  int	reduced_p_size;

  DM        fda2;	// Data Structure for vectors with 2 variables
  Vec	    lUstar;
  Cmpnts2   **komega_plane;

   // for sens
  SNES snes;
  Mat  J,PJ,J2;//Jacobian
  Vec  RFC;
  ISColoring iscoloring;
  MatFDColoring matfdcoloring;
  

} UserCtx;

typedef struct {
  UserCtx *user;
  PetscInt thislevel;
  DM       packer;
  // DMMG        *dmmg;
} MGCtx;

typedef struct {
  PetscInt mglevels;
  PetscInt thislevel;

  PetscBool  isc, jsc, ksc;
  MGCtx *mgctx;
  //packer
  DM       packer;
  SNES     snespacker;
} UserMG;

typedef struct {
  PetscReal     P;    //Press on the surface elmt
  PetscInt      n_P; //number of Press Pts on the elmt
  PetscReal     Tow_ws, Tow_wt; //wall shear stress of the elmt
  PetscReal     Tow_wn; // normal stress 

  PetscInt      Clsnbpt_i,Clsnbpt_j,Clsnbpt_k; //Closest Near Bndry Pt to each surf elmt 
  PetscInt      icell,jcell,kcell;
  PetscInt      FoundAroundcell;
  PetscInt      Need3rdPoint;
  //PetscInt      Aroundcellnum;
} SurfElmtInfo;

typedef struct {
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
  PetscReal    F_x_old,F_y_old,F_z_old; //Forces & Area
  PetscReal    F_x_real,F_y_real,F_z_real; //Forces & Area
  PetscReal    M_x,M_y,M_z; // Moments
  PetscReal    M_x_old,M_y_old,M_z_old; //Forces & Area
  PetscReal    M_x_real,M_y_real,M_z_real; //Forces & Area
  PetscReal    M_x_rm2,M_y_rm2,M_z_rm2; //Forces & Area
  PetscReal    M_x_rm3,M_y_rm3,M_z_rm3; //Forces & Area
  PetscReal    x_c,y_c,z_c; // center of rotation(mass)
  PetscReal    a_c[3];  //initial center of rotataion
  PetscReal    Mdpdn_x, Mdpdn_y,Mdpdn_z;
  PetscReal    Mdpdn_x_old, Mdpdn_y_old,Mdpdn_z_old;
  PetscReal    Power; //power

  PetscReal    clone;
  PetscInt     pbc[3];
  PetscReal    I_inv[3][3];
  PetscReal    L_n[3],L_o[3],L_r[3]; //Angular momentum
  PetscReal    alpha[3];
  PetscReal    acc[3];
  PetscReal    R[3][3],q[4],q_r[4];
  // Aitkin's iteration
  PetscReal    dS[6],dS_o[6],atk,atk_o;

  // for force calculation
  SurfElmtInfo  *elmtinfo;
  IBMInfo       *fsi_intp;

  //Max & Min of ibm domain where forces are calculated
  PetscReal    Max_xbc,Min_xbc;
  PetscReal    Max_ybc,Min_ybc;
  PetscReal    Max_zbc,Min_zbc;

  // CV bndry
  PetscInt     CV_ys,CV_ye,CV_zs,CV_ze;
} FSInfo;

typedef struct {
  PetscInt   n_time, n_midp, n_subit; //number of time steps and number of control points
  PetscReal  *x_midp, *y_midp, *z_midp; //control points x/y_midp[time][n_midp]
  PetscReal  *x_com, *y_com, *head_ang; //position of center of mass and head angle
  PetscReal  *s1,*s2,*s3; /* spline coefficients for spatial interpolation
			     S(xx) = Y_mipd[i] + s1[i] * w + s2[i] * w**2 + s3[i] * w**3
			     where w = xx - x[i]
			  */
  PetscReal  *st1,*st2,*st3; /* spline coefficients for the temporal interpolation    
				S(xx) = Y_mipd[i] + s1[i] * w + s2[i] * w**2 + s3[i] * w**3
			extern int testfilter_ik, testfilter_1d;
	where w = xx - x[i]
			     */
  Mat        Mphi;
  PetscReal  xmin,xmax,ymin,ymax,zmin,zmax;//min max of ibm
} Cstart;



#define COEF_TIME_ACCURACY 1.5
//PetscReal COEF_TIME_ACCURACY=1.5;
extern int i_periodic, j_periodic, k_periodic;
extern int ii_periodic, jj_periodic, kk_periodic;
extern int les, wallfunction, rans;
extern int i_homo_filter, j_homo_filter, k_homo_filter;
extern int testfilter_ik, testfilter_1d;
extern int clark;

#if defined (__CPLUSPLUS__) || defined (__cplusplus)
PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm);

void initlist(List *ilist);
void insertnode(List *ilist, PetscInt Node);
void destroy(List *ilist);
void InitIBMList(IBMList *ilist);
void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp);
void DestroyIBMList(IBMList *ilist);

PetscInt intsect_triangle(PetscReal orig[3], PetscReal dir[3],
			  PetscReal vert0[3], PetscReal vert1[3], PetscReal vert2[3],
			  PetscReal *t, PetscReal *u, PetscReal *v);
PetscInt ISPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3,
			   PetscReal nfx, PetscReal nfy, PetscReal nfz);
PetscInt ISLineTriangleIntp(Cmpnts p1, Cmpnts p2, IBMNodes *ibm, PetscInt ln_v);
PetscErrorCode triangle_intp2(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			      IBMInfo *ibminfo);
PetscErrorCode triangle_intp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			     IBMInfo *ibminfo);
PetscErrorCode triangle_intpp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, 
			      IBMInfo *ibminfo);
PetscErrorCode Dis_P_Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, 
			  PetscReal *d);

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	     Cmpnts *Ub, double nx, double ny, double nz);
void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	       Cmpnts *Ub, double nx, double ny, double nz);
void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, 
		    PetscReal *ustar, double nx, double ny, double nz);

PetscErrorCode ibm_read_tecplot(IBMNodes *ibm, PetscInt nt, PetscBool flg);
PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_read_Ansys(IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_read_Icem(IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_read_Ansys_vol(IBMVNodes *ibmv, PetscInt ibi);
PetscErrorCode ibmv_cent_of_mass(IBMVNodes *ibmv, FSInfo *FSinfo,PetscInt ibi);
PetscErrorCode ibm_surface_out(IBMNodes *ibm, PetscInt ti,PetscInt ibi);
PetscErrorCode ibm_surface_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti);
PetscErrorCode ibm_volume_VTKOut(IBMVNodes *ibmv, PetscInt ibi, PetscInt ti);
PetscErrorCode calc_ibm_normal(IBMNodes *ibm);
PetscErrorCode calc_ibm_velocity(IBMNodes *ibm, PetscReal delti);
PetscErrorCode calc_ibm_volumeFlux(IBMNodes *ibm, PetscReal delti, PetscReal *VolumeFlux);
PetscErrorCode ibm_scale(IBMNodes *ibm, PetscReal scale[]);
PetscErrorCode ibm_placement(IBMNodes *ibm, FSInfo *fsi,
			     PetscInt Nibm_y, PetscInt Nibm_z,
			     PetscReal dYibm, PetscReal dZibm, 
			     PetscInt ibi);

PetscErrorCode ibm_read_fish(IBMNodes *ibm,PetscReal cl);
PetscErrorCode fish_init(PetscReal *delti);
PetscErrorCode fish_swim(IBMNodes *ibm, PetscReal time
			 ,PetscReal delti);
PetscErrorCode wing_motion(IBMNodes *ibm, FSInfo *fsi, PetscInt ti, PetscReal *dt);
PetscErrorCode read_midpoint(Cstart *cstart);
PetscErrorCode form_cstart(Cstart *cstart,IBMNodes *ibm, FSInfo *fsi,
			   PetscReal dt, PetscInt ibi, PetscInt nonIntertial);

PetscErrorCode cop_init(IBMNodes *ibm);
PetscErrorCode cop_init_time(PetscInt regime, PetscReal *delti);
PetscErrorCode cop_swim(IBMNodes *ibm, PetscReal time
			,PetscReal delti);
PetscErrorCode Initalize_Projecting(IBMNodes * ibm );
PetscErrorCode FsiInitialize(PetscInt n_elmt, FSInfo *fsi,PetscInt ibi);
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ibi, PetscInt ti);
PetscErrorCode FSI_DATA_Output(FSInfo *FSinf, PetscInt ti);
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, 
				 PetscReal dt, PetscInt ibi);
PetscErrorCode Elmt_Move_FSI_TRANS(FSInfo *FSinfo, IBMNodes *ibm
				   , PetscInt ibi);

PetscErrorCode Calc_forces_SI(FSInfo *FSinfo,UserCtx *user,
			      IBMNodes *ibm,PetscInt ti, 
			      PetscInt ibi, PetscInt bi);
PetscErrorCode Calc_FSI_pos_intg(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt);
PetscErrorCode Calc_FSI_Ang_intg(FSInfo *FSinfo,IBMNodes *ibm, 
				 PetscReal dt,PetscInt itrSC,
				 PetscInt ibi, UserCtx *user, PetscReal S_init);
PetscErrorCode Calc_FSI_Ang(FSInfo *FSinfo,IBMNodes *ibm, 
			    PetscReal dt, PetscReal dtime,
			    PetscInt ibi, UserCtx *user);
PetscErrorCode Calc_FSI_vel_particle(FSInfo *FSinfo,IBMNodes *ibm,PetscReal dt);
PetscErrorCode calc_quarternion(FSInfo *FSinfo,PetscReal dt,PetscInt ibi);
PetscErrorCode CollisionDetectionOfCylinders(FSInfo *fsi,IBMNodes *ibm, 
					     PetscInt NumberOfBodies);
PetscErrorCode CollisionDetectionOfParticles(UserCtx *user,IBMNodes *ibm,FSInfo *fsi);
PetscErrorCode Forced_Motion(FSInfo *fsi,
			     PetscReal dt);
PetscErrorCode SwingCylinder(FSInfo *fsi, IBMNodes *ibm);

PetscErrorCode LV_beat(IBMNodes *ibm, PetscReal time
		       ,PetscReal delti, PetscInt ibi);

PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, 
				   PetscInt ibi);
PetscErrorCode ibm_search_advanced_rev(UserCtx *user, IBMNodes *ibm, 
				       PetscInt ibi);

PetscErrorCode ibm_interpolation_advanced(UserCtx *user, 
					  IBMNodes *ibm,			  
					  PetscInt ibi,
					  PetscInt Add_dUndt);

PetscErrorCode SetInitialGuessToOne(UserCtx *user);
PetscErrorCode InflowFlux(UserCtx *user);
PetscErrorCode OutflowFlux(UserCtx *user);
PetscErrorCode FormBCS(UserCtx *user);
PetscErrorCode Block_Interface_U(UserCtx *user);
PetscErrorCode Block_Blank_Correction_adv(UserCtx *user, //Vec lUcor, 					  
					  PetscInt flg) ;
PetscErrorCode Block_Interface_Correction(UserCtx *user);
PetscErrorCode InletRead(UserCtx *user);
PetscErrorCode fluxin(UserCtx *user);

PetscErrorCode FormFunction1(UserCtx *user, Vec Rhs);
PetscErrorCode CalcFrameVelocity(UserCtx *user,  Cmpnts u_c, 
				 Cmpnts omega_c, Cmpnts a_c);

PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, 
			    FSInfo *fsi, Cstart *cstart,
			    PetscInt itr_sc,
			    PetscInt tistart, 
			    PetscBool *DoSCLoop);
PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, 
			   FSInfo *fsi);

void Calculate_dxdydz(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz);
void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, double *dudc, double *dvdc, double *dwdc, double *dude, double *dvde, double *dwde, double *dudz, double *dvdz, double *dwdz);
void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***K, PetscReal ***nvert, double *dkdc, double *dkde, double *dkdz);
void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
					double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
					double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );
void Compute_dscalar_dxyz ( double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
							double dkdc, double dkde, double dkdz, double *dk_dx, double *dk_dy, double *dk_dz);

void K_Omega_IC(UserCtx *user);
void K_Omega_Set_Constant(UserCtx *user);
void Solve_K_Omega(UserCtx *user);
PetscErrorCode Solve_K_Omega_DMMG(UserCtx *user);
void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat);
void Compute_eddy_viscosity_LES(UserCtx *user);
PetscErrorCode KE_Output(UserCtx *user);

PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImplicitMomentumSolver1(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Implicit_DMMG(UserMG *usermg);
PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Implicit_FD_Jacobian(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Initialize_Hypre_Var(PetscInt nb);
void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonSolver_MG(UserMG *usermg);
PetscErrorCode UpdatePressure(UserCtx *user);
PetscErrorCode Projection(UserCtx *user);
PetscErrorCode Divergence(UserCtx *user);


PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);

PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user);
PetscErrorCode GridDivergence(UserCtx *user);
PetscErrorCode Contra2Cart(UserCtx *user);

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user, PetscInt bi);
PetscErrorCode Viscosity(UserCtx *user,FSInfo *FSinfo,IBMNodes *ibm);
PetscErrorCode MG_Finalize(UserMG *usermg);
PetscErrorCode Bed_Change(IBMNodes *ibm, PetscReal delti, PetscInt tstep);
PetscErrorCode rotate_ibm(IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode rotate_force(FSInfo *fsi, PetscInt ibi, PetscInt bi);
PetscErrorCode rotate_quarternion(IBMNodes *ibm, FSInfo *fsi, PetscReal dt);
PetscErrorCode FormAreaMoment(UserCtx *user, Cmpnts a_c);
PetscErrorCode MomentAreaDivergence(UserCtx *user);
PetscErrorCode MomentumJet(UserCtx *user, PetscInt Kplane);
PetscErrorCode Particle_Duplicate_Indentify (UserCtx *user,IBMNodes *ibm, FSInfo *fsi, PetscInt *X, PetscInt *Z);
PetscErrorCode Particle_Clone (IBMNodes *ibm1, FSInfo *fsi1,IBMNodes *ibm2, FSInfo *fsi2,PetscInt X, PetscInt Z, PetscReal dt,PetscInt ibi);
PetscErrorCode ibm_duplicate(IBMNodes *ibm,IBMNodes *ibm0,PetscReal r_x, PetscReal r_y, PetscReal r_z);

PetscErrorCode ComputelNface(UserCtx *ctx);
PetscInt lidxLocal(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscErrorCode CheckRowDiagonal_packer(Mat J,UserCtx *ctx,PetscInt i,PetscInt j,PetscInt k,PetscInt direction);
PetscErrorCode CheckRow_packer(Mat J,UserCtx *ctx,PetscInt i,PetscInt j,PetscInt k,PetscInt direction);
PetscErrorCode ComputelNface(UserCtx *ctx);
PetscInt lidxLocal(PetscInt i, PetscInt j, PetscInt k, UserCtx *user);
PetscErrorCode CheckCol_packer(Mat J,UserCtx *ctx,PetscInt m,PetscInt n,PetscInt p);
PetscErrorCode  Implicit_SNES(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode  Implicit_SNES_Packer(UserMG *usermg,IBMNodes *ibm, FSInfo *fsi);


void Calculate_normal(struct Cmpntsgdb csi, struct Cmpntsgdb eta, struct Cmpntsgdb zet, double ni[3], 
		      double nj[3], double nk[3]);

/////////////////////////////////////data
PetscInt point_cell_advanced(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
			     IBMNodes *ibm, PetscInt ncx, PetscInt ncy,
			     PetscInt ncz, PetscReal dcx, PetscReal dcy,
			     PetscReal xbp_min, PetscReal ybp_min,
			     PetscReal zbp_max, List *cell_trg,
			     PetscInt flg);


PetscErrorCode Interpolation_matrix(UserMG *usermg);

#endif
