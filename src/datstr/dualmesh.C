/*
 * dualmesh.C
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#include "../header/hpfem.h"

SolRec::SolRec(double *doublekeyrangein, int size, int prime, double XR[], double YR[],
    int ifrestart) :
		HashTable(doublekeyrangein, size, prime, XR, YR, ifrestart) {

	double zero_sol[3]={0.,0.,0.};
	double kact_z=0.;
	solution_zero=new Solution(zero_sol,kact_z);

}

Solution* SolRec::get_zero_solution(){
	return solution_zero;
}

//DualMesh::DualMesh(double *dxy, int power, double *XRange, double *YRange) {
//
//	if (power < 0)
//		power = 0;
//
//	dx = dxy[0] / pow(2.0, power);
//	dy = dxy[1] / pow(2.0, power);
//	//printf("dx=%g dy=%g  XRange={%g,%g} YRange={%g,%g}\n",dx,dy,XRange[0],XRange[1],YRange[0],YRange[1]);
//
//	xminmax[0] = XRange[0];
//	xminmax[1] = XRange[1];
//	yminmax[0] = YRange[0];
//	yminmax[1] = YRange[1];
//
//	Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
//	Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer
//
//	dualcell = new DualCell**[Ny];
//	for (int i = 0; i < Ny; ++i)
//		dualcell[i] = new DualCell*[Nx];
////	DualCell dualcell[Ny][Nx];
//
//	double position[2];
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j) {
//			ind_to_pos(i, j, position);
//			unsigned key[2] = { i, j };
//			dualcell[i][j] = new DualCell(key, position);
//		}
//
//	double zero_sol[3] = { 0., 0., 0. };
//
//	zero_solution = new Solution(zero_sol,0.);
//
//}
//
//void DualMesh::update_sol(Element* Curr_El, Solution* solution) {
//
//	double *coord = Curr_El->get_coord();
//	double *dxy = Curr_El->get_dx();
//
//	double xstart = coord[0] - 0.5 * dxy[0], xstop = coord[0] + 0.5 * dxy[0], ystart = coord[1]
//			- 0.5 * dxy[1], ystop = coord[1] + 0.5 * dxy[1];
//
//	int ixstart = (int) ((xstart - xminmax[0]) / dx + 0.5);
//	int ixstop = (int) ((xstop - xminmax[0]) / dx + 0.5);
//	int iystart = (int) ((ystart - yminmax[0]) / dy + 0.5);
//	int iystop = (int) ((ystop - yminmax[0]) / dy + 0.5);
//
//	if (ixstart < 0)
//		ixstart = 0;
//	if (ixstop == ixstart) {
//		ixstart = (int) ((xstart - xminmax[0]) / dx);
//		ixstop = ixstart + 1;
//	}
//	if (ixstop > Nx)
//		ixstop = Nx;
//
//	if (iystart < 0)
//		iystart = 0;
//	if (iystop == iystart) {
//		iystart = (int) ((ystart - yminmax[0]) / dy);
//		iystop = iystart + 1;
//	}
//	if (iystop > Ny)
//		iystop = Ny;
//
//	for (int iy = iystart; iy < iystop; ++iy)
//		for (int ix = ixstart; ix < ixstop; ++ix)
//			dualcell[iy][ix]->put_solution(solution);
//
//}
//
//void DualMesh::allocCellMem() {
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			dualcell[i][j]->allocMem();
//}
//
//int DualMesh::get_Nx() {
//	return Nx;
//}
//
//int DualMesh::get_Ny() {
//	return Ny;
//}
//
//void DualMesh::ind_to_pos(int indy, int indx, double* position) {
//
//	position[0] = xminmax[0] + dx * (indx + .5);
//	position[1] = yminmax[0] + dy * (indy + .5);
//
//}
//
//void DualMesh::pos_to_ind(double y, double x, int* ind) {
//
//	ind[0] = (int) ((y - yminmax[0]) / dy);
//	ind[1] = (int) ((x - xminmax[0]) / dx);
//
//}
//
//DualCell* DualMesh::get_dualcell(int i, int j) {
//
//	if (i < 0 || j < 0 || i > (Ny - 1) || j > (Nx - 1))
//		return NULL;
//
//	return dualcell[i][j];
//
//}
//
//double DualMesh::get_dx() {
//	return dx;
//}
//
//double DualMesh::get_dy() {
//	return dy;
//}
//
//void DualMesh::calc_flux() {
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->update_flux(this, 0);
//
//}
//
//void DualMesh::calc_slopes() {
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->calc_slopes(this);
//
//}
//
//void DualMesh::initialize_dual_flow(MatProps* matprops_ptr) {
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->calc_topo_data(this, matprops_ptr);
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->calc_gravity(matprops_ptr);
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->calc_d_gravity(this);
//
//	for (int i = 0; i < Ny; ++i)
//		for (int j = 0; j < Nx; ++j)
//			(this->get_dualcell(i, j))->calc_slopes(this);
//
//}
//
//const double* DualMesh::get_xminmax() const {
//	return xminmax;
//}
//
//const double* DualMesh::get_yminmax() const {
//	return yminmax;
//}
//
//void DualMesh::set_sol_zero(Element* Curr_El){
//
//	update_sol( Curr_El, zero_solution);
//}
//
//DualMesh::~DualMesh() {
//
//	for (int i = 0; i < Ny; ++i)
//		delete[] dualcell[i];
//
//	delete[] dualcell;
//}
