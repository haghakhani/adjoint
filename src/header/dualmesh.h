/*
 * dualmesh.h
 *
 *  Created on: May 25, 2015
 *      Author: haghakha
 */

#ifndef DUALMESH_H
#define DUALMESH_H

#include "element2.h"
#include "jacobian.h"

class SolRec: public HashTable {
private:
	Solution* solution_zero;

public:
	SolRec(double *doublekeyrangein, int size, int prime, double  XR[], double YR[], int ifrestart);
	Solution* get_zero_solution();

//	wrtie_sol_to_disk

};

//private:
//	//! number of cells in the x direction on the map
//	int Nx;
//
//	//! number of cells in the y direction on the map
//	int Ny;
//
//	//! length of a cell in the x direction
//	double dx;
//
//	//! length of a cell in the y direction
//	double dy;
//
//	//! min and max x coordinate on the map
//	double xminmax[2];
//
//	//! min and max y coordinate on the map
//	double yminmax[2];
//
//	//!
//	DualCell ***dualcell;
//
//	Solution* zero_solution;
//
//public:
//
//	DualMesh(double *dxy, int power, double *XRange, double *YRange);
//
//	void update_sol(Element* Curr_El, Solution* solution);
//
//	void allocCellMem();
//
//	int get_Nx();
//
//	int get_Ny();
//
//	double get_dx();
//
//	double get_dy();
//
//	DualCell* get_dualcell(int i, int j);
//
//	void pos_to_ind(double y, double x, int* ind);
//
//	void ind_to_pos(int indy, int indx, double* position);
//
//	void calc_flux();
//
//	void calc_slopes();
//
//	void initialize_dual_flow(MatProps* matprops_ptr);
//
//	const double* get_xminmax() const;
//
//	const double* get_yminmax() const;
//
//	void calc_func_sens();
//
//	void set_sol_zero(Element* Curr_El);
//
//
//	~DualMesh();
//};

#endif /* DUALMESH_H */
