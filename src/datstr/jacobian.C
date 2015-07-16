/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: jacobian.C  2014-04-10 10:33:10 haghakha $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

Solution::Solution(double* curr_sol, double kactxy) {

	for (int i = 0; i < NUM_STATE_VARS; ++i)
		states[i] = curr_sol[i];

	kact = kactxy;

}

double* Solution::get_solution() {
	return states;
}

double Solution::get_kact() {
	return kact;
}

Solution::~Solution() {

}

Jacobian::Jacobian(unsigned* key, double* position) {
	for (int i = 0; i < KEYLENGTH; ++i) {
		Jacobian::key[i] = key[i];
		Jacobian::position[i] = position[i];
	}

}

double* Jacobian::get_position() {
	return position;
}

void Jacobian::put_solution(Solution* solution, int iter) {

	//because insert is faster than emplace
	solContainer[iter] = solution;

	return;
}

Solution* Jacobian::get_solution(int iter) {

	return solContainer[iter];
}

unsigned* Jacobian::get_key() {
	return key;
}

Jacobian::~Jacobian() {

	solContainer.clear();

}

//DualCell::DualCell(unsigned* key, double* position) :
//		Jacobian(key, position) {
//
//	curr_adjoint = NULL;
//	prev_adjoint = NULL;
//	func_sens = NULL;
//	state_vars = NULL;
//	prev_state_vars = NULL;
//	flux = NULL;
//	gravity = NULL;
//	d_gravity = NULL;
//	d_state_vars = NULL;
//	zeta = NULL;
//	curvature = NULL;
//
//	kact = 0.;
//	elevation = 0.;
//
//}
//
//double* DualCell::get_state_vars() {
//	return state_vars;
//}
//
//double* DualCell::get_prev_state_vars() {
//	return prev_state_vars;
//}
//
//double* DualCell::get_curr_adjoint() {
//	return curr_adjoint;
//}
//
//double* DualCell::get_prev_adjoint() {
//	return prev_adjoint;
//}
//
//void DualCell::rev_state_vars(int iter) {
//
//	for (int i = 0; i < NUM_STATE_VARS; i++)
//		state_vars[i] = *((solvector.at(iter))->get_solution() + i);
//
//	for (int i = 0; i < NUM_STATE_VARS; i++)
//		prev_state_vars[i] = *((solvector.at(iter - 1))->get_solution() + i);
//
//	kact = (solvector.at(iter))->get_kact();
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		prev_adjoint[i] = curr_adjoint[i];
//
//}
//
//double DualCell::get_kact() {
//	return kact;
//}
//
//void DualCell::xdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft) {
//
//	//sfs is abreviation of states, fluxes and speeds
//
//	double a = 0., Vel = 0.;
//
//	if (prev_state_vars[0] < GEOFLOW_TINY || lgft)
//
//		for (int i = 0; i < NUM_STATE_VARS; ++i)
//			for (int j = 0; j < NUM_STATE_VARS; ++j)
//				sfs[i][j] = 0.0;
//
//	else {
//
//		for (int i = 0; i < NUM_STATE_VARS; i++)
//			sfs[0][i] = prev_state_vars[i];
//
//		//velocity in x-dir
//		Vel = sfs[0][1] / sfs[0][0];
//
//		// sound-speed : a^2 = k_ap*h*g(3)
//		a = sqrt(kact * sfs[0][0] * gravity[2]);
//
//		//fluxes
//		sfs[1][0] = sfs[0][1];
//		sfs[1][1] = sfs[0][1] * Vel + 0.5 * a * a * sfs[0][0];
//		sfs[1][2] = sfs[0][2] * Vel;
//
//		//wave speeds
//		sfs[2][0] = Vel - a;
//		sfs[2][1] = Vel;
//		sfs[2][2] = Vel + a;
//
//	}
//
//	for (int i = 0; i < NUM_STATE_VARS; i++)
//		for (int j = 0; j < NUM_STATE_VARS; j++)
//			if (isnan(sfs[i][j]))
//				cout << "flux is NAN" << endl;
//
//}
//
//void DualCell::ydir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft) {
//
//	//sfs is abreviation of states, fluxes and speeds
//
//	double a = 0., Vel = 0.;
//
//	if (prev_state_vars[0] < GEOFLOW_TINY || lgft)
//
//		for (int i = 0; i < NUM_STATE_VARS; ++i)
//			for (int j = 0; j < NUM_STATE_VARS; ++j)
//				sfs[i][j] = 0.0;
//
//	else {
//
//		for (int i = 0; i < NUM_STATE_VARS; i++)
//			sfs[0][i] = prev_state_vars[i];
//
//		//velocity in y-dir
//		Vel = sfs[0][2] / sfs[0][0];
//
//		// sound-speed : a^2 = k_ap*h*g(3)
//		a = sqrt(kact * sfs[0][0] * gravity[2]);
//
//		//fluxes
//		sfs[1][0] = sfs[0][2];
//		sfs[1][1] = sfs[0][1] * Vel;
//		sfs[1][2] = sfs[0][2] * Vel + 0.5 * a * a * sfs[0][0];
//
//		//wave speeds
//		sfs[2][0] = Vel - a;
//		sfs[2][1] = Vel;
//		sfs[2][2] = Vel + a;
//
//	}
//
//	for (int i = 0; i < NUM_STATE_VARS; i++)
//		for (int j = 0; j < NUM_STATE_VARS; j++)
//			if (isnan(sfs[i][j]))
//				cout << "flux is NAN" << endl;
//
//}
//
//void DualCell::zdir_flux(double sfs[NUM_STATE_VARS][NUM_STATE_VARS], int lgft, int side) {
//
//	if (side % 2 == 0)
//		xdir_flux(sfs, lgft);
//	else
//		ydir_flux(sfs, lgft);
//}
//
//void DualCell::update_flux_x(DualMesh* dualmesh, int* lgft) {
//
//	update_flux(dualmesh, lgft, 0);
//
//}
//
//void DualCell::update_flux_y(DualMesh* dualmesh, int* lgft) {
//
//	update_flux(dualmesh, lgft, 1);
//
//}
//
//void DualCell::update_flux(DualMesh* dualmesh, int* lgft, int side) {
//
//	int Nx = dualmesh->get_Nx();
//	int Ny = dualmesh->get_Ny();
//
//	double sfsl[NUM_STATE_VARS][NUM_STATE_VARS], sfsr[NUM_STATE_VARS][NUM_STATE_VARS];
//
//	DualCell* right_cell;
//	int dummy_lgft[2] = { 0, 0 };
//
//	if (!lgft)
//		lgft = dummy_lgft;
//
//	if (key[1] == Nx - 1 && side == 0) {
//
//		DualCell* left_cell = dualmesh->get_dualcell(key[0], key[1] - 1);
//		for (int i = 0; i < 3; ++i)
//			flux[side * NUM_STATE_VARS + i] = *(left_cell->get_xflux() + i);
//
//	} else if (key[0] == Ny - 1 && side == 1) {
//
//		DualCell* left_cell = dualmesh->get_dualcell(key[0] - 1, key[1]);
//		for (int i = 0; i < 3; ++i)
//			flux[side * NUM_STATE_VARS + i] = *(left_cell->get_yflux() + i);
//
//	} else {
//
//		zdir_flux(sfsl, lgft[0], side);
//
//		if (side == 0)
//			right_cell = dualmesh->get_dualcell(key[0], key[1] + 1);
//		else
//			right_cell = dualmesh->get_dualcell(key[0] + 1, key[1]);
//
//		right_cell->zdir_flux(sfsr, lgft[1], side);
//
//		double* flux;
//
//		if (side == 0)
//			flux = get_xflux();
//		else
//			flux = get_yflux();
//
//		riemannflux_cell(sfsl, sfsr, flux);
//
//	}
//
//}
//
//void DualCell::update_flux(DualMesh* dualmesh, int* lgft) {
//
//	for (int side = 0; side < 2; ++side)
//
//		update_flux(dualmesh, lgft, side);
//
//}
//
//void DualCell::calc_gravity(MatProps* matprops_ptr) {
//
//	double max_slope = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1]);
//	double max_angle = atan(max_slope);
//
//	double down_slope_gravity = 9.8 * sin(max_angle);
//	if (dabs(down_slope_gravity) > GEOFLOW_TINY) {
//		gravity[0] = -down_slope_gravity * zeta[0] / max_slope;
//		gravity[1] = -down_slope_gravity * zeta[1] / max_slope;
//		gravity[2] = 9.8 * cos(max_angle);
//	} else {
//		gravity[0] = 0;
//		gravity[1] = 0;
//		gravity[2] = 9.8;
//	}
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		gravity[i] = gravity[i] / matprops_ptr->GRAVITY_SCALE;
//
//}
//
//void DualCell::calc_d_gravity(DualMesh* dualmesh) {
//
//	/* x direction */
//	DualCell* cellp = dualmesh->get_dualcell(key[0], key[1] + 1);
//	DualCell* cellm = dualmesh->get_dualcell(key[0], key[1] - 1);
//	double dx = dualmesh->get_dx();
//
//	if (cellp != NULL && cellm != NULL) {
//
//		double dp, dm;
//		dp = (cellp->gravity[2] - gravity[2]) / dx;
//		dm = (gravity[2] - cellm->gravity[2]) / dx;
//
//		d_gravity[0] = .5 * (dp + dm);  // weighted average
//
//	} else if (cellm != NULL) {
//
//		d_gravity[0] = (gravity[2] - cellm->gravity[2]) / dx;
//
//	} else if (cellp != NULL) {
//
//		d_gravity[0] = (cellp->gravity[2] - gravity[2]) / dx;
//
//	} else
//		//no neighbors on either side -- assume that the ground is flat
//		d_gravity[0] = 0;
//
//	/* y direction */
//	cellp = dualmesh->get_dualcell(key[0] + 1, key[1]);
//	cellm = dualmesh->get_dualcell(key[0] - 1, key[1]);
//
//	if (cellp != NULL && cellm != NULL) {
//
//		double dp, dm;
//
//		dp = (cellp->gravity[2] - gravity[2]) / dx;
//		dm = (gravity[2] - cellm->gravity[2]) / dx;
//
//		d_gravity[1] = .5 * (dp + dm);  // weighted average
//
//	} else if (cellm != NULL) {
//
//		d_gravity[1] = (gravity[2] - cellm->gravity[2]) / dx;
//
//	} else if (cellp != NULL) {
//
//		d_gravity[1] = (cellp->gravity[2] - gravity[2]) / dx;
//
//	} else
//		//no neighbors on either side -- assume that the ground is flat
//		d_gravity[1] = 0;
//
//}
//
//void DualCell::calc_slopes(DualMesh* dualmesh) {
//
//	int Ny = dualmesh->get_Ny();
//	int Nx = dualmesh->get_Nx();
//	double dx = dualmesh->get_dx();
//	double dy = dualmesh->get_dy();
//
//	if (key[0] == 0 || key[0] == Ny - 1 || key[1] == 0 || key[1] == Nx - 1) {
//		for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
//			d_state_vars[i] = 0;
//		return;
//	} else {
//		/* x direction */
//		DualCell* cellp = dualmesh->get_dualcell(key[0], key[1] + 1);
//		DualCell* cellm = dualmesh->get_dualcell(key[0], key[1] - 1);
//		for (int i = 0; i < NUM_STATE_VARS; ++i) {
//
//			double dp = (*(cellp->get_prev_state_vars() + i) - prev_state_vars[i]) / dx;
//
//			double dm = (prev_state_vars[i] - *(cellm->get_prev_state_vars() + i)) / dx;
//
//			double dc = .5 * (dp + dm);
//
//			d_state_vars[i] = .5 * (c_sgn(dp) + c_sgn(dm)) * c_dmin1(dabs(dp), dabs(dm), dabs(dc));
//		}
//
//		/* y direction */
//		cellp = dualmesh->get_dualcell(key[0] + 1, key[1]);
//		cellm = dualmesh->get_dualcell(key[0] - 1, key[1]);
//
//		for (int i = 0; i < NUM_STATE_VARS; ++i) {
//
//			double dp = (*(cellp->get_prev_state_vars() + i) - prev_state_vars[i]) / dy;
//
//			double dm = (prev_state_vars[i] - *(cellm->get_prev_state_vars() + i)) / dy;
//
//			double dc = .5 * (dp + dm);
//
//			d_state_vars[i + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm))
//			    * c_dmin1(dabs(dp), dabs(dm), dabs(dc));
//
//		}
//
//	}
//
//}
//
//void DualCell::allocMem() //in forward run we just save the solution and in backward run we compute the jacobian
//{
//	new_jacobianMat();
//
//	curr_adjoint = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		curr_adjoint[i] = 0.;
//
//	prev_adjoint = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		prev_adjoint[i] = 0.;
//
//	func_sens = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		func_sens[i] = 0.;
//
//	state_vars = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		state_vars[i] = 0.;
//
//	prev_state_vars = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		prev_state_vars[i] = 0.;
//
//	flux = new double[DIMENSION * NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		flux[i] = 0.;
//
//	gravity = new double[NUM_STATE_VARS];
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		gravity[i] = 0.;
//
//	d_gravity = new double[DIMENSION];
//	for (int i = 0; i < DIMENSION; ++i)
//		d_gravity[i] = 0.;
//
//	d_state_vars = new double[DIMENSION * NUM_STATE_VARS];
//	for (int i = 0; i < DIMENSION * NUM_STATE_VARS; ++i)
//		d_state_vars[i] = 0.;
//
//	zeta = new double[DIMENSION];
//	for (int i = 0; i < DIMENSION; ++i)
//		zeta[i] = 0.;
//
//	curvature = new double[DIMENSION];
//	for (int i = 0; i < DIMENSION; ++i)
//		curvature[i] = 0.;
//
//	return;
//}
//
//double* DualCell::get_curvature() {
//	return curvature;
//}
//
//double* DualCell::get_d_gravity() {
//	return d_gravity;
//}
//
//double* DualCell::get_d_state_vars() {
//	return d_state_vars;
//}
//
//double DualCell::get_elevation() {
//	return elevation;
//}
//
//double* DualCell::get_gravity() {
//	return gravity;
//}
//
//double* DualCell::get_xflux() {
//	return flux;
//}
//
//double* DualCell::get_yflux() {
//	return (flux + NUM_STATE_VARS);
//}
//
//double* DualCell::get_zeta() {
//	return zeta;
//}
//
//void DualCell::calc_topo_data(DualMesh* dualmesh, MatProps* matprops_ptr) {
//
//	double dx = dualmesh->get_dx();
//	double dy = dualmesh->get_dy();
//
//	double resolution = (dx + dy) * (matprops_ptr->LENGTH_SCALE) / 2.0;
//	double xcoord = position[0] * (matprops_ptr->LENGTH_SCALE);
//	double ycoord = position[1] * (matprops_ptr->LENGTH_SCALE);
//	int i = Get_elevation(resolution, xcoord, ycoord, &elevation);
//	elevation = elevation / matprops_ptr->LENGTH_SCALE;
//	i = Get_slope(resolution, xcoord, ycoord, zeta, (zeta + 1));
//	i = Get_curvature(resolution, xcoord, ycoord, curvature, (curvature + 1));
//	curvature[0] = curvature[0] * (matprops_ptr->LENGTH_SCALE);
//	curvature[1] = curvature[1] * (matprops_ptr->LENGTH_SCALE);
//
//}
//
//double* DualCell::get_funcsens() {
//	return func_sens;
//}
//
//void DualCell::print_cell_info(int iter) {
//
//	FILE *fp;
//	char filename[256];
//	sprintf(filename, "cell_NY%3d_NX%3d_iter%04d", key[0], key[1], iter);
//	fp = fopen(filename, "a");
//
//	fprintf(fp, "cell position: x=%f , y=%f \n", position[0], position[1]);
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "state_vars[%d]:      %e ", i, state_vars[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "prev_state_vars[%d]: %e ", i, prev_state_vars[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "func_sens[%d]:       %e ", i, func_sens[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "curr_adjoint[%d]:    %e ", i, curr_adjoint[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "prev_adjoint[%d]:    %e ", i, prev_adjoint[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS; ++i)
//		fprintf(fp, "gravity[%d]:         %e ", i, gravity[i]);
//	fprintf(fp, "\n");
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
//		fprintf(fp, "flux[%d]: %e ", i, flux[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
//		fprintf(fp, "d_state_vars[%d]: %e ", i, d_state_vars[i]);
//	fprintf(fp, "\n");
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < DIMENSION; ++i)
//		fprintf(fp, "d_gravity[%d]: %e ", i, d_gravity[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < DIMENSION; ++i)
//		fprintf(fp, "zeta[%d]:      %e ", i, zeta[i]);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < DIMENSION; ++i)
//		fprintf(fp, "curvature[%d]: %e ", i, curvature[i]);
//	fprintf(fp, "\n");
//	fprintf(fp, "\n");
//
//	fprintf(fp, "elevation= %e \n", elevation);
//	fprintf(fp, "kact=      %e \n", kact);
//	fprintf(fp, "\n");
//
//	for (int i = 0; i < 5; i++) {
//		fprintf(fp, "Matrix=  %d,\n", i);
//
//		for (int j = 0; j < NUM_STATE_VARS; j++) {
//			for (int k = 0; k < NUM_STATE_VARS; k++) {
//				fprintf(fp, "%10.8f  ", jacobianMat[i][j][k]);
//				if (dabs(jacobianMat[i][j][k]) > 10.)
//					fprintf(fp, "Jedi begir mano \n");
//			}
//			fprintf(fp, "\n");
//		}
//	}
//
//	fprintf(fp, "====================================\n");
//
//	fclose(fp);
//
//}
//
//void DualCell::print_cell_neighb_info(DualMesh*dualmesh, int iter) {
//
//	for (int effel = 0; effel < 5; ++effel)
//		if (effel == 0)
//			print_cell_info(iter);
//
//		else {
//
//			int a, b;
//			set_ab(&a, &b, effel);
//
//			DualCell* neigh_cell = dualmesh->get_dualcell(key[0] + a, key[1] + b);
//			if (neigh_cell != NULL)
//				neigh_cell->print_cell_info(iter);
//			else
//				cout << "for effel  " << effel << " is NULL " << endl;
//
//		}
//
//}
//
//DualCell::~DualCell() {
//
//	delete[] curr_adjoint;
//	delete[] prev_adjoint;
//	delete[] func_sens;
//	delete[] state_vars;
//	delete[] prev_state_vars;
//	delete[] flux;
//	delete[] gravity;
//	delete[] d_gravity;
//	delete[] d_state_vars;
//	delete[] zeta;
//	delete[] curvature;
//
//}
//
//void riemannflux_cell(double sfsl[NUM_STATE_VARS][NUM_STATE_VARS],
//    double sfsr[NUM_STATE_VARS][NUM_STATE_VARS], double* flux) {
//
//	if ((sfsl[0][0] == 0.) && (sfsr[0][0] == 0.))
//		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
//			flux[ivar] = 0.;
//	else {
//
//		double sl, sr;
//		if (sfsl[0][0] == 0.) {
//			sl = min(0., 2.0 * sfsr[2][0] - sfsr[2][1]);
//			sr = max(0., 2.0 * sfsr[2][2] - sfsr[2][1]);
//		} else if (sfsr[0][0] == 0.) {
//			sl = min(0., 2.0 * sfsl[2][0] - sfsl[2][1]);
//			sr = max(0., 2.0 * sfsl[2][2] - sfsl[2][1]);
//		} else {
//			sl = min(0., min(sfsl[2][0], sfsr[2][0]));
//			sr = max(0., max(sfsl[2][2], sfsr[2][2]));
//		}
//
//		if (sl >= 0.0)
//			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
//				flux[ivar] = sfsl[1][ivar];
//
//		else if (sr <= 0.0)
//			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
//				flux[ivar] = sfsr[1][ivar];
//
//		else
//			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
//				flux[ivar] = (sr * sfsl[1][ivar] - sl * sfsr[1][ivar]
//				    + sl * sr * (sfsr[0][ivar] - sfsl[0][ivar])) / (sr - sl);
//	}
//}
//
