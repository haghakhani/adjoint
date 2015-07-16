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
 *******************************************************************
 */
//Jan 30, 2015
//haghakha
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define KEY0   3941335040
#define KEY1   0
#define ITER   9
#define EFFELL 0
#define J      0
#define JACIND 0

//#define DEBUG

void reset_resflag(ResFlag resflag[EFF_ELL]);

void calc_jacobian(MeshCTX* meshctx, PropCTX* propctx, PertElemInfo* eleminfo,
    double const increment) {

	int myid, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;

	int neigh_flag;
	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El = NULL;

	int iter = timeprops_ptr->iter;
	double tiny = GEOFLOW_TINY;

	//this array holds ResFlag for element itself and its neighbors
	ResFlag resflag[EFF_ELL];
	reset_resflag(resflag);

#ifdef DEBUGFILE
	ofstream myfile;
	myfile.open("debug.txt",ios::app);
#endif

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					int boundary = 0;

					//this part handles if the Curr_El is a boundary element
					//cout<<"I am running do not worry"<<endl;
					for (int neighnum = 0; neighnum < 4; neighnum++)
						if (*(Curr_El->get_neigh_proc() + neighnum) == INIT) {
							boundary = 1;
							// this for the element that are on the boundary
							Curr_El->set_jacobian();
							break;
						}
					if (!boundary) {

						double state_vars[NUM_STATE_VARS], gravity[NUM_STATE_VARS];
						double d_gravity[DIMENSION], curvature[DIMENSION];
						double d_state_vars_old[NUM_STATE_VARS * DIMENSION];
						double *prev_state_vars = Curr_El->get_prev_state_vars();
						Curr_El->calc_stop_crit(matprops_ptr); //this function updates bedfric properties
						double bedfrict = Curr_El->get_effect_bedfrict();
						double dx[DIMENSION] = { *(Curr_El->get_dx()), *(Curr_El->get_dx() + 1) };
						double kactxy[DIMENSION];
						double orgSrcSgn[DIMENSION];
						int check_stop[DIMENSION] = { 0, 0 };

						for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
							gravity[ind] = *(Curr_El->get_gravity() + ind);

						for (int ind = 0; ind < DIMENSION; ++ind) {

							d_gravity[ind] = *(Curr_El->get_d_gravity() + ind);
							curvature[ind] = *(Curr_El->get_curvature() + ind);

						}

						for (int ind = 0; ind < NUM_STATE_VARS * DIMENSION; ++ind)
							d_state_vars_old[ind] = *(Curr_El->get_d_state_vars() + ind);

						if (timeprops_ptr->iter < 51)
							matprops_ptr->frict_tiny = 0.1;
						else
							matprops_ptr->frict_tiny = 0.000000001;

						orgSourceSgn(Curr_El, matprops_ptr->frict_tiny, orgSrcSgn);

						double fluxold[4][NUM_STATE_VARS];
						record_flux(El_Table, NodeTable, Curr_El->pass_key(), matprops_ptr, myid, fluxold);

						double dt = timeprops_ptr->dt.at(iter - 1);	//at final time step we do not need the computation of adjoint and we always compute it for the previouse time so we need iter.
						double dtdx = dt / dx[0];
						double dtdy = dt / dx[1];

						double org_res[NUM_STATE_VARS];

						int gggflag = 0;

						if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1 && iter == ITER)
							gggflag = 1;

						residual(org_res, Curr_El->get_state_vars(), prev_state_vars, fluxold[0], //4
						    fluxold[1], fluxold[2], fluxold[3], dtdx, dtdy, dt, d_state_vars_old, //7
						    (d_state_vars_old + NUM_STATE_VARS), curvature, //2
						    matprops_ptr->intfrict, //1
						    bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()), //4
						    matprops_ptr->frict_tiny, orgSrcSgn, 0./*=increment*/, //3
						    matprops_ptr->epsilon, check_stop); //2

						for (int ind = 0; ind < NUM_STATE_VARS; ++ind)
							state_vars[ind] = *(Curr_El->get_state_vars() + ind);

						for (int effelement = 0; effelement < EFF_ELL; effelement++) { //0 for the element itself, and the rest id for neighbour elements

#ifdef DEBUG
							int gggflag = 0;

							if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
									&& iter == ITER && effelement == EFFELL)
							gggflag = 1;
#endif

							if (effelement == 0 && prev_state_vars[0] == 0.)

								//this is a void element so the residual vector does not change by changing it's values
								Curr_El->set_jacobianMat_zero(effelement);

							else if (effelement > 4
							    && compare_key((Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH),
							        (Curr_El->get_neighbors() + (effelement - 5) * KEYLENGTH)))

								Curr_El->set_jacobianMat_zero(effelement);

							else if (effelement != 0 && void_neigh_elem(El_Table, Curr_El, effelement))

								//this is a void neighbor element so the residual of the curr_el does not depend on this neighbor
								Curr_El->set_jacobianMat_zero(effelement);

							else {

								for (int j = 0; j < NUM_STATE_VARS; j++) { //there is a problem here: I do not need to compute for first the component of adjoint

									double vec_res[NUM_STATE_VARS];
									double total_res[NUM_STATE_VARS] = { 0., 0., 0. };

									int scheme = 0;
									for (; scheme < 2; scheme++) {

										//this flag shows that the pileheight before adding the increment is below or above the GEOFLOW_TINY.
										//if it is below GEOFLOW_TINY then there is no need to update fluxes and kactxy
										int updateflux, srcflag;
										reset_resflag(resflag);

										// here we modify increment to one time compute forward and one time compute backward difference if it is necessary
										double signe = pow(-1., scheme);
										double incr = signe * increment;
										increment_state(El_Table, Curr_El, incr, effelement, j, &updateflux, &srcflag,
										    resflag);

										calc_flux_slope_kact(El_Table, NodeTable, Curr_El, matprops_ptr, myid,
										    effelement, updateflux, srcflag, resflag);

										//Attention make sure that NUM_STATE_VARS are selected correctly
										//Actually we just need 3, but there is an excessive for first adjoint
										//const int state_num=NUM_STATE_VARS-2;
										double flux[4][NUM_STATE_VARS];
										record_flux(El_Table, NodeTable, Curr_El->pass_key(), matprops_ptr, myid, flux);

#ifdef DEBUG
										if (*(Curr_El->pass_key()) == KEY0 && *(Curr_El->pass_key() + 1) == KEY1
												&& effelement == EFFELL && iter == ITER && j == J)
										flux_debug(Curr_El, fluxold[0], fluxold[2], fluxold[1], fluxold[3], flux[0],
												flux[2], flux[1], flux[3], effelement, j, iter, dt);
#endif

										double *d_state_vars = Curr_El->get_d_state_vars();

										//here we compute the residuals
										residual(vec_res, state_vars, prev_state_vars, flux[0],									//4
										    flux[1], flux[2], flux[3], dtdx, dtdy, dt, d_state_vars,								//7
										    (d_state_vars + NUM_STATE_VARS), curvature,									//2
										    matprops_ptr->intfrict,									//1
										    bedfrict, gravity, d_gravity, *(Curr_El->get_kactxy()),									//4
										    matprops_ptr->frict_tiny, orgSrcSgn, incr,									//3
										    matprops_ptr->epsilon, check_stop, srcflag, 0);									//2

#ifdef DEBUGFILE
										myfile << "Elem Key[0]= " << *(Curr_El->pass_key()) << "  Key[1]= "
										<< *(Curr_El->pass_key()) << " iter= " << iter << " eff_el= " << effelement
										<< " j= " << j << " residual[0]= " << vec_res[0] << " residual[1]= "
										<< vec_res[1] << " residual[2]= " << vec_res[2] << endl;
#endif

										//we have to return everything back
//										restore(El_Table, NodeTable, Curr_El, effelement, j, incr, fluxold,
//										    d_state_vars_old);
										restore(El_Table, NodeTable, Curr_El, matprops_ptr, effelement, j, myid, incr);

										for (int ind = 0; ind < NUM_STATE_VARS; ind++)
											total_res[ind] += signe * vec_res[ind];

//											if (scheme == 0 && fabs(vec_res[0] / incr) < 5.
//													&& fabs(vec_res[1] / incr) < 5.
//													&& fabs(vec_res[2] / incr) < 5.)
//												//this means that forward difference is enough, and we do not need to compute central difference
										break;

									}

									double jacincr = increment; //= scheme == 0 ? increment : 2. * increment;

									Curr_El->set_jacobian(effelement, total_res, j,
									// following term is necessary to consider the scheme that whether it is forward difference or central difference
									    jacincr); //sets the propper components of the Jacobian for this element

								}
							}
						}
					}
				}
				currentPtr = currentPtr->next;
			}
		}
	}

#ifdef DEBUGFILE
	myfile.close();
#endif

}

int jac_mat_index(int effelement, int xp) {

	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
// effelement - 1 shows which neighbor we are processing
	if (effelement == 0)
		return 0;
	else if (effelement - 1 == xp)
		return 1;
	else if (effelement - 1 == yp)
		return 2;
	else if (effelement - 1 == xm)
		return 3;
	else if (effelement - 1 == ym)
		return 4;
	else
		exit(1);
}

//this function returns 1 if the neighbor element is void
//int void_neigh_cell(DualMesh* dualmesh, DualCell* cell, int effelement) {
//
//	int a, b;
//	set_ab(&a, &b, effelement);
//
//	int i = *(cell->get_key()), j = *(cell->get_key() + 1);
//
//	DualCell* neigh_cell = dualmesh->get_dualcell(i + a, j + b); //basically we are checking all neighbor elements, and start from xp neighbor
//
//	assert(neigh_cell);
//
//	if (*(neigh_cell->get_prev_state_vars()) == 0.)
//		return 1;
//
//	return 0;
//}

//this function returns 1 if the neighbor element is void
int void_neigh_elem(HashTable* El_Table, Element* Curr_El, int effelement) {

	Element* neigh_elem = (Element*) (El_Table->lookup(
	    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));

	assert(neigh_elem);

	if (*(neigh_elem->get_prev_state_vars()) == 0.)
		return 1;

	return 0;
}

void restore(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El, MatProps* matprops_ptr,
    int effelement, int j, int myid, double increment) {

	Element* neigh_elem;
	double *prev_state_vars = Curr_El->get_prev_state_vars();

	double dummydt = 0., outflow = 0.;
	int order_flag = 1;

	ResFlag resflag;
	resflag.callflag = 1;
	resflag.lgft = 0;

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	if (effelement == 0) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		prev_state_vars[j] -= increment; //changing the state varibales at the element itself
		Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
		    &outflow, resflag, resflag); //chenge of the state_vars causes the all around fluxes change, this update xp ,yp

		if ((*(Curr_El->get_neigh_proc() + xm)) != INIT) { //we have to make sure that there exit an element in xm side

			Element* elem_xm = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + xm * KEYLENGTH));
			elem_xm->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag, resflag); //this update the flux on share edge with xm
		}
		if ((*(Curr_El->get_neigh_proc() + ym)) != INIT) { //we have to make sure that there exit an element in ym side

			Element* elem_ym = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + ym * KEYLENGTH));
			elem_ym->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag, resflag); //this update the flux on share edge with ym
		}

	} else if ((*(Curr_El->get_neigh_proc() + (effelement - 1))) != INIT) {

		neigh_elem = (Element*) (El_Table->lookup(
		    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
		*(neigh_elem->get_prev_state_vars() + j) -= increment;

		if ((effelement - 1) == xp || (effelement - 1) == yp)
			Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag, resflag); //if we change the state variables in xp or yp, just the flux at this element has to be updated
		else
			neigh_elem->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag, resflag); //otherwise the flux at the corresponding element has to be updated

	}
	Curr_El->get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
	return;
}

void restore(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El, int effelement, int j,
    double increment, double fluxold[4][NUM_STATE_VARS],
    double d_state_vars_old[DIMENSION * NUM_STATE_VARS]) {

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	if (effelement == 0) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		double *prev_state_vars = Curr_El->get_prev_state_vars();
		//changing the state varibales at the element itself
		prev_state_vars[j] -= increment;

		// in begining of calc_jacobian we made sure that the Curr_El is not on the boundary
		Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

		Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

		Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

		Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

		for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
			nxp->flux[ivar] = fluxold[0][ivar];
			nyp->flux[ivar] = fluxold[1][ivar];
			nxm->flux[ivar] = fluxold[2][ivar];
			nym->flux[ivar] = fluxold[3][ivar];
		}

	} else if ((*(Curr_El->get_neigh_proc() + (effelement - 1))) != INIT) {

		Element* neigh_elem = (Element*) (El_Table->lookup(
		    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
		*(neigh_elem->get_prev_state_vars() + j) -= increment;

		if ((effelement - 1) == xp || (effelement - 1) == yp) {

			Node* nxp = (Node*) NodeTable->lookup(Curr_El->getNode() + (xp + 4) * 2);

			Node* nyp = (Node*) NodeTable->lookup(Curr_El->getNode() + (yp + 4) * 2);

			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
				nxp->flux[ivar] = fluxold[0][ivar];
				nyp->flux[ivar] = fluxold[1][ivar];
			}

		} else {

			Node* nxm = (Node*) NodeTable->lookup(Curr_El->getNode() + (xm + 4) * 2);

			Node* nym = (Node*) NodeTable->lookup(Curr_El->getNode() + (ym + 4) * 2);

			for (int ivar = 0; ivar < NUM_STATE_VARS; ivar++) {
				nxm->flux[ivar] = fluxold[2][ivar];
				nym->flux[ivar] = fluxold[3][ivar];
			}
		}
	}
	double* d_state_vars = Curr_El->get_d_state_vars();

	for (int i = 0; i < NUM_STATE_VARS * DIMENSION; ++i)
		d_state_vars[i] = d_state_vars_old[i];

}

//void calc_flux_slope_kact(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El,
//    MatProps* matprops_ptr, int myid, int effelement, int updateflux, int srcflag,
//    ResFlag resflag[EFF_ELL]) {
//
//	ResFlag dummyresflag;
//	dummyresflag.callflag = 1;
//	dummyresflag.lgft = 0;
//
//	Element* neigh_elem;
//
//	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
//	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
//
//	Curr_El->get_slopes_prev(El_Table, NodeTable, matprops_ptr->gamma);	//we also have to update the d_state_vars for the current element
//
////	double *d_state_vars = Curr_El->get_d_state_vars();
////										if (srcflag && effelement == 0) {
////
////											gmfggetcoef_(Curr_El->get_prev_state_vars(), d_state_vars,
////													(d_state_vars + NUM_STATE_VARS), dx,
////													&(matprops_ptr->bedfrict[Curr_El->get_material()]),
////													&(matprops_ptr->intfrict), &kactxy[0], &kactxy[1],
////													&tiny, &(matprops_ptr->epsilon));
////
////											Curr_El->put_kactxy(kactxy);
////											Curr_El->calc_stop_crit(matprops_ptr);
////										}
//
//	if (effelement == 0 && updateflux) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars
//
//		Curr_El->calc_fluxes(El_Table, NodeTable, myid, resflag[0], dummyresflag);
//		//change of the state_vars causes the all around fluxes change, this update xp ,yp
//		// earlier in this file, we made sure that this element is not a boundary element
//		// so here we do not require to check the neighbor elements to make sure they are not located on the boundary
//
//		Element* elem_xm = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + xm * KEYLENGTH));
//		assert(elem_xm);
//		elem_xm->calc_xflux(El_Table, NodeTable, myid, resflag[xm + 1], resflag[0]);
//
//		Element* elem_ym = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + ym * KEYLENGTH));
//		assert(elem_ym);
//		elem_ym->calc_yflux(El_Table, NodeTable, myid, resflag[ym + 1], resflag[0]);
//
//	} else if (effelement != 0 && updateflux) {
//
//		if ((effelement - 1) % 4 == xp)
//			Curr_El->calc_xflux(El_Table, NodeTable, myid, resflag[0], resflag[xp + 1]);
//		//if we change the state variables in xp or yp, just the flux at this element has to be updated
//
//		else if ((effelement - 1) % 4 == yp)
//			Curr_El->calc_yflux(El_Table, NodeTable, myid, resflag[0], resflag[ym + 1]);
//
//		else if ((effelement - 1) % 4 == xm) {
//
//			neigh_elem = (Element*) (El_Table->lookup(
//			    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
//			assert(neigh_elem);
//
//			//otherwise the flux at the corresponding element has to be updated
//			neigh_elem->calc_xflux(El_Table, NodeTable, myid, resflag[xm + 1], resflag[0]);
//
//		} else {
//
//			neigh_elem = (Element*) (El_Table->lookup(
//			    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
//			assert(neigh_elem);
//			//otherwise the flux at the corresponding element has to be updated
//			neigh_elem->calc_yflux(El_Table, NodeTable, myid, resflag[ym + 1], resflag[0]);
//
//		}
//
//	}
//
//}

void increment_state(HashTable* El_Table, Element* Curr_El, double increment, int effelement, int j,
    int* updateflux, int* srcflag, ResFlag resflag[EFF_ELL]) {

	*updateflux = 1;
	*srcflag = 1;
	double *prev_state_vars = Curr_El->get_prev_state_vars();

	if (effelement == 0) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		if (j == 0 && prev_state_vars[j] < GEOFLOW_TINY) {
			*updateflux = 0;
			*srcflag = 0;
			resflag[effelement].lgft = 1;
		}

		prev_state_vars[j] += increment; //changing the state varibales at the element itself

	} else {

		Element* neigh_elem = (Element*) (El_Table->lookup(
		    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
		assert(neigh_elem);

		if (j == 0 && *(neigh_elem->get_prev_state_vars() + j) < GEOFLOW_TINY) {
			*updateflux = 0;
			resflag[effelement].lgft = 1;
		}

		*(neigh_elem->get_prev_state_vars() + j) += increment;

	}
	return;
}

void reset_resflag(ResFlag resflag[EFF_ELL]) {

	for (int i = 0; i < EFF_ELL; i++) {
		resflag[i].callflag = 1;
		resflag[i].lgft = 0;
	}

	return;
}

void calc_flux_slope_kact(HashTable* El_Table, HashTable* NodeTable, Element* Curr_El,
    MatProps* matprops_ptr, int myid, int effelement, int updateflux, int srcflag,
    ResFlag resflag[5]) {

	double dummydt = 0., outflow = 0.;
	int order_flag = 1;

	ResFlag dummyresflag;
	dummyresflag.callflag = 1;
	dummyresflag.lgft = 0;

	Element* neigh_elem;

	int xp = Curr_El->get_positive_x_side(); //finding the direction of element
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	Curr_El->get_slopes_prev(El_Table, NodeTable, matprops_ptr->gamma);	//we also have to update the d_state_vars for the current element

	double *d_state_vars = Curr_El->get_d_state_vars();

//										if (srcflag && effelement == 0) {
//
//											gmfggetcoef_(Curr_El->get_prev_state_vars(), d_state_vars,
//													(d_state_vars + NUM_STATE_VARS), dx,
//													&(matprops_ptr->bedfrict[Curr_El->get_material()]),
//													&(matprops_ptr->intfrict), &kactxy[0], &kactxy[1],
//													&tiny, &(matprops_ptr->epsilon));
//
//											Curr_El->put_kactxy(kactxy);
//											Curr_El->calc_stop_crit(matprops_ptr);
//										}

	if (effelement == 0 && updateflux) { //this part of code add an increment to the state variables to find the Jacobian, but the problem is since it is called after correct, the increment shoud be added to the prev_state_vars

		Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
		    &outflow, resflag[0], dummyresflag); //change of the state_vars causes the all around fluxes change, this update xp ,yp
		// earlier in this file, we made sure that this element is not a boundary element
		// so here we do not require to check the neighbor elements to make sure they are not located on the boundary

		Element* elem_xm = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + xm * KEYLENGTH));
		assert(elem_xm);
		elem_xm->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
		    &outflow, resflag[3], resflag[0]); //this update the flux on share edge with xm

		Element* elem_ym = (Element*) (El_Table->lookup(Curr_El->get_neighbors() + ym * KEYLENGTH));
		assert(elem_ym);
		elem_ym->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
		    &outflow, resflag[4], resflag[0]); //this update the flux on share edge with ym

	} else if (effelement != 0 && updateflux) {

		if ((effelement - 1) == xp)

			Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag[0], resflag[1]); //if we change the state variables in xp or yp, just the flux at this element has to be updated
		else if ((effelement - 1) == yp)
			Curr_El->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag[0], resflag[2]);

		else if ((effelement - 1) == xm) {

			neigh_elem = (Element*) (El_Table->lookup(
			    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
			assert(neigh_elem);

			neigh_elem->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag[3], resflag[0]); //otherwise the flux at the corresponding element has to be updated
//otherwise the flux at the corresponding element has to be updated

		} else {

			neigh_elem = (Element*) (El_Table->lookup(
			    Curr_El->get_neighbors() + (effelement - 1) * KEYLENGTH));
			assert(neigh_elem);

			neigh_elem->calc_edge_states(El_Table, NodeTable, matprops_ptr, myid, dummydt, &order_flag,
			    &outflow, resflag[4], resflag[0]); //otherwise the flux at the corresponding element has to be updated
//otherwise the flux at the corresponding element has to be updated
		}

	}

	return;

}

