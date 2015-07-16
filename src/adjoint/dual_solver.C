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
 * $Id: calc_jacobian.C 164 2013-02-27 15:27:22Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "../header/hpfem.h"

#define DEBUG1

#define KEY0   3920807148
#define KEY1   1321528399
#define ITER   10
#define ITER    187
#define J       0

void dual_solver(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx, PertElemInfo* eleminfo) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	const int rescomp = 1;
	const double increment = INCREMENT;
	const int maxiter = timeprops_ptr->iter;

	unsigned keyy[2] = { 3796806314, 2863311530 };
//	if (checkElement(El_Table, NULL, keyy))
//		cout << "I found the suspecious element \n";

	allocJacoMat(El_Table);

	reset_adaption_flag(El_Table);

	double functional = 0.0, dt;

	cout << "computing ADJOINT time step " << maxiter << endl;

	calc_adjoint(meshctx, propctx);

//	uinform_refine(meshctx, propctx, numprocs, myid);
//
//	error_compute(meshctx, propctx, maxiter, myid, numprocs);
//
//	double UNREFINE_TARGET = .01;	//dummy value is not used in the function
//	unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
//	    rescomp);

	int tecflag = 2;
	meshplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

	tecflag = 1;
	for (int iter = maxiter; iter > 0; --iter) {

		timeprops_ptr->iter = iter;
		cout << "computing ADJOINT time step " << iter - 1 << endl;
		dt = timeprops_ptr->dt.at(iter - 1);
		timeprops_ptr->adjiter++;

		// we need this even for  iter = maxiter because after refine and unrefine
		// the state variables are not same as forward run
		// reverse_states(El_Table, solrec, iter, &refinelist, &unrefinelist);

		setup_dual_flow(solrec, meshctx, propctx, iter);

//		print_Elem_Table(El_Table,timeprops_ptr->iter,1);

//		cout << "num_elem= " << num_nonzero_elem(El_Table) << "\n";

		allocJacoMat(El_Table);

		timeprops_ptr->adjoint_time(iter - 1);

		//this array holds ResFlag for element itself and its neighbors
		ResFlag resflag;
		resflag.callflag = 1;
		resflag.lgft = 0;

		calc_flux(meshctx, propctx, myid, resflag);

//		if (checkElement(El_Table, NodeTable, NULL, keyy))
//			cout << "I found the suspecious element \n";

		slopes(El_Table, NodeTable, matprops_ptr, 1);

		compute_functional(El_Table, &functional, timeprops_ptr);

		eleminfo->update_dual_func(functional);

//		calc_d_gravity(El_Table);

		calc_jacobian(meshctx, propctx, eleminfo, INCREMENT);

//		check_state_vars_with_record(El_Table, solrec, iter);

//		print_jacobian(El_Table, iter);

		calc_adjoint(meshctx, propctx);

		if (eleminfo->iter == iter - 1)
			fill_pertelem_info(El_Table, eleminfo);

		//for first adjoint iteration there is no need to compute Jacobian and adjoint can be computed from the functional
		//sensitivity w.r.t to parameters

//		uinform_refine(meshctx, propctx, numprocs, myid);
//
//		error_compute(meshctx, propctx, iter, myid, numprocs);
//
//		// in dual weighted error estimation if solver performs n step, we'll have n+1
//		// solution and n+1 adjoint solution, but we'll have just n residual and as a
//		// result n error estimate. The point is that at initial step (0'th step),
//		// we know the solution from initial condition  so the error of 0th step is zero,
//		// and we have to compute the error for other time steps.
//
//		double UNREFINE_TARGET = .01;	//dummy value is not used in the function
//		unrefine(El_Table, NodeTable, UNREFINE_TARGET, myid, numprocs, timeprops_ptr, matprops_ptr,
//		    rescomp);

		if (timeprops_ptr->adjiter/*timeprops_ptr->ifadjoint_out() /*|| adjiter == 1*/)
			meshplotter(El_Table, NodeTable, matprops_ptr, timeprops_ptr, mapname_ptr, 0., tecflag);

	}

	return;
}

int num_nonzero_elem(HashTable *El_Table) {
	int num = 0;			//myid
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

//void initSolRec(HashTable* El_Table, HashTable* NodeTable, DualMesh *dualmesh, double dt,
//    int myid) {
//
//	HashEntryPtr* buck = El_Table->getbucketptr();
//	HashEntryPtr currentPtr;
//	Element* Curr_El;
//	int num = 0;
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//
//					Solution *solution = new Solution(Curr_El->get_state_vars(), *(Curr_El->get_kactxy()));
//
//					dualmesh->update_sol(Curr_El, solution);
//
//				}
//				currentPtr = currentPtr->next;
//			}
//		}
//	}
//	return;
//}

void record_solution(MeshCTX* meshctx, PropCTX* propctx, SolRec* solrec) {

	HashTable* El_Table = meshctx->el_table;

	TimeProps* timeptr = propctx->timeprops;

	HashEntryPtr* buck = El_Table->getbucketptr();
	HashEntryPtr currentPtr;
	Element* Curr_El;

	if (timeptr->iter == 1) {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) {
						Jacobian *jacobian = new Jacobian(Curr_El->pass_key(), Curr_El->get_coord());

						if (*(Curr_El->get_state_vars()) > 0.) {
							Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
							    *(Curr_El->get_kactxy()));
							jacobian->put_solution(solution, timeptr->iter - 1);

						} else
							jacobian->put_solution(solrec->get_zero_solution(), timeptr->iter - 1);

						solrec->add(jacobian->get_key(), jacobian);

					}
					currentPtr = currentPtr->next;
				}
			}
		}

	} else {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);
					if (Curr_El->get_adapted_flag() > 0) {

						int aa = 0, bb = 1;
						unsigned keyy[2] = { 3781669179, 330382100 };
						if (*(Curr_El->pass_key()) == keyy[0] && *(Curr_El->pass_key() + 1) == keyy[1]
						    && timeptr->iter == 11)
							bb = aa;

						Jacobian *jacobian = (Jacobian *) solrec->lookup(Curr_El->pass_key());
						if (jacobian) {
							if (*(Curr_El->get_state_vars()) > 0.) {
								Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
								    *(Curr_El->get_kactxy()));
								jacobian->put_solution(solution, timeptr->iter - 1);

							} else
								jacobian->put_solution(solrec->get_zero_solution(), timeptr->iter - 1);

						} else {
							jacobian = new Jacobian(Curr_El->pass_key(), Curr_El->get_coord());
							solrec->add(jacobian->get_key(), jacobian);

							if (*(Curr_El->get_state_vars()) > 0.) {
								Solution *solution = new Solution(Curr_El->get_prev_state_vars(),
								    *(Curr_El->get_kactxy()));
								jacobian->put_solution(solution, timeptr->iter - 1);

							} else
								jacobian->put_solution(solrec->get_zero_solution(), timeptr->iter - 1);

						}
					}
					currentPtr = currentPtr->next;
				}
			}
		}
	}
}

//void initSolRec(HashTable* El_Table, HashTable* NodeTable, DualMesh *dualmesh, double dt,
//    int myid) {
//
//	HashEntryPtr* buck = El_Table->getbucketptr();
//	HashEntryPtr currentPtr;
//	Element* Curr_El;
//	int num = 0;
//
//	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) { // this part allocate memory and initialize jacobian matrices inside the corresponding Jacobian
//		if (*(buck + i)) {
//			currentPtr = *(buck + i);
//			while (currentPtr) {
//				Curr_El = (Element*) (currentPtr->value);
//				if (Curr_El->get_adapted_flag() > 0) {
//
//					Solution *solution = new Solution(Curr_El->get_state_vars(), *(Curr_El->get_kactxy()));
//
//					dualmesh->update_sol(Curr_El, solution);
//
//				}
//				currentPtr = currentPtr->next;
//			}
//		}
//	}
//	return;
//}

void allocJacoMat(HashTable *El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->new_jacobianMat();
				currentPtr = currentPtr->next;
			}
		}

}

double tiny_sgn(double num, double tiny) {
	if (dabs(num) < tiny)
		return 0.;
	else if (num > tiny)
		return 1.;
	else
		return -1.;
}

int num_nonzero_elem(HashTable *El_Table, int type) {
	int num = 0;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() == type)
					num++;
				currentPtr = currentPtr->next;
			}
		}

	return (num);
}

void reverse_states(HashTable* El_Table, HashTable* solrec, int iter, ElemPtrList* refinelist,
    ElemPtrList* unrefinelist) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

#ifdef DEBUG
	if (checkElement(El_Table))
	exit(22);
	for (int i = 0; i < nonz1; i++) {
		dbgvec[i] = 0;
		pass[i] = 0;
	}
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				currentPtr = currentPtr->next;
				if (Curr_El->get_adapted_flag() > 0) {

					int index = Curr_El->get_sol_rec_ind();
					dbgvec[index] += 1;
					pass[index] = 1;

				}
			}
		}
	}
	for (int i = 0; i < nonz1; i++) {
		if (dbgvec[i] != 1) {
			cout << "these elements have problem:   " << i << endl
			<< "the value is  " << dbgvec[i] << endl;
			exit(EXIT_FAILURE);

		} else if (pass[i] != 1) {
			cout << "this index has not been passed  " << i << endl;
			exit(EXIT_FAILURE);

		}
	}

	delete[] dbgvec;

	cout << "number of elements after unrefinement"
	<< num_nonzero_elem(El_Table) << endl;
//		getchar();
	if (checkElement(El_Table))
	exit(22);
#endif

	int reg = 0, ref = 0, unref = 0;
	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					Curr_El->rev_state_vars(solrec, El_Table, iter, &reg, &unref, &ref, refinelist,
					    unrefinelist);

				}
				currentPtr = currentPtr->next;
			}
		}
	}

//	cout << "reg= " << reg << "  ref=" << ref << "  unref= " << unref << endl;
}

void print_jacobian(HashTable* El_Table, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++) {
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->print_jacobian(iter);

				currentPtr = currentPtr->next;
			}
		}
	}

}

void compute_functional(HashTable* El_Table, double* functional, TimeProps* timeprops_ptr) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();
	double const *dx;
	double const *state_vars;
	double const *prev_state_vars;
	double dt;

	dt = timeprops_ptr->dt.at(timeprops_ptr->iter - 1);

//	printf("iter=%4d  dt=%8f \n", timeprops_ptr->iter, dt);

//we do not have make it zero here, because we want to compute the integration over the time and space
//*functional = 0.0;

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {

					dx = Curr_El->get_coord();
					state_vars = Curr_El->get_state_vars();
					prev_state_vars = Curr_El->get_prev_state_vars();

					// we used trapezoidal integration on time
					// flag is for time step 0

					*functional += .5
					    * (state_vars[0] * state_vars[0] + prev_state_vars[0] * prev_state_vars[0]) * dx[0]
					    * dx[1] * dt;

				}
				currentPtr = currentPtr->next;
			}
		}

}

void orgSourceSgn(Element* Curr_El, double frictiny, double* orgSgn) {

	double* d_state_vars_x = Curr_El->get_d_state_vars();
	double* d_state_vars_y = d_state_vars_x + NUM_STATE_VARS;
	double* prev_state_vars = Curr_El->get_prev_state_vars();
	double h_inv;
	double tmp = 0.0;
	double velocity[2];
	for (int i = 0; i < 2; i++)
		orgSgn[i] = 0.0;

	if (prev_state_vars[0] > GEOFLOW_TINY) {

		velocity[0] = prev_state_vars[1] / prev_state_vars[0];
		velocity[1] = prev_state_vars[2] / prev_state_vars[0];

	} else {
		for (int k = 0; k < DIMENSION; k++)
			velocity[k] = 0.;
	}

	if (prev_state_vars[0] > 0.0)
		h_inv = 1. / prev_state_vars[0];

	tmp = h_inv * (d_state_vars_y[1] - velocity[0] * d_state_vars_y[0]);
	orgSgn[0] = tiny_sgn(tmp, frictiny);

	tmp = h_inv * (d_state_vars_x[2] - velocity[1] * d_state_vars_x[0]);
	orgSgn[1] = tiny_sgn(tmp, frictiny);

}

void refinement_report(HashTable* El_Table) {

	int newbuffer = 0, buffer = 0, newson = 0, newfather = 0, norecadapt = 0, tobedeleted = 0,
	    oldfather = 0, oldson = 0;

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				switch (Curr_El->get_adapted_flag()) {
					case 5:
						newbuffer++;
						break;
					case 4:
						buffer++;
						break;
					case 3:
						newson++;
						break;
					case 2:
						newfather++;
						break;
					case 1:
						norecadapt++;
						break;
					case 0:
						tobedeleted++;
						break;
					case -6:
						oldfather++;
						break;
					case -7:
						oldson++;
						break;
					default:
						cout << "this case is irregular \n";
				}
				currentPtr = currentPtr->next;
			}
		}

	cout << " new buffer: " << newbuffer << "\n buffer:     " << buffer << "\n newson:     " << newson
	    << "\n newfather:  " << newfather << "\n norecadapt: " << norecadapt << "\n tobedeleted:"
	    << tobedeleted << "\n oldfather:  " << oldfather << "\n oldson:     " << oldson << "\n";

}

void dual_refine_unrefine(MeshCTX* meshctx, PropCTX* propctx, ElemPtrList* refinelist,
    ElemPtrList* unrefinelist) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	ElemPtrList NewFatherList, OtherProcUpdate;

	int rescomp = 1;

//	unsigned keyy[2] = { 3410598297, 2576980374 };

	delete_unused_elements_nodes(El_Table, NodeTable, myid);

	double target = 0.1;

//	cout << "1 \n";
//	refinement_report(El_Table);

	if (refinelist->get_num_elem()) {
		// start refinement
		for (int i = 0; i < refinelist->get_num_elem(); ++i) {
			refine(refinelist->get(i), El_Table, NodeTable, matprops_ptr, rescomp);
			(refinelist->get(i))->put_adapted_flag(OLDFATHER);
			(refinelist->get(i))->put_refined_flag(1);
		}

//		cout << "2 \n";
//		refinement_report(El_Table);

		refine_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) refinelist, timeprops_ptr);

//		cout << "3 \n";
//		refinement_report(El_Table);

		move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);

//		cout << "4 \n";
//		refinement_report(El_Table);

		int refdel = 0;
		int hash_size = El_Table->get_no_of_buckets();
		for (int i = 0; i < hash_size; i++) {
			HashEntryPtr entryp = *(El_Table->getbucketptr() + i);
			while (entryp) {
				Element* EmTemp = (Element*) (entryp->value);
				entryp = entryp->next;

				if (EmTemp->get_adapted_flag() == TOBEDELETED) {
					El_Table->remove(EmTemp->pass_key());
					refdel++;

				}
			}
		}
//		cout << "5 \n";
//		refinement_report(El_Table);
//		cout << "number of deleted elem after ref " << refdel << "  number of ref list  "
//		    << refinelist->get_num_elem() << endl;
	}

	if (unrefinelist->get_num_elem()) {
		int newfthelem = 0;
		Element* brothers[4];

//		cout << "6 \n";
//		refinement_report(El_Table);

		int counter = 0, unrefined = 0;

		do {

			NewFatherList.trashlist();

			for (int i = 0; i < unrefinelist->get_num_elem(); ++i) {
				Element* Curr_El = (unrefinelist->get(i));
				if ((Curr_El->get_which_son() == 0) && (Curr_El->get_adapted_flag() != OLDSON))

					Curr_El->find_brothers(El_Table, NodeTable, target, myid, matprops_ptr, &NewFatherList,
					    &OtherProcUpdate, rescomp);

			}

			unrefined += NewFatherList.get_num_elem();

//			cout << "7 \n";
//			refinement_report(El_Table);

			unrefine_neigh_update(El_Table, NodeTable, myid, (void*) &NewFatherList);

//			cout << "8 \n";
//			refinement_report(El_Table);

			unrefine_interp_neigh_update(El_Table, NodeTable, numprocs, myid, (void*) &OtherProcUpdate);

			for (int k = 0; k < NewFatherList.get_num_elem(); k++)
				delete_oldsons(El_Table, NodeTable, myid, NewFatherList.get(k));

			reset_adaption_flag(El_Table);
//			reset_newson_adaption_flag(El_Table);

//			cout<<"this is the counter  "<<counter++<<endl;

		} while (unrefined != (unrefinelist->get_num_elem() / 4));

//		cout << "9 \n";
//		refinement_report(El_Table);
//
//		cout << "number of created elem after unref " << unrefined << "  number of unref list  "
//		    << unrefinelist->get_num_elem() << endl;

		move_data(numprocs, myid, El_Table, NodeTable, timeprops_ptr);
	}

//	calc_d_gravity(El_Table);

	refinelist->trashlist();
	unrefinelist->trashlist();
	reset_adaption_flag(El_Table);

}

void reset_adaption_flag(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
					Curr_El->put_refined_flag(0);
				}

				currentPtr = currentPtr->next;
			}
		}

}

void reset_newson_adaption_flag(HashTable* El_Table) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() == NEWSON) {
					Curr_El->put_adapted_flag(NOTRECADAPTED);
					Curr_El->put_refined_flag(0);
				}

				currentPtr = currentPtr->next;
			}
		}

}

void setup_dual_flow(SolRec* solrec, MeshCTX* meshctx, PropCTX* propctx, int iter) {

	HashTable* El_Table = meshctx->el_table;
	HashTable* NodeTable = meshctx->nd_table;

	TimeProps* timeprops_ptr = propctx->timeprops;
	MapNames* mapname_ptr = propctx->mapnames;
	MatProps* matprops_ptr = propctx->matprops;
	int myid = propctx->myid, numprocs = propctx->numproc;

	ElemPtrList refinelist, unrefinelist;
	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	if (iter % 5 == 4) {

		for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
			if (*(buck + i)) {
				currentPtr = *(buck + i);
				while (currentPtr) {
					Curr_El = (Element*) (currentPtr->value);

					if (Curr_El->get_adapted_flag() > 0)
						Curr_El->check_refine_unrefine(solrec, El_Table, iter, &refinelist, &unrefinelist);

					currentPtr = currentPtr->next;
				}
			}

		dual_refine_unrefine(meshctx, propctx, &refinelist, &unrefinelist);

//			setup_geoflow(El_Table, NodeTable, myid, numprocs, matprops_ptr, timeprops_ptr);

		calc_d_gravity(El_Table);

	}

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					Curr_El->update_state(solrec, El_Table, iter);

				currentPtr = currentPtr->next;
			}
		}

}

void check_state_vars_with_record(HashTable* El_Table, HashTable* solrec, int iter) {

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				if (Curr_El->get_adapted_flag() > 0)
					if (Curr_El->check_state(solrec, El_Table, iter))
						cout << "the idea did not work \n";

				currentPtr = currentPtr->next;
			}
		}

}

void print_Elem_Table(HashTable* El_Table,int iter,int place){

	ofstream myfile;
	char filename [50];
	sprintf(filename,"El_Table_%d_%08d",place,iter);

	myfile.open(filename,ios::app);

	HashEntryPtr currentPtr;
	Element *Curr_El;
	HashEntryPtr *buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			currentPtr = *(buck + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);

				myfile<<*(Curr_El->pass_key())<<" "<<*(Curr_El->pass_key()+1)<<" ";
				for (int i=0;i<8*2;++i)
					myfile<<*(Curr_El->get_neighbors()+i)<<" ";
				for (int i=0;i<8;++i)
					myfile<<*(Curr_El->get_neigh_gen()+i)<<" ";
				myfile<<endl;

				currentPtr = currentPtr->next;
			}
		}

	myfile.close();

}
