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
 * $Id: get_coef_and_eigen.C,v 1.4 2004/08/11 15:58:46 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define DO_EROSION
#define DEBUG

#define KEY0   3941335040
#define KEY1   0
#define ITER   9
#define X      1.97853625e-01
#define Y      2.30663500e-01

#include "../header/hpfem.h"

void correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, Element *EmTemp, double *forceint, double *forcebed,
    double *eroded, double *deposited) {
//	Element *EmTemp = (Element *) EmTemp_in;
	double *dx = EmTemp->get_dx();
	double dtdx = dt / dx[0];
	double dtdy = dt / dx[1];
	double kactxy[DIMENSION];

	double tiny = GEOFLOW_TINY;
	int xp = EmTemp->get_positive_x_side();
	int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;

	int ivar, i, j, k;
	double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
	double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

	Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode() + (xp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxp[ivar] = nxp->flux[ivar];

	Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode() + (yp + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxyp[ivar] = nyp->flux[ivar];

	Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode() + (xm + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxxm[ivar] = nxm->flux[ivar];

	Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode() + (ym + 4) * 2);
	for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
		fluxym[ivar] = nym->flux[ivar];

#ifdef DO_EROSION
	int do_erosion = 1;
#else
	int do_erosion=0;
#endif

#ifdef STOPCRIT_CHANGE_SOURCE
	int IF_STOPPED=EmTemp->get_stoppedflags();
#else
	int IF_STOPPED = !(!EmTemp->get_stoppedflags());
#endif

	double *state_vars = EmTemp->get_state_vars();
	double *prev_state_vars = EmTemp->get_prev_state_vars();
	double *d_state_vars = EmTemp->get_d_state_vars();
	double *gravity = EmTemp->get_gravity();
	double *d_gravity = EmTemp->get_d_gravity();
	double *zeta = EmTemp->get_zeta();
	double *curvature = EmTemp->get_curvature();
	double bedfrict = EmTemp->get_effect_bedfrict();
	double *Influx = EmTemp->get_influx();
	double terminal_vel = matprops_ptr->v_terminal;
	double navslip_coef = matprops_ptr->navslip_coef;

	double Vsolid[DIMENSION];

	if (state_vars[0] > GEOFLOW_TINY) {
		for (i = 0; i < DIMENSION; i++)
			kactxy[i] = *(EmTemp->get_effect_kactxy() + i);

		// fluid velocities
		Vsolid[0] = state_vars[1] / state_vars[0];
		Vsolid[1] = state_vars[2] / state_vars[0];

	} else {
		for (i = 0; i < DIMENSION; i++) {
			kactxy[i] = matprops_ptr->epsilon;
			Vsolid[i] = 0.;
		}
		bedfrict = matprops_ptr->bedfrict[EmTemp->get_material()];
	}

	double V_avg[DIMENSION];
	V_avg[0] = Vsolid[0];
	V_avg[1] = Vsolid[1];
	EmTemp->convect_dryline(V_avg, dt); //this is necessary

	int debuging, ggg = 0;
	if (*(EmTemp->pass_key()) == KEY0 && *(EmTemp->pass_key() + 1) == KEY1 && timeprops->iter == ITER)
//	if (dabs(*(EmTemp->get_coord()) - X) < INCREMENT&& dabs(*(EmTemp->get_coord()+1) - Y)<INCREMENT
//	&& timeprops->iter == ITER)
		debuging = ggg = 1;

	double dragforce[2] = { 0., 0. };
	correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt,
	    d_state_vars, (d_state_vars + NUM_STATE_VARS), &(zeta[0]), &(zeta[1]), curvature,
	    &(matprops_ptr->intfrict), &bedfrict, gravity, d_gravity, kactxy, &(matprops_ptr->frict_tiny),
	    forceint, forcebed, dragforce, &do_erosion, eroded, Vsolid, &terminal_vel,
	    &(matprops_ptr->epsilon), &IF_STOPPED, Influx);

#ifdef DEBUG
	if (*(EmTemp->pass_key()) == KEY0 && *(EmTemp->pass_key() + 1) == KEY1 && timeprops->iter == ITER) {
		int state_num = NUM_STATE_VARS;
		FILE *fp;
		fp = fopen("debugfile", "w");

		fprintf(fp,
		    "In corrector time step %d with dt=%f dtdx=%f dtdy=%f kactx=%f , kacty=%f , x=%6f, y=%6f \n state vars are: \n",
		    timeprops->iter, dt, dtdx, dtdy, *(EmTemp->get_kactxy()), *(EmTemp->get_kactxy() + 1),
		    *(EmTemp->get_coord()), *(EmTemp->get_coord() + 1));
		for (i = 0; i < state_num; i++)
			fprintf(fp, "%10e ", state_vars[i]);
		fprintf(fp, "\n prev state vars are: \n");
		for (i = 0; i < state_num; i++)
			fprintf(fp, "%10e ", prev_state_vars[i]);
		fprintf(fp, "\n xp: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", fluxxp[state]);
		fprintf(fp, "\n xm: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", fluxxm[state]);
		fprintf(fp, "\n yp: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", fluxyp[state]);
		fprintf(fp, "\n ym: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", fluxym[state]);
		fprintf(fp, "\n dUdx: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", d_state_vars[state]);
		fprintf(fp, "\n dUdy: \n");
		for (int state = 0; state < state_num; state++)
			fprintf(fp, "%10e ", d_state_vars[state + NUM_STATE_VARS]);
		fprintf(fp, "\n");

		fclose(fp);
	}
#endif

	EmTemp->put_drag(dragforce);
	*forceint *= dx[0] * dx[1];
	*forcebed *= dx[0] * dx[1];
	*eroded *= dx[0] * dx[1];

	bool print_vars = false;
	for (i = 0; i < NUM_STATE_VARS; i++)
		if (state_vars[i] > 1e3) // || dabs(state_vars[5])>5 )
			print_vars = false;

	if (print_vars) {
		double tempU[NUM_STATE_VARS];
		for (i = 0; i < NUM_STATE_VARS; i++)
			tempU[i] = prev_state_vars[i] - dtdx * (fluxxp[i] - fluxxm[i])
			    - dtdy * (fluxyp[i] - fluxym[i]);
		printf("ElemKey: %d   ElemKey2: %d\n ", *EmTemp->pass_key(), *(EmTemp->pass_key() + 1));
		printf("Kactxy = %10.5e, %10.5e\n", kactxy[0], kactxy[1]);
		printf("BedFrict: %10.5e: IntFrict: %10.5e\n", bedfrict, matprops_ptr->intfrict);

		printf("state_vars: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", state_vars[i]);
		printf("\n");

		printf("prev_state_vars: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", prev_state_vars[i]);
		printf("\n");

		printf("Ustore: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("%10.5e, ", tempU[i]);
		printf("\n");

		printf("fluxes: \n");
		for (i = 0; i < NUM_STATE_VARS; i++)
			printf("fluxxp:%10.5e, fluxxm:%10.5e, fluxyp:%10.5e, fluxym:%10.5e \n ", fluxxp[i], fluxxm[i],
			    fluxyp[i], fluxym[i]);

		exit(1);
	}

	if (EmTemp->get_stoppedflags() == 2)
		*deposited = state_vars[0] * dx[0] * dx[1];
	else
		*deposited = 0.0;

	if (EmTemp->get_stoppedflags())
		*eroded = 0.0;

	EmTemp->calc_shortspeed(1.0 / dt);
	return;
}
