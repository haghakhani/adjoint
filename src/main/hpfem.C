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
 * $Id: hpfem.C 211 2009-06-16 20:02:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#if HAVE_HDF5
#include "../header/GMFG_hdfapi.h"
#endif

//#define LOAD_BAL_DEBUG  //turns on a whole mess of mpi barriers so it makes run time more sensitive to load imbalances i.e. more sensitive to the load balance weights, it just makes it easier to adjust the constants.
//#define PERFTEST
#define TARGETPROC  -1
#define TARGETPROCA -1

int REFINE_LEVEL = 3;

int main(int argc, char *argv[]) {

	//give the information of the element that you want to check, and also the time step
	PertElemInfo* eleminfo;

	MPI_Init(&argc, &argv);

	//test flag indicates that the code is in test mode
	for (int test = 0; test < 1; test++) {

		int i; //-- counters

		double functional = 0.;

		HashTable* BT_Node_Ptr;
		HashTable* BT_Elem_Ptr;

		SolRec* solrec;

		//-- MPI
		int myid, master, numprocs;
		int namelen;
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		MPI_Status status;

		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		MPI_Get_processor_name(processor_name, &namelen);

		double start, end;

		start = MPI_Wtime();

		/* create new MPI datastructures for class objects */
		MPI_New_Datatype();

		/* read original data from serial preprocessing
		 code and then initialize element
		 stiffness routines info */
		int material_count = 0;
		double epsilon = 1., intfrictang = 1, *bedfrictang = NULL, gamma = 1;
		double frict_tiny = .1, mu = 1.0e-03, rho = 1600;
		char **matnames = NULL;
		int xdmerr;

		StatProps statprops;
		MatProps matprops(material_count, matnames, intfrictang, bedfrictang, mu, rho, epsilon, gamma,
		    frict_tiny, 1.0, 1.0, 1.0);
		TimeProps timeprops;
		timeprops.starttime = time(NULL);

		MapNames mapnames;
		PileProps pileprops;
		FluxProps fluxprops;
		OutLine outline;
		DISCHARGE discharge;

		int adaptflag;
		int rescomp = 0;
		int adjflag = 0;
		double end_time = 10000.0;

		int OUTPUT = 0; //Hossein generated this flag to turn off writing output in forward run
		/*
		 * viz_flag is used to determine which viz output to use
		 * nonzero 1st bit of viz_flag means output tecplotxxxx.plt
		 * nonzero 2nd bit of viz_flag means output mshplotxxxx.plt
		 * nonzero 4th bit of viz_flag means output xdmfxxxxx.xmf
		 * nonzero 5th bit of viz_flag means output grass_sites files

		 order_flag == 1 means use first order method
		 order_flag == 2 means use second order method
		 */
		int viz_flag = 0, order_flag, savefileflag = 1; //savefileflag will be flipped so first savefile will end in 0
		int Init_Node_Num, Init_Elem_Num, srctype;
		double v_star; // v/v_slump
		double nz_star; /* temporary... used for negligible velocity as stopping
		 criteria paper... plan to include in v_star implicitly
		 later */

		Read_data(myid, &matprops, &pileprops, &statprops, &timeprops, &fluxprops, &adaptflag,
		    &viz_flag, &order_flag, &mapnames, &discharge, &outline, &srctype);

		if (!loadrun(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, &matprops, &timeprops, &mapnames,
		    &adaptflag, &order_flag, &statprops, &discharge, &outline)) {
			Read_grid(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, &matprops, &outline, &solrec);

			setup_geoflow(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops);

			move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

			AssertMeshErrorFree(BT_Elem_Ptr, BT_Node_Ptr, numprocs, myid, -1.0);

			//initialize pile height and if appropriate perform initial adaptation
			init_piles(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, adaptflag, &matprops, &timeprops,
			    &mapnames, &pileprops, &fluxprops, &statprops);

			setup_geoflow(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops);
		}

		if (myid == 0) {
			for (int imat = 1; imat <= matprops.material_count; imat++)
				printf("bed friction angle for \"%s\" is %g\n", matprops.matnames[imat],
				    matprops.bedfrict[imat] * 180.0 / PI);

			printf("internal friction angle is %g, epsilon is %g \n method order = %i\n",
			    matprops.intfrict * 180.0 / PI, matprops.epsilon, order_flag);
			printf("REFINE_LEVEL=%d\n", REFINE_LEVEL);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		calc_stats(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops, &statprops, &discharge, 0.0);

		output_discharge(&matprops, &timeprops, &discharge, myid);

		move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

		MeshCTX meshctx;
		meshctx.el_table = BT_Elem_Ptr;
		meshctx.nd_table = BT_Node_Ptr;

		PropCTX propctx;
		propctx.timeprops = &timeprops;
		propctx.matprops = &matprops;
		propctx.mapnames = &mapnames;
		propctx.numproc = numprocs;
		propctx.myid = myid;

//		record_solution(&meshctx, &propctx, solrec);

		if (myid == 0)
			output_summary(&timeprops, &statprops, savefileflag);

		if (viz_flag & 1)
			tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar,
			    adjflag);

		if (viz_flag & 2)
			meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar);

#ifdef HAVE_HDF5
		if(viz_flag&8)
		xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_NEW);
#endif

		if (viz_flag & 16) {
			if (myid == 0)
				grass_sites_header_output(&timeprops);
			grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops);
		}

		if (!test)
			find_test_elem(BT_Elem_Ptr, &eleminfo, timeprops.iter);

		if (test) {
			perturbU(BT_Elem_Ptr, eleminfo, timeprops.iter);
		}

		/*
		 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		 Time Stepping Loop

		 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		 */
		long element_counter = 0; // for performance count elements/timestep/proc
		int ifstop = 0;
		double max_momentum = 100;  //nondimensional

		/* ifend(0.5*statprops.vmean) is a hack, the original intent (when we were
		 intending to use vstar as a stopping criteria) whas to have the
		 calculation when vstar dropped back down below 1, instead we're
		 using the ifend() function to stop the simulation when the volume
		 averaged velocity falls back down below 2 meters... this hack is only
		 for the colima hazard map runs, otherwise pass ifend() a constant
		 valued */

		while (!(timeprops.ifend(0)) && !ifstop) //(timeprops.ifend(0.5*statprops.vmean)) && !ifstop)
		{

			/*
			 *  mesh adaption routines
			 */
			double TARGET = .05;
			//double UNREFINE_TARGET = .005;
			double UNREFINE_TARGET = .01;
			int h_count = 0;
			if (timeprops.iter < 50)
				matprops.frict_tiny = 0.1;
			else
				matprops.frict_tiny = 0.000000001;

			//check for changes in topography and update if necessary
			//may want to put an "if(timeprops.iter %20==0)" (20 is arbitrary) here
			if (timeprops.iter == 200) {
				update_topo(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops, &mapnames);
			}

			if ((adaptflag != 0) && (timeprops.iter % 5 == 4)) {
				AssertMeshErrorFree(BT_Elem_Ptr, BT_Node_Ptr, numprocs, myid, -2.0);

//				unsigned keyy[2] = { 635356396, 1321528399 };
//				if (checkElement(BT_Elem_Ptr, NULL, keyy))
//					cout << "I found the suspecious element \n";
//				refinement_report(BT_Elem_Ptr);

				H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_count, TARGET, &matprops, &fluxprops, &timeprops, 5);

//				if (checkElement(BT_Elem_Ptr, NULL, keyy))
//					cout << "I found the suspecious element \n";
//				refinement_report(BT_Elem_Ptr);

				move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

				unrefine(BT_Elem_Ptr, BT_Node_Ptr, UNREFINE_TARGET, myid, numprocs, &timeprops, &matprops,
				    rescomp);

//				if (checkElement(BT_Elem_Ptr, NULL, keyy))
//					cout << "I found the suspecious element \n";
//				refinement_report(BT_Elem_Ptr);

				MPI_Barrier(MPI_COMM_WORLD);      //for debug

				move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work

				if ((numprocs > 1) && (timeprops.iter % 10 == 9)) {

					repartition2(BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

					move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
				}
				move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

				calc_d_gravity(BT_Elem_Ptr);
			}

			step(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops, &pileprops, &fluxprops,
			    &statprops, &order_flag, &outline, &discharge, adaptflag);

//			print_Elem_Table(BT_Elem_Ptr,timeprops.iter,0);

//			unsigned keyy[2] = { 3796806314, 2863311530 };
//			if (checkElement(BT_Elem_Ptr, BT_Node_Ptr,NULL, keyy))
//				cout << "I found the suspecious element \n";
//			cout << "num_elem= " << num_nonzero_elem(BT_Elem_Ptr) << "\n";

			record_solution(&meshctx, &propctx, solrec);

			compute_functional(BT_Elem_Ptr, &functional, &timeprops);

			if (test)
				eleminfo->update_pert_func(functional);
			else
				eleminfo->update_forw_func(functional);

			/*
			 * output results to file
			 */
			//if (OUTPUT) {
			if (timeprops.ifoutput() && OUTPUT) {
				move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

				output_discharge(&matprops, &timeprops, &discharge, myid);

				if (myid == 0) {
					output_summary(&timeprops, &statprops, savefileflag);
				}

				if (viz_flag & 1)
					tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar,
					    adjflag);

				if (viz_flag & 2)
					meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar);

#ifdef HAVE_HDF5
				if(viz_flag&8)
				xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_OLD);
#endif

				if (viz_flag & 16) {
					if (myid == 0)
						grass_sites_header_output(&timeprops);
					grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops);
				}
			}
			//}

#ifdef PERFTEST
			int countedvalue=timeprops.iter%2+1;
			int e_buckets=BT_Elem_Ptr->get_no_of_buckets();
			HashEntry* entryp;
			for(i=0; i<e_buckets; i++)
			{
				entryp = *(BT_Elem_Ptr->getbucketptr() + i);
				while(entryp)
				{
					Element * EmTemp = (Element*)entryp->value;
					assert(EmTemp);
					assert(EmTemp->get_counted()!=countedvalue);

					if((EmTemp->get_adapted_flag()>=NOTRECADAPTED)&&
							(EmTemp->get_adapted_flag()<=BUFFER)
					) {
						//if this element doesn't belong on this processor don't involve
						element_counter++;
						EmTemp->put_counted(countedvalue);
					}
					entryp = entryp->next;
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}

		MPI_Barrier(MPI_COMM_WORLD);

		move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);
		MPI_Barrier(MPI_COMM_WORLD);

		output_discharge(&matprops, &timeprops, &discharge, myid);
		MPI_Barrier(MPI_COMM_WORLD);

		if (myid == 0)
			output_summary(&timeprops, &statprops, savefileflag);

		//printf("hpfem.C 1: xcen=%g\n",statprops.xcen);

		if (viz_flag & 1)
			tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar,
			    adjflag);
		//printf("hpfem.C 2: xcen=%g\n",statprops.xcen);
		MPI_Barrier(MPI_COMM_WORLD);

		if (viz_flag & 2)
			meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar);
		MPI_Barrier(MPI_COMM_WORLD);

#ifdef HAVE_HDF5
		if(viz_flag&8)
		xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_CLOSE);
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if (viz_flag & 16) {
			if (myid == 0)
				grass_sites_header_output(&timeprops);
			grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// write out ending warning, maybe flow hasn't finished moving
		sim_end_warning(BT_Elem_Ptr, &matprops, &timeprops, statprops.vstar);
		MPI_Barrier(MPI_COMM_WORLD);

		if (test) {
			eleminfo->compute_sensitivity();
			eleminfo->dumpinfo();
			eleminfo->dump_functional(&timeprops);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		// write out ending warning, maybe flow hasn't finished moving
		sim_end_warning(BT_Elem_Ptr, &matprops, &timeprops, statprops.vstar);
		MPI_Barrier(MPI_COMM_WORLD);

		if (!test)
			dual_solver(solrec, &meshctx, &propctx, eleminfo);

		//write out the final pile statistics (and run time)
		if (myid == 0)
			out_final_stats(&timeprops, &statprops);

		MPI_Barrier(MPI_COMM_WORLD);

		//write out stochastic simulation statistics
		//if(statprops.lhs.runid>=0)
		if (myid == 0)
			output_stoch_stats(&matprops, &statprops);
		MPI_Barrier(MPI_COMM_WORLD);

		//output maximum flow depth a.k.a. flow outline
		OutLine outline2;
		double dxy[2];
		dxy[0] = outline.dx;
		dxy[1] = outline.dy;
		outline2.init2(dxy, outline.xminmax, outline.yminmax);
		int NxNyout = outline.Nx * outline.Ny;
		MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout,
		MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (myid == 0)
			outline2.output(&matprops, &statprops);

#ifdef PERFTEST  
		long m = element_counter, ii;

		MPI_Allreduce ( &element_counter, &ii, 1,
				MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

		end=MPI_Wtime();
		char perffilename[256];
		sprintf(perffilename,"perform%04d.%04d",numprocs,myid);
		FILE *fpperf=fopen(perffilename,"w");
		fprintf(fpperf,"%d Finished -- used %ld elements of %ld total in %e seconds, %e\n",myid,m,ii,end-start, ii/(end-start));
		fclose(fpperf);
#endif

		Delete_Table(BT_Elem_Ptr, BT_Node_Ptr);
//		if (!test)
//			delete dualmesh;

	}
	delete eleminfo;
	MPI_Finalize();
	return (0);

}
