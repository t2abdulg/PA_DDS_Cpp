#include <stdio.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <limits>
#include <fstream>
#include <direct.h>
#include <iomanip>

#include "solution.h"

long RND_COUNT = 1000000; // numebr of random numbers the code reads or writes from/to files
long rnd_seeder_uniform_index = -1; //this index increases by one, each time that a random number is read from the random numners vector
long rnd_seeder_normal_index = -1;
std::vector<float> rnd_seeder_uniform, rnd_seeder_normal; // vectors in which I store random numbers read from text files




#include "my_functions.h"

using namespace std;

// use this to generate new random numbers and store them in text files
void random_writer()
{
	srand(1);

	fstream myfile;
	myfile.open("uniform.txt", ios::out);
	for (int i = 0; i < RND_COUNT; i++)
	{
		myfile << ((float)rand() / RAND_MAX) << "\n";
	}
	myfile.close();

	myfile.open("normal.txt", ios::out);
	for (int i = 0; i < RND_COUNT; i++)
	{
		myfile << norm_dist(0.0, 1.0) << "\n";
	}
	myfile.close();
}

//use this to read random number from files and store them in arrays
void random_reader(){
	//READING FROM THE FILES

	fstream myfile;
	string line;

	myfile.open("uniform.txt", ios::in);
	//getline(myfile, line);
	for (int i = 0; i < RND_COUNT; i++){
		myfile >> line;
		rnd_seeder_uniform.push_back(std::stof(line));
	}
	myfile.close();

	myfile.open("normal.txt", ios::in);
	for (int i = 0; i < RND_COUNT; i++){
		myfile >> line;
		rnd_seeder_normal.push_back(std::stof(line));
	}
	myfile.close();
}



void DDS(
	string working_folder,
	string project_name,
	int trial, int its,
	int num_dec,
	int num_objs,
	bool use_manual_random_feed,
	vector<float> S_min,
	vector<float> S_max,
	int ZDT,
	int dominance_flag,
	int maxiter,
	int Select_metric,
	float fraction1,
	bool track_best_solution,
	bool resumability,
	bool this_is_a_resume
	)
{




	string pareto_filename = working_folder;
	pareto_filename.append("//");
	pareto_filename.append(project_name);
	pareto_filename.append("_nondom_sol_");
	if (this_is_a_resume){
		pareto_filename.append("RESUME");
	}
	else{
		pareto_filename.append(std::to_string(trial));
	}
	pareto_filename.append(".txt");




	string solutions_filename = working_folder;
	solutions_filename.append("//");
	solutions_filename.append(project_name);
	solutions_filename.append("_nondom_pts_");
	if (this_is_a_resume){
		solutions_filename.append("RESUME");
	}
	else{
		solutions_filename.append(std::to_string(trial));
	}
	solutions_filename.append(".txt");




	string track_filename = working_folder;
	track_filename.append("//");
	track_filename.append(project_name);
	track_filename.append("_track_best_");
	track_filename.append(std::to_string(trial));
	track_filename.append(".txt");
	fstream track_best_file;
	if (track_best_solution){ track_best_file.open(track_filename, ios::out); }



	//there is only one backup for whole trials. It means once a trial is finished, then the backup is overwritten by the next trial.
	string backup_filename = working_folder;
	backup_filename.append("//");
	backup_filename.append(project_name);
	backup_filename.append("_backup");
	backup_filename.append(".txt");






	vector <solution> Archive;// the main Archive


	clock_t begin = clock();






	int istart = 1;


	if (!this_is_a_resume){




		std::cout << "Trial " << trial << "\n";
		std::cout << "Creating initial solutions. \n";


		//creating initials
		for (int i = 1; i < its + 1; i++){

			solution stest(num_dec, num_objs);
			for (int j = 0; j < num_dec; j++){
				if (use_manual_random_feed){
					//Where you see %%%%%%%%%%%%%%%%%, in the next lines there will be random numers calling
					rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					stest.dv[j] = S_min[j] + (S_max[j] - S_min[j])*rnd_seeder_uniform[rnd_seeder_uniform_index];
				}
				else{
					stest.dv[j] = S_min[j] + (S_max[j] - S_min[j])*((float)rand() / (RAND_MAX));
				}
			}
			for (int j = 0; j < num_objs; j++){
				stest.f[j] = F(stest, j, num_dec, num_objs, ZDT);
			}

			if (i == 1){
				Archive.push_back(stest);
			}
			else{

				Archive = Update_Archive(Archive, stest, num_dec, num_objs, dominance_flag);

				/*if (dominance_flag != -1){
				Archive.push_back(stest);
				}*/
			}

		}
		//finished creating initials




	}
	else{//this is a resume run
		//need to read resuability settings
		//there should be a DDS_inp.txt
		//there should be a track_best file
		//there should be a 
		//assign values to maxiter, istart, and fill Archive

		//read the archive and put it in the Archive.
		

		fstream myfile_res;
		string line;

		myfile_res.open(backup_filename, ios::in);
		myfile_res >> line;	istart = std::stoi(line) + 1;
		myfile_res >> line; dominance_flag = std::stoi(line);
		myfile_res >> line; int Archive_size = std::stoi(line);

		//reading that Archive

		for (int i = 0; i < Archive_size; i++){
			solution stemp(num_dec, num_objs);
			for (int j = 0; j < num_objs; j++){
				myfile_res >> line; stemp.f[j] = std::stof(line);
			}
			Archive.push_back(stemp);
		}

		for (int i = 0; i < Archive_size; i++){
			for (int j = 0; j < num_dec; j++){
				myfile_res >> line; Archive[i].dv[j] = std::stof(line);
			}
		}
		myfile_res.close();
	}






	int iLeft = maxiter - its;

	//Calculating Selection Metric z:
	Archive = Calc_Z(Archive, num_dec, num_objs, use_manual_random_feed, Select_metric);
	//solution x_cur = SelectFrom(Archive);

	solution sbest(num_dec, num_objs);// = Archive[Archive.size() - 1];

	std::cout << "Main loop iterations: \n";


	//MAIN LOOP
	for (int i = istart; i < iLeft + 1; i++){
		std::cout << " " << i;

		if (dominance_flag == -1){
			sbest = SelectFrom(Archive, use_manual_random_feed);
		}
		else{
			sbest = Archive[Archive.size() - 1];
		}

		// %% DDS

		float Pn = 1 - log10(i) / log10(iLeft);
		int dvn_count = 0;

		vector <float> randnums;
		for (int j = 0; j < num_dec; j++){
			if (use_manual_random_feed){
				rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				randnums.push_back(rnd_seeder_uniform[rnd_seeder_uniform_index]);
			}
			else{
				randnums.push_back(((float)rand() / (RAND_MAX)));
			}
		}
		solution stest = sbest;
		for (int j = 0; j < num_dec; j++){
			if (randnums[j] < Pn){
				dvn_count += 1;
				float new_value = neigh_value_continuous(sbest.dv[j], S_min[j], S_max[j], fraction1, use_manual_random_feed);
				stest.dv[j] = new_value;// % change relevant dec var value in stest

			}
		}
		if (dvn_count == 0){
			int dec_var;
			if (use_manual_random_feed){
				rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				dec_var = ceil(num_dec * rnd_seeder_uniform[rnd_seeder_uniform_index]) - 1;
			}
			else{
				dec_var = ceil(num_dec * ((float)rand() / (RAND_MAX))) - 1;
				if (dec_var < 0){ dec_var = 0; }
			}



			float new_value = neigh_value_continuous(sbest.dv[dec_var], S_min[dec_var], S_max[dec_var], fraction1, use_manual_random_feed);
			stest.dv[dec_var] = new_value;
			//stest.dv[0] = new_value;

			if (new_value == 1){
				int asjkdslifhsl = 34234;
			}

		}

		//solution stest = perturb_sol(sbest, Pn, S_min, S_max, fraction1);

		for (int j = 0; j < num_objs; j++){
			stest.f[j] = F(stest, j, num_dec, num_objs, ZDT);
		}




		//checking to see if x_curr dominates x_new
		bool sbest_dominates_stest = false;

		if (dominion_status(stest, sbest, num_objs) == 2){
			sbest_dominates_stest = true;
			dominance_flag = -1;
		}
		else{

			//check to see if it is a duplicate, if yes, flag = 0, do not Update_Archive, and do not add

			bool stest_is_duplicate = false;
			for (int j = 0; j < Archive.size(); j++){

				//if (equal(stest.f.begin(), stest.f.end(), Archive[j].f.begin())){

				/*float sum1 = 0; float sum2 = 0;
				for (int jj = 0; jj < num_objs; jj++){
				sum1 += stest.f[jj] ;
				sum2 += Archive[j].f[jj] ;
				}*/

				if (stest.f == Archive[j].f){
					//if (fabs(sum1 - sum2) < 0.00001){
					stest_is_duplicate = true;
					dominance_flag = 0;
					break;
				}
			}



			if (stest_is_duplicate == false){ Archive = Update_Archive(Archive, stest, num_dec, num_objs, dominance_flag); }

			stest_is_duplicate = false;

			if (dominance_flag != -1){

				Archive = Calc_Z(Archive, num_dec, num_objs, use_manual_random_feed, Select_metric);
			}

		}


		if (track_best_solution){
			//recording the best solution

			for (int itrack = 0; itrack < num_objs; itrack++)
			{
				track_best_file << sbest.f[itrack] << "\t";
			}

			for (int itrack = 0; itrack < num_dec; itrack++)
			{
				track_best_file << sbest.dv[itrack] << "\t";
			}
			track_best_file << "\n";


		}


		if (!this_is_a_resume){
			if (resumability){
				//recording the whole information

				fstream backup_file;
				//remove(backup_filename.c_str());
				backup_file.open(backup_filename, ios::out);


				backup_file << i << "\n";
				backup_file << dominance_flag << "\n";
				backup_file << Archive.size() << "\n";


				for (int ires = 0; ires < Archive.size(); ires++){
					for (int jres = 0; jres < num_objs - 1; jres++){
						backup_file << Archive[ires].f[jres] << "\t";
					}
					backup_file << Archive[ires].f[num_objs - 1] << "\n";
				}


				for (int ires = 0; ires < Archive.size(); ires++)
				{
					for (int jres = 0; jres < num_dec - 1; jres++){
						backup_file << Archive[ires].dv[jres] << "\t";
					}
					backup_file << Archive[ires].dv[num_dec - 1] << "\n";
				}

				backup_file.close();

			}
		}









	}

	clock_t end = clock();

	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "\n Process time in sec = " << elapsed_secs << endl;

	dim = 0;
	std::sort(Archive.begin(), Archive.begin() + Archive.size(), sort_by_f);

	if (track_best_solution){ track_best_file.close(); }


	//reporting the results
	fstream myfile;
	myfile.open(pareto_filename, ios::out);
	for (int i = 0; i < Archive.size(); i++)
	{
		for (int j = 0; j < num_objs - 1; j++){
			myfile << Archive[i].f[j] << "\t";
		}
		myfile << Archive[i].f[num_objs - 1] << "\n";
	}
	myfile.close();


	myfile.open(solutions_filename, ios::out);
	for (int i = 0; i < Archive.size(); i++)
	{
		for (int j = 0; j < num_dec - 1; j++){
			myfile << Archive[i].dv[j] << "\t";
		}
		myfile << Archive[i].dv[num_dec - 1] << "\n";
	}
	myfile.close();


	std::cout << Archive.size() << " points achieved. Output file built. \n Done.~~~~~~~~~~~~~~~~~~~~~~~ \n \n";


}





void main(void)
{


	//read the input file
	fstream myfile_inp;
	string line;

	myfile_inp.open("DDS_inp.txt", ios::in);

	myfile_inp >> line; myfile_inp >> line;

	myfile_inp >> line;	const std::string project_name = line; myfile_inp >> line;//project name
	myfile_inp >> line;	int trials = std::stoi(line); myfile_inp >> line;//number of optimization trials to run
	myfile_inp >> line;	int maxiter = std::stoi(line); myfile_inp >> line;//maximum number of objective function evaluations per optimization trial
	myfile_inp >> line;	bool use_manual_random_feed = std::stoi(line); myfile_inp >> line;
	myfile_inp >> line;	int random_seed = std::stoi(line); myfile_inp >> line;//random seeder for C++ random number generators
	myfile_inp >> line;	float fraction1 = std::stof(line); myfile_inp >> line;
	myfile_inp >> line;	int Select_metric = std::stoi(line); myfile_inp >> line;//selection metric --> 0: Random,    1: Crowding distance,     2: Hypervolume Contribution (ESTIMATE)     3: Hypervolume Contribution (EXACT)

	myfile_inp >> line;	int num_dec = std::stoi(line); myfile_inp >> line;//number of decesion variables (depends on the objective functions and the problem we are solving)
	myfile_inp >> line;	int num_objs = std::stoi(line); myfile_inp >> line;//number of objective functions
	myfile_inp >> line;	int ZDT = std::stoi(line); myfile_inp >> line;//number of objective function for bio-dimension problems , we have ZDT1, ZDT2, ZDT3, and ZDT4
	myfile_inp >> line;	bool track_best_solution = std::stoi(line); myfile_inp >> line;//
	myfile_inp >> line;	bool resumability = std::stoi(line); myfile_inp >> line;//
	myfile_inp >> line;	bool this_is_a_resume = std::stoi(line); myfile_inp >> line;//


	myfile_inp.close();


	int dominance_flag = 0;

	if (this_is_a_resume){ trials = 1; }

	string selectionM;
	switch (Select_metric)
	{
	case 0: {			  selectionM = "RND";	}
		break;
	case 1:	{			  selectionM = "CD";	}
		break;
	case 2:	{			  selectionM = "HVA";	}
		break;
	case 3:	{			  selectionM = "HVC";	}
		break;
	}


	string working_folder = project_name;
	working_folder.append("_iter");
	working_folder.append(std::to_string(maxiter));
	working_folder.append("_");
	working_folder.append(selectionM);

	_mkdir(working_folder.c_str());

	myfile_inp.close();

	if (use_manual_random_feed){
		random_writer();
		random_reader();
	}
	else{
		srand(random_seed);
	}

	std::cout << "Random numbers loaded. \n";

	int its = max(5.0, 0.005 * maxiter); //number of initial solutions


	//read dv bounds
	vector<float> S_min, S_max; //low and high bounds


	fstream myfile_bounds;
	myfile_bounds.open("dv.txt", ios::in);
	myfile_bounds >> line; myfile_bounds >> line; myfile_bounds >> line; myfile_bounds >> line; myfile_bounds >> line;
	for (int i = 0; i < num_dec; i++){ //assign custom bounds regarding the objective function selected, for all ZDTs it is the same
		myfile_bounds >> line;
		myfile_bounds >> line;
		S_min.push_back(std::stof(line));
		myfile_bounds >> line;
		S_max.push_back(std::stof(line));
		myfile_bounds >> line;
	}
	myfile_bounds.close();



	//copy DDS_inp.txt and dv.txt to the project folder
	string destt = working_folder; destt.append("//DDS_inp.txt");
	copyFile("DDS_inp.txt", destt.c_str());

	destt = working_folder; destt.append("//dv.txt");
	copyFile("dv.txt", destt.c_str());



	//if this is a resume then there will be only one trial and the filenames will be changed to *_RESUME.txt

	for (int trial = 1; trial < trials + 1; trial++){
		DDS(working_folder, project_name, trial, its, num_dec, num_objs, use_manual_random_feed, S_min, S_max, ZDT, dominance_flag, maxiter, Select_metric, fraction1, track_best_solution, resumability, this_is_a_resume);
	}




	char aaa; std::cin >> aaa;
}


