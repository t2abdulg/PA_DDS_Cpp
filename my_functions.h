#include <math.h>
#include <vector>
#include <iostream>  
#include <algorithm>  
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <bitset>


int dim;
int dimension;
double dSqrtDataNumber;
double volume;


bool sort_by_z(const solution &lhs, const solution &rhs) { return lhs.z < rhs.z; }

bool sort_by_f(const solution &lhs, const solution &rhs) { return lhs.f[dim] < rhs.f[dim]; }

void stream(double regLow[], double regUp[], const vector<double*>& cubs, int lev, double cov);
bool cmp(double* a, double* b);



void copyFile(const char *SRC, const char* DEST)
{
	std::ifstream src(SRC, std::ios::binary);
	std::ofstream dest(DEST, std::ios::binary);
	dest << src.rdbuf();
	src.close(); dest.close();
}

//I read random numbers from files so I do not use this function acctualy
float norm_dist(float mu, float sigma)
{
	float num = 0;
	while (1)
	{
		float v1 = 2.0 * ((float)rand() / RAND_MAX) - 1.0;
		float v2 = 2.0 * ((float)rand() / RAND_MAX) - 1.0;
		float w = v1 * v1 + v2 * v2;

		if (w <= 1)
		{
			float y = ::sqrt(-2.0 * ::log(w) / (float)w) * sigma;
			num = v1 * y + mu;
			break;
		}


	}

	return num;
}



int dominion_status(solution &x1, solution &x2, int &num_objs){

	int ds = 1;

	for (int i = 0; i < num_objs; i++){
		if (x1.f[i] > x2.f[i]){

			goto try_2;
		}
	}

	return 1;

try_2:

	ds = 2;
	for (int i = 0; i < num_objs; i++){
		if (x1.f[i] < x2.f[i]){

			return 0;
		}
	}

	return 2;

}





//the objective functions
float F(solution &x, int &fun_num,int &num_dec, int &num_objs, int &ZDT)
{
	float sum = 0;
	
	if (num_objs == 2){
		switch (fun_num)
		{
		case 0:
		{
				 sum = x.dv[0];//ZDT
		}
			break;
		case 1:
		{
				  
				  switch (ZDT)
				  {
				  case 1:
				  {
							for (int i = 2; i < x.dv.size() + 1; i++){
								sum += x.dv[i - 1];
							}
							sum = 1 + 9.0 * (float)sum / (float)(x.dv.size() - 1);
							sum = sum * (1.0 - sqrt(x.dv[0] / sum));
				  }
					  break;
				  case 2:
				  {
							for (int i = 2; i < x.dv.size() + 1; i++){
								sum += x.dv[i - 1];
							}
							sum = 1 + 9.0 * (float)sum / (float)(x.dv.size() - 1);
							sum = sum * (1.0 - ::powf(x.dv[0] / sum, 2));

				  }
					  break;
				  case 3:
				  {
							for (int i = 2; i < x.dv.size() + 1; i++){
								sum += x.dv[i - 1];
							}
							sum = 1 + 9.0 * (float)sum / (float)(x.dv.size() - 1);
							sum = sum * (1 - ::powf(x.dv[0] / sum, 0.5) - (x.dv[0] / sum) * sin(10 * 3.1415 * x.dv[0]));

				  }
					  break;
				  case 4:
				  {
							for (int i = 2; i < x.dv.size() + 1; i++){
								sum += x.dv[i - 1];
							}
							sum = 1 + 10.0 * (float)sum / (float)(x.dv.size() - 1);
							sum = sum * (1 - ::powf(x.dv[0] / sum, 0.5) - (x.dv[0] / sum) * sin(10 * 3.1415 * x.dv[0]));
				  }
					  break;
				  }

		}
			break;
		}
	}

	else if (num_objs == 3){

		float g_X = 0;

		//DTLZ2
		for (int i = 3; i < 13; i++){
			g_X += ::powf(x.dv[i - 1] - 0.5, 2.0);
		}

		switch (fun_num)
		{
		case 0:
		{
				  sum = ::cosf(x.dv[0] * 3.1415 / 2.0)*::cosf(x.dv[1] * 3.1415 / 2.0)*(1 + g_X);
		}
			break;
		case 1:
		{
				  sum = ::cosf(x.dv[0] * 3.1415 / 2.0)*::sinf(x.dv[1] * 3.1415 / 2.0)*(1 + g_X);
		}
			break;
		case 2:
		{
				  sum = ::sinf(x.dv[0] * 3.1415 / 2.0)*(1 + g_X);
		}
			break;
		}

	}
	
	return sum;
}

float neigh_value_continuous(float &s, float &s_min, float &s_max, float &r, bool &use_manual_random_feed){

	float s_range = s_max - s_min;

	float snew;
	if (use_manual_random_feed){
		rnd_seeder_normal_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		snew = s + rnd_seeder_normal[rnd_seeder_normal_index] * r * s_range;
	}
	else{
		snew = s + norm_dist(0, 1) * r * s_range;
	}
	//float snew = s + 1.0 * r * s_range;

	float P_Abs_or_Ref;
	if (use_manual_random_feed){
		rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		P_Abs_or_Ref = rnd_seeder_uniform[rnd_seeder_uniform_index];
	}
	else{
		P_Abs_or_Ref = ((float)rand() / (RAND_MAX));
	}

	if (snew < s_min) {


		if (P_Abs_or_Ref <= 0.5){
			snew = s_min + (s_min - snew);
		}
		else{
			snew = s_min;
		}
		if (snew > s_max){ snew = s_min; }
	}
	else if (snew > s_max){


		if (P_Abs_or_Ref <= 0.5){
			snew = s_max - (snew - s_max);
		}
		else{
			snew = s_max;
		}

		if (snew < s_min){
			snew = s_max;
		}

	}
	return snew;
}

solution SelectFrom(vector<solution> &archive, bool &use_manual_random_feed){

	//sort(archive.begin(), archive.begin() + archive.size(), sort_by_z);

	vector<float> z_cum(archive.size());
	z_cum[0] = 0;

	for (int i = 0; i < archive.size(); i++){
		if (i == 0){
			z_cum[i] = archive[i].z;
		}
		else{
			z_cum[i] = z_cum[i - 1] + archive[i].z;
		}

	}

	float t;
	if (use_manual_random_feed){
		rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		t = rnd_seeder_uniform[rnd_seeder_uniform_index] * z_cum[archive.size() - 1];
	}
	else{
		t = ((float)rand() / (RAND_MAX)) * z_cum[archive.size() - 1];
	}

	int ii = 0;

	for (int i = 0; i < archive.size(); i++){
		if (z_cum[i] >= t){
			ii = i;
			break;
		}
	}

	return archive[ii];
}

vector<solution> Update_Archive(vector<solution> &archive, solution &candidate, int &num_dec, int &num_objs, int &dominance_flag){
	vector<int> list_of_items_to_be_removed;

	dominance_flag = 0;

	for (int i = 0; i < archive.size(); i++){

		int dominion_status_result = dominion_status(candidate, archive[i],  num_objs);

		if (dominion_status_result == 2){ //candidate is dominated
			dominance_flag = -1;
			//x_new_is_archived = false;
			return archive;
		}
		else if (dominion_status_result == 1){
			list_of_items_to_be_removed.push_back(i);
		}
	}

	//if you are here, it means that the candidate is not dominated and evenmore maybe it has dominated some other solutions.

	if (list_of_items_to_be_removed.size() > 0){

		dominance_flag = 1;

		for (int i = list_of_items_to_be_removed.size() - 1; i > -1; i--){
			archive.erase(archive.begin() + list_of_items_to_be_removed[i]);
		}
	}
	archive.push_back(candidate);

	return archive;
}

double  HV(int data_n, int dim_n, double ref[], vector<double*> &points)
{

	int i, j;
	int dataNumber = data_n;

	dimension = dim_n;

	vector<double*> pointsInitial(dataNumber);
	for (int n = 0; n < dataNumber; n++) {
		pointsInitial[n] = new double[dimension];
		for (int i = 0; i < dimension; i++) {

			pointsInitial[n][i] = points[n][i];
		}
	}
	double* refPoint = new double[dimension];

	for (i = 0; i < dimension; i++) {
		//fileRef >> word;
		refPoint[i] = ref[i];
	}

	// initialize volume
	volume = 0.0;
	// sqrt of dataNumber
	dSqrtDataNumber = sqrt((double)dataNumber);

	// initialize region
	double* regionLow = new double[dimension - 1];
	double* regionUp = new double[dimension - 1];
	for (j = 0; j < dimension - 1; j++)  {
		// determine minimal j coordinate
		double min = 10000000.0;
		for (i = 0; i < dataNumber; i++) {
			if (pointsInitial[i][j] < min) {
				min = pointsInitial[i][j];
			}
		}
		regionLow[j] = min;
		regionUp[j] = refPoint[j];
	}

	// sort pointList according to d-th dimension
	sort(pointsInitial.begin(), pointsInitial.end(), cmp);
	// call stream initially
	stream(regionLow, regionUp, pointsInitial, 0, refPoint[dimension - 1]);

	// print hypervolume
	return volume;
	fflush(stdout);
}

vector<solution> Calc_Z(vector<solution> &archive, int &num_dec, int &num_objs, bool &use_manual_random_feed, int &Select_metric){

	switch (Select_metric)
	{
	case 0://RND
	{
			   for (int i = 0; i < archive.size(); i++){
				   archive[i].z = 1.0;
			   }
	}
		break;
	case 1://CD
	{
			   for (int i = 0; i < archive.size(); i++){
				   archive[i].z = 0;
			   }
			   for (int obj = 0; obj < archive[0].f.size(); obj++){

				   dim = obj;
				   sort(archive.begin(), archive.begin() + archive.size(), sort_by_f);

				   for (int i = 1; i < archive.size() - 1; i++){
					   archive[i].z += abs(archive[i - 1].f[obj] - archive[i + 1].f[obj]) / abs(archive[0].f[obj] - archive[archive.size() - 1].f[obj]);
				   }

				   archive[0].z = archive[1].z;

				   archive[archive.size() - 1].z = archive[archive.size() - 2].z;

			   }
	}
		break;
	case 2://HVC_ESTIMATE
	{
			   for (int i = 0; i < archive.size(); i++){
				   archive[i].z = 0;
			   }


			   //define boundries

			   //float* f_low_bound = new float[m];
			   //float* f_high_bound = new float[m];

			   vector<float>f_low_bound(num_objs);
			   vector<float>f_high_bound(num_objs);

			   //float range_ceo = 0.5;

			   for (int i = 0; i < num_objs; i++){
				   dim = i;
				   sort(archive.begin(), archive.begin() + archive.size(), sort_by_f);

				   f_low_bound[i] = archive[0].f[i];
				   f_high_bound[i] = archive[archive.size() - 1].f[i];

				   //float range_dim_i = abs(f_high_bound[i] - f_low_bound[i]);

				   //f_low_bound[i] -= range_ceo * range_dim_i;
				   //f_high_bound[i] += range_ceo * range_dim_i;

				   /*f_low_bound[i] = 0;
				   f_high_bound[i] = 3;*/

			   }

			   //prepare monte carlo archive
			   long dots_num = 100;
			   vector <solution> mc_points;

			   for (int i = 0; i < dots_num; i++){

				   solution dot(num_dec, num_objs);

				   for (int j = 0; j < num_objs; j++){
					   if (use_manual_random_feed){
						   rnd_seeder_uniform_index += 1;//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						   dot.f[j] = f_low_bound[j] + (f_high_bound[j] - f_low_bound[j])*rnd_seeder_uniform[rnd_seeder_uniform_index];
					   }
					   else{
						   dot.f[j] = f_low_bound[j] + (f_high_bound[j] - f_low_bound[j])*((float)rand() / (RAND_MAX));
					   }
				   }

				   mc_points.push_back(dot);
			   }


			   for (int i = 0; i < dots_num; i++){
				   int jj = archive.size(); bool any_good = false;

				   for (int j = 0; j < archive.size(); j++){

					   if (dominion_status(archive[j], mc_points[i],  num_objs) == 1){
						   jj = j;
						   any_good = true;
						   break;
					   }
				   }

				   if (any_good){
					   for (int k = jj + 1; k < archive.size(); k++){
						   if (dominion_status(archive[k], mc_points[i], num_objs) == 1){
							   //this dot means nothing, goto the next dot.
							   any_good = false;
							   break;
						   }
					   }
				   }

				   //this dot is only dominated by archive[jj]
				   if (any_good){
					   archive[jj].z += 1;
					   any_good = false;
				   }


			   }


			   //normalize z

			   float best_z = 0;

			   for (int i = 0; i < archive.size(); i++){
				   archive[i].z = (archive[i].z / (float)dots_num);
				   if (archive[i].z > best_z){
					   best_z = archive[i].z;
				   }

			   }






			   //for (int i = 0; i < m; i++){

			   //dim = i;

			   //sort(archive.begin(), archive.begin() + archive.size(), sort_by_z);

			   //archive[0].z = best_z;
			   //archive[archive.size() - 1].z = best_z;


			   for (int i = 0; i < archive.size(); i++){
				   if (archive[i].z == 0){
					   archive[i].z = 0.5*best_z;

				   }
			   }


			   //archive[0].z = archive[1].z;
			   //archive[archive.size() - 1].z = archive[archive.size() - 2].z;


			   //destructing 
			   /*for (int j = 0; j < f_low_bound.size(); j++){
				   delete f_low_bound[j];
				   }*/

			   f_low_bound.clear();
			   f_high_bound.clear();








	}
		break;
	case 3://HVC_EXACT
	{

			   int dataNumber = archive.size();
			   int dimension = num_objs;

			   double* refPoint = new double[dimension];
			   for (int i = 0; i < dimension; i++){
				   dim = i;
				   sort(archive.begin(), archive.begin() + archive.size(), sort_by_f);

				   refPoint[i] = 1.00001*(archive[archive.size() - 1].f[i]);
			   }






			   vector<double*> pointsInitial(dataNumber); // will be destructed.
			   for (int i = 0; i < dataNumber; i++) {
				   pointsInitial[i] = new double[dimension];


				   for (int j = 0; j < num_objs; j++){
					   pointsInitial[i][j] = archive[i].f[j];
				   }

			   }

			   float HyperVolume = HV(dataNumber, dimension, refPoint, pointsInitial);

			   //destructing pointsInitial
			   for (int j = 0; j < pointsInitial.size(); j++){
				   delete[] pointsInitial[j];
			   }


			   float best_z = 0;

			   for (int i = 0; i < dataNumber; i++) {



				   vector<int> included_points;
				   for (int j = 0; j < dataNumber; j++){
					   included_points.push_back(j);
				   }

				   included_points.erase(included_points.begin() + i);



				   vector<double*> pointsInitial_sub(dataNumber - 1); //will be removed



				   for (int j = 0; j < dataNumber - 1; j++){
					   pointsInitial_sub[j] = new double[dimension];

					   for (int k = 0; k < num_objs; k++){
						   pointsInitial_sub[j][k] = archive[included_points[j]].f[k];
					   }

				   }

				   archive[i].z = HyperVolume - HV(dataNumber - 1, dimension, refPoint, pointsInitial_sub);

				   if (archive[i].z>best_z){ best_z = archive[i].z; }


				   //destructing pointsInitial_sub
				   for (int j = 0; j < pointsInitial_sub.size(); j++){
					   delete[] pointsInitial_sub[j];
				   }



			   }

			   //destructing refPoint
			   //			   delete [] refPoint;


			   // taking care of the edges
			   for (int i = 0; i < dimension; i++){
				   dim = i;
				   sort(archive.begin(), archive.begin() + archive.size(), sort_by_f);


				   archive[0].z = best_z;
				   archive[archive.size() - 1].z = best_z;
			   }



	}
	}
	return archive;
}

#pragma region Region_1


inline bool cmp(double* a, double* b) {
	return (a[dimension - 1] < b[dimension - 1]);
}
inline bool covers(const double* cub, const double regLow[]) {
	int i_hvc;
	for (i_hvc = 0; i_hvc<dimension - 1; i_hvc++) {
		if (cub[i_hvc] > regLow[i_hvc]) {
			return false;
		}
	}
	return true;
}
inline bool partCovers(const double* cub, const double regUp[])
{
	int i_hvc;
	for (i_hvc = 0; i_hvc < dimension - 1; i_hvc++)
	{
		if (cub[i_hvc] >= regUp[i_hvc])
		{
			return false;
		}
	}
	return true;
}
inline int containsBoundary(const double* cub, const double regLow[], const int split) {
	// condition only checked for split>0
	if (regLow[split] >= cub[split]){
		// boundary in dimension split not contained in region, thus
		// boundary is no candidate for the splitting line
		return -1;
	}
	else {
		int j_hvc;
		for (j_hvc = 0; j_hvc < split; j_hvc++) { // check boundaries
			if (regLow[j_hvc] < cub[j_hvc]) {
				// boundary contained in region
				return 1;
			}
		}
	}
	// no boundary contained in region
	return 0;
}
inline double getMeasure(const double regLow[], const double regUp[]) {
	double vol_hvc;
	int i_hvc;
	vol_hvc = 1.0;
	for (i_hvc = 0; i_hvc < dimension - 1; i_hvc++) {
		vol_hvc *= (regUp[i_hvc] - regLow[i_hvc]);
	}
	return vol_hvc;
}
inline int isPile(const double* cub, const double regLow[], const double regUp[]) {
	int pile_hvc;
	int k_hvc;

	pile_hvc = dimension;
	// check all dimensions of the node
	for (k_hvc = 0; k_hvc<dimension - 1; k_hvc++) {
		// k-boundary of the node's region contained in the cuboid? 
		if (cub[k_hvc] > regLow[k_hvc]) {
			if (pile_hvc != dimension) {
				// second dimension occured that is not completely covered
				// ==> cuboid is no pile
				return -1;
			}
			pile_hvc = k_hvc;
		}
	}
	// if pile == this.dimension then
	// cuboid completely covers region
	// case is not possible since covering cuboids have been removed before

	// region in only one dimenison not completly covered 
	// ==> cuboid is a pile 
	return pile_hvc;
}
inline double computeTrellis(const double regLow[], const double regUp[], const double trellis[]) {

	int i_hvc, j_hvc;
	double vol_hvc;
	int numberSummands_hvc;
	double summand_hvc;
	bitset<16> bitvector_hvc;

	vol_hvc = 0.0;
	summand_hvc = 0.0;
	numberSummands_hvc = 0;

	// calculate number of summands
	bitset<16> nSummands;
	for (i_hvc = 0; i_hvc < dimension - 1; i_hvc++) {
		nSummands[i_hvc] = 1;
	}
	numberSummands_hvc = nSummands.to_ulong();

	double* valueTrellis = new double[dimension - 1];
	double* valueRegion = new double[dimension - 1];
	for (i_hvc = 0; i_hvc < dimension - 1; i_hvc++) {
		valueTrellis[i_hvc] = trellis[i_hvc] - regUp[i_hvc];
	}
	for (i_hvc = 0; i_hvc < dimension - 1; i_hvc++) {
		valueRegion[i_hvc] = regUp[i_hvc] - regLow[i_hvc];
	}


	double* dTemp = new double[numberSummands_hvc / 2 + 1];

	// sum
	for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) {

		// set bitvector length to fixed value 16
		// TODO Warning: dimension-1 <= 16 is assumed
		bitvector_hvc = (long)i_hvc;

		// construct summand
		// 0: take factor from region
		// 1: take factor from cuboid
		summand_hvc = 1.0;
		for (j_hvc = 0; j_hvc < dimension - 2; j_hvc++) {
			if (bitvector_hvc[j_hvc]) {
				summand_hvc *= valueTrellis[j_hvc];
			}
			else {
				summand_hvc *= valueRegion[j_hvc];
			}
		}
		summand_hvc *= valueRegion[dimension - 2];

		// determine sign of summand
		vol_hvc -= summand_hvc;
		dTemp[i_hvc] = -summand_hvc;

		// add summand to sum
		// sign = (int) pow((double)-1, (double)counterOnes+1);
		// vol += (sign * summand); 
	}


	bitvector_hvc = (long)i_hvc;
	summand_hvc = 1.0;
	for (j_hvc = 0; j_hvc < dimension - 1; j_hvc++) {
		if (bitvector_hvc[j_hvc]) {
			summand_hvc *= valueTrellis[j_hvc];
		}
		else {
			summand_hvc *= valueRegion[j_hvc];
		}
	}
	vol_hvc -= summand_hvc;

	for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) {
		summand_hvc = dTemp[i_hvc];
		summand_hvc *= regUp[dimension - 2] - trellis[dimension - 2];
		summand_hvc /= valueRegion[dimension - 2];
		vol_hvc -= summand_hvc;
	}
	//she has already set it{
	delete[] valueTrellis;
	delete[] valueRegion;
	//she has already set it}

	return vol_hvc;
}
// return median of the list of boundaries considered as a set
// TODO linear implementation
inline double getMedian(vector<double>& bounds) {
	// do not filter duplicates
	unsigned int i_hvc;
	if (bounds.size() == 1) {
		return bounds[0];
	}
	else if (bounds.size() == 2) {
		return bounds[1];
	}
	vector<double>::iterator median;
	median = bounds.begin();
	for (i_hvc = 0; i_hvc <= bounds.size() / 2; i_hvc++){
		median++;
	}
	partial_sort(bounds.begin(), median, bounds.end());
	return bounds[bounds.size() / 2];
}
// recursive calculation of hypervolume
inline void stream(double regionLow[], double regionUp[], const vector<double*>& points, int split, double cover) {

	//--- init --------------------------------------------------------------//

	double coverOld_hvc;
	coverOld_hvc = cover;
	unsigned int coverIndex_hvc = 0;
	int c_hvc;

	//--- cover -------------------------------------------------------------//

	// identify first covering cuboid
	double dMeasure = getMeasure(regionLow, regionUp);
	while (cover == coverOld_hvc && coverIndex_hvc < points.size()) {
		if (covers(points[coverIndex_hvc], regionLow)) {
			// new cover value
			cover = points[coverIndex_hvc][dimension - 1];
			volume += dMeasure * (coverOld_hvc - cover);
		}
		else coverIndex_hvc++;
	}

	/* coverIndex shall be the index of the first point in points which
	* is ignored in the remaining process
	*
	* It may occur that that some points in front of coverIndex have the same
	* d-th coordinate as the point at coverIndex. This points must be discarded
	* and therefore the following for-loop checks for this points and reduces
	* coverIndex if necessary.
	*/
	for (c_hvc = coverIndex_hvc; c_hvc > 0; c_hvc--) {
		if (points[c_hvc - 1][dimension - 1] == cover) {
			coverIndex_hvc--;
		}
	}

	// abort if points is empty
	if (coverIndex_hvc == 0) {
		return;
	}
	// Note: in the remainder points is only considered to index coverIndex



	//--- allPiles  ---------------------------------------------------------//

	bool allPiles = true;
	unsigned int iii;

	int* piles_hvc = new int[coverIndex_hvc];
	for (iii = 0; iii < coverIndex_hvc; iii++) {
		piles_hvc[iii] = isPile(points[iii], regionLow, regionUp);
		if (piles_hvc[iii] == -1) {
			allPiles = false;

			//she has already set it{
			delete[] piles_hvc;
			//she has already set it{

			break;
		}
	}

	/*
	* tre llis[i] contains the values of the minimal i-coordinate of
	* the i-piles.
	* If there is no i-pile the default value is the upper bpund of the region.
	* The 1-dimensional KMP of the i-piles is: reg[1][i] - tre llis[i]
	*
	*/

	if (allPiles) { // sweep

		// initialize trellis with region's upper bound
		double* trellis = new double[dimension - 1];
		for (c_hvc = 0; c_hvc < dimension - 1; c_hvc++) {
			trellis[c_hvc] = regionUp[c_hvc];
		}

		double current = 0.0;
		double next = 0.0;
		iii = 0;
		do { // while(next != coverNew)
			current = points[iii][dimension - 1];
			do { // while(next == current)
				if (points[iii][piles_hvc[iii]] < trellis[piles_hvc[iii]]) {
					trellis[piles_hvc[iii]] = points[iii][piles_hvc[iii]];
				}
				iii++; // index of next point
				if (iii < coverIndex_hvc) {
					next = points[iii][dimension - 1];
				}
				else {
					next = cover;
				}

			} while (next == current);
			volume += computeTrellis(regionLow, regionUp, trellis) * (next - current);
		} while (next != cover);
	}
	//--- split -------------------------------------------------------------//
	// inner node of partition tree
	else{
		double bound = -1.0;
		vector<double> boundaries;
		vector<double> noBoundaries;

		do {
			for (iii = 0; iii<coverIndex_hvc; iii++) {
				int contained = containsBoundary(points[iii], regionLow, split);
				if (contained == 1) {
					boundaries.push_back(points[iii][split]);
				}
				else if (contained == 0) {
					noBoundaries.push_back(points[iii][split]);
				}
			}

			if (boundaries.size() >  0) {
				bound = getMedian(boundaries);
				//bound = getRandom(boundaries);
			}
			else if (noBoundaries.size() > dSqrtDataNumber) {
				bound = getMedian(noBoundaries);
				//bound = getRandom(noBoundaries);
			}
			else {
				split++;
			}
		} while (bound == -1.0);



		//I think here I can destruct boundaries vector







		double dLast;
		vector<double*> pointsChild;
		pointsChild.reserve(coverIndex_hvc);

		// left child
		// reduce maxPoint
		dLast = regionUp[split];
		regionUp[split] = bound;
		for (iii = 0; iii < coverIndex_hvc; iii++) {
			if (partCovers(points[iii], regionUp)) {
				pointsChild.push_back(points[iii]);
			}
		}
		if (!pointsChild.empty()) {
			stream(regionLow, regionUp, pointsChild, split, cover);
		}

		// right child
		// increase minPoint
		pointsChild.clear();
		regionUp[split] = dLast;
		dLast = regionLow[split];
		regionLow[split] = bound;
		for (iii = 0; iii < coverIndex_hvc; iii++) {
			if (partCovers(points[iii], regionUp)) {
				pointsChild.push_back(points[iii]);
			}
		}
		if (!pointsChild.empty()) {
			stream(regionLow, regionUp, pointsChild, split, cover);
		}
		regionLow[split] = dLast;

	}// end inner node

} // end stream



#pragma endregion HV_Subroutins