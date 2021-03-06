#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
//#include <sys/time.h>
#include <cstdlib>
#include <float.h>
#include <random>
#include <stdint.h>
#include <algorithm>
#include <cstdio>
#include <climits>
#include <vector>

#define _USE_MATH_DEFINES
using namespace std;

int NE = 75;
int NI = 25;
double PEE = 0.15;//probability of postsynaptic connections
double PIE = 0.5;
double PEI = 0.5;
double PII = 0.4;
double kickE = 7000.0;//external drive rate
double kickI = 7000.0;
double HitEE = 1000.0 / 1.4;//delay time rate
double HitIE = 1000.0 / 1.2;
double HitI = 1000.0 / 4.5;
int E_spike = 0;
int I_spike = 0;
int S_e = 24;
int S_i = 48;
double a_ee = 0.5;
double a_ie = 0.5;
double a_ei = 0.79;
double a_ii = 0.21;

vector<double> P_BE_Ex(NE + 1);
vector<double> P_GE_Ex(NE + 1);
vector<double> P_BI_Ex(NI + 1);
vector<double> P_GI_Ex(NI + 1);
vector<double> P_BE_E(NE + 1);
vector<double> P_BI_E(NI + 1);
vector<double> P_GE_E(NE + 1);
vector<double> P_GI_E(NI + 1);
vector<double> P_GE_I(NE + 1);
vector<double> P_GI_I(NI + 1);




template <class T>
class Vector : public vector<T> {
private:
	T sum;
public:

	double get_sum()
	{
		return sum;
	}

	void maintain()
	{
		T tmp_sum = 0;
		for (auto&& it : (*this))
			tmp_sum += it;
		//       cout<<tmp_sum<<endl;
		sum = tmp_sum;
	}

	void switch_element(const int index, T element)
		//switch [index] by element and maintain sum
	{
		sum = sum - (*this)[index] + element;
		(*this)[index] = element;
	}
	//The first three function is only for Clock[]
};



int find_index(Vector<double>& array, mt19937& mt, uniform_real_distribution<double>& u)
//find a random element in a positive array. The probability of chosen is proportional to the value.
{
	double sum = array.get_sum();
	double tmp = u(mt);
	int index = 0;
	int size = array.size();
	double tmp_double = array[0] / sum;
	//    cout<<tmp_double<<" "<<tmp<<endl;
	while (tmp_double < tmp && index < size - 1)
	{
		index++;
		tmp_double += array[index] / sum;
	}
	return index;
}


void update(vector<double>& time_spike, vector<int>& num_spike, vector<double>& record_time_point, vector<int>& record_N_GE, vector<int>& record_N_GI, vector<int>& record_HE,
	vector<int>& record_HI, Vector<double>& Clock, int N_GE, int N_GI, int HE, int HI, const double terminate_time, mt19937& mt, uniform_real_distribution<double>& u)
{



	double current_time = 0.0;
	double record_time = 0.0;
	int count = 0;
	while (current_time < terminate_time)
	{
		if (current_time - record_time > 0.005)
		{
			record_time = current_time;
			record_time_point.push_back(current_time);
			record_HI.push_back(HI);
			record_HE.push_back(HE);

			record_N_GE.push_back(N_GE);
			record_N_GI.push_back(N_GI);

		}
		current_time += -log(1 - u(mt)) / Clock.get_sum();
		int index = find_index(Clock, mt, u);
		int whichHit;

		//            cout<<"time "<<current_time <<" index "<<index<<endl;
		switch (index)
		{
		case 0:
			count++;
			// external E 
			if (u(mt)<1-N_GE*1.0/NE)
			{
				if (u(mt) < P_BE_Ex[N_GE])
				{
					N_GE += 1;
				}
			}
			else
			{
				
				if (u(mt) < P_GE_Ex[N_GE])
				{
					E_spike++;
					N_GE -= 1;
					HE += S_e;
					time_spike.push_back(current_time);
					num_spike.push_back(0); 
					Clock.switch_element(2, HitEE * HE * a_ee);
					Clock.switch_element(3, HitIE * HE * a_ie);
				}
			}
			break;
		case 1:
			// external I
			if (u(mt) < 1 - N_GI*1.0 / NI)
			{
				if (u(mt) < P_BI_Ex[N_GI])
				{
					N_GI += 1;
				}
			}
			else
			{
				if (u(mt) < P_GI_Ex[N_GI])
				{
					I_spike++;
					N_GI -= 1;
					HI += S_i;
					time_spike.push_back(current_time);
					num_spike.push_back(1);
					Clock.switch_element(4, HitI * HI * a_ei);
					Clock.switch_element(5, HitI * HI * a_ii);
				}
			}
			break;
		case 2:
			if (u(mt) < 1 - N_GE*1.0 / NE)
			{
				if (u(mt) < P_BE_E[N_GE])
				{
					N_GE += 1;
				}
			}
			else
			{
				if (u(mt) < P_GE_E[N_GE])
				{
					E_spike++;
					N_GE -= 1;
					HE += S_e;
					time_spike.push_back(current_time);
					num_spike.push_back(0);
				}
			}
			HE--;
			Clock.switch_element(2, HitEE * HE * a_ee);
			Clock.switch_element(3, HitIE * HE * a_ie);
			break;
		case 3:
			if (u(mt) < 1 - N_GI*1.0 / NI)
			{
				if (u(mt) < P_BI_E[N_GI])
				{
					N_GI += 1;
				}
			}
			else
			{
				if (u(mt) < P_GI_E[N_GI])
				{
					I_spike++;
					N_GI -= 1;
					HI += S_i;
					time_spike.push_back(current_time);
					num_spike.push_back(1);
					Clock.switch_element(4, HitI* HI* a_ei);
					Clock.switch_element(5, HitI* HI* a_ii);
				}
			}
			HE--;
			Clock.switch_element(2, HitEE * HE * a_ee);
			Clock.switch_element(3, HitIE * HE * a_ie);
			break;
		case 4:
			// I pending spikes
			
			if (u(mt) < N_GE*1.0 / NE)
			{
				if (u(mt) < P_GE_I[N_GE])
				{
					N_GE -= 1;
				}
			}
			HI--;
			Clock.switch_element(4, HitI * HI * a_ei);
			Clock.switch_element(5, HitI * HI * a_ii);
			break;
		case 5:
			if (u(mt) < N_GI*1.0 / NI)
			{
				if (u(mt) < P_GI_I[N_GI])
				{
					N_GI -= 1;
				}
			}
			HI--;
			Clock.switch_element(4, HitI * HI * a_ei);
			Clock.switch_element(5, HitI * HI * a_ii);
			break;
		}
	}
}



int main()
{
	ifstream P("Probability.txt");
	for (int i = 0; i <= NE; i++)
	{
		P >> P_BE_Ex[i];
	}

	for (int i = 0; i <= NE; i++)
	{
		P >> P_GE_Ex[i];
	}
	for (int i = 0;i <= NE;i++)
	{
		P >> P_BE_E[i];
	}
	for (int i = 0;i <= NE;i++)
	{
		P >> P_GE_E[i];
	}
	for (int i = 0;i <= NE;i++)
	{
		P >> P_GE_I[i];
	}
	cout << endl;
	for (int i = 0;i <= NI;i++)
	{
		P >> P_BI_Ex[i];
	}
	cout << endl;
	for (int i = 0;i <= NI;i++)
	{
		P >> P_GI_Ex[i];
	}
	cout << endl;
	for (int i = 0;i <= NI;i++)
	{
		P >> P_BI_E[i];
	}
	cout << endl;
	for (int i = 0;i <= NI;i++)
	{
		P >> P_GI_E[i];
	}
	cout << endl;
	for (int i = 0;i <= NI;i++)
	{
		P >> P_GI_I[i];
	}
	//  struct timeval t1, t2;
	//  gettimeofday(&t1, NULL);
	ofstream myfile;
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> u(0, 1);
	myfile.open("spike_info_coarse_largesize.txt");

	int N_GE = 0;
	int N_GI = 0;
	int HE = 0;
	int HI = 0;


	Vector<double> Clock;
	Clock.reserve(6);
	// 0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII

	Clock.push_back(NE * kickE);
	Clock.push_back(NI * kickI);
	for (int i = 2; i < 6; i++)
		Clock.push_back(0);
	Clock.maintain();
	double terminate_time = 1;
	vector<double> time_spike;
	time_spike.reserve(100000);
	vector<int> num_spike;// 0 represents E spike and 1 represents I spike.
	num_spike.reserve(100000);
	vector<double> record_time_point;
	record_time_point.reserve(250000);
	vector<int> record_N_GE;
	record_N_GE.reserve(250000);
	vector<int> record_N_GI;
	record_N_GI.reserve(250000);
	vector<int> record_HE;
	record_HE.reserve(250000);
	vector<int> record_HI;
	record_HI.reserve(250000);



	update(time_spike, num_spike, record_time_point, record_N_GE, record_N_GI, record_HE, record_HI, Clock, N_GE, N_GI, HE, HI, terminate_time, mt, u);

	cout << "E spike rate= " << (double)E_spike / (terminate_time * NE) << endl;
	cout << "I spike rate = " << (double)I_spike / (terminate_time * NI) << endl;
	int spike_count = time_spike.size();
	for (int i = 0; i < spike_count; i++)
	{
		myfile << time_spike[i] << "  " << num_spike[i] << endl;
	}
	myfile.close();
	/*myfile.open("trajectory_info_sample_largesize.txt");
	for (int i = 0; i < record_HI.size(); i++)
	{
		myfile << record_time_point[i] << " " << record_N_GE[i] << " " << record_N_GI[i] << " " << record_HE[i] << " " << record_HI[i] << endl;
	}
	myfile.close();*/
	//myfile.open("membrane_potential_sample.txt");
	//for (int i = 0; i < record_HI.size(); i++)
	//{
	//    for (int j = 0; j < NE; j++)
	//    {
	//        myfile << V_e_distribution[i * NE + j] << " ";
	//    }
	//    for (int j = 0; j < NI; j++)
	//    {
	//        myfile << V_i_distribution[i * NI + j] << " ";
	//    }
	//    myfile << endl;
	//}
	//myfile.close();




  //  gettimeofday(&t2, NULL);
  //  double delta = ((t2.tv_sec - t1.tv_sec) * 1000000u +
  //      t2.tv_usec - t1.tv_usec) / 1.e6;

   // cout << "total CPU time = " << delta << endl;

}







