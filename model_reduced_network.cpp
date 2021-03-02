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

int NE = 300;
int NI = 100;
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

vector<double> P_BE_Ex;
P_BE_Ex.reserve(NE + 1);
vector<double> P_GE_Ex;
P_GE_Ex.reserve(NE + 1);
vector<double> P_BI_Ex;
P_BI_Ex.reserve(NI + 1);
vector<double> P_GI_Ex;
P_GI_Ex.reserve(NI + 1);
vector<double> P_BE_E;
P_BE_E.reserve(NE + 1);
vector<double> P_BI_E;
P_BI_E.reserve(NI + 1);
vector<double> P_GE_E;
P_GE_E.reserve(NE + 1);
vector<double> P_GI_E;
P_GI_E.reserve(NI + 1);
vector<double> P_GE_I;
P_GE_I.reserve(NE + 1);
vector<double> P_GI_I;
P_GI_I.reserve(NI + 1);



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

	bool remove_elt(unsigned index)
		//remove element [index]
	{
		if (this->empty() || index < 0 || index >= this->size()) {
			cout << "Empty vector" << endl;
			return false;
		}
		if (index == this->size() - 1) {
			this->pop_back();
		}
		else {
			T elt = (*this)[(*this).size() - 1];
			(*this)[index] = elt;
			(*this).pop_back();
		}
		return true;
	}

	int select(mt19937& mt, uniform_real_distribution<double>& u)
	{
		if (this->size() == 0)
		{
			cout << "Size cannot be zero " << endl;
			return -1;
		}
		else
		{
			int index = floor(u(mt) * this->size());
			T tmp = (*this)[index];
			remove_elt(index);
			return tmp;
		}
	}

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

void spikeE(const int whichHit, Vector<double>& Clock, vector<int>& VE, Vector<int>& HEE, Vector<int>& HIE, mt19937& mt, uniform_real_distribution<double>& u)
{
	E_spike++;
	VE[whichHit] = 0;
	for (int i = 0; i < NE; i++)
	{
		if (u(mt) < PEE)
		{
			HEE.push_back(i);
		}
	}
	for (int i = 0; i < NI; i++)
	{
		if (u(mt) < PIE)
		{
			HIE.push_back(i);
		}
	}
	Clock.switch_element(2, HitEE * HEE.size());
	Clock.switch_element(3, HitIE * HIE.size());

}
void spikeI(const int whichHit, Vector<double>& Clock, vector<int>& VI, Vector<int>& HEI, Vector<int>& HII, mt19937& mt, uniform_real_distribution<double>& u)
{
	I_spike++;
	VI[whichHit] = 0;
	for (int i = 0; i < NE; i++)
	{
		if (u(mt) < PEI )
		{
			HEI.push_back(i);
		}
	}
	for (int i = 0; i < NI; i++)
	{
		if (u(mt) < PII )
		{
			HII.push_back(i);
		}
	}
	Clock.switch_element(4, HitI * HEI.size());
	Clock.switch_element(5, HitI * HII.size());
}

void update(vector<double>& time_spike, vector<int>& num_spike, vector<double>& record_time_point, vector<int>& record_N_GE, vector<int>& record_N_GI, vector<int>& record_HE,
	vector<int>& record_HI, Vector<double>& Clock, Vector<int>& VE, Vector<int>& VI, Vector<int>& HEE, Vector<int>& HEI,
	Vector<int>& HIE, Vector<int>& HII, const double terminate_time, mt19937& mt, uniform_real_distribution<double>& u)
{
	VE.maintain();
	VI.maintain();

	int N_GE = VE.get_sum();
	int N_GI = VI.get_sum();
	double current_time = 0.0;
	double record_time = 0.0;
	int count = 0;
	while (current_time < terminate_time)
	{
		if (current_time - record_time > 0.0005)
		{
			record_time = current_time;
			record_time_point.push_back(current_time);
			record_HI.push_back(HII.size() + HEI.size());
			record_HE.push_back(HEE.size() + HIE.size());

			record_N_GE.push_back(N_GE);
			record_N_GI.push_back(N_GI);

		}
		current_time += -log(1 - u(mt)) / Clock.get_sum();
		int index = find_index(Clock, mt, u);
		int whichHit;
		count++;

		//            cout<<"time "<<current_time <<" index "<<index<<endl;
		switch (index)
		{
		case 0:
			// external E 
			whichHit = floor(u(mt) * NE);
			if (VE[whichHit] == 0)
			{
				if (u(mt) < P_BE_Ex[N_GE])
				{
					VE[whichHit] == 1;
					N_GE += 1;
				}
			}
			else
			{
				if (u(mt) < P_GE_Ex[N_GE])
				{
					spikeE(whichHit, Clock, VE, HEE, HIE, mt, u);
					time_spike.push_back(current_time);
					num_spike.push_back(whichHit);
				}
			}
			break;
		case 1:
			// external I
			whichHit = floor(u(mt) * NI);
			if (VI[whichHit] == 0)
			{
				if (u(mt) < P_BI_Ex[N_GI])
				{
					VE[whichHit] == 1;
					N_GI += 1;
				}
			}
			else
			{
				if (u(mt) < P_GI_Ex[N_GI])
				{
					spikeI(whichHit, Clock, VI, HEI, HII, mt, u);
					time_spike.push_back(current_time);
					num_spike.push_back(whichHit);
				}
			}
			break;
		case 2:

			whichHit = HEE.select(mt, u);
			//                cout<<"ID = "<<whichHit<<" V = "<<VE[whichHit]<<endl<<" status = "<<awakeE[whichHit]<<endl;

			
			if (VE[whichHit] == 0)
			{
				if (u(mt) < P_BE_E[N_GE])
				{
					VE[whichHit] == 1;
					N_GE += 1;
				}
			}
			else
			{
				if (u(mt) < P_GE_E[N_GE])
				{
					spikeE(whichHit, Clock, VE, HEE, HIE, mt, u);
					time_spike.push_back(current_time);
					num_spike.push_back(whichHit);
				}
			}

			
			

			//                cout<<" after "<<VE[whichHit]<<endl;
			Clock.switch_element(2, HitEE * HEE.size());
			break;
		case 3:
			whichHit = HIE.select(mt, u);
			//                cout<<HEE.size()<<endl;
			
			if (VI[whichHit] == 0)
			{
				if (u(mt) < P_BI_E[N_GI])
				{
					VE[whichHit] == 1;
					N_GI += 1;
				}
			}
			else
			{
				if (u(mt) < P_GI_E[N_GI])
				{
					spikeI(whichHit, Clock, VI, HEI, HII, mt, u);
					time_spike.push_back(current_time);
					num_spike.push_back(whichHit);
				}
			}
			
			Clock.switch_element(3, HitIE * HIE.size());
			break;
		case 4:
			// I pending spikes
			whichHit = HEI.select(mt, u);
			if (VE[whichHit] == 1)
			{
				if (u(mt) < P_GE_I[N_GE])
				{
					VE[whichHit] == 0;
					N_GE -= 1;
				}
			}
			Clock.switch_element(4, HitI * HEI.size());
			break;
		case 5:
			whichHit = HII.select(mt, u);
			if (VI[whichHit] == 1)
			{
				if (u(mt) < P_GI_I[N_GI])
				{
					VE[whichHit] == 0;
					N_GI -= 1;
				}
			}
			Clock.switch_element(5, HitI * HII.size());
			break;
		}
	}
}



int main()
{
	//  struct timeval t1, t2;
	//  gettimeofday(&t1, NULL);
	ofstream myfile;
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> u(0, 1);
	myfile.open("spike_info_reduced_largesize.txt");
	Vector<int> VE;//membrane potential: 0 represents base, 1 represents gate.
	Vector<int> VI;
	VE.reserve(NE);
	VI.reserve(NI);
	for (int i = 0; i < NE;i++)
		VE.push_back(0);
	for (int i = 0; i < NI;i++)
		VI.push_back(0);
	Vector<int> HEE, HEI, HII, HIE;
	//HEE, HEI ... pending postsynaptic kick
	//Eref, Iref: indices of neurons at state R
	HEE.reserve(100000);
	HEI.reserve(100000);
	HII.reserve(100000);
	HIE.reserve(100000);

	Vector<double> Clock;
	Clock.reserve(6);
	// 0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII

	Clock.push_back(NE * kickE);
	Clock.push_back(NI * kickI);
	for (int i = 2; i < 6; i++)
		Clock.push_back(0);
	Clock.maintain();
	double terminate_time = 40;
	vector<double> time_spike;
	time_spike.reserve(100000);
	vector<int> num_spike;
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



	update(time_spike, num_spike, record_time_point, record_N_GE, record_N_GI, record_HE, record_HI, Clock, VE, VI, HEE, HEI, HIE, HII, terminate_time, mt, u);

	cout << "E spike rate= " << (double)E_spike / (terminate_time * NE) << endl;
	cout << "I spike rate = " << (double)I_spike / (terminate_time * NI) << endl;
	int spike_count = time_spike.size();
	for (int i = 0; i < spike_count; i++)
	{
		myfile << time_spike[i] << "  " << num_spike[i] << endl;
	}
	myfile.close();
	myfile.open("trajectory_info_sample_largesize.txt");
	for (int i = 0; i < record_HI.size(); i++)
	{
		myfile << record_time_point[i] << " " << record_N_GE[i] << " " << record_N_GI[i] << " " << record_HE[i] << " " << record_HI[i] << endl;
	}
	myfile.close();
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







