#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
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

//int NE = 300;
//int NI = 100;
//int Level = 100;//membrane potential
//double PEE = 0.15;//probability of postsynaptic connections
//double PIE = 0.5;
//double PEI = 0.5;
//double PII = 0.4;
//double SEE = 5;//strength of postsynaptic connection
//double SIE = 3;
//double SEI = -2;
//double SII = -2;
//double kickE = 2000.0;//external drive rate
//double kickI = 2000.0;
//double Ref = 250.0;//time rate at state R
//double HitE = 1000.0 / 2.0;//delay time rate
//double HitI = 1000.0 / 4.0;//
//int Reverse = -66;//reverse potential
//int E_spike = 0;
//int I_spike = 0;
//int gate = 74;

double factor = 1.5;
int NE = 300;
int NI = 70;
int NS = 30;
double SEE = 5;//strength of postsynaptic connection
double SIE = 2;
double SSE = 2;
double SEI = -4.91;
double SII = -4.91;
double SSI = -4.91;
double SES = -4.91/factor;
double SIS = -4.91;
double SSS = 0;

int Level = 100;//membrane potential
double PEE = 0.15;//probability of postsynaptic connections
double PIE = 0.5;
double PSE = 0.5;
double PEI = 0.5;
double PII = 0.4;
double PSI = 0.4;
double PES = 0.5*factor;
double PIS = 0.4;
double PSS = 0.4;

double kickE = 7000.0;//external drive rate
double kickI = 7000.0;
double kickS = 700;
double Ref = 250.0;//time rate at state R
double HitEE = 1000.0 / 1.4;//delay time rate
double HitIE = 1000.0 / 1.2;
double HitSE = 1000.0 / 1.2;
double HitI = 1000.0 / 4.5;
double HitES = 1000.0 / 20;
double HitIS = 1000.0 / 4.5;
double HitSS = 1000.0 / 0.1;

double LeakE = 10*1000.0 / 20;
double LeakI = 10*1000.0 / 16.7;
double LeakS = 10*1000.0 / 20;
double gLeak = -1+pow((1/2.71828),0.1);  
int Reverse = -66;//reverse potential

int E_spike = 0;
int I_spike = 0;
int S_spike = 0;

int real2int(const double x, mt19937& mt, uniform_real_distribution<double>& u)
//choose a random integer n if x is not an integer, such that the expectation of n is equal to x
{
    int xf = floor(x);
    double q = x - (double)xf;
    double y = 0;
    if (u(mt) < q)
        y = xf + 1;
    else
        y = xf;
    return y;
}

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

    int find_min() {
        if (!this->empty()) {
            auto min = (*this)[0];
            auto index = 0;
            auto iter = 0;
            for (auto&& item : (*this)) {
                if (item < min) {
                    min = item;
                    index = iter;
                }
                iter++;
            }
            return index;
        }
        return -1;
    }

    void print_vector() {
        for (auto&& it : (*this)) {
            cout << it << " ";
        }
        cout << endl;
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



void spikeE(const int whichHit, Vector<double>& Clock, vector<int>& VE, Vector<int>& HEE, Vector<int>& HIE, Vector<int>& HSE, Vector<int>& Eref,
    vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS, mt19937& mt, uniform_real_distribution<double>& u)
{
    E_spike++;
    VE[whichHit] = 0;
    awakeE[whichHit] = 0;
    Eref.push_back(whichHit);
    Clock.switch_element(15, Ref * Eref.size());
    for (int i = 0; i < NE; i++)
    {
        if (u(mt) < PEE && awakeE[i])
        {
            HEE.push_back(i);
        }
    }
    for (int i = 0; i < NI; i++)
    {
        if (u(mt) < PIE && awakeI[i])
        {
            HIE.push_back(i);
        }
    }
    for (int i = 0; i < NS; i++)
    {
        if (u(mt) < PSE && awakeI[i])
        {
            HSE.push_back(i);
        }
    }
    Clock.switch_element(6, HitEE * HEE.size());
    Clock.switch_element(7, HitIE * HIE.size());
    Clock.switch_element(8, HitSE * HSE.size());

}

void spikeI(const int whichHit, Vector<double>& Clock, vector<int>& VI, Vector<int>& HEI, Vector<int>& HII, Vector<int>& HSI, Vector<int>& Iref,
    vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS, mt19937& mt, uniform_real_distribution<double>& u)
{
    I_spike++;
    VI[whichHit] = 0;
    awakeI[whichHit] = 0;
    Iref.push_back(whichHit);
    Clock.switch_element(16, Ref * Iref.size());
    for (int i = 0; i < NE; i++)
    {
        if (u(mt) < PEI && awakeE[i])
        {
            HEI.push_back(i);
        }
    }
    for (int i = 0; i < NI; i++)
    {
        if (u(mt) < PII && awakeI[i])
        {
            HII.push_back(i);
        }
    }
    for (int i = 0; i < NS; i++)
    {
        if (u(mt) < PSI && awakeS[i])
        {
            HSI.push_back(i);
        }
    }
    Clock.switch_element(9, HitI * HEI.size());
    Clock.switch_element(10, HitI * HII.size());
    Clock.switch_element(11, HitI * HSI.size());
}

void spikeS(const int whichHit, Vector<double>& Clock, vector<int>& VS, Vector<int>& HES, Vector<int>& HIS, Vector<int>& HSS, Vector<int>& Sref,
    vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS, mt19937& mt, uniform_real_distribution<double>& u)
{
    S_spike++;
    VS[whichHit] = 0;
    awakeS[whichHit] = 0;
    Sref.push_back(whichHit);
    Clock.switch_element(17, Ref * Sref.size());
    for (int i = 0; i < NE; i++)
    {
        if (u(mt) < PES && awakeE[i])
        {
            HES.push_back(i);
        }
    }
    for (int i = 0; i < NI; i++)
    {
        if (u(mt) < PIS && awakeI[i])
        {
            HIS.push_back(i);
        }
    }
    for (int i = 0; i < NS; i++)
    {
        if (u(mt) < PSS && awakeS[i])
        {
            HSS.push_back(i);
        }
    }
    Clock.switch_element(12, HitES * HES.size());
    Clock.switch_element(13, HitIS * HIS.size());
    Clock.switch_element(14, HitSS * HSS.size());
}




void update(vector<double>& time_spike, vector<int>& num_spike, vector<double>& record_time_point, vector<int>& N_GE, vector<int>& N_GI, vector<int>& total_HE,
    vector<int>& total_HI, vector<int>& V_e_distribution, vector<int>& V_i_distribution, vector<int>& V_s_distribution, Vector<double>& Clock, vector<int>& VE,
    vector<int>& VI, vector<int>& VS, Vector<int>& HEE, Vector<int>& HEI, Vector<int>& HES, Vector<int>& HIE, Vector<int>& HII, Vector<int>& HIS,
    Vector<int>& HSE, Vector<int>& HSI, Vector<int>& HSS,
    Vector<int>& Eref, Vector<int>& Iref, Vector<int>& Sref, vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS,
    const double terminate_time, mt19937& mt, uniform_real_distribution<double>& u)
{
    double current_time = 0.0;
    double record_time = 0.0;
    double Leak = 0;
    int count = 0;
    int whichHit;
    while (current_time < terminate_time)
    {
        if (current_time - record_time > 0.0005)
        {
            record_time = current_time;
            record_time_point.push_back(current_time);
            total_HI.push_back(HII.size() + HEI.size());
            total_HE.push_back(HEE.size() + HIE.size());
            for (int i = 0; i < NE; i++)
            {
                V_e_distribution.push_back(VE[i] + (1 - awakeE[i]) * (Reverse - 1));
            }
            for (int i = 0; i < NI; i++)
            {
                V_i_distribution.push_back(VI[i] + (1 - awakeI[i]) * (Reverse - 1));
            }
            for (int i = 0; i < NS; i++)
            {
                V_s_distribution.push_back(VS[i] + (1 - awakeS[i]) * (Reverse - 1));
            }
        }
        current_time += -log(1 - u(mt)) / Clock.get_sum();
        int index = find_index(Clock, mt, u);
        count++;

        //            cout<<"time "<<current_time <<" index "<<index<<endl;
        switch (index)
        {
        case 0:
            whichHit = floor(u(mt) * NE);
            if (awakeE[whichHit])
            {
                VE[whichHit]++;
                if (VE[whichHit] >= Level)
                {
                    spikeE(whichHit, Clock, VE, HEE, HIE, HSE, Eref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit);
                }
            }
            break;
        case 1:
            whichHit = floor(u(mt) * NI);
            if (awakeI[whichHit])
            {
                VI[whichHit]++;
                if (VI[whichHit] >= Level)
                {
                    spikeI(whichHit, Clock, VI, HEI, HII, HSI, Iref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit + NE);
                }
            }
            break;
        case 2:
            whichHit = floor(u(mt) * NS);
            if (awakeS[whichHit])
            {
                VS[whichHit]++;
                if (VS[whichHit] >= Level)
                {
                    spikeS(whichHit, Clock, VS, HES, HIS, HSS, Sref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit + NE + NI);
                }
            }
            break;
            //0, Edrive, 1, Idrive, 2, Sdrive, 3, LeakE, 4, LeakI, 5, LeakS, 6-14, H.., 15, Eref, 16, Iref, 17, Sref
        case 3:
            whichHit = floor(u(mt) * NE);
            //                cout<<"ID = "<<whichHit<<" V = "<<VE[whichHit]<<endl<<" status = "<<awakeE[whichHit]<<endl;
            if (awakeE[whichHit])
            {
                Leak = gLeak * VE[whichHit];
                VE[whichHit] += real2int(Leak, mt, u);
            }
            //                cout<<" after "<<VE[whichHit]<<endl;
            break;
        case 4:
            whichHit = floor(u(mt) * NI);
            //                cout<<HEE.size()<<endl;
            if (awakeI[whichHit])
            {
                Leak = gLeak * VI[whichHit];
                VI[whichHit] += real2int(Leak, mt, u);
            }
            break;
        case 5:
            whichHit = floor(u(mt) * NS);
            if (awakeS[whichHit])
            {
                Leak = gLeak * VS[whichHit];
                VS[whichHit] += real2int(Leak, mt, u);
            }
            break;
        case 6:
            whichHit = HEE.select(mt, u);
            if (awakeE[whichHit])
            {
                VE[whichHit] += real2int(SEE, mt, u);
                if (VE[whichHit] >= Level)
                {
                    spikeE(whichHit, Clock, VE, HEE, HIE, HSE, Eref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit);
                }
            }
            Clock.switch_element(6, HitEE * HEE.size());
            break;
        case 7:
            whichHit = HIE.select(mt, u);
            if (awakeI[whichHit])
            {
                VI[whichHit] += real2int(SIE, mt, u);
                if (VI[whichHit] >= Level)
                {
                    spikeI(whichHit, Clock, VI, HEI, HII, HSI, Iref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit + NE);
                }
            }
            Clock.switch_element(7, HitIE * HIE.size());
            break;
        case 8:
            whichHit = HSE.select(mt, u);
            if (awakeS[whichHit])
            {
                VS[whichHit] += real2int(SSE, mt, u);
                if (VS[whichHit] >= Level)
                {
                    spikeS(whichHit, Clock, VS, HES, HIS, HSS, Sref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit + NE + NI);
                }
            }
            Clock.switch_element(8, HitSE * HSE.size());
            break;
        case 9:
            whichHit = HEI.select(mt, u);
            if (awakeE[whichHit])
            {
                VE[whichHit] += real2int(SEI * (VE[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VE[whichHit] < Reverse)
                    VE[whichHit] = Reverse;
            }
            Clock.switch_element(9, HitI * HEI.size());
            break;
        case 10:
            whichHit = HII.select(mt, u);
            if (awakeI[whichHit])
            {
                VI[whichHit] += real2int(SII * (VI[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VI[whichHit] < Reverse)
                    VI[whichHit] = Reverse;
            }
            Clock.switch_element(10, HitI * HII.size());
            break;
        case 11:
            whichHit = HSI.select(mt, u);
            if (awakeS[whichHit])
            {
                VS[whichHit] += real2int(SSI * (VS[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VS[whichHit] < Reverse)
                    VS[whichHit] = Reverse;
            }
            Clock.switch_element(11, HitI * HSI.size());
            break;
        case 12:
            whichHit = HES.select(mt, u);
            if (awakeE[whichHit])
            {
                VE[whichHit] += real2int(SES * (VE[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VE[whichHit] < Reverse)
                    VE[whichHit] = Reverse;
            }
            Clock.switch_element(12, HitES * HES.size());
            break;
        case 13:
            whichHit = HIS.select(mt, u);
            if (awakeI[whichHit])
            {
                VI[whichHit] += real2int(SIS * (VI[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VI[whichHit] < Reverse)
                    VI[whichHit] = Reverse;
            }
            Clock.switch_element(13, HitIS * HIS.size());
            break;
        case 14:
            whichHit = HSS.select(mt, u);
            if (awakeS[whichHit])
            {
                VS[whichHit] += real2int(SSS * (VS[whichHit] - Reverse) / (Level - Reverse), mt, u);
                if (VS[whichHit] < Reverse)
                    VS[whichHit] = Reverse;
            }
            Clock.switch_element(14, HitSS * HSS.size());
            break;
        case 15:
            whichHit = Eref.select(mt, u);
            awakeE[whichHit] = 1;
            Clock.switch_element(15, Ref * Eref.size());
            break;
        case 16:
            whichHit = Iref.select(mt, u);
            awakeI[whichHit] = 1;
            Clock.switch_element(16, Ref * Iref.size());
            break;
        case 17:
            whichHit = Sref.select(mt, u);
            awakeS[whichHit] = 1;
            Clock.switch_element(17, Ref * Sref.size());
            break;
        }
    }
}

int main()
{
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    ofstream myfile;

    mt19937 mt(time(NULL));
    uniform_real_distribution<double> u(0, 1);

    vector<int> VE(NE);//membrane potential
    vector<int> VI(NI);
    vector<int> VS(NS);
    for (auto i : VE)
        i = 0;
    for (auto i : VI)
        i = 0;
    for (auto i : VS)
        i = 0;
    Vector<int> HEE, HEI, HES, HII, HIE, HIS, HSE, HSI, HSS, Eref, Iref, Sref;
    //HEE, HEI ... pending postsynaptic kick
    //Eref, Iref: indices of neurons at state R
    HEE.reserve(100000);
    HEI.reserve(100000);
    HES.reserve(100000);
    HIE.reserve(100000);
    HII.reserve(100000);
    HIS.reserve(100000);
    HSE.reserve(100000);
    HSI.reserve(100000);
    HSS.reserve(100000);

    Eref.reserve(NE);
    Iref.reserve(NI);
    Sref.reserve(NS);
    vector<int> awakeE(NE);
    for (auto& i : awakeE)
        i = 1;
    vector<int> awakeI(NI);
    for (auto& i : awakeI)
        i = 1;
    vector<int> awakeS(NS);
    for (auto& i : awakeS)
        i = 1;
    Vector<double> Clock;
    Clock.reserve(18);
    //0, Edrive, 1, Idrive, 2, Sdrive, 3, LeakE, 4, LeakI, 5, LeakS, 6-14, H.., 15, Eref, 16, Iref, 17, Sref
    //6-14: [EE,IE,SE,EI,II,SI,ES,IS,SS]

    Clock.push_back(NE * kickE);
    Clock.push_back(NI * kickI);
    Clock.push_back(NS * kickS);
    Clock.push_back(NE * LeakE);
    Clock.push_back(NI * LeakI);
    Clock.push_back(NS * LeakS);

    for (int i = 6; i < 18; i++)
        Clock.push_back(0);
    Clock.maintain();
    cout << "start" << endl;
    double terminate_time = 10.0;
    vector<double> time_spike;
    time_spike.reserve(100000);
    vector<int> num_spike;
    num_spike.reserve(100000);
    vector<double> record_time_point;
    record_time_point.reserve(250000);


    vector<int> N_GE;
    N_GE.reserve(250000);
    vector<int> N_GI;
    N_GI.reserve(250000);


    vector<int> total_HE;
    total_HE.reserve(250000);
    vector<int> total_HI;
    total_HI.reserve(250000);


    vector<int> V_e_distribution;
    V_e_distribution.reserve(8000000);
    vector<int> V_i_distribution;
    V_i_distribution.reserve(8000000);
    vector<int> V_s_distribution;
    V_s_distribution.reserve(8000000);





    update(time_spike, num_spike, record_time_point, N_GE, N_GI, total_HE, total_HI, V_e_distribution, V_i_distribution, V_s_distribution, Clock, VE, VI, VS,
        HEE, HEI, HES, HIE, HII, HIS, HSE, HSI, HSS, Eref, Iref, Sref, awakeE, awakeI, awakeS, terminate_time, mt, u);

    cout << "E spike rate= " << (double)E_spike / (terminate_time * NE) << endl;
    cout << "I spike rate = " << (double)I_spike / (terminate_time * NI) << endl;
    cout << "S spike rate = " << (double)S_spike / (terminate_time * NS) << endl;
    int spike_count = time_spike.size();

    myfile.open("spike_info_ref.txt");
    for (int i = 0; i < spike_count; i++)
    {
        myfile << time_spike[i] << "  " << num_spike[i] << endl;
    }
    myfile.close();

    myfile.open("trajectory_info_ref.txt");
    for (int i = 0; i < total_HI.size(); i++)
    {
        myfile << record_time_point[i] << " " << N_GE[i] << " " << N_GI[i] << " " << total_HE[i] << " " << total_HI[i] << endl;
    }
    myfile.close();
    myfile.open("membrane_potential_sample_ref.txt");
    for (int i = 0; i < total_HI.size(); i++)
    {
        for (int j = 0; j < NE; j++)
        {
            myfile << V_e_distribution[i * NE + j] << " ";
        }
        for (int j = 0; j < NI; j++)
        {
            myfile << V_i_distribution[i * NI + j] << " ";
        }
        for (int j = 0; j < NS; j++)
        {
            myfile << V_s_distribution[i * NS + j] << " ";
        }
        myfile << endl;
    }
    myfile.close();

    gettimeofday(&t2, NULL);
    double delta = ((t2.tv_sec - t1.tv_sec) * 1000000u +
        t2.tv_usec - t1.tv_usec) / 1.e6;

    cout << "total CPU time = " << delta << endl;

}







