#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>

#include <string>
#include <cstdlib>
#include <float.h>
#include <random>
#include <stdint.h>
#include <algorithm>
#include <cstdio>
#include <climits>
#include <vector>
#include <direct.h>
#include <io.h>

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

double factor;
int NE;
int NI;
int NS;
double SEE;//strength of postsynaptic connection
double SIE;
double SSE;
double SEEN;//strength of postsynaptic connection
double SIEN;
double SSEN;
double PEEN;
double PIEN;
double PSEN;
double SEI;
double SII;
double SSI;
double SES;
double SIS;
double SSS;

int Level;//membrane potential
double PEE;//probability of postsynaptic connections
double PIE;
double PSE;
double PEI;
double PII;
double PSI;
double PES;
double PIS;
double PSS;

double kickE;//external drive rate
double kickI;
double kickS;
double Ref;//time rate at state R
double HitEE;//delay time rate
double HitIE;
double HitSE;
double HitI;
double HitES;
double HitIS;
double HitSS;
double HitEEN;
double HitIEN;
double HitSEN;

double factor2;
double LeakE;
double LeakI;
double LeakS;
double gLeak;
double terminate_time;
int Reverse;//reverse potential

// double factor = 1.5;
// int NE = 300;
// int NI = 70;
// int NS = 30;
// double SEE = 5;//strength of postsynaptic connection
// double SIE = 2;
// double SSE = 2;
// double SEI = -4.91;
// double SII = -4.91;
// double SSI = -4.91;
// double SES = -4.91/factor;
// double SIS = -4.91;
// double SSS = 0;

// int Level = 100;//membrane potential
// double PEE = 0.15;//probability of postsynaptic connections
// double PIE = 0.5;
// double PSE = 0.5;
// double PEI = 0.5;
// double PII = 0.4;
// double PSI = 0.4;
// double PES = 0.5*factor;
// double PIS = 0.4;
// double PSS = 0.4;

// double kickE = 7000.0;//external drive rate
// double kickI = 7000.0;
// double kickS = 700;
// double Ref = 250.0;//time rate at state R
// double HitEE = 1000.0 / 1.4;//delay time rate
// double HitIE = 1000.0 / 1.2;
// double HitSE = 1000.0 / 1.2;
// double HitI = 1000.0 / 4.5;
// double HitES = 1000.0 / 20;
// double HitIS = 1000.0 / 4.5;
// double HitSS = 1000.0 / 0.1;

// double LeakE = 10*1000.0 / 20;
// double LeakI = 10*1000.0 / 16.7;
// double LeakS = 10*1000.0 / 20;
// double gLeak = -1+pow((1/2.71828),0.1);  
// int Reverse = -66;//reverse potential

int E_spike = 0;
int I_spike = 0;
int S_spike = 0;

string StringToNumstr(string s)
{
    int i = 0;
    while ((s[i] < '0') || (s[i] > '9')) i++;
    s = s.substr(i, s.length() - i);
    return s;
}

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



void spikeE(const int whichHit, Vector<double>& Clock, vector<int>& VE, Vector<int>& HEE, Vector<int>& HIE, Vector<int>& HSE, Vector<int>& HEEN, Vector<int>& HIEN,
    Vector<int>& HSEN, Vector<int>& Eref, vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS, mt19937& mt, uniform_real_distribution<double>& u)
{
    E_spike++;
    VE[whichHit] = 0;
    awakeE[whichHit] = 0;
    Eref.push_back(whichHit);
    Clock.switch_element(18, Ref * Eref.size());
    for (int i = 0; i < NE; i++)
    {
        if (u(mt) < PEE && awakeE[i])
        {
            if (u(mt) < PEEN)
                HEEN.push_back(i);
            else
                HEE.push_back(i);
        }
    }
    for (int i = 0; i < NI; i++)
    {
        if (u(mt) < PIE && awakeI[i])
        {
            if (u(mt) < PIEN)
                HIEN.push_back(i);
            else
                HIE.push_back(i);
        }
    }
    for (int i = 0; i < NS; i++)
    {
        if (u(mt) < PSE && awakeI[i])
        {
            if (u(mt) < PSEN)
                HSEN.push_back(i);
            else
                HSE.push_back(i);
        }
    }
    Clock.switch_element(6, HitEE * HEE.size());
    Clock.switch_element(7, HitIE * HIE.size());
    Clock.switch_element(8, HitSE * HSE.size());
    Clock.switch_element(15, HitEEN * HEEN.size());
    Clock.switch_element(16, HitIEN * HIEN.size());
    Clock.switch_element(17, HitSEN * HSEN.size());

}

void spikeI(const int whichHit, Vector<double>& Clock, vector<int>& VI, Vector<int>& HEI, Vector<int>& HII, Vector<int>& HSI, Vector<int>& Iref,
    vector<int>& awakeE, vector<int>& awakeI, vector<int>& awakeS, mt19937& mt, uniform_real_distribution<double>& u)
{
    I_spike++;
    VI[whichHit] = 0;
    awakeI[whichHit] = 0;
    Iref.push_back(whichHit);
    Clock.switch_element(19, Ref * Iref.size());
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
    Clock.switch_element(20, Ref * Sref.size());
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




void update(vector<double>& time_spike, vector<int>& num_spike, vector<double>& record_time_point,
    vector<int>& total_HEE, vector<int>& total_HIE, vector<int>& total_HSE,
    vector<int>& total_HEEN, vector<int>& total_HIEN, vector<int>& total_HSEN,
    vector<int>& total_HEI, vector<int>& total_HII, vector<int>& total_HSI,
    vector<int>& total_HES, vector<int>& total_HIS, vector<int>& total_HSS,
    vector<int>& V_e_distribution, vector<int>& V_i_distribution, vector<int>& V_s_distribution,
    Vector<double>& Clock, vector<int>& VE, vector<int>& VI, vector<int>& VS,
    Vector<int>& HEE, Vector<int>& HEI, Vector<int>& HES,
    Vector<int>& HEEN, Vector<int>& HIEN, Vector<int>& HSEN,
    Vector<int>& HIE, Vector<int>& HII, Vector<int>& HIS,
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
            total_HIE.push_back(HIE.size());
            total_HEE.push_back(HEE.size());
            total_HSE.push_back(HSE.size());
            total_HIEN.push_back(HIEN.size());
            total_HEEN.push_back(HEEN.size());
            total_HSEN.push_back(HSEN.size());
            total_HII.push_back(HII.size());
            total_HEI.push_back(HEI.size());
            total_HSI.push_back(HSI.size());
            total_HIS.push_back(HIS.size());
            total_HES.push_back(HES.size());
            total_HSS.push_back(HSS.size());
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
                    spikeE(whichHit, Clock, VE, HEE, HIE, HSE, HEEN, HIEN, HSEN, Eref, awakeE, awakeI, awakeS, mt, u);
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
                    spikeE(whichHit, Clock, VE, HEE, HIE, HSE, HEEN, HIEN, HSEN, Eref, awakeE, awakeI, awakeS, mt, u);
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
            whichHit = HEEN.select(mt, u);
            if (awakeE[whichHit])
            {
                VE[whichHit] += real2int(SEEN, mt, u);
                if (VE[whichHit] >= Level)
                {
                    spikeE(whichHit, Clock, VE, HEE, HIE, HSE, HEEN, HIEN, HSEN, Eref, awakeE, awakeI, awakeS, mt, u);
                    time_spike.push_back(current_time);
                    num_spike.push_back(whichHit);
                }
            }
            Clock.switch_element(15, HitEEN * HEEN.size());
            break;
        case 16:
            whichHit = HIEN.select(mt, u);
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
            Clock.switch_element(16, HitIEN * HIEN.size());
            break;
        case 17:
            whichHit = HSEN.select(mt, u);
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
            Clock.switch_element(17, HitSEN * HSEN.size());
            break;
        case 18:
            whichHit = Eref.select(mt, u);
            awakeE[whichHit] = 1;
            Clock.switch_element(18, Ref * Eref.size());
            break;
        case 19:
            whichHit = Iref.select(mt, u);
            awakeI[whichHit] = 1;
            Clock.switch_element(19, Ref * Iref.size());
            break;
        case 20:
            whichHit = Sref.select(mt, u);
            awakeS[whichHit] = 1;
            Clock.switch_element(20, Ref * Sref.size());
            break;
        }
    }
}

int main()
{
    ifstream inf;
    inf.open("..//model_SLN_full_params.txt");
    string s;
    getline(inf, s);
    factor = stod(StringToNumstr(s));
    getline(inf, s);
    NE = stoi(StringToNumstr(s));
    getline(inf, s);
    NI = stoi(StringToNumstr(s));
    getline(inf, s);
    NS = stoi(StringToNumstr(s));
    getline(inf, s);
    SEE = stod(StringToNumstr(s));//strength of postsynaptic connection
    getline(inf, s);
    SIE = stod(StringToNumstr(s));
    getline(inf, s);
    SSE = stod(StringToNumstr(s));
    getline(inf, s);
    SEEN = stod(StringToNumstr(s));//strength of postsynaptic connection
    getline(inf, s);
    SIEN = stod(StringToNumstr(s));
    getline(inf, s);
    SSEN = stod(StringToNumstr(s));
    getline(inf, s);
    PEEN = stod(StringToNumstr(s));//strength of postsynaptic connection
    getline(inf, s);
    PIEN = stod(StringToNumstr(s));
    getline(inf, s);
    PSEN = stod(StringToNumstr(s));
    getline(inf, s);
    SEI = -stod(StringToNumstr(s));
    getline(inf, s);
    SII = -stod(StringToNumstr(s));
    getline(inf, s);
    SSI = -stod(StringToNumstr(s));
    getline(inf, s);
    SES = -stod(StringToNumstr(s)) / factor;
    getline(inf, s);
    SIS = -stod(StringToNumstr(s));
    getline(inf, s);
    SSS = stod(StringToNumstr(s));

    getline(inf, s);
    Level = stoi(StringToNumstr(s));//membrane potential
    getline(inf, s);
    PEE = stod(StringToNumstr(s));//probability of postsynaptic connections
    getline(inf, s);
    PIE = stod(StringToNumstr(s));
    getline(inf, s);
    PSE = stod(StringToNumstr(s));
    getline(inf, s);
    PEI = stod(StringToNumstr(s));
    getline(inf, s);
    PII = stod(StringToNumstr(s));
    getline(inf, s);
    PSI = stod(StringToNumstr(s));
    getline(inf, s);
    PES = stod(StringToNumstr(s)) / factor;
    getline(inf, s);
    PIS = stod(StringToNumstr(s));
    getline(inf, s);
    PSS = stod(StringToNumstr(s));

    getline(inf, s);
    kickE = stod(StringToNumstr(s));//external drive rate
    getline(inf, s);
    kickI = stod(StringToNumstr(s));
    getline(inf, s);
    kickS = stod(StringToNumstr(s));
    getline(inf, s);
    Ref = stod(StringToNumstr(s));//time rate at state R
    getline(inf, s);
    HitEE = 1000.0 / stod(StringToNumstr(s));//delay time rate
    getline(inf, s);
    HitIE = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitSE = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitI = 1000.0 / stod(StringToNumstr(s));//
    getline(inf, s);
    HitES = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitIS = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitSS = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitEEN = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitIEN = 1000.0 / stod(StringToNumstr(s));
    getline(inf, s);
    HitSEN = 1000.0 / stod(StringToNumstr(s));

    getline(inf, s);
    factor2 = stod(StringToNumstr(s));
    getline(inf, s);
    LeakE = 1000.0 * factor2 / stod(StringToNumstr(s));
    getline(inf, s);
    LeakI = 1000.0 * factor2 / stod(StringToNumstr(s));
    getline(inf, s);
    LeakS = 1000.0 * factor2 / stod(StringToNumstr(s));
    getline(inf, s);
    gLeak = -1 + pow((1 / 2.71828), 1 / factor2);
    Reverse = -stoi(StringToNumstr(s));//reverse potential
    getline(inf, s);
    terminate_time = stod(StringToNumstr(s));

    // std::cout << NE << endl;
    // std::cout << NI << endl;
    // std::cout << NS << endl;
    // std::cout << SEE << endl;
    // std::cout << SIE << endl;
    // std::cout << SSE << endl;
    // std::cout << SEI << endl;
    // std::cout << SII << endl;
    // std::cout << SSI << endl;
    // std::cout << SES << endl;
    // std::cout << SIS << endl;
    // std::cout << SSS << endl;
    // std::cout << SEEN << endl;
    // std::cout << SIEN << endl;
    // std::cout << SSEN << endl;
    // std::cout << Level << endl;
    // std::cout << PEE << endl;
    // std::cout << PIE << endl;
    // std::cout << PEI << endl;
    // std::cout << PII << endl;
    // std::cout << kickE << endl;
    // std::cout << kickI << endl;
    // std::cout << Ref << endl;
    // std::cout << HitEE << endl;
    // std::cout << HitIE << endl;
    // std::cout << HitSE << endl;
    // std::cout << HitI << endl;
    // std::cout << HitES << endl;
    // std::cout << HitIS << endl;
    // std::cout << HitSS << endl;
    // std::cout << HitEEN << endl;
    // std::cout << HitIEN << endl;
    // std::cout << HitSEN << endl;
    // std::cout << Reverse << endl;
    // std::cout << terminate_time << endl;
    // std::cin.get();
   

    //struct timeval t1, t2;
    //gettimeofday(&t1, NULL);
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
    Vector<int> HEE, HEI, HES, HII, HIE, HIS, HSE, HSI, HSS, HEEN, HIEN, HSEN, Eref, Iref, Sref;
    //HEE, HEI ... pending postsynaptic kick
    //Eref, Iref: indices of neurons at state R
    HEE.reserve(100000);
    HEI.reserve(100000);
    HES.reserve(100000);
    HEEN.reserve(100000);
    HIEN.reserve(100000);
    HSEN.reserve(100000);
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
    Clock.reserve(21);
    //0, Edrive, 1, Idrive, 2, Sdrive, 3, LeakE, 4, LeakI, 5, LeakS, 6-14, H.., 18, Eref, 19, Iref, 20, Sref
    //6-17: [EE,IE,SE,EI,II,SI,ES,IS,SS,EEN,IEN,SEN]

    Clock.push_back(NE * kickE);
    Clock.push_back(NI * kickI);
    Clock.push_back(NS * kickS);
    Clock.push_back(NE * LeakE);
    Clock.push_back(NI * LeakI);
    Clock.push_back(NS * LeakS);

    for (int i = 6; i < 21; i++)
        Clock.push_back(0);
    Clock.maintain();
    cout << "start" << endl;
    vector<double> time_spike;
    time_spike.reserve(100000);
    vector<int> num_spike;
    num_spike.reserve(100000);
    vector<double> record_time_point;
    record_time_point.reserve(250000);


    vector<int> total_HEE;
    total_HEE.reserve(250000);
    vector<int> total_HIE;
    total_HIE.reserve(250000);
    vector<int> total_HSE;
    total_HSE.reserve(250000);
    vector<int> total_HEEN;
    total_HEEN.reserve(250000);
    vector<int> total_HIEN;
    total_HIEN.reserve(250000);
    vector<int> total_HSEN;
    total_HSEN.reserve(250000);
    vector<int> total_HEI;
    total_HEI.reserve(250000);
    vector<int> total_HII;
    total_HII.reserve(250000);
    vector<int> total_HSI;
    total_HSI.reserve(250000);
    vector<int> total_HES;
    total_HES.reserve(250000);
    vector<int> total_HIS;
    total_HIS.reserve(250000);
    vector<int> total_HSS;
    total_HSS.reserve(250000);


    vector<int> V_e_distribution;
    V_e_distribution.reserve(8000000);
    vector<int> V_i_distribution;
    V_i_distribution.reserve(8000000);
    vector<int> V_s_distribution;
    V_s_distribution.reserve(8000000);





    update(time_spike, num_spike, record_time_point,
        total_HEE, total_HIE, total_HSE, total_HEEN, total_HIEN, total_HSEN, total_HEI, total_HII, total_HSI, total_HES, total_HIS, total_HSS,
        V_e_distribution, V_i_distribution, V_s_distribution,
        Clock, VE, VI, VS,
        HEE, HEI, HES, HEEN, HIEN, HSEN, HIE, HII, HIS, HSE, HSI, HSS,
        Eref, Iref, Sref, awakeE, awakeI, awakeS, terminate_time, mt, u);

    cout << "E spike rate= " << (double)E_spike / (terminate_time * NE) << endl;
    cout << "I spike rate = " << (double)I_spike / (terminate_time * NI) << endl;
    cout << "S spike rate = " << (double)S_spike / (terminate_time * NS) << endl;
    int spike_count = time_spike.size();
    std::string save_path = "..//outputs";
    if (_access(save_path.c_str(), 0) == -1)
        _mkdir(save_path.c_str());
    save_path = "..//outputs//model_SLN_full";
    if (_access(save_path.c_str(), 0) == -1)
        _mkdir(save_path.c_str());

    myfile.open("..//outputs//model_SLN_full//spike.txt");
    for (int i = 0; i < spike_count; i++)
    {
        myfile << time_spike[i] << "  " << num_spike[i] << endl;
    }
    myfile.close();

    myfile.open("..//outputs//model_SLN_full//H.txt");
    for (int i = 0; i < total_HIE.size(); i++)
    {
        myfile << record_time_point[i] << " " << total_HEE[i] << " " << total_HIE[i] << " " << total_HSE[i] << " " << total_HEEN[i] << " " << total_HIEN[i] << " " << total_HSEN[i] << " " << total_HEI[i] << " " << total_HII[i] << " " << total_HSI[i] << " " << total_HES[i] << " " << total_HIS[i] << " " << total_HSS[i] << endl;
    }
    myfile.close();
    myfile.open("..//outputs//model_SLN_full//V.txt");
    for (int i = 0; i < total_HIE.size(); i++)
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

    //gettimeofday(&t2, NULL);
    //double delta = ((t2.tv_sec - t1.tv_sec) * 1000000u +
    //    t2.tv_usec - t1.tv_usec) / 1.e6;

    //cout << "total CPU time = " << delta << endl;
}







