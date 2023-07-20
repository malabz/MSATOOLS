#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <tuple>
#include <cmath>
#include "../SuffixTree/SuffixTree.hpp"
#include "../Utils/Utils.hpp"

/*#include <windows.h>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")
using namespace std;
void showMemoryInfo(void)
{
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
    cout << "内存使用:" << pmc.WorkingSetSize / 1000 << "K/" << pmc.PeakWorkingSetSize / 1000 << "K + " << pmc.PagefileUsage / 1000 << "K/" << pmc.PeakPagefileUsage / 1000 << "K" << endl;
}
#undef max
#undef min*/
float extend = 2;

using RandomAccessIterator = std::vector<unsigned char>::const_iterator;
using gap_vector_type = std::vector<size_t>;
using char_vector_type = std::vector<unsigned char>;
using triple = std::array<size_t, 3>;
using quadra = std::array<size_t, 4>;
using insert = std::vector<std::tuple<int, int>>;
using in = std::tuple<int, int>;

namespace Kband
{ 
    int match = 1;
    int mismatch = -2;
    int d = 3;
    int e = 1;
    int my_INT_MIN;
    const int thresh0 = 1201;

    unsigned char* A = new unsigned char[thresh0];
    unsigned char* B = new unsigned char[thresh0];
    int (*pm)[thresh0] = new int[3][thresh0];
    int (*pm2)[thresh0] = new int[3][thresh0];
    int (*pmt1)[thresh0] = pm;
    int (*pmt2)[thresh0] = pm2;
    int (*pmt)[thresh0] = pm;
    unsigned char (*bt)[thresh0] = new unsigned char[thresh0][thresh0];
    std::vector<unsigned char> seq_A(2*thresh0);
    std::vector<unsigned char> seq_B(2*thresh0);

    inline int score(unsigned char xi, unsigned char yi)
    {
        if (xi == yi)
            return match;
        else
            return mismatch;
    }
    inline bool InsiderStrip(int i, int j, int k, int diff = 0)
    {
        return ((-k <= (j - i)) && ((j - i) <= (k + diff)));
    }
    inline int index(int i, int l)
    {
        return (i + l) % l;
    }
    inline int maxi(int a, int b)
    {
        if (a > b)return a;
        else return b;
    }
    inline int maxi(int a, int b, int c)
    {
        int max; if (a > b)max = a; else max = b;
        if (max > c)return max; else return c;
    }
    void Init(int m, int k, int diff)
    {
        for (int i = 0; i < (m + 1); i++)
        {
            for (int j = 0; j < (diff + 2 * k + 1); j++)
                bt[i][j] = '\0';
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < (diff + 2 * k + 1); j++)
                pm[i][j] = my_INT_MIN;
        }
        pm[0][k] = 0;
        bt[0][k] = '\16';
        for (int j = 0; j < (diff + k + 1); j++)
        {
            pm[1][j + k] = -d - e * (j - 1);
            bt[0][j + k] = (char)8;
        }
        for (int i = 0; i < (k + 1); i++)
            bt[i][k - i] = '\3';
    }
    void InitTwo( int ii, int k, int diff)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < (diff + 2 * k + 1); j++)
                pm2[i][j] = my_INT_MIN;
        }
        if (ii < k + 1)
            pm2[2][index(k - ii,diff + 2 * k + 1)] = -d - e * (ii - 1);
    }
    int ChooseWay(int p0, int p1, int p2, bool state = true)
    {
        if (p0 >= p1)
        {
            if (p0 >= p2)
                return state ? 16 : 0;
            else
                return state ? 48 : 2;
        }
        else if (p1 >= p2)
            return state ? 32 : 1;
        else
            return state ? 48 : 2;
    }
    inline int parse(int b, int s)
    {
        //b = (int)b;
        b = (b >> (4 - s * 2));
        return (b & 3) - 1;
    }

    void PSA_AGP_Kband_more(RandomAccessIterator lhs_first, RandomAccessIterator lhs_last,
            RandomAccessIterator rhs_first, RandomAccessIterator rhs_last, int a_begin, int b_begin, utils::more_block& more_gap, int cmatch = 1, int cmismatch = -2, int cd = 3, int ce = 1)
    {
        utils::m_block * Block = new utils::m_block();
        more_gap.push_back(Block);
        Block->start1 = a_begin;
        Block->start2 = b_begin;
        
        match = cmatch;
        mismatch = cmismatch;
        d = cd;
        e = ce;
        pm = pmt1;
        pm2 = pmt2;
        int i = 0, j, z, diff, k = 1, m, n, b_j, l, old, bt1, bt2, bt3, newi, l2, channel;
        bool state_ex = false;//交换标识
        unsigned char* tmp;
        int a_len = std::distance(lhs_first, lhs_last);
        int b_len = std::distance(rhs_first, rhs_last);
        
        Block->end1 = a_begin + a_len;
        Block->end2 = b_begin + b_len;
        insert& a_gap = Block->gap1;
        insert& b_gap = Block->gap2;

        i = 0;
        while (lhs_first < lhs_last) A[i++] = *(lhs_first++);
        i = b_len - 1;
        while (rhs_first < rhs_last) Kband::B[i--] = 5 - (int)*(rhs_first++);
        /*for (i = 0; i < a_len; i++)
        {
            std::cout << (int)A[i];
        }std::cout << "\n";
        for (i = 0; i < b_len; i++)
        {
            std::cout << (int)B[i];
        }std::cout << "\n";*/
        if (a_len > b_len)  //保证  B长，A短
        {
            tmp = A; A = B; B = tmp;
            diff = a_len - b_len;
            i = a_len; a_len = b_len; b_len = i;
            state_ex = true; //交换标识
        }
        else
            diff = b_len - a_len;
        m = a_len, n = b_len;
        my_INT_MIN = old = -d * n;

        while (k <= m)
        {// init
            //Init(bt, pm, m, k, diff);
            pm = pmt1;
            pm2 = pmt2;
            for (i = 0; i < (m + 1); i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    bt[i][j] = '\0';
            }
            for (i = 0; i < 3; i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    pm[i][j] = my_INT_MIN;
            }
            pm[0][k] = 0;
            bt[0][k] = '\16';
            for (j = 0; j < (diff + k + 1); j++)
            {
                pm[1][j + k] = -d - e * (j - 1);
                bt[0][j + k] = (char)8;
            }
            for (i = 0; i < (k + 1); i++)
                bt[i][k - i] = '\3';
            l = diff + 2 * k + 1;
            //end-init
            for (i = 1; i < (m + 1); i++)
            {
                //InitTwo(pm2, i, k, diff);
                for (int q = 0; q < 3; q++)
                {
                    for (j = 0; j < (diff + 2 * k + 1); j++)
                        pm2[q][j] = my_INT_MIN;
                }
                if (i < k + 1)
                    pm2[2][index(k - i, diff + 2 * k + 1)] = -d - e * (i - 1);
                l2 = diff + 2 * k + 1;
                //end-init
                for (int z = -k; z < (diff + k + 1); z++)
                {
                    j = z;
                    if ((1 <= (j + i)) && ((j + i) <= n))
                    {
                        j = j + k;
                        bt1 = bt2 = bt3 = 0;
                        bt1 = ChooseWay(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]);
                        pm2[0][index(j, l2)] = maxi(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]) + score(A[i - 1], B[j + i - k - 1]);
                        if (InsiderStrip(i, j + i - k - 1, k, diff))// x : B[j] ~_
                        {
                            pm2[1][index(j, l2)] = maxi(pm2[0][index(j - 1, l2)] - d, pm2[1][index(j - 1, l2)] - e);
                            if ((pm2[0][index(j - 1, l2)] - d) > (pm2[1][index(j - 1, l2)] - e)) bt2 = 4;
                            else bt2 = 8;
                        }
                        if (InsiderStrip(i - 1, j + i - k, k, diff))// y : A[i] ~_
                        {
                            pm2[2][index(j, l2)] = maxi(pm[0][index(j + 1, l)] - d, pm[2][index(j + 1, l)] - e);
                            if ((pm[0][index(j + 1, l)] - d) > (pm[2][index(j + 1, l)] - e)) bt3 = 1;
                            else bt3 = 3;
                        }
                        bt[i][index(j, l)] = (char)(bt1 + bt2 + bt3);
                    }
                }
                pmt = pm;
                pm = pm2;
                pm2 = pmt;
            }
            newi = maxi(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k]);
            if (old == newi || (k * 2) > m) break;
            else { old = newi; k *= 2; }
        }
        channel = ChooseWay(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k], false);

        //traceback
        i = m;
        b_j = n;
        j = diff + k;

        seq_A.clear();
        seq_B.clear();

        while (i > 0 || j > k)
        {
            if (channel == 0)
            {
                channel = parse(bt[i][j], 0);
                seq_A.push_back(A[--i]);
                seq_B.push_back(B[--b_j]);
            }
            else if (channel == 1)
            {
                channel = parse(bt[i][j], 1);
                seq_A.push_back('\7');
                seq_B.push_back(B[b_j - 1]);
                b_j--;
                j--;
            }
            else if (channel == 2)
            {
                channel = parse(bt[i][j], 2);
                seq_A.push_back(A[i - 1]);
                seq_B.push_back('\7');
                i--;
                j++;
            }
            else
            {
                std::cout << "channel error!\n";
                exit(-1);
            }
        }

        int j1 = 0, j2 = 0, num1 = 0, num2 = 0, match_num = 0;
        if (state_ex)
        {
            for (i = seq_A.size() - 1; i >= 0; i--)
            {
                if (seq_A[i] == '\7') { num1++; }
                else { if (num1 != 0)b_gap.push_back(in(j1, num1)); num1 = 0; j1++; }
                if (seq_B[i] == '\7') { num2++; }
                else
                {
                    if (num2 != 0)a_gap.push_back(in(j2, num2)); num2 = 0; j2++;
                    if (seq_A[i] == seq_B[i]) match_num++;
                }
            }
            if (num1 != 0)b_gap.push_back(in(j1, num1));
            if (num2 != 0)a_gap.push_back(in(j2, num2));
        }
        else
        {
            for (i = seq_A.size() - 1; i >= 0; i--)
            {
                if (seq_A[i] == '\7') { num1++; }
                else { if (num1 != 0)a_gap.push_back(in(j1, num1)); num1 = 0; j1++; }
                if (seq_B[i] == '\7') { num2++; }
                else
                {
                    if (num2 != 0)b_gap.push_back(in(j2, num2)); num2 = 0; j2++;
                    if (seq_A[i] == seq_B[i]) match_num++;
                }
            }
            if (num1 != 0)a_gap.push_back(in(j1, num1));
            if (num2 != 0)b_gap.push_back(in(j2, num2));
        }
    }


    std::tuple<insert, insert>
    PSA_AGP_Kband3(RandomAccessIterator lhs_first, RandomAccessIterator lhs_last,
        RandomAccessIterator rhs_first, RandomAccessIterator rhs_last, int a_begin, int b_begin, insert& SNP, int cmatch = 1, int cmismatch = -2, int cd = 3, int ce = 1)
    {
        match = cmatch;
        mismatch = cmismatch;
        d = cd;
        e = ce;
        pm = pmt1;
        pm2 = pmt2;
        int i = 0, j, z, diff, k = 1, m, n, b_j, l, old, bt1, bt2, bt3, newi, l2, channel;
        bool state_ex = false;//交换标识
        unsigned char* tmp;
        int a_len = std::distance(lhs_first, lhs_last);
        int b_len = std::distance(rhs_first, rhs_last);

        insert a_gap;
        insert b_gap;
        if (a_len == 0)                             //a为0
        {
            a_gap.push_back(in(a_begin, b_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }
        else if (b_len == 0)                        //b为0
        {
            b_gap.push_back(in(b_begin, a_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }

        while (lhs_first < lhs_last) A[i++] = *(lhs_first++);
        i = 0;
        while (rhs_first < rhs_last) B[i++] = *(rhs_first++);
        /*for (i = 0; i < a_len; i++)
        {
            std::cout << (int)A[i];
        }std::cout << "\n";
        for (i = 0; i < b_len; i++)
        {
            std::cout << (int)B[i];
        }std::cout << "\n";*/
        if (a_len > b_len)  //保证  B长，A短
        {
            tmp = A; A = B; B = tmp;
            diff = a_len - b_len;
            i = a_len; a_len = b_len; b_len = i;
            state_ex = true; //交换标识
        }
        else
            diff = b_len - a_len;
        m = a_len, n = b_len;
        my_INT_MIN = old = -d * n;

        while (k <= m)
        {// init
            //Init(bt, pm, m, k, diff);
            pm = pmt1;
            pm2 = pmt2;
            for (i = 0; i < (m + 1); i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    bt[i][j] = '\0';
            }
            for (i = 0; i < 3; i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    pm[i][j] = my_INT_MIN;
            }
            pm[0][k] = 0;
            bt[0][k] = '\16';
            for (j = 0; j < (diff + k + 1); j++)
            {
                pm[1][j + k] = -d - e * (j - 1);
                bt[0][j + k] = (char)8;
            }
            for (i = 0; i < (k + 1); i++)
                bt[i][k - i] = '\3';
            l = diff + 2 * k + 1;
            //end-init
            for (i = 1; i < (m + 1); i++)
            {
                //InitTwo(pm2, i, k, diff);
                for (int q = 0; q < 3; q++)
                {
                    for (j = 0; j < (diff + 2 * k + 1); j++)
                        pm2[q][j] = my_INT_MIN;
                }
                if (i < k + 1)
                    pm2[2][index(k - i, diff + 2 * k + 1)] = -d - e * (i - 1);
                l2 = diff + 2 * k + 1;
                //end-init
                for (int z = -k; z < (diff + k + 1); z++)
                {
                    j = z;
                    if ((1 <= (j + i)) && ((j + i) <= n))
                    {
                        j = j + k;
                        bt1 = bt2 = bt3 = 0;
                        bt1 = ChooseWay(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]);
                        pm2[0][index(j, l2)] = maxi(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]) + score(A[i - 1], B[j + i - k - 1]);
                        if (InsiderStrip(i, j + i - k - 1, k, diff))// x : B[j] ~_
                        {
                            pm2[1][index(j, l2)] = maxi(pm2[0][index(j - 1, l2)] - d, pm2[1][index(j - 1, l2)] - e);
                            if ((pm2[0][index(j - 1, l2)] - d) > (pm2[1][index(j - 1, l2)] - e)) bt2 = 4;
                            else bt2 = 8;
                        }
                        if (InsiderStrip(i - 1, j + i - k, k, diff))// y : A[i] ~_
                        {
                            pm2[2][index(j, l2)] = maxi(pm[0][index(j + 1, l)] - d, pm[2][index(j + 1, l)] - e);
                            if ((pm[0][index(j + 1, l)] - d) > (pm[2][index(j + 1, l)] - e)) bt3 = 1;
                            else bt3 = 3;
                        }
                        bt[i][index(j, l)] = (char)(bt1 + bt2 + bt3);
                    }
                }
                pmt = pm;
                pm = pm2;
                pm2 = pmt;
            }
            newi = maxi(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k]);
            if (old == newi || (k * 2) > m) break;
            else { old = newi; k *= 2; }
        }
        channel = ChooseWay(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k], false);

        //traceback
        i = m;
        b_j = n;
        j = diff + k;

        seq_A.clear();
        seq_B.clear();

        while (i > 0 || j > k)
        {
            if (channel == 0)
            {
                channel = parse(bt[i][j], 0);
                seq_A.push_back(A[--i]);
                seq_B.push_back(B[--b_j]);
            }
            else if (channel == 1)
            {
                channel = parse(bt[i][j], 1);
                seq_A.push_back('\7');
                seq_B.push_back(B[b_j - 1]);
                b_j--;
                j--;
            }
            else if (channel == 2)
            {
                channel = parse(bt[i][j], 2);
                seq_A.push_back(A[i - 1]);
                seq_B.push_back('\7');
                i--;
                j++;
            }
            else
            {
                std::cout << "channel error!\n";
                exit(-1);
            }
        }

        int j1 = 0, j2 = 0, num1 = 0, num2 = 0, match_num = 0;
        if (state_ex)
        {
            for (i = seq_A.size() - 1; i >= 0; i--)
            {
                if (seq_A[i] == '\7') { num1++; }
                else { if (num1 != 0)b_gap.push_back(in(b_begin + j1, num1)); num1 = 0; j1++; }
                if (seq_B[i] == '\7') { num2++; }
                else
                {
                    if (num2 != 0)a_gap.push_back(in(a_begin + j2, num2)); num2 = 0; j2++;
                    if (seq_A[i] == seq_B[i]) match_num++;
                }
            }
            if (num1 != 0)b_gap.push_back(in(b_begin + j1, num1));
            if (num2 != 0)a_gap.push_back(in(a_begin + j2, num2));

            if (b_len <= 10 && (b_len - a_len) < 4)
            {
                if (seq_A.back() != '\7' && seq_B.back() != '\7' && (seq_A.back() != seq_B.back()))
                    SNP.push_back(std::make_tuple(a_begin, b_begin));
                if (seq_A[0] != '\7' && seq_B[0] != '\7' && (seq_A[0] != seq_B[0]))
                    SNP.push_back(std::make_tuple(a_begin + b_len - 1, b_begin + a_len - 1));
            }
            else if (b_len > 10 && (b_len - a_len) < 10 && match_num > (int)(a_len * 0.8))
            {
                //std::cout << b_len << " " << a_len << " " << match_num << "\n";
                j1 = j2 = 0;
                for (i = seq_A.size() - 1; i >= 0; i--)
                {
                    if (seq_A[i] != '\7' && seq_B[i] != '\7')
                    {
                        if (seq_A[i] != seq_B[i])
                            SNP.push_back(std::make_tuple(a_begin + j2, b_begin + j1));
                        j1++; j2++;
                    }
                    else
                    {
                        if (seq_A[i] != '\7') j1++;
                        if (seq_B[i] != '\7') j2++;
                    }
                }
            }

        }
        else
        {
            for (i = seq_A.size() - 1; i >= 0; i--)
            {
                if (seq_A[i] == '\7') { num1++; }
                else { if (num1 != 0)a_gap.push_back(in(a_begin + j1, num1)); num1 = 0; j1++; }
                if (seq_B[i] == '\7') { num2++; }
                else
                {
                    if (num2 != 0)b_gap.push_back(in(b_begin + j2, num2)); num2 = 0; j2++;
                    if (seq_A[i] == seq_B[i]) match_num++;
                }
            }
            if (num1 != 0)a_gap.push_back(in(a_begin + j1, num1));
            if (num2 != 0)b_gap.push_back(in(b_begin + j2, num2));

            if (a_len == 1 && b_len == 1)
                SNP.push_back(std::make_tuple(a_begin, b_begin));
            else if (b_len <= 10 && (b_len - a_len) < 4)
            {
                if (seq_A.back() != '\7' && seq_B.back() != '\7' && (seq_A.back() != seq_B.back()))
                    SNP.push_back(std::make_tuple(a_begin, b_begin));
                if (seq_A[0] != '\7' && seq_B[0] != '\7' && (seq_A[0] != seq_B[0]))
                    SNP.push_back(std::make_tuple(a_begin + a_len - 1, b_begin + b_len - 1));
            }
            else if ((b_len > 10) && ((b_len - a_len) < 10) && (match_num > (int)(a_len * 0.8)))
            {
                j1 = j2 = 0;
                //std::cout << b_len << " " << a_len << " " << match_num << "\n";
                for (i = seq_A.size() - 1; i >= 0; i--)
                {
                    if ((seq_A[i] != '\7') && (seq_B[i] != '\7'))
                    {
                        if (seq_A[i] != seq_B[i])
                            SNP.push_back(std::make_tuple(a_begin + j1, b_begin + j2));
                        j1++; j2++;
                    }
                    else
                    {
                        if (seq_A[i] != '\7') j1++;
                        if (seq_B[i] != '\7') j2++;
                    }
                }
            }
        }


        //EmptySet();
        return std::make_tuple(std::move(a_gap), std::move(b_gap));
        
    }


    void PSA_AGP_Kband_SW(std::vector<unsigned char>& B_ans, std::vector<unsigned char>& A_ans, int _b_len, int _a_len, int &ans_b, int &ans_a)
    {
        pm = pmt1;
        pm2 = pmt2;
        int i = 0, j, z, diff, k = 1, m, n, b_j, l, old, bt1, bt2, bt3, newi, l2, channel;
        ans_a=0, ans_b=0;
        bool state_ex = false;//交换标识
        unsigned char* tmp;
        int a_len = _a_len;
        int b_len = _b_len;
        
        if (a_len > b_len)  //保证  B长，A短
        {
            tmp = A; A = B; B = tmp;
            diff = a_len - b_len;
            i = a_len; a_len = b_len; b_len = i;
            state_ex = true; //交换标识
        }
        else
            diff = b_len - a_len;
        m = a_len, n = b_len;
        my_INT_MIN = old = -d * n;

        while (k <= m)
        {
            pm = pmt1;
            pm2 = pmt2;
            for (i = 0; i < (m + 1); i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    bt[i][j] = '\0';
            }
            for (i = 0; i < 3; i++)
            {
                for (int j = 0; j < (diff + 2 * k + 1); j++)
                    pm[i][j] = my_INT_MIN;
            }
            pm[0][k] = 0;
            bt[0][k] = '\16';
            for (j = 0; j < (diff + k + 1); j++)
            {
                pm[1][j + k] = -d - e * (j - 1);
                bt[0][j + k] = (char)8;
            }
            for (i = 0; i < (k + 1); i++)
                bt[i][k - i] = '\3';
            l = diff + 2 * k + 1;
            //end-init
            for (i = 1; i < (m + 1); i++)
            {
                //InitTwo(pm2, i, k, diff);
                for (int q = 0; q < 3; q++)
                {
                    for (j = 0; j < (diff + 2 * k + 1); j++)
                        pm2[q][j] = my_INT_MIN;
                }
                if (i < k + 1)
                    pm2[2][index(k - i, diff + 2 * k + 1)] = -d - e * (i - 1);
                l2 = diff + 2 * k + 1;
                //end-init
                for (int z = -k; z < (diff + k + 1); z++)
                {
                    j = z;
                    if ((1 <= (j + i)) && ((j + i) <= n))
                    {
                        j = j + k;
                        bt1 = bt2 = bt3 = 0;
                        bt1 = ChooseWay(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]);
                        pm2[0][index(j, l2)] = maxi(pm[0][index(j, l)], pm[1][index(j, l)], pm[2][index(j, l)]) + score(A[i - 1], B[j + i - k - 1]);
                        if (InsiderStrip(i, j + i - k - 1, k, diff))// x : B[j] ~_
                        {
                            pm2[1][index(j, l2)] = maxi(pm2[0][index(j - 1, l2)] - d, pm2[1][index(j - 1, l2)] - e);
                            if ((pm2[0][index(j - 1, l2)] - d) > (pm2[1][index(j - 1, l2)] - e)) bt2 = 4;
                            else bt2 = 8;
                        }
                        if (InsiderStrip(i - 1, j + i - k, k, diff))// y : A[i] ~_
                        {
                            pm2[2][index(j, l2)] = maxi(pm[0][index(j + 1, l)] - d, pm[2][index(j + 1, l)] - e);
                            if ((pm[0][index(j + 1, l)] - d) > (pm[2][index(j + 1, l)] - e)) bt3 = 1;
                            else bt3 = 3;
                        }
                        bt[i][index(j, l)] = (char)(bt1 + bt2 + bt3);
                    }
                }
                pmt = pm;
                pm = pm2;
                pm2 = pmt;
            }
            
            newi = maxi(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k]);
            if (old == newi || (k * 2) > m) break;
            else { old = newi; k *= 2; }
        }
        channel = ChooseWay(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k], false);

        i = m;
        b_j = n;
        j = diff + k;

        seq_A.clear();
        seq_B.clear();

        while (i > 0 || j > k)
        {
            if (channel == 0)
            {
                channel = parse(bt[i][j], 0);
                seq_A.push_back(A[--i]);
                seq_B.push_back(B[--b_j]);
            }
            else if (channel == 1)
            {
                channel = parse(bt[i][j], 1);
                seq_A.push_back('\7');
                seq_B.push_back(B[b_j - 1]);
                b_j--;
                j--;
            }
            else if (channel == 2)
            {
                channel = parse(bt[i][j], 2);
                seq_A.push_back(A[i - 1]);
                seq_B.push_back('\7');
                i--;
                j++;
            }
            else
            {
                std::cout << "channel error!\n";
                exit(-1);
            }
        }

        if (state_ex)
        {
            A_ans = seq_B;
            B_ans = seq_A;
        }
        else
        {
            A_ans = seq_A;
            B_ans = seq_B;
        }
        bool gap1 = false;   //false代表首个gap
        int score_two=0, max_score=0, max_i= seq_A.size() - 1;
        for (i = seq_A.size() - 1; i >= 0; i--)
        {
            if ((seq_A[i] != '\7') && (seq_B[i] != '\7'))//都不是空格
            {
                if (gap1) gap1 = false; //下一次为首个空格
                if (seq_A[i] == seq_B[i]) score_two += match; //匹配
                else score_two += mismatch;   // 不匹配
            }
            else if ((seq_A[i] != '\7') || (seq_B[i] != '\7'))
            {
                if (gap1 == false)//首个空格罚分
                {
                    score_two -= d; gap1 = true;//下次非首个空格
                }
                else score_two -= e;//非首个空格罚分
            }
            if (score_two > max_score)
            {
                max_score = score_two;
                max_i = i;
            }
        }
        //std::cout << max_score << " " << max_i << " SCORE\n";
        for (i = seq_A.size() - 1; i >= max_i; i--)
        {
            if (seq_A[i] != '\7')ans_a++;
            if (seq_B[i] != '\7')ans_b++;
        }
        //std::cout << ans_a << " " << ans_b << " \n";
        /*
        for (i = seq_A.size() - 1; i >= 0; i--)
            std::cout << (int)seq_A[i];
        std::cout << "\n";
        for (i = seq_B.size() - 1; i >= 0; i--)
            std::cout << (int)seq_B[i];
        std::cout << "\n";
        */
    }
}

namespace BWT_MINI //星比对命名空间
{
    //BWT套件！！！
    //第一步，从传入全部的同源区段，[[A.index，B.index，length]...]，每个B.index选出一个A.index，相对距离最近
    std::vector<triple> _MultiReg(const std::vector<triple>& common_substrings) //最优路径
    {
        std::vector<triple> optimal_common_substrings;
        if (common_substrings.empty()) return optimal_common_substrings;
        int start = common_substrings[0][0];  //初始化
        int end = common_substrings[0][0] + common_substrings[0][2];//
        int i;
        float a_length, tmp, pre_tmp, b_length = common_substrings.rbegin()[0][1] + common_substrings.rbegin()[0][2]; //B的长度
        for (i = 1; i < common_substrings.size(); i++) // 找出所有后缀号列表中，最小后缀号做  start起点，最大后缀号 + length_i做end终点
        {
            if (common_substrings[i][0] < start) start = common_substrings[i][0];
            if ((common_substrings[i][0] + common_substrings[i][2]) > end) end = common_substrings[i][0] + common_substrings[i][2];
        }
        a_length = end - start;
        i = 0;
        optimal_common_substrings.push_back(common_substrings[i]);//first
        pre_tmp = common_substrings[i][0] / a_length - common_substrings[i][1] / b_length;//加入第一个元素
        pre_tmp = fabs(pre_tmp);
        for (i = i + 1; i < common_substrings.size(); i++)
        {
            if (optimal_common_substrings.rbegin()[0][1] != common_substrings[i][1]) //b.index 不同
            {
                pre_tmp = common_substrings[i][0] / a_length - common_substrings[i][1] / b_length;
                pre_tmp = fabs(pre_tmp);
                optimal_common_substrings.push_back(common_substrings[i]);
            }
            else
            {
                tmp = common_substrings[i][0] / a_length - common_substrings[i][1] / b_length;
                tmp = fabs(tmp);
                if (tmp < pre_tmp)
                {
                    optimal_common_substrings.rbegin()[0][0] = common_substrings[i][0];
                    pre_tmp = tmp;
                }
            }
        }
        return optimal_common_substrings;
    }
    //回溯
    std::vector<int> _trace_back(const std::vector<triple>& common_substrings, int* p)
    {
        std::vector<int> ansi;
        int* tmp = p;
        for (int i = 0; i < common_substrings.size(); i++)
        {
            if (*p < tmp[i])
                p = &tmp[i];
        }
        int j, i = p - tmp;
        ansi.push_back(i);
        p = tmp;
        while (i > 0)
        {
            j = i - 1;
            if (p[i] == common_substrings[i][2]) break;
            while (j >= 0)
            {

                if ((common_substrings[i][0] >= (common_substrings[j][0] + common_substrings[j][2])) && (p[i] == (p[j] + common_substrings[i][2])))
                {
                    ansi.push_back(j);
                    i = j;
                    break;
                }
                j--;
            }
        }
        reverse(ansi.begin(), ansi.end());//反转
        return ansi;
    }
    //第二步，依据动态规划，选出合适的不重叠的同源区段
    std::vector<triple> _optimal_path(const std::vector<triple>& common_substrings) //最优路径
    {
        std::vector<triple> optimal_common_substrings = _MultiReg(common_substrings);//调用第一步结
        int m = optimal_common_substrings.size();
        if (m <= 1) return optimal_common_substrings;
        int* p = new int[m];
        for (int i = 0; i < m; i++) p[i] = optimal_common_substrings[i][2];
        for (int i = 1; i < m; i++)
            for (int j = 0; j < i; j++)
                if (optimal_common_substrings[i][0] >= (optimal_common_substrings[j][0] + optimal_common_substrings[j][2]))
                    p[i] = (p[i] > (p[j] + optimal_common_substrings[i][2])) ? p[i] : (p[j] + optimal_common_substrings[i][2]);
        std::vector<int> ansi = _trace_back(optimal_common_substrings, p);

        std::vector<triple> ans_common_substrings;
        for (int i = 0; i < ansi.size(); i++) ans_common_substrings.push_back(optimal_common_substrings[ansi[i]]);
        std::vector<triple>().swap(optimal_common_substrings);//清空一块内存
        std::vector<int>().swap(ansi);
        delete[] p;
        return ans_common_substrings;
    }
/*
    std::vector<triple> pre_MultiReg(const std::vector<triple>& common_substrings) //最优路径
    {
        
        std::vector<triple> optimal_common_substrings;
        if (common_substrings.empty()) return optimal_common_substrings;
        int i;
        int a_length, tmp, pre_tmp; //B的长度

        i = 0;
        optimal_common_substrings.push_back(common_substrings[i]);//first
        pre_tmp = common_substrings[i][0];//加入第一个元素
        for (i = i + 1; i < common_substrings.size(); i++)
        {
            if (optimal_common_substrings.rbegin()[0][1] != common_substrings[i][1]) //b.index 不同
            {
                pre_tmp = common_substrings[i][0];
                optimal_common_substrings.push_back(common_substrings[i]);
            }
            else
            {
                if (common_substrings[i][0] < pre_tmp)
                {
                    optimal_common_substrings.rbegin()[0][0] = common_substrings[i][0];
                    pre_tmp = common_substrings[i][0];
                }
            }
        }

        return optimal_common_substrings;

    }
    std::vector<triple> pre_optimal_path(const std::vector<triple>& common_substrings) //最优路径
    {
        std::vector<triple> optimal_common_substrings = pre_MultiReg(common_substrings);//调用第一步结果

        int m = optimal_common_substrings.size();
        if (m <= 1) return optimal_common_substrings;
        int* p = new int[m];
        for (int i = 0; i < m; i++) p[i] = optimal_common_substrings[i][2];
        for (int i = 1; i < m; i++)
            for (int j = 0; j < i; j++)
                if (optimal_common_substrings[i][0] >= (optimal_common_substrings[j][0] + optimal_common_substrings[j][2]))
                    p[i] = (p[i] > (p[j] + optimal_common_substrings[i][2])) ? p[i] : (p[j] + optimal_common_substrings[i][2]);
        std::vector<int> ansi = _trace_back(optimal_common_substrings, p);

        std::vector<triple> ans_common_substrings;
        for (int i = 0; i < ansi.size(); i++) ans_common_substrings.push_back(optimal_common_substrings[ansi[i]]);
        std::vector<triple>().swap(optimal_common_substrings);//清空一块内存
        std::vector<int>().swap(ansi);
        delete[] p;
        return ans_common_substrings;
    }
*/
    std::tuple<insert, insert>
        BWT_mini(RandomAccessIterator lhs_first, RandomAccessIterator lhs_last,
            RandomAccessIterator rhs_first, RandomAccessIterator rhs_last, int a_begin, int b_begin, int thresh, insert& SNP, int d = 3, int e = 1)
    {
        int a_len = std::distance(lhs_first, lhs_last);
        int b_len = std::distance(rhs_first, rhs_last);
        std::tuple<insert, insert> before;
        std::tuple<insert, insert> after;
        int a, b, c, dd;
        //std::cout << a_len << " " << b_len << "--len\n";
        insert a_gap;
        insert b_gap;
        if (a_len == 0)                             //a为0
        {
            a_gap.push_back(in(a_begin, b_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }
        else if (b_len == 0)                        //b为0
        {
            b_gap.push_back(in(b_begin, a_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }
        /*else if (((a_len > b_len ? (a_len / b_len) : (b_len / a_len)) > 10) && (std::min(a_len, b_len) < 10000))
        {
            if (a_len > b_len)  //长串在前
            {
                int* fun = get_bwt_LCS(lhs_first, lhs_last, rhs_first, rhs_last);//
                int max = fun[0], i1 = fun[1], i2 = fun[2];
                std::cout << max << " " << i1 << " " << i2 << "--1max\n";
                if ((int)(i2 * extend) > i1)//前i1 不够长
                {
                    a = 0; b = i1; c = 0; dd = i2;
                    before = PSA_AGP_Kband3(lhs_first + a, lhs_first + b, rhs_first + c, rhs_first + dd, a_begin + a, b_begin + c);
                }
                else
                {
                    a = i1 - (int)(i2 * extend); b = i1; c = 0; dd = i2;
                    before = PSA_AGP_Kband3(lhs_first + a, lhs_first + b, rhs_first + c, rhs_first + dd, a_begin + a, b_begin + c);
                    if (std::get<1>(before).size() > 0 && std::get<0>(std::get<1>(before)[0]) == b_begin)
                        std::get<1>(std::get<1>(before)[0]) += a;
                    else
                        std::get<1>(before).push_back(in(b_begin, a));
                }


                if ((int)((b_len - i2 - max) * extend) > (a_len - i1 - max))//后i1 不够长
                {
                    a = i1 + max; b = a_len; c = i2 + max; dd = b_len;
                    after = PSA_AGP_Kband3(lhs_first + a, lhs_first + b, rhs_first + c, rhs_first + dd, a_begin + a, b_begin + c);
                }
                else
                {
                    a = i1 + max; b = a + (int)((b_len - i2 - max) * extend); c = i2 + max; dd = b_len;
                    after = PSA_AGP_Kband3(lhs_first + a, lhs_first + b, rhs_first + c, rhs_first + dd, a_begin + a, b_begin + c);
                    if (std::get<1>(after).size() > 0 && std::get<0>(std::get<1>(after).back()) == b_begin + dd)
                        std::get<1>(std::get<1>(after).back()) += (a_len - b);
                    else
                        std::get<1>(after).push_back(in((b_begin + dd), (a_len - b)));
                }
                std::get<0>(before).insert(std::get<0>(before).end(), std::get<0>(after).begin(), std::get<0>(after).end());
                std::get<1>(before).insert(std::get<1>(before).end(), std::get<1>(after).begin(), std::get<1>(after).end());
                insert().swap(std::get<0>(after));
                insert().swap(std::get<1>(after));
                std::cout << max << " " << i1 << " " << i2 << "--end\n";
                return std::move(before);
            }
            else
            {
                a = a_len; a_len = b_len; b_len = a;
                a = a_begin; a_begin = b_begin; b_begin = a;

                int* fun = get_bwt_LCS(rhs_first, rhs_last, lhs_first, lhs_last);//
                int max = fun[0], i1 = fun[1], i2 = fun[2];
                std::cout << max << " " << i1 << " " << i2 << "--2max\n";
                if ((int)(i2 * extend) > i1)//前i1 不够长
                {
                    a = 0; b = i1; c = 0; dd = i2;
                    before = PSA_AGP_Kband3(rhs_first + a, rhs_first + b, lhs_first + c, lhs_first + dd, a_begin + a, b_begin + c);
                }
                else
                {
                    a = i1 - (int)(i2 * extend); b = i1; c = 0; dd = i2;
                    before = PSA_AGP_Kband3(rhs_first + a, rhs_first + b, lhs_first + c, lhs_first + dd, a_begin + a, b_begin + c);
                    if (std::get<1>(before).size() > 0 && std::get<0>(std::get<1>(before)[0]) == b_begin)
                        std::get<1>(std::get<1>(before)[0]) += a;
                    else
                        std::get<1>(before).push_back(in(b_begin, a));
                }


                if ((int)((b_len - i2 - max) * extend) > (a_len - i1 - max))//后i1 不够长
                {
                    a = i1 + max; b = a_len; c = i2 + max; dd = b_len;
                    after = PSA_AGP_Kband3(rhs_first + a, rhs_first + b, lhs_first + c, lhs_first + dd, a_begin + a, b_begin + c);
                }
                else
                {
                    a = i1 + max; b = a + (int)((b_len - i2 - max) * extend); c = i2 + max; dd = b_len;
                    after = PSA_AGP_Kband3(rhs_first + a, rhs_first + b, lhs_first + c, lhs_first + dd, a_begin + a, b_begin + c);
                    if (std::get<1>(after).size() > 0 && std::get<0>(std::get<1>(after).back()) == b_begin + dd)
                        std::get<1>(std::get<1>(after).back()) += (a_len - b);
                    else
                        std::get<1>(after).push_back(in((b_begin + dd), (a_len - b)));
                }
                std::get<0>(before).insert(std::get<0>(before).end(), std::get<0>(after).begin(), std::get<0>(after).end());
                std::get<1>(before).insert(std::get<1>(before).end(), std::get<1>(after).begin(), std::get<1>(after).end());
                insert().swap(std::get<0>(after));
                insert().swap(std::get<1>(after));
                std::cout << max << " " << i1 << " " << i2 << "--end\n";
                return make_tuple(std::move(std::get<1>(before)), std::move(std::get<0>(before)));
                
            }
        }
        */
        else if ((a_len < thresh) && (b_len < thresh)) //A，B都小于阈值
        {
            /*
            auto [seq_A, seq_B] = PSA_AGP_Kband2(lhs_first, lhs_last, rhs_first, rhs_last);
            int i, j = 0;
            int num = 0;
            for (i = 0; i < seq_A.size(); ++i)
            {
                if (seq_A[i] == '\7') { num++; }
                else { if (num != 0)a_gap.push_back(in(a_begin + j, num)); num = 0; j++; }
            }
            if (num != 0)a_gap.push_back(in(a_begin + j, num));
            j = 0; num = 0;
            for (i = 0; i < seq_B.size(); ++i)
            {
                if (seq_B[i] == '\7') { num++; }
                else { if (num != 0)b_gap.push_back(in(b_begin + j, num)); num = 0; j++; }
            }
            if (num != 0)b_gap.push_back(in(b_begin + j, num));
            char_vector_type().swap(seq_A);
            char_vector_type().swap(seq_B);
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
            */
            return Kband::PSA_AGP_Kband3(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, SNP);
        }
        else if (b_len <= a_len)  //保持A长
        {
            suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> mini_bwt(lhs_first, lhs_last, nucleic_acid_pseudo::GAP);//实例化后缀树对象st，比对同源区域
            std::vector<triple> common_substrings = _optimal_path(mini_bwt.get_common_substrings(rhs_first, rhs_last, 15, false)); //common_substrings [A.index,B.index,length]
            if (common_substrings.size() == 0)
                common_substrings = _optimal_path(mini_bwt.get_common_substrings(rhs_first, rhs_last, 1, false));
            std::vector<quadra> intervals;
            intervals.reserve(common_substrings.size() + 1);
            if (common_substrings.empty())
            {
                intervals.push_back(quadra({ 0, (size_t)a_len, 0, (size_t)b_len }));
            }
            else
            {
                if (common_substrings[0][0] || common_substrings[0][1])
                    intervals.push_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

                for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                    if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                        common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                        intervals.push_back(quadra
                        ({
                            common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                            common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                            }));

                if (common_substrings.back()[0] + common_substrings.back()[2] != (size_t)a_len ||
                    common_substrings.back()[1] + common_substrings.back()[2] != (size_t)b_len)
                    intervals.push_back(quadra
                    ({
                        common_substrings.back()[0] + common_substrings.back()[2], (size_t)a_len,
                        common_substrings.back()[1] + common_substrings.back()[2], (size_t)b_len
                        }));
            }
            //三元组同源序列
            std::vector<triple>().swap(common_substrings);
            //四元组 mis片段
            for (size_t j = 0; j != intervals.size(); ++j)
            {
                const size_t centre_begin = intervals[j][0];
                const size_t centre_end = intervals[j][1];
                const size_t sequence_begin = intervals[j][2];
                const size_t sequence_end = intervals[j][3];
                auto [lhs_gaps, rhs_gaps] = BWT_mini(lhs_first + centre_begin, lhs_first + centre_end, rhs_first + sequence_begin, rhs_first + sequence_end, a_begin + centre_begin, b_begin + sequence_begin, thresh, SNP); //分治到比thresh小，然后k-band
                a_gap.insert(a_gap.end(), lhs_gaps.begin(), lhs_gaps.end());
                b_gap.insert(b_gap.end(), rhs_gaps.begin(), rhs_gaps.end());
                insert().swap(lhs_gaps);
                insert().swap(rhs_gaps);
            }
            std::vector<quadra>().swap(intervals);
            //EmptySet();
            return std::make_tuple(std::move(a_gap), std::move(b_gap));//返回000005000  gap-vecter
        }
        else
        {
            suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> mini_bwt(rhs_first, rhs_last, nucleic_acid_pseudo::GAP);//实例化后缀树对象st，比对同源区域
            std::vector<triple> common_substrings = _optimal_path(mini_bwt.get_common_substrings(lhs_first, lhs_last, 15, false)); //common_substrings [A.index,B.index,length]
            if (common_substrings.size() == 0)
                common_substrings = _optimal_path(mini_bwt.get_common_substrings(lhs_first, lhs_last, 1, false));
            std::vector<quadra> intervals;
            intervals.reserve(common_substrings.size() + 1);
            if (common_substrings.empty())
            {
                intervals.push_back(quadra({ 0, (size_t)b_len, 0, (size_t)a_len }));
            }
            else
            {
                if (common_substrings[0][0] || common_substrings[0][1])
                    intervals.push_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

                for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                    if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                        common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                        intervals.push_back(quadra
                        ({
                            common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                            common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                            }));

                if (common_substrings.back()[0] + common_substrings.back()[2] != (size_t)b_len ||
                    common_substrings.back()[1] + common_substrings.back()[2] != (size_t)a_len)
                    intervals.push_back(quadra
                    ({
                        common_substrings.back()[0] + common_substrings.back()[2], (size_t)b_len,
                        common_substrings.back()[1] + common_substrings.back()[2], (size_t)a_len
                        }));
            }
            //三元组同源序列
            std::vector<triple>().swap(common_substrings);
            //四元组 mis片段

            for (size_t j = 0; j != intervals.size(); ++j)
            {
                const size_t centre_begin = intervals[j][0];
                const size_t centre_end = intervals[j][1];
                const size_t sequence_begin = intervals[j][2];
                const size_t sequence_end = intervals[j][3];
                auto [lhs_gaps, rhs_gaps] = BWT_mini(lhs_first + sequence_begin, lhs_first + sequence_end, rhs_first + centre_begin, rhs_first + centre_end, a_begin + sequence_begin, b_begin + centre_begin, thresh, SNP); //分治到比thresh小，然后k-band
                a_gap.insert(a_gap.end(), lhs_gaps.begin(), lhs_gaps.end());
                b_gap.insert(b_gap.end(), rhs_gaps.begin(), rhs_gaps.end());
                insert().swap(lhs_gaps);
                insert().swap(rhs_gaps);
            }

            std::vector<quadra>().swap(intervals);
            EmptySet();
            return std::make_tuple(std::move(a_gap), std::move(b_gap));//返回000005000  gap-vecter
        }

    }

    std::tuple<insert, insert>
        BWT_mini_stop(RandomAccessIterator lhs_first, RandomAccessIterator lhs_last,
            RandomAccessIterator rhs_first, RandomAccessIterator rhs_last, int a_begin, int b_begin, int thresh, insert& SNP, int d = 3, int e = 1)
    {
        int a_len = std::distance(lhs_first, lhs_last);
        int b_len = std::distance(rhs_first, rhs_last);
        std::tuple<insert, insert> before;
        std::tuple<insert, insert> after;

        insert a_gap;
        insert b_gap;
        if (a_len == 0)                             //a为0
        {
            a_gap.push_back(in(a_begin, b_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }
        else if (b_len == 0)                        //b为0
        {
            b_gap.push_back(in(b_begin, a_len));
            return std::make_tuple(std::move(a_gap), std::move(b_gap));
        }
        else if (b_len <= a_len)  //保持A长
        {
            suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> mini_bwt(lhs_first, lhs_last, nucleic_acid_pseudo::GAP);//实例化后缀树对象st，比对同源区域
            std::vector<triple> common_substrings = _optimal_path(mini_bwt.get_common_substrings(rhs_first, rhs_last, 15, false)); //common_substrings [A.index,B.index,length]
            if (common_substrings.size() == 0)
                common_substrings = _optimal_path(mini_bwt.get_common_substrings(rhs_first, rhs_last, 1, false));
            std::vector<quadra> intervals;
            intervals.reserve(common_substrings.size() + 1);
            if (common_substrings.empty())
            {
                intervals.push_back(quadra({ 0, (size_t)a_len, 0, (size_t)b_len }));
            }
            else
            {
                if (common_substrings[0][0] || common_substrings[0][1])
                    intervals.push_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

                for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                    if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                        common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                        intervals.push_back(quadra
                        ({
                            common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                            common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                            }));

                if (common_substrings.back()[0] + common_substrings.back()[2] != (size_t)a_len ||
                    common_substrings.back()[1] + common_substrings.back()[2] != (size_t)b_len)
                    intervals.push_back(quadra
                    ({
                        common_substrings.back()[0] + common_substrings.back()[2], (size_t)a_len,
                        common_substrings.back()[1] + common_substrings.back()[2], (size_t)b_len
                        }));
            }
            //三元组同源序列
            std::vector<triple>().swap(common_substrings);
            //四元组 mis片段
            for (size_t j = 0; j != intervals.size(); ++j)
            {
                if ((intervals[j][1] - intervals[j][0]) > (intervals[j][3] - intervals[j][2]))
                    b_gap.push_back(in(b_begin + intervals[j][3], intervals[j][1] - intervals[j][0] - intervals[j][3] + intervals[j][2]));
                else
                    a_gap.push_back(in(a_begin + intervals[j][1], intervals[j][3] - intervals[j][2] - intervals[j][1] + intervals[j][0]));
            }
            std::vector<quadra>().swap(intervals);
            //EmptySet();
            return std::make_tuple(std::move(a_gap), std::move(b_gap));//返回000005000  gap-vecter
        }
        else
        {
            suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> mini_bwt(rhs_first, rhs_last, nucleic_acid_pseudo::GAP);//实例化后缀树对象st，比对同源区域
            std::vector<triple> common_substrings = _optimal_path(mini_bwt.get_common_substrings(lhs_first, lhs_last, 15, false)); //common_substrings [A.index,B.index,length]
            if (common_substrings.size() == 0)
                common_substrings = _optimal_path(mini_bwt.get_common_substrings(lhs_first, lhs_last, 1, false));
            std::vector<quadra> intervals;
            intervals.reserve(common_substrings.size() + 1);
            if (common_substrings.empty())
            {
                intervals.push_back(quadra({ 0, (size_t)b_len, 0, (size_t)a_len }));
            }
            else
            {
                if (common_substrings[0][0] || common_substrings[0][1])
                    intervals.push_back(quadra({ 0, common_substrings[0][0], 0, common_substrings[0][1] }));

                for (size_t j = 0, end_index = common_substrings.size() - 1; j != end_index; ++j)
                    if (common_substrings[j][0] + common_substrings[j][2] != common_substrings[j + 1][0] ||
                        common_substrings[j][1] + common_substrings[j][2] != common_substrings[j + 1][1])
                        intervals.push_back(quadra
                        ({
                            common_substrings[j][0] + common_substrings[j][2], common_substrings[j + 1][0],
                            common_substrings[j][1] + common_substrings[j][2], common_substrings[j + 1][1]
                            }));

                if (common_substrings.back()[0] + common_substrings.back()[2] != (size_t)b_len ||
                    common_substrings.back()[1] + common_substrings.back()[2] != (size_t)a_len)
                    intervals.push_back(quadra
                    ({
                        common_substrings.back()[0] + common_substrings.back()[2], (size_t)b_len,
                        common_substrings.back()[1] + common_substrings.back()[2], (size_t)a_len
                        }));
            }
            //三元组同源序列
            std::vector<triple>().swap(common_substrings);
            //四元组 mis片段

            for (size_t j = 0; j != intervals.size(); ++j)
            {
                if ((intervals[j][1] - intervals[j][0]) > (intervals[j][3] - intervals[j][2]))
                    a_gap.push_back(in(a_begin + intervals[j][3], intervals[j][1] - intervals[j][0] - intervals[j][3] + intervals[j][2]));
                else
                    b_gap.push_back(in(b_begin + intervals[j][1], intervals[j][3] - intervals[j][2] - intervals[j][1] + intervals[j][0]));
            }

            std::vector<quadra>().swap(intervals);
            EmptySet();
            return std::make_tuple(std::move(a_gap), std::move(b_gap));//返回000005000  gap-vecter
        }

    }
}


//顶层调用，递归 返回匹配好的串vector
std::tuple<insert, insert>
new_main_Align(RandomAccessIterator lhs_first, RandomAccessIterator lhs_last,
    RandomAccessIterator rhs_first, RandomAccessIterator rhs_last, int a_begin, int b_begin, int thresh, insert& SNP_i, bool _mg_tag, utils::more_block& more_gap, int d = 3, int e = 1)
{
    float thresh_ni = 0.01;
    int a_len = std::distance(lhs_first, lhs_last);
    int b_len = std::distance(rhs_first, rhs_last);
    insert a_gaps;
    insert b_gaps;
    int gap_num = 0, i = 0;
    //极端长度差序列硬比。
    /*if (a_len != 0 && b_len != 0)
    {
        if (a_len > b_len)
        {
            if (a_len / b_len > 100 || (a_len / b_len > 10 && b_len < thresh))
            {
                return BWT_MINI::BWT_mini_stop(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, thresh, SNP_i);
                b_gaps.push_back(in(b_begin + b_len, a_len - b_len));
                return std::make_tuple(std::move(a_gaps), std::move(b_gaps));
            }
        }
        else
        {
            if (b_len / a_len > 100 || (b_len / a_len > 10 && a_len < thresh))
            {
                return BWT_MINI::BWT_mini_stop(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, thresh, SNP_i);
                a_gaps.push_back(in(a_begin + a_len, b_len - a_len));
                return std::make_tuple(std::move(a_gaps), std::move(b_gaps));
            }
        }
    }
    */
    if ((a_len > 100000) || (b_len > 100000))
    {
        std::cout <<"long: " << a_len << " " << b_len << " \n";
    }
    //&& a_len*1.0/b_len < 10 && b_len*1.0/a_len < 10
    if (a_len == 0)                             //a为0
    {
        a_gaps.push_back(in(a_begin, b_len));
        return std::make_tuple(std::move(a_gaps), std::move(b_gaps));
    }
    else if (b_len == 0)                        //b为0
    {
        b_gaps.push_back(in(b_begin, a_len));
        return std::make_tuple(std::move(a_gaps), std::move(b_gaps));
    }
    else if (_mg_tag && a_len > 100 && b_len > 100 && a_len * 1.0 / b_len < 10 && b_len * 1.0 / a_len < 10)
    {
        if (a_len >= b_len)
        {
            if ((a_len < thresh) && (b_len < thresh)) //A，B都小于阈值
            {
                auto[agap, bgap] =  Kband::PSA_AGP_Kband3(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, SNP_i);
                gap_num = 0;
                for (int k = 0; k < agap.size(); k++)
                    gap_num += std::get<1>(agap[k]);
                for (int k = 0; k < bgap.size(); k++)
                    gap_num += std::get<1>(bgap[k]);
                if ((gap_num - a_len + b_len) / 2.0 / b_len > thresh_ni)
                {
                    Kband::PSA_AGP_Kband_more(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, more_gap);
                }
                return std::make_tuple(std::move(agap), std::move(bgap));
            }
            else
            {
                auto [agap, bgap] = BWT_MINI::BWT_mini(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, thresh, SNP_i); //MINI bwt变换求
                gap_num = 0;
                for (int k = 0; k < agap.size(); k++)
                    gap_num += std::get<1>(agap[k]);
                for (int k = 0; k < bgap.size(); k++)
                    gap_num += std::get<1>(bgap[k]);
                if ((gap_num - a_len + b_len) / 2.0 / b_len > thresh_ni)
                {
                    insert SNP;
                    std::vector<unsigned char> A(a_len);
                    std::vector<unsigned char> B(b_len);
                    utils::m_block* Block = new utils::m_block();
                    more_gap.push_back(Block);
                    Block->start1 = a_begin;
                    Block->start2 = b_begin;
                    Block->end1 = a_begin + a_len;
                    Block->end2 = b_begin + b_len;
                    
                    i = 0;
                    while (lhs_first < lhs_last) A[i++] = (int)*(lhs_first++);
                    i = b_len-1;
                    while (rhs_first < rhs_last) B[i--] = 5 - (int)*(rhs_first++);
                    auto [agap_, bgap_] = BWT_MINI::BWT_mini(A.cbegin(), A.cend(), B.cbegin(), B.cend(), 0, 0, thresh, SNP); //MINI bwt变换求
                    Block->gap1.assign(agap_.cbegin(), agap_.cend());
                    Block->gap2.assign(bgap_.cbegin(), bgap_.cend());
                    insert().swap(agap_);
                    insert().swap(bgap_);
                    std::vector<unsigned char>().swap(A);
                    std::vector<unsigned char>().swap(B);
                    std::vector<std::tuple<int, int>>().swap(SNP);
                }
                return std::make_tuple(std::move(agap), std::move(bgap));
            }
        }
        else
        {
            if ((a_len < thresh) && (b_len < thresh)) //A，B都小于阈值
            {
                auto [agap, bgap] = Kband::PSA_AGP_Kband3(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, SNP_i);
                gap_num = 0;
                for (int k = 0; k < agap.size(); k++)
                    gap_num += std::get<1>(agap[k]);
                for (int k = 0; k < bgap.size(); k++)
                    gap_num += std::get<1>(bgap[k]);
                if ((gap_num - b_len + a_len) / 2.0 / a_len > thresh_ni)
                {
                    Kband::PSA_AGP_Kband_more(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, more_gap);
                }    
                return std::make_tuple(std::move(agap), std::move(bgap));
            }
            else
            {
                auto [agap, bgap] = BWT_MINI::BWT_mini(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, thresh, SNP_i); //MINI bwt变换求
                gap_num = 0;
                for (int k = 0; k < agap.size(); k++)
                    gap_num += std::get<1>(agap[k]);
                for (int k = 0; k < bgap.size(); k++)
                    gap_num += std::get<1>(bgap[k]);
                if ((gap_num - b_len + a_len) / 2.0 / a_len > thresh_ni)
                {
                    insert SNP;
                    std::vector<unsigned char> A(a_len);
                    std::vector<unsigned char> B(b_len);
                    utils::m_block* Block = new utils::m_block();
                    more_gap.push_back(Block);
                    Block->start1 = a_begin;
                    Block->start2 = b_begin;
                    Block->end1 = a_begin + a_len;
                    Block->end2 = b_begin + b_len;

                    i = 0;
                    while (lhs_first < lhs_last) A[i++] = (int)*(lhs_first++);
                    i = b_len - 1;
                    while (rhs_first < rhs_last) B[i--] = 5 - (int)*(rhs_first++);
                    auto [agap_, bgap_] = BWT_MINI::BWT_mini(A.cbegin(), A.cend(), B.cbegin(), B.cend(), 0, 0, thresh, SNP); //MINI bwt变换求
                    Block->gap1.assign(agap_.cbegin(), agap_.cend());
                    Block->gap2.assign(bgap_.cbegin(), bgap_.cend());
                    insert().swap(agap_);
                    insert().swap(bgap_);
                    std::vector<unsigned char>().swap(A);
                    std::vector<unsigned char>().swap(B);
                    std::vector<std::tuple<int, int>>().swap(SNP);
                }
                return std::make_tuple(std::move(agap), std::move(bgap));
            }
        }
    }
    else if ((a_len < thresh) && (b_len < thresh)) //A，B都小于阈值
    {
        /*
        auto [seq_A, seq_B] = PSA_AGP_Kband2(lhs_first, lhs_last, rhs_first, rhs_last);
        int i, j = 0;
        int num = 0;
        for (i = 0; i < seq_A.size(); ++i)
        {
            if (seq_A[i] == '\7') { num++; }
            else { if (num != 0)a_gaps.push_back(in(a_begin + j, num)); num = 0; j++; }
        }
        if (num != 0)a_gaps.push_back(in(a_begin + j, num));
        j = 0; num = 0;
        for (i = 0; i < seq_B.size(); ++i)
        {
            if (seq_B[i] == '\7') { num++; }
            else { if (num != 0)b_gaps.push_back(in(b_begin + j, num)); num = 0; j++; }
        }
        if (num != 0)b_gaps.push_back(in(b_begin + j, num));
        char_vector_type().swap(seq_A);
        char_vector_type().swap(seq_B);
        return std::make_tuple(std::move(a_gaps), std::move(b_gaps));
        */

        return Kband:: PSA_AGP_Kband3(lhs_first, lhs_last, rhs_first, rhs_last,a_begin,b_begin, SNP_i);
    }
    else
    {
        return BWT_MINI::BWT_mini(lhs_first, lhs_last, rhs_first, rhs_last, a_begin, b_begin, thresh, SNP_i); //MINI bwt变换求
    }
}