#include <iostream>  
#include <string> 
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits.h>
#include <list>
#include <chrono>
#include<algorithm>
#include <cstring> 


#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#include <process.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <malloc.h>
#include <unistd.h>
#include <sys/resource.h>
#include <pthread.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L;        /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L;        /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif
#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;            /* Unsupported. */
#endif
}

inline void GetMemoryUsage()
{
    int mem = getPeakRSS() / 1024 / 1024;
    std::cout << "memory peak usage: " << mem << "MB" << std::endl;
}


using namespace std;
struct Insertion2   //插入
{
    size_t index; //索引位置
    size_t n_num; //从该位置开始有number个n
    size_t gap_num; //从该位置开始有number个gap
};
struct Insertion   //插入
{
    size_t index; //索引位置
    size_t number; //从该位置开始有number个gap
};

void pre_to_REMOVE_N(std::istream& is, std::ostream& os, std::ostream& ns) //读入数据
{
    size_t seq_num = UINT_MAX; //序列条数
    size_t index;
    size_t num;
    size_t setwidth = 0;
    bool tag;
    std::string namei;  //每条序列的name
    std::string each_line;
    std::string each_sequence;
    std::vector<Insertion> n_gap;

    while (seq_num)
    {
        seq_num /= 10;
        setwidth++;
    }//求位数
    ns << UINT_MAX << "\n";//占位用

    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);

            if (flag)
            {
                //each_sequence 
                os << "> " << namei << "\n";//写出os
                index = 0; num = 0; tag = true;
                for (int j = 0; j < each_sequence.size(); j++)
                {
                    if ((int)each_sequence[j] == 78 || (int)each_sequence[j] == 110)
                    {
                        if (tag)//首个N
                            tag = false;
                        num++;
                    }
                    else
                    {
                        os << each_sequence[j];//写出os
                        if (num != 0)
                        {
                            n_gap.push_back(Insertion({ index, num }));
                            tag = true;
                            num = 0;
                        }
                        index++;
                    }
                }if (num != 0) n_gap.push_back(Insertion({ index, num }));
                os << "\n";//写出os
                //写出ns
                ns << "> " << namei << "\n" << n_gap.size() << " " << index << "\n";
                for (int i = 0; i < n_gap.size(); i++)
                {
                    ns << n_gap[i].index << " " << n_gap[i].number << " ";
                    ns << "\n";
                }

                std::vector<Insertion>().swap(n_gap);
                each_sequence.clear();
            }

            flag = true;
            namei.assign(each_line);
            seq_num++;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    //each_sequence 
    os << "> " << namei << "\n";//写出os
    index = 0; num = 0; tag = true;
    for (int j = 0; j < each_sequence.size(); j++)
    {
        if ((int)each_sequence[j] == 78 || (int)each_sequence[j] == 110)
        {
            if (tag)//首个N
                tag = false;
            num++;
        }
        else
        {
            os << each_sequence[j];//写出os
            if (num != 0)
            {
                n_gap.push_back(Insertion({ index, num }));
                tag = true;
                num = 0;
            }
            index++;
        }
    }
    if (num != 0) n_gap.push_back(Insertion({ index, num }));

    os << "\n";//写出os
    //写出ns
    ns << "> " << namei << "\n" << n_gap.size() << " " << index << "\n";
    for (int i = 0; i < n_gap.size(); i++)
    {
        ns << n_gap[i].index << " " << n_gap[i].number << " ";
        ns << "\n";
    }

    std::vector<Insertion>().swap(n_gap);
    each_sequence.clear();
    ns.seekp(0, ios::beg);
    ns << std::setw(setwidth) << std::left << seq_num;
}

inline void insert_others(Insertion2 That, Insertion2 This, std::vector<std::vector<Insertion2>>& more_insertions, int k, int sequence_number, int* ii)
{
    for (int m = 0; m < sequence_number; m++)
    {
        if (m == k)
        {
            while ((ii[m] < more_insertions[m].size()) && (more_insertions[m][ii[m]].index < That.index)) ii[m]++;
            if (ii[m] == more_insertions[m].size())
                more_insertions[m].push_back(That);
            else if (more_insertions[m][ii[m]].index == That.index)
            {
                more_insertions[m][ii[m]].n_num += That.n_num;
                more_insertions[m][ii[m]].gap_num += That.gap_num;
            }
            else
                more_insertions[m].insert(more_insertions[m].begin() + (ii[m]++), That);
        }
        else
        {
            while ((ii[m] < more_insertions[m].size()) && (more_insertions[m][ii[m]].index < This.index)) ii[m]++;
            if (ii[m] == more_insertions[m].size())
                more_insertions[m].push_back(This);
            else if (more_insertions[m][ii[m]].index == This.index)
            {
                more_insertions[m][ii[m]].n_num += This.n_num;
                more_insertions[m][ii[m]].gap_num += This.gap_num;
            }
            else
                more_insertions[m].insert(more_insertions[m].begin() + (ii[m]++), This);
        }
    }
}

void insertion_gap_N(std::vector<std::list<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions)
{
    //新版插入，带N
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* ii = new int[sequence_number]();
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;

    //变insertions
    for (int k = 0; k < sequence_number; k++)
    {
        i = 0; j = 0; g_num = 0;
        for (int kk = 0; kk < sequence_number; kk++) ii[kk] = 0;
        std::vector<Insertion2>().swap(i_all_insertions);//清空 i_all_insertions
        while (i < insertions[k].size() && j < N_insertions[k].size())
        {
            if (insertions[k][i].index == N_insertions[k][j].index)//相等变number
            {
                g_num += insertions[k][i].number;
                if (N_insertions[k][j].number > insertions[k][i].number)
                {
                    i_all_insertions.push_back(Insertion2({ insertions[k][i].index ,insertions[k][i].number, 0 }));  //gap 个 N
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
                }
                else
                {
                    i_all_insertions.push_back(Insertion2({ insertions[k][i].index ,N_insertions[k][j].number ,insertions[k][i].number - N_insertions[k][j].number }));
                }
                i++;
                j++;
            }
            else if (insertions[k][i].index > N_insertions[k][j].index) //不等插新
            {
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
                j++;
            }
            else
            {
                g_num += insertions[k][i].number;
                i_all_insertions.push_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
                i++;
            }
        }
        while (j < N_insertions[k].size()) //后续多余的N，继续插入
        {
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size()) //后续多余的gap，继续插入
        {
            g_num += insertions[k][i].number;
            i_all_insertions.push_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.push_back(i_all_insertions);
    }

    //消除  一列多个N导致de多个插空
    int min_size, min_i;
    std::vector<Insertion> multi;
    std::vector<Insertion> tmp;
    for (int k = 0; k < sequence_number; k++)
        if (k == 0 || more_insertions[k].size() < min_size)
        {
            min_size = more_insertions[k].size();
            min_i = k;
        }
    i = 0; j = 0;
    if (min_i == 0) diff = 1;
    else diff = 0;
    while (i < more_insertions[diff].size() && j < more_insertions[min_i].size())
    {
        if (more_insertions[diff][i].index == more_insertions[min_i][j].index)
        {
            if (more_insertions[diff][i].gap_num == 0 || more_insertions[min_i][j].gap_num == 0);
            else if (more_insertions[diff][i].gap_num > more_insertions[min_i][j].gap_num)
                multi.push_back(Insertion({ more_insertions[diff][i].index ,more_insertions[min_i][j].gap_num }));
            else multi.push_back(Insertion({ more_insertions[diff][i].index ,more_insertions[diff][i].gap_num }));
            i++; j++;
        }
        else if (more_insertions[diff][i].index > more_insertions[min_i][j].index) j++;
        else i++;
    }

    if (min_i == 0)
        for (int k = 2; k < sequence_number; k++)
        {
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.push_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.push_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }
    else
        for (int k = 1; k < sequence_number; k++)
        {
            if (k == min_i) continue;
            i = 0; j = 0;
            std::vector<Insertion>().swap(tmp);
            while (i < more_insertions[k].size() && j < multi.size())
            {
                if (more_insertions[k][i].index == multi[j].index)
                {
                    if (more_insertions[k][i].gap_num == 0);
                    else if (more_insertions[k][i].gap_num > multi[j].number)
                        tmp.push_back(Insertion({ more_insertions[k][i].index ,multi[j].number }));
                    else tmp.push_back(Insertion({ more_insertions[k][i].index ,more_insertions[k][i].gap_num }));
                    i++; j++;
                }
                else if (more_insertions[k][i].index > multi[j].index) j++;
                else i++;
            }
            multi = tmp;
        }

    //插入到原串
    for (i = 0; i < sequence_number; i++)
    {
        auto iter = begin(sequences[i]);
        pre = all_insertions[i][all_insertions[i].size() - 1].index;
        std::advance(iter, pre);
        for (j = all_insertions[i].size() - 1; j >= 0; j--)
        {
            diff = all_insertions[i][j].index - pre;
            std::advance(iter, diff);
            sequences[i].insert(iter, all_insertions[i][j].n_num, 'N'); //与vector插入方式相反
            sequences[i].insert(iter, all_insertions[i][j].gap_num, '-');

            pre = all_insertions[i][j].index + all_insertions[i][j].gap_num + all_insertions[i][j].n_num;
        }

        if (more_insertions[i].size() > 0)
        {
            iter = begin(sequences[i]); mi = multi.size() - 1;
            pre = more_insertions[i][more_insertions[i].size() - 1].index;
            std::advance(iter, pre);
            for (j = more_insertions[i].size() - 1; j >= 0; j--)
            {
                diff = more_insertions[i][j].index - pre;
                std::advance(iter, diff);
                sequences[i].insert(iter, more_insertions[i][j].n_num, 'N'); //与vector插入方式相反
                if ((mi >= 0) && (more_insertions[i][j].index == multi[mi].index))
                {
                    sequences[i].insert(iter, (more_insertions[i][j].gap_num - multi[mi].number), '-');
                    pre = more_insertions[i][j].index + more_insertions[i][j].gap_num + more_insertions[i][j].n_num - multi[mi].number;
                    mi--;
                }
                else
                {
                    sequences[i].insert(iter, more_insertions[i][j].gap_num, '-');
                    pre = more_insertions[i][j].index + more_insertions[i][j].gap_num + more_insertions[i][j].n_num;
                }
            }
        }
    }
    //释放空间
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;

}

void post_to_INSERT_N(std::istream& is, std::istream& ns, std::ostream& os)
{
    size_t seq_num = 0; //序列条数
    size_t index;
    size_t num;
    size_t setwidth = 0;
    size_t len_N, s_index;
    int i;
    bool tag;
    std::string namei;  //每条序列的name
    std::string each_line;
    std::string each_sequence;

    ns >> seq_num;
    std::vector <std::string> names(seq_num);
    std::vector <size_t> length(seq_num);
    std::vector<std::vector<Insertion>> N_insertions(seq_num);
    std::vector<std::vector<Insertion>> insertions(seq_num);
    std::vector<std::list<unsigned char>> sequences(seq_num);


    //N
    for (i = 0; i < seq_num; i++)
    {
        for (; std::getline(ns, each_line); )
        {
            if (each_line[0] == '>')
            {
                each_line.erase(0, 1);
                if ((int)(*each_line.rbegin()) == 13)
                    each_line.pop_back();
                each_line.erase(0, each_line.find_first_not_of(" "));
                each_line.erase(each_line.find_last_not_of(" ") + 1);
                names[i].assign(each_line);
                ns >> len_N >> length[i];
                if (len_N != 0)N_insertions[i].resize(len_N);
                for (int j = 0; j < len_N; j++)//名字
                    ns >> N_insertions[i][j].index >> N_insertions[i][j].number;
                break;
            }
        }
    }
    for (i = 0; i < seq_num; i++) 
        sequences[i].resize(length[i]);

     //GAP
    each_line.clear();
    each_sequence.clear();
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);

            if (flag)
            {
                //each_sequence 
                i = find(names.begin(), names.end(), namei) - names.begin();
                if (i == seq_num)
                {
                    std::cout << "Error: file name mismatch!\n";
                    exit(-1);
                }
                index = 0; num = 0; tag = true;
                auto itor = sequences[i].begin();
                for (int j = 0; j < each_sequence.size(); j++)
                {
                    if ((int)each_sequence[j] == 45)
                    {
                        if (tag)//首个-
                            tag = false;
                        num++;
                    }
                    else
                    {
                        if (std::distance(sequences[i].begin(), itor) == length[i])
                        {
                            std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                            exit(-1);
                        }
                        *itor++ = (unsigned char)each_sequence[j];
                        if (num != 0)
                        {
                            insertions[i].push_back(Insertion({ index, num }));
                            tag = true;
                            num = 0;
                        }
                        index++;
                    }
                }
                if (num != 0) insertions[i].push_back(Insertion({ index, num }));
                if (std::distance(sequences[i].begin(), itor) != length[i])
                {
                    std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                    exit(-1);
                }
                each_sequence.clear();
            }

            flag = true;
            namei.assign(each_line);
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    //last_sequence 
    i = find(names.begin(), names.end(), namei) - names.begin();
    if (i == seq_num)
    {
        std::cout << "Error: file name mismatch!\n";
        exit(-1);
    }
    index = 0; num = 0; tag = true;
    auto itor = sequences[i].begin();
    for (int j = 0; j < each_sequence.size(); j++)
    {
        if ((int)each_sequence[j] == 45)
        {
            if (tag)//首个-
                tag = false;
            num++;
        }
        else
        {
            if (std::distance(sequences[i].begin(), itor) == length[i])
            {
                std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                exit(-1);
            }
            *itor++ = (unsigned char)each_sequence[j];
            if (num != 0)
            {
                insertions[i].push_back(Insertion({ index, num }));
                tag = true;
                num = 0;
            }
            index++;
        }
    }
    if (num != 0) insertions[i].push_back(Insertion({ index, num }));
    if (std::distance(sequences[i].begin(), itor) != length[i])
    {
        std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
        exit(-1);
    }
    each_sequence.clear();

    //插入
    insertion_gap_N(sequences, insertions, N_insertions);

    //写出
    for (i = 0; i < seq_num; i++)
    {
        os << "> " << names[i] << "\n";
        for (auto iter = sequences[i].begin(); iter != sequences[i].end(); iter++)
            os << (char)*iter;
        os << "\n";
    }
}

void direct_to_INSERT_N(std::istream& is, std::istream& ios, std::ostream& os)
{
    size_t seq_num = 0; //序列条数
    size_t index;
    size_t num;
    int i;
    bool tag;
    std::string namei;  //每条序列的name
    std::string each_line;
    std::string each_sequence;

    std::vector<std::vector<Insertion>> N_insertions;
    std::vector<Insertion> n_gap;
    std::vector <std::string> names;
    std::vector <size_t> length;

    //N
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);

            if (flag)
            {
                index = 0; num = 0; tag = true;
                for (int j = 0; j < each_sequence.size(); j++)
                {
                    if ((int)each_sequence[j] == 78 || (int)each_sequence[j] == 110)
                    {
                        if (tag)//首个N
                            tag = false;
                        num++;
                    }
                    else
                    {
                        if (num != 0)
                        {
                            n_gap.push_back(Insertion({ index, num }));
                            tag = true;
                            num = 0;
                        }
                        index++;
                    }
                }if (num != 0) n_gap.push_back(Insertion({ index, num }));
                N_insertions.push_back(n_gap);
                length.push_back(index);
                std::vector<Insertion>().swap(n_gap);
                each_sequence.clear();
            }

            flag = true;
            namei.assign(each_line);
            names.push_back(namei);
            seq_num++;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }
    index = 0; num = 0; tag = true;
    for (int j = 0; j < each_sequence.size(); j++)
    {
        if ((int)each_sequence[j] == 78 || (int)each_sequence[j] == 110)
        {
            if (tag)//首个N
                tag = false;
            num++;
        }
        else
        {
            if (num != 0)
            {
                n_gap.push_back(Insertion({ index, num }));
                tag = true;
                num = 0;
            }
            index++;
        }
    }
    if (num != 0) n_gap.push_back(Insertion({ index, num }));
    length.push_back(index);
    N_insertions.push_back(n_gap);
    each_sequence.clear();
    
    std::vector<std::vector<Insertion>> insertions(seq_num);
    std::vector<std::list<unsigned char>> sequences(seq_num);
    for (i = 0; i < seq_num; i++)
        sequences[i].resize(length[i]);

    
    //GAP
    each_line.clear();
    each_sequence.clear();
    for (bool flag = false; std::getline(ios, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);

            if (flag)
            {
                //each_sequence 
                i = find(names.begin(), names.end(), namei) - names.begin();
                if (i == seq_num)
                {
                    std::cout << "Error: file name mismatch!\n";
                    exit(-1);
                }
                index = 0; num = 0; tag = true;
                auto itor = sequences[i].begin();
                for (int j = 0; j < each_sequence.size(); j++)
                {
                    if ((int)each_sequence[j] == 45)
                    {
                        if (tag)//首个-
                            tag = false;
                        num++;
                    }
                    else
                    {
                        if (std::distance(sequences[i].begin(), itor) == length[i])
                        {
                            std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                            exit(-1);
                        }
                        *itor++ = (unsigned char)each_sequence[j];
                        if (num != 0)
                        {
                            insertions[i].push_back(Insertion({ index, num }));
                            tag = true;
                            num = 0;
                        }
                        index++;
                    }
                }
                if (num != 0) insertions[i].push_back(Insertion({ index, num }));
                if (std::distance(sequences[i].begin(), itor) != length[i])
                {
                    std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                    exit(-1);
                }
                each_sequence.clear();
            }

            flag = true;
            namei.assign(each_line);
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    //last_sequence 
    i = find(names.begin(), names.end(), namei) - names.begin();
    if (i == seq_num)
    {
        std::cout << "Error: file name mismatch!\n";
        exit(-1);
    }
    index = 0; num = 0; tag = true;
    auto itor = sequences[i].begin();
    for (int j = 0; j < each_sequence.size(); j++)
    {
        if ((int)each_sequence[j] == 45)
        {
            if (tag)//首个-
                tag = false;
            num++;
        }
        else
        {
            if (std::distance(sequences[i].begin(), itor) == length[i])
            {
                std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
                exit(-1);
            }
            *itor++ = (unsigned char)each_sequence[j];
            if (num != 0)
            {
                insertions[i].push_back(Insertion({ index, num }));
                tag = true;
                num = 0;
            }
            index++;
        }
    }
    if (num != 0) insertions[i].push_back(Insertion({ index, num }));
    if (std::distance(sequences[i].begin(), itor) != length[i])
    {
        std::cout << "Error: The \"" << names[i] << "\" sequence lengths of the two files do not match!\n";
        exit(-1);
    }
    each_sequence.clear();


    //插入
    insertion_gap_N(sequences, insertions, N_insertions);

    for (i = 0; i < seq_num; i++)
    {
        os << "> " << names[i] << "\n";
        for (auto iter = sequences[i].begin(); iter != sequences[i].end(); iter++)
            os << (char)*iter;
        os << "\n";
    }
}

void cout_help()
{
    std::cout << "-----There are three functions-----\n";
    std::cout << "1 Remove N and record. \n";
    std::cout << "        main.exe -1 in_original.fasta out_removeN.fasta out_recordN.tmp\n";
    std::cout << "\n";
    std::cout << "2 Insert N back into the result file. \n";
    std::cout << "        main.exe -2 in_ans_withoutN.fasta out_ans_withN.fasta in_recordN.tmp\n";
    std::cout << "\n";
    std::cout << "3 Generate the final result file with N from the source file and the result file. \n";
    std::cout << "        main.exe -3 in_original.fasta in_ans_withoutN.fasta out_ans_withN.fasta\n";
    std::cout << "\n";
}

int main(int argc, char* argv[])
{
    if (argc == 2)
    {
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0)
            cout_help();
        else std::cout << "Error parameters!\n";
        exit(0);
    }
    if (argc != 5)
    {
        std::cout << "Incorrect number( "<< argc <<" ) of parameters!\n";
        exit(-1);
    }
    if(argv[1][0]!='-')
    {
        std::cout << "Error parameters!\n";
        std::cout << argv[1] << "\n";
        exit(-1);
    }
    int tag = atoi(argv[1]+1);
    std::string in_file_name(argv[2]);
    std::string out_file_name(argv[3]);
    std::string N_file_name(argv[4]);

    const auto start_point = std::chrono::high_resolution_clock::now(); //记录起始时间
    if (tag == 1)//1i,2o
    {
        std::cout << "***** Start to perform function 1 *****\n";
        std::ifstream is(in_file_name, std::ios::binary | std::ios::in); //判断输入路径合法否
        if (!is)
        {
            std::cout << "cannot access file " << in_file_name << '\n';
            exit(0);
        }
        std::ofstream os(out_file_name, std::ios::binary | std::ios::out);
        if (!os)
        {
            std::cout << "cannot write file " << out_file_name << '\n';
            exit(0);
        }
        std::ofstream ns(N_file_name, std::ios::binary | std::ios::out);
        if (!ns)
        {
            std::cout << "cannot write file " << N_file_name << '\n';
            exit(0);
        }
        pre_to_REMOVE_N(is, os, ns);
        is.close();
        os.close();
        ns.close();
    }
    else if (tag == 2)//2i,1o
    {
        std::cout << "***** Start to perform function 2 *****\n";
        std::ifstream is(in_file_name, std::ios::binary | std::ios::in); //IN
        if (!is)
        {
            std::cout << "cannot access file " << in_file_name << '\n';
            exit(0);
        }
        std::ifstream ns(N_file_name, std::ios::binary | std::ios::in);//withoutN
        if (!ns)
        {
            std::cout << "cannot access file " << N_file_name << '\n';
            exit(0);
        }
        std::ofstream os(out_file_name, std::ios::binary | std::ios::out); //withN
        if (!os)
        {
            std::cout << "cannot write file " << out_file_name << '\n';
            exit(0);
        }
        post_to_INSERT_N(is, ns, os);
        is.close();
        ns.close();
        os.close();
    }
    else if (tag == 3)//2i,1o
    {
        std::cout << "***** Start to perform function 3 *****\n";
        std::ifstream is(in_file_name, std::ios::binary | std::ios::in); //判断输入路径合法否
        if (!is)
        {
            std::cout << "cannot access file " << in_file_name << '\n';
            exit(0);
        }
        std::ifstream ios(out_file_name, std::ios::binary | std::ios::in);
        if (!ios)
        {
            std::cout << "cannot access file " << out_file_name << '\n';
            exit(0);
        }
        std::ofstream os(N_file_name, std::ios::binary | std::ios::out);
        if (!os)
        {
            std::cout << "cannot write file " << N_file_name << '\n';
            exit(0);
        }
        direct_to_INSERT_N(is, ios, os);
        is.close();
        ios.close();
        os.close();
    }
    else
    {
        std::cout << "error number!\n";
        exit(-1);
    }
    std::cout << "***** Program execution completed *****\n";
    std::cout << "time consumes    : " << (int)chrono::duration <double, milli>((std::chrono::high_resolution_clock::now() - start_point)).count() << "ms \n"; //输出耗费时间
    GetMemoryUsage();
    return 0;
}
