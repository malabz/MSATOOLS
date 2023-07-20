#include "Utils.hpp"
#include "Pseudo.hpp"
#include "Fasta.hpp"
#include "Arguments.hpp"

#include <stdio.h>
#include <regex>
#include <iostream>
#include <limits>
#include <iomanip>
#include <list>
#include <fstream>

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

#if defined(_WIN32)
void getFiles_win(std::string path, std::vector<std::string>& files)
{
    //文件句柄
    intptr_t hFile = 0;
    //文件信息
    struct _finddata_t fileinfo;
    std::string p;
    if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
    {
        do
        {
            //如果是目录,迭代之
            //如果不是,加入列表
            if ((fileinfo.attrib & _A_SUBDIR))
            {
                if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
                    getFiles_win(p.assign(path).append("\\").append(fileinfo.name), files);
            }
            else
            {
                files.push_back(p.assign(path).append(fileinfo.name));
            }
        } while (_findnext(hFile, &fileinfo) == 0);
        _findclose(hFile);
    }
}
#elif defined(__unix__) || defined(__unix) || defined(unix)
void getFiles_linux(std::string path, std::vector<std::string>& filenames)
{
    DIR* pDir;
    struct dirent* ptr;
    if (!(pDir = opendir(path.c_str()))) {
        std::cout << "Folder doesn't Exist!" << std::endl;
        return;
    }
    while ((ptr = readdir(pDir)) != 0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0) {
            filenames.push_back(path + ptr->d_name);
        }
    }
    closedir(pDir);
}
#endif

int NSCORE = 0;
int HOXD70[6][6] = { {},{0,91,-114,-31,-123,NSCORE},{0,-114,100,-125,-31,NSCORE},{0,-31,-125,100,-114,NSCORE},
    {0,-123,-31,-114,91,NSCORE},{0,NSCORE,NSCORE,NSCORE,NSCORE,NSCORE} };//评分矩阵
int cs[8] = { 0,91,100,100,91,0,0,0 };
int d = 400, e = 30; //首个gap和后续gap
int stop_g = 5;

std::string utils::remove_white_spaces(const std::string& str) //去除空格
{
    static const std::regex white_spaces("\\s+");
    return std::regex_replace(str, white_spaces, "");
}

std::vector<unsigned char> utils::to_pseudo(const std::string& str)//预处理，char->int
{
    std::vector<unsigned char> pseu;
    pseu.reserve((size_t)(str.size()+1));

    for (auto i : str) pseu.push_back(to_pseudo(i));
    return pseu;
}

std::string utils::from_pseudo(const std::vector<unsigned char>& pseu)//后处理，int->char
{
    static constexpr char map[nucleic_acid_pseudo::NUMBER]{ '-', 'c', 'g', 'a', 't', 'n' };

    std::string str;
    str.reserve(pseu.size());

    for (auto i : pseu) str.push_back(map[i]);
    return str;
}

unsigned char* _get_map() //统一DNA/RNA 及大小写
{
    using namespace nucleic_acid_pseudo;

    static unsigned char map[std::numeric_limits<unsigned char>::max()];
    memset(map, N, sizeof(map));

    // map['-'] = GAP; // we could not process sequences with '-'
    map['c'] = map['C'] = C;
    map['g'] = map['G'] = G;
    map['a'] = map['A'] = A;
    map['t'] = map['T'] = map['u'] = map['U'] = T;

    return map;
}

static const unsigned char* _map = _get_map();

unsigned char utils::to_pseudo(char ch)  //预处理，char->int
{
    return _map[ch];
}

std::vector<std::vector<unsigned char>> utils::read_to_pseudo(std::istream& is, std::vector<std::string>& name,
    std::vector<std::vector<utils::Insertion>>& N_gap) //读入数据
{
    std::vector<std::vector<unsigned char>> sequences;

    std::string each_line;
    std::string each_sequence;
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            //each_line.erase(std::remove(each_line.begin(), each_line.end(), ' '), each_line.end());
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);
            replace(each_line.begin(), each_line.end(), ' ', '-');
            name.push_back(each_line);
            
            if (flag)
            {
                sequences.push_back(to_pseudo(each_sequence));
                each_sequence.clear();
            }
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    sequences.push_back(to_pseudo(each_sequence));

    /*for (int i = 0; i < sequences.size(); i++)
    {
        for (int j = 0; j < sequences[i].size(); j++)
            std::cout << (int)sequences[i][j] << " ";
        std::cout << "\n";
    }*/
    //去N并记录
    std::vector<utils::Insertion > n_gap;
    size_t index = 0, num = 0, ai = 0;
    bool tag;
    for (int i = 0; i < sequences.size(); i++)
    {
        index = 0; num = 0; tag = true;
        for (int j = 0; j < sequences[i].size(); j++)
        {
            if ((int)sequences[i][j] == 5)
            {
                if (tag)//首个N
                    tag = false;
                num++;
            }
            else
            {
                if (num != 0)
                {
                    n_gap.push_back(utils::Insertion({ index, num }));
                    tag = true;
                    num = 0;
                }
                index++;
            }
        }
        if (num != 0) n_gap.push_back(utils::Insertion({ index, num }));
        N_gap.push_back(n_gap);
        std::vector<utils::Insertion>().swap(n_gap);
        sequences[i].erase(std::remove(sequences[i].begin(), sequences[i].end(), '\5'), sequences[i].end());
    }

    /*for (int i = 0; i < sequences.size(); i++)
    {
        for (int j = 0; j < sequences[i].size(); j++)
            std::cout << (int)sequences[i][j]<< " ";
        std::cout << "\n";
    }*/
    return sequences;  //返回sequences序列数据
}

void utils::read_to_pseudo(std::vector<std::vector<unsigned char>> &sequences,
    std::istream& is, std::vector<std::size_t>& all_Length, std::vector<std::size_t> &Length, std::vector<std::string>& name,
    std::vector<std::vector<utils::Insertion>>& N_gap) //读入数据
{
    std::string each_line;
    std::string each_sequence;
    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            //each_line.erase(std::remove(each_line.begin(), each_line.end(), ' '), each_line.end());
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);
            //replace(each_line.begin(), each_line.end(), ' ', '-');
            name.push_back(each_line);

            if (flag)
            {
                sequences.push_back(to_pseudo(each_sequence));
                each_sequence.clear();
            }
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }
    sequences.push_back(to_pseudo(each_sequence));
    //去N并记录
    std::vector<utils::Insertion> n_gap;
    size_t index = 0, num = 0, ai = 0;
    bool tag;
    //N-gap
    for (int i = 0; i < sequences.size(); i++)
    {
        index = 0; num = 0; tag = true;
        for (int j = 0; j < sequences[i].size(); j++)
        {
            if ((int)sequences[i][j] == 5)
            {
                if (tag)//首个N
                    tag = false;
                num++;
            }
            else
            {
                if (num != 0)
                {
                    n_gap.push_back(utils::Insertion({ index, num }));
                    tag = true;
                    num = 0;
                }
                index++;
            }
        }
        if (num != 0) 
        { 
            n_gap.push_back(utils::Insertion({ index, num }));
        }
        N_gap.push_back(n_gap);
        std::vector<utils::Insertion>().swap(n_gap);
        all_Length.push_back(sequences[i].size());
        sequences[i].erase(std::remove(sequences[i].begin(), sequences[i].end(), '\5'), sequences[i].end());
        Length.push_back(sequences[i].size());
    }
}

bool utils::next_reads(utils::reads & read, std::istream& is) //读入数据
{
    std::string each_line;
    std::string each_sequence;

    while (std::getline(is, each_line))
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;
            
        if (each_line[0] == '@')
        {
            each_line.erase(0, 1);
            if ((int)(*each_line.rbegin()) == 13)
                each_line.pop_back();
            each_line.erase(0, each_line.find_first_not_of(" "));
            each_line.erase(each_line.find_last_not_of(" ") + 1);
            replace(each_line.begin(), each_line.end(), ' ', '-');
            read.name = each_line;
            break;
        }
        else
            exit(-2);
    }
    if (each_line.size() != 0) 
    {
        std::getline(is, each_sequence);
        if ((int)(*each_sequence.rbegin()) == 13)
            each_sequence.pop_back();
        is.seekg(each_sequence.size() + 3, std::ios::cur);
    }

    if (each_sequence.size() == 0)
        return false;
    read.all_length = each_sequence.size();
    read.read = to_pseudo(each_sequence);
    //去N并记录
    std::vector<utils::Insertion>().swap(read.n_gap);
    size_t index = 0, num = 0, ai = 0;
    bool tag = true;
    //N-gap
    for (int j = 0; j < read.read.size(); j++)
    {
        if ((int)read.read[j] == 5)
        {
            if (tag)//首个N
                tag = false;
            num++;
        }
        else
        {
            if (num != 0)
            {
                read.n_gap.push_back(utils::Insertion({ index, num }));
                tag = true;
                num = 0;
            }
            index++;
        }
    }
    if (num != 0) read.n_gap.push_back(utils::Insertion({ index, num }));
    read.read.erase(std::remove(read.read.begin(), read.read.end(), '\5'), read.read.end());
    read.length = read.read.size();
    
    read.read_ni.resize(read.length);
    for (int j = 0, k = read.length - 1; j < read.length; j++, k--)
        read.read_ni[k] = 5 - (int)read.read[j];
    return true;
}


//按拉长fasta输出,按原文件插入old_chao
void utils::insert_and_write(std::ostream& os, std::istream& is, const std::vector<std::vector<Insertion>>& insertions) //写回比对结果
{
    const size_t sequence_number = insertions.size();

    std::string each_line;
    std::string each_sequence;
    std::string each_sequence_aligned;
    for (unsigned count = 0, length, flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
            continue;

        if (each_line[0] == '>')
        {
            if (flag)
            {
                if (count == 0)
                {
                    length = each_sequence.size();
                    for (auto insertion : insertions[0])
                        length += insertion.number;

                    each_sequence_aligned.reserve(length);
                    each_sequence.reserve(length);
                }

                utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
                    insertions[count].cbegin(), insertions[count].cend(), std::back_inserter(each_sequence_aligned), '-');

                if (arguments::output_matrix)
                    os << each_sequence_aligned;
                else
                    Fasta::cut_and_write(os, each_sequence_aligned);
                os << '\n';

                each_sequence.clear();
                each_sequence_aligned.clear();
                ++count;
            }

            if (arguments::output_matrix == false)
                os << each_line << '\n';
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
    }

    utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
        insertions.back().cbegin(), insertions.back().cend(), std::back_inserter(each_sequence_aligned), '-');

    if (arguments::output_matrix)
        os << each_sequence_aligned;
    else
        Fasta::cut_and_write(os, each_sequence_aligned);
    os << '\n';
}

void insert_others(utils::Insertion2 That, utils::Insertion2 This, std::vector<std::vector<utils::Insertion2>>& more_insertions,
    int k, int sequence_number, int* ii)
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

//对于小内存，开不了大链表，按拉长fasta输出,按文件写回来插入new_zhou,完成后直接结束
void utils::insert_and_write_file(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
    std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions,
    std::vector<std::string>& name, std::vector<bool>& sign)
    //写回比对结果
{
    //新版插入，带N
    std::string each_sequence, each_line;
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    int i = 0, j = 0, pre = 0, diff, mi, k;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);
    char chars[7] = { 'A','C','G','T','N','N','-' };
    size_t g_num;

    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
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
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
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
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
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
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size()) //后续多余的gap，继续插入
        {
            g_num += insertions[k][i].number;
            i_all_insertions.push_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.push_back(i_all_insertions);
        /*for (int m = 0; m < sequence_number; m++)
        {
            for (int d = 0; d < insertions[m].size(); d++)
            {
                std::cout << insertions[m][d].index << " " << insertions[m][d].number << "\n";
            }
            std::cout << "\n";
        }
            std::cout << k<<" \n";*/
    }

    /*for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < insertions[k].size(); d++)
        {
            std::cout << insertions[k][d].index << " " << insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < N_insertions[k].size(); d++)
        {
            std::cout << N_insertions[k][d].index << " " << N_insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < all_insertions[k].size(); d++)
        {
            std::cout << all_insertions[k][d].index << " " << all_insertions[k][d].n_num<<" "<< all_insertions[k][d].gap_num << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < more_insertions[k].size(); d++)
        {
            std::cout << more_insertions[k][d].index << " " << more_insertions[k][d].n_num << " " << more_insertions[k][d].gap_num << "\n";
        }
        std::cout << k<<"\n";
    }*/

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
        if (sign[i]) os << "> " << name[i] << " + " << "\n";
        else os << "> " << name[i] << " - " << "\n";
        std::ofstream tmpo(arguments::tmp_file_name, std::ios::binary | std::ios::out);
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmpo << chars[(int)(sequences[i][k++]) - 1];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmpo << 'N';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmpo << '-';
        }while (k < sequences[i].size())tmpo << chars[(int)(sequences[i][k++]) - 1];
        tmpo.close();
        //read
        std::ifstream tmpi(arguments::tmp_file_name, std::ios::binary | std::ios::in);
        while (std::getline(tmpi, each_line))
        {
            if (each_line.size() == 0 || (each_line.size() == 1 && (int)each_line[0] == 13)) //跳过空行
                continue;
            each_sequence += each_line;
            if ((int)(*each_line.rbegin()) == 13)
                each_sequence.pop_back();
        }
        tmpi.close();
        sequences[i].assign(each_sequence.begin(), each_sequence.end());
        each_sequence.clear();
        each_line.clear();
        //second insert
        mi = 0;
        k = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (k < more_insertions[i][j].index)
                os << sequences[i][k++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                os << 'N';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    os << '-';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    os << '-';
        }while (k < sequences[i].size()) os << sequences[i][k++];
        std::vector<unsigned char>().swap(sequences[i]);
        os << "\n";
    }os << "\n";
    remove(arguments::tmp_file_name.data());
    //释放空间
    /*std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;*/
}

//对于小内存，不可以开大链表，按重写入vector来插入N和gap
int* utils::vector_insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences, 
    std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions) //写回比对结果
{
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    std::vector<unsigned char> tmp_vector;
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;
    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            len_sequences[i] += N_insertions[i][j].number;
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
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
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
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
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
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size()) //后续多余的gap，继续插入
        {
            g_num += insertions[k][i].number;
            i_all_insertions.push_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.push_back(i_all_insertions);
        /*for (int m = 0; m < sequence_number; m++)
        {
            for (int d = 0; d < insertions[m].size(); d++)
            {
                std::cout << insertions[m][d].index << " " << insertions[m][d].number << "\n";
            }
            std::cout << "\n";
        }
            std::cout << k<<" \n";*/
    }

    /*for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < insertions[k].size(); d++)
        {
            std::cout << insertions[k][d].index << " " << insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < N_insertions[k].size(); d++)
        {
            std::cout << N_insertions[k][d].index << " " << N_insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < all_insertions[k].size(); d++)
        {
            std::cout << all_insertions[k][d].index << " " << all_insertions[k][d].n_num<<" "<< all_insertions[k][d].gap_num << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < more_insertions[k].size(); d++)
        {
            std::cout << more_insertions[k][d].index << " " << more_insertions[k][d].n_num << " " << more_insertions[k][d].gap_num << "\n";
        }
        std::cout << k<<"\n";
    }*/

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
    int more = 0, all_size = 0, ti = 0, k = 0;
    for (i = 0; i < sequence_number; i++)
    {
        if (i == 0)
        {
            more = sequences[i].size();
            for (j = 0; j < all_insertions[i].size(); j++)
                more += (all_insertions[i][j].n_num + all_insertions[i][j].gap_num);
            all_size = more;
            mi = 0;
            for (j = 0; j < more_insertions[i].size(); j++)
            {
                all_size += more_insertions[i][j].n_num;
                if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
                    all_size -= multi[mi++].number;
                all_size += more_insertions[i][j].gap_num;
            }
            tmp_vector.resize(more);
        }

        ti = 0;
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmp_vector[ti++] = sequences[i][k++];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmp_vector[ti++] = '\5';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmp_vector[ti++] = '\7';
        }while (k < sequences[i].size())tmp_vector[ti++] = sequences[i][k++];

        sequences[i].resize(all_size);
        mi = 0;
        k = 0;
        ti = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (ti < more_insertions[i][j].index)
                sequences[i][k++] = tmp_vector[ti++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                sequences[i][k++] = '\5';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    sequences[i][k++] = '\7';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    sequences[i][k++] = '\7';
        }while (ti < more) sequences[i][k++] = tmp_vector[ti++];

    }

    //释放空间
    std::vector<unsigned char>().swap(tmp_vector);
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;
    
    return len_sequences;
}


//对于大内存，可以开大链表，插入N和gap到内存
int* utils::insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences, 
    std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions)
{
    //新版插入，带N
    const size_t sequence_number = insertions.size();
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;
    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
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
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
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
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
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
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
            j++;
        }
        while (i < insertions[k].size()) //后续多余的gap，继续插入
        {
            g_num += insertions[k][i].number;
            i_all_insertions.push_back(Insertion2({ insertions[k][i].index, 0,insertions[k][i].number }));
            i++;
        }
        all_insertions.push_back(i_all_insertions);
        /*for (int m = 0; m < sequence_number; m++)
        {
            for (int d = 0; d < insertions[m].size(); d++)
            {
                std::cout << insertions[m][d].index << " " << insertions[m][d].number << "\n";
            }
            std::cout << "\n";
        }
            std::cout << k<<" \n";*/
    }

    /*for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < insertions[k].size(); d++)
        {
            std::cout << insertions[k][d].index << " " << insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < N_insertions[k].size(); d++)
        {
            std::cout << N_insertions[k][d].index << " " << N_insertions[k][d].number << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < all_insertions[k].size(); d++)
        {
            std::cout << all_insertions[k][d].index << " " << all_insertions[k][d].n_num<<" "<< all_insertions[k][d].gap_num << "\n";
        }
        std::cout << k << "\n";
    }
    for (int k = 0; k < sequence_number; k++)
    {
        for (int d = 0; d < more_insertions[k].size(); d++)
        {
            std::cout << more_insertions[k][d].index << " " << more_insertions[k][d].n_num << " " << more_insertions[k][d].gap_num << "\n";
        }
        std::cout << k<<"\n";
    }*/

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
        if (all_insertions[i].size() > 1000)
        {
            std::list<unsigned char> tmp(sequences[i].begin(), sequences[i].end());  //vector转list
            auto iter = begin(tmp);
            pre = all_insertions[i][all_insertions[i].size() - 1].index;
            std::advance(iter, pre);
            for (j = all_insertions[i].size() - 1; j >= 0; j--)
            {
                diff = all_insertions[i][j].index - pre;
                std::advance(iter, diff);
                tmp.insert(iter, all_insertions[i][j].n_num, '\5'); //与vector插入方式相反
                tmp.insert(iter, all_insertions[i][j].gap_num, '\7');
                pre = all_insertions[i][j].index + all_insertions[i][j].gap_num + all_insertions[i][j].n_num;
            }
            if (more_insertions[i].size() > 0)
            {
                iter = begin(tmp); mi = multi.size() - 1;
                pre = more_insertions[i][more_insertions[i].size() - 1].index;
                std::advance(iter, pre);
                for (j = more_insertions[i].size() - 1; j >= 0; j--)
                {
                    diff = more_insertions[i][j].index - pre;
                    std::advance(iter, diff);

                    tmp.insert(iter, more_insertions[i][j].n_num, '\5'); //与vector插入方式相反
                    if ((mi >= 0) && (more_insertions[i][j].index == multi[mi].index))
                    {
                        tmp.insert(iter, (more_insertions[i][j].gap_num - multi[mi].number), '\7');
                        pre = more_insertions[i][j].index + more_insertions[i][j].gap_num + more_insertions[i][j].n_num - multi[mi].number;
                        mi--;
                    }
                    else
                    {
                        tmp.insert(iter, more_insertions[i][j].gap_num, '\7');
                        pre = more_insertions[i][j].index + more_insertions[i][j].gap_num + more_insertions[i][j].n_num;
                    }
                }
            }
            sequences[i].assign(tmp.begin(), tmp.end());
        }
        else
        {
            for (j = all_insertions[i].size() - 1; j >= 0; j--)
            {
                sequences[i].insert(sequences[i].begin() + all_insertions[i][j].index, all_insertions[i][j].gap_num, '\7');//先插空，空在N后
                sequences[i].insert(sequences[i].begin() + all_insertions[i][j].index, all_insertions[i][j].n_num, '\5'); //后插N，N在空前
            }
            mi = multi.size() - 1;
            for (j = more_insertions[i].size() - 1; j >= 0; j--)
            {

                if (mi >= 0 && more_insertions[i][j].index == multi[mi].index)
                {
                    sequences[i].insert(sequences[i].begin() + more_insertions[i][j].index, more_insertions[i][j].gap_num - multi[mi].number, '\7');//先插空，空在N后
                    mi--;
                }
                else sequences[i].insert(sequences[i].begin() + more_insertions[i][j].index, more_insertions[i][j].gap_num, '\7');//先插空，空在N后
                sequences[i].insert(sequences[i].begin() + more_insertions[i][j].index, more_insertions[i][j].n_num, '\5'); //后插N，N在空前
            }
        }
    }


    /*
    char chars[7] = {'A','C','G','T','N','N','-'};
    for (i = 0; i < sequence_number; i++)
    {
        for (int k = 0; k < sequences[i].size(); k++) std::cout << chars[(int)(sequences[i][k]) - 1];
        std::cout << "\n";
    }
    std::cout << "?3\n";
    exit(2);*/
    //释放空间
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            len_sequences[i] += N_insertions[i][j].number;
    return len_sequences;
}

//按拉长fasta输出,按内存原串插入
void utils::insert_and_write_fasta(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
    std::vector<std::vector<Insertion>>& insertions, std::vector<std::vector<Insertion>>& N_insertions,
    std::vector<std::string>& name, std::vector<bool>& sign)
    //写回比对结果
{
    std::cout << "0memory usage: " << getPeakRSS() << " B" << std::endl;//输出内存耗费
       //插入到原串
    const size_t sequence_number = insertions.size();
    std::cout << "insertions " << sequence_number << " " << insertions[0].size() << " " << insertions[1].size() << "\n";
    int score = 0, length = 0, name_len = 0;

    char chars[7] = { 'A','C','G','T','N','N','-' };
    for (int i = 0; i < name.size(); i++)
        if (name_len < name[i].size())
            name_len = name[i].size();
    const auto align_start1 = std::chrono::high_resolution_clock::now(); //记录插入起始时间
/*
    //新版插入，插N插GAP
#if defined(_WIN32)//内存小
    len_sequences = vector_insertion_gap_N(sequences, insertions, N_insertions);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    len_sequences = insertion_gap_N(sequences, insertions, N_insertions);
#endif
    //原版插入，不插N
*/
//插入到原串
//*************insert_begin***********
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    int* len_sequences = new int[sequence_number];
    int* ii = new int[sequence_number]();
    int* score_two = NULL;
    std::vector<unsigned char> tmp_vector;
    int i = 0, j = 0, pre = 0, diff, mi;
    std::vector<std::vector<Insertion2>> more_insertions(sequence_number);

    size_t g_num;
    for (int k = 0; k < sequence_number; k++)  len_sequences[k] = sequences[k].size();
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            len_sequences[i] += N_insertions[i][j].number;
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
                    insert_others(Insertion2({ insertions[k][i].index + g_num ,N_insertions[k][j].number - insertions[k][i].number, 0 }), utils::Insertion2({ insertions[k][i].index + g_num ,0, N_insertions[k][j].number - insertions[k][i].number }), more_insertions, k, sequence_number, ii);
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
                insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
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
            insert_others(Insertion2({ N_insertions[k][j].index + g_num ,N_insertions[k][j].number, 0 }), utils::Insertion2({ N_insertions[k][j].index + g_num ,0, N_insertions[k][j].number }), more_insertions, k, sequence_number, ii);
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
    //释放空间
    for (i = 0; i < N_insertions[i].size(); i++)
    {
        std::vector<Insertion>().swap(N_insertions[i]);
        std::vector<Insertion>().swap(insertions[i]);
    }
    std::vector<std::vector<Insertion>>().swap(N_insertions);
    std::vector<std::vector<Insertion>>().swap(insertions);
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
    int more = 0, all_size = 0, ti = 0, k = 0;
    std::ofstream tmpo(arguments::score_file, std::ios::binary | std::ios::out);
    tmpo << "I" << "\t" << "Match_num" << "\t" << "Len_sequence" << "\t" << "M/L" << "\t" << "Name" << "\n";
    tmpo << 0 << "\t" << "--" << "\t" << len_sequences[0] << "\t" << "---" << "\t" << name[0] << "\n";
    for (i = 0; i < sequence_number; i++)
    {
        if (i == 0)
        {
            more = sequences[i].size();
            for (j = 0; j < all_insertions[i].size(); j++)
                more += (all_insertions[i][j].n_num + all_insertions[i][j].gap_num);
            all_size = more;
            mi = 0;
            for (j = 0; j < more_insertions[i].size(); j++)
            {
                all_size += more_insertions[i][j].n_num;
                if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
                    all_size -= multi[mi++].number;
                all_size += more_insertions[i][j].gap_num;
            }
            tmp_vector.resize(more);
        }

        ti = 0;
        k = 0;
        for (j = 0; j < all_insertions[i].size(); j++)
        {
            while (k < all_insertions[i][j].index)
                tmp_vector[ti++] = sequences[i][k++];
            for (int p = 0; p < all_insertions[i][j].n_num; p++)
                tmp_vector[ti++] = '\5';
            for (int p = 0; p < all_insertions[i][j].gap_num; p++)
                tmp_vector[ti++] = '\7';
        }
        while (k < sequences[i].size())tmp_vector[ti++] = sequences[i][k++];

        sequences[i].resize(all_size);
        mi = 0;
        k = 0;
        ti = 0;
        for (j = 0; j < more_insertions[i].size(); j++)
        {
            while (ti < more_insertions[i][j].index)
                sequences[i][k++] = tmp_vector[ti++];
            for (int p = 0; p < more_insertions[i][j].n_num; p++)
                sequences[i][k++] = '\5';
            if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            {
                for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                    sequences[i][k++] = '\7';
                mi++;
            }
            else
                for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                    sequences[i][k++] = '\7';
        }while (ti < more) sequences[i][k++] = tmp_vector[ti++];
        os << "> " << name[i];
        if (sign[i]) os << "(+)\n";
        else os << "(-)\n";
        for (k = 0; k < sequences[i].size(); k++) os << chars[(int)(sequences[i][k]) - 1];
        os << "\n";
        if (i != 0)
        {
            score_two = Compare_two(sequences[0], sequences[i], 0, 0, 1, 0);
            tmpo << i << "\t" << score_two[0] << "\t" << len_sequences[i] << "\t" << score_two[0] * 1.0 / len_sequences[i] << "\t" << name[i] << "\n";
            std::vector<unsigned char>().swap(sequences[i]);
        }
    }
    tmpo.close();
    std::vector<unsigned char>().swap(sequences[0]);
    std::vector<unsigned char>().swap(tmp_vector);
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    std::vector<Insertion>().swap(tmp);
    delete[] ii;

    //*************insert_end*************

    const auto align_start2 = std::chrono::high_resolution_clock::now(); //记录插入结束时间
    std::cout << align_start2 - align_start1 << "  insert_and_write\n";

}

//按完全匹配输出maf
void utils::insert_and_write_maf(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, int thresh) //写回比对结果
{
    const size_t sequence_number = insertions.size();
    std::cout << "insertions " << sequence_number << " " << insertions[0].size() << " " << insertions[1].size() << "\n";
    int* si = new int[sequence_number];
    int* len_sequences;
    int* score_two = NULL;
    for (int i = 0; i < sequence_number; i++) si[i] = 0;
    int score = 0, length = 0, ai = 0;
    int name_len = 0, pre = 0, diff;
    char chars[5] = { 'A','C','G','T','N' };
    for (int i = 0; i < name.size(); i++)
        if (name_len < name[i].size())
            name_len = name[i].size();
    const auto align_start1 = std::chrono::high_resolution_clock::now(); //记录比对起始时间

    //插入
    len_sequences = vector_insertion_gap_N(sequences, insertions, N_insertions);

    const auto align_start2 = std::chrono::high_resolution_clock::now(); //记录比对起始时间
    std::cout << align_start2 - align_start1 << "  insert\n";

    //求分数
    std::ofstream tmpo(arguments::score_file, std::ios::binary | std::ios::out);
    tmpo << "I" << "\t" << "Match_num" << "\t" << "Len_sequence" << "\t" << "M/L" << "\t" << "Name" << "\n";
    tmpo << 0 << "\t" << "--" << "\t" << len_sequences[0] << "\t" << "---" << "\t" << name[0] << "\n";
    for (int i = 0; i < sequence_number; i++)
    {
        if (i != 0)
        {
            score_two = Compare_two(sequences[0], sequences[i], 0, 0, 1, 0);
            tmpo << i << "\t" << score_two[0] << "\t" << len_sequences[i] << "\t" << score_two[0] * 1.0 / len_sequences[i] << "\t" << name[i] << "\n";
        }
    }
    tmpo.close();

    ai = 0;
    bool tag = false;
    unsigned char first;
    std::string zf;
    int i = 0, start = 0;

    std::cout << "all_length = " << sequences[0].size() << "\n";
    while (ai < sequences[0].size())
    {
        first = sequences[0][ai];
        for (i = 0; i < sequence_number; i++) if ((int)sequences[i][ai] != 7) si[i]++;
        if ((int)first == 7)
            tag = false;
        else
        {
            for (i = 1; i < sequence_number; i++) if (sequences[i][ai] != first) { tag = false; break; }
            if (i == sequence_number)
            {
                tag = true;
                if (length == 0)
                    start = ai;
                length++;
            }
        }

        if (tag == false)
        {
            if (length >= thresh)
            {
                os << "a score=" << length << "\n";
                for (i = 0; i < sequence_number; i++)
                {
                    if (sign[i])zf = " + "; else zf = " - ";
                    if ((int)sequences[i][ai] != 7)
                        os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << si[i] - length - 1 << " " << length << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                    else
                        os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << si[i] - length << " " << length << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                    for (int k = start; k < start + length; k++) os << chars[(int)(sequences[i][k]) - 1];
                    os << "\n";
                }
                os << "\n";
            }
            length = 0;
        }
        ai++;
    }
}


//比对两条序列得分
int* utils::Compare_two(std::vector<unsigned char>& s1, std::vector<unsigned char>& s2, int d, int e, int m, int mis)//, int d = 3, int e = 1, int m = 1, int mis = -2
{
    int N = 0;
    int length = s1.size();
    int* score_two = new int[2];  //score_two[0]得分,score_two[1]长度
    score_two[1] = length;
    if (length != s2.size())//长度不同，抛出异常
    {
        std::cerr << length << " != " << s2.size() << " !the length of s1 and s2 is wrong!\n";
        exit(1);
    }
    score_two[0] = 0.0;// score_two[0]得分
    bool gap1 = false;   //false代表首个gap
    for (int i = 0; i < length; i++)
    {
        if (((int)s1[i] != 7) && ((int)s2[i] != 7))//都不是空格
        {
            if (gap1) gap1 = false; //下一次为首个空格
            if (((int)s1[i] == 5) || ((int)s2[i] == 5)) score_two[0] += N; //匹配N
            else if (s1[i] == s2[i]) score_two[0] += m; //匹配
            else score_two[0] += mis;   // 不匹配
        }
        else if (((int)s1[i] != 7) || ((int)s2[i] != 7))
        {
            if (gap1 == false)//首个空格罚分
            {
                score_two[0] -= d; gap1 = true;//下次非首个空格
            }
            else score_two[0] -= e;//非首个空格罚分
        }
    }
    //score_two[2] = score_two[0]*1.0/length;
    //std::cout << "score: " << (int)score_two[0] << " / " << (int)score_two[1] << " = " << score_two[2] << "\n";
    return score_two;   //返回比对得分
}

//自定义输出依赖函数
typedef struct
{
    int index;
    int value;
}sort_st;
inline bool compare(sort_st a, sort_st b)
{
    return a.value < b.value; //升序排列，如果改为return a.value<b.value，则为降序
}
std::vector<int> index_100(std::vector<int>& nscore, int score100)
{
    int i;
    float sum = 0;
    std::vector<int> ans;
    //std::cout << "? " << score100 << " " << nscore[0] << " " << nscore[1] << " \n";
    for (i = 0; i < nscore.size(); i++)
        if (nscore[i] >= score100)
            break;
    if (i == nscore.size()) return ans;

    std::vector <sort_st> sort_array(nscore.size());
    for (int i = 0; i < nscore.size(); ++i) {
        sort_array[i].index = i;
        sort_array[i].value = nscore[i];
        sum += nscore[i];
    }
    sort(sort_array.begin(), sort_array.end(), compare);
    while (1)
    {
        if (sum / sort_array.size() >= score100)
            break;
        else
        {
            sum -= sort_array[0].value;
            sort_array.erase(sort_array.begin());
        }
    }
    for (int i = 0; i < sort_array.size(); ++i) ans.push_back(sort_array[i].index);
    sort(ans.begin(), ans.end());
    return ans;

}
//按用户自定义输出maf
/*void utils::insert_and_write_100_maf(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, int thresh1, int thresh2, int thresh3)//长度，条数，分数
{

    const size_t sequence_number = insertions.size();
    std::cout << "insertions " << sequence_number << " " << insertions[0].size() << " " << insertions[1].size() << "\n";
    int* si = new int[sequence_number];
    int* pre_si = new int[sequence_number];
    int* start = new int[sequence_number];
    bool* gap_tag = new bool[sequence_number];
    bool new_tag = true;
    int* len_sequences;
    int* gap_num = new int[sequence_number];//当前block gap数量
    int* pre_gap_num = new int[sequence_number];//前面几个block gap数量
    for (int i = 0; i < sequence_number; i++) { si[i] = gap_num[i] = pre_gap_num[i] = 0; gap_tag[i] = true; }
    int  length = 0, ai = 0;
    float score = 0;
    float ave_score = 0;
    int name_len = 0, pre = 0, diff;
    char chars[7] = { 'A','C','G','T','N','N','-' };
    for (int i = 0; i < name.size(); i++)
        if (name_len < name[i].size())
            name_len = name[i].size();
    const auto align_start1 = std::chrono::high_resolution_clock::now(); //插入前

    //插入
#if defined(_WIN32)//内存小
    len_sequences = vector_insertion_gap_N(sequences, insertions, N_insertions);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    len_sequences = insertion_gap_N(sequences, insertions, N_insertions);
#endif  

    const auto align_start2 = std::chrono::high_resolution_clock::now(); //插入后
    std::cout << align_start2 - align_start1 << "  insert\n";

    //求分数
    int* score_two = Compare_two(sequences[0], sequences[1]);
    const auto align_start3 = std::chrono::high_resolution_clock::now(); //求分后
    std::cout << align_start3 - align_start2 << "  score\n";
    std::cout << "score: " << score_two[0] << " / " << score_two[1] << " = " << score_two[0] * 1.0 / score_two[1] << "\n";

    std::string zf;
    int i = 0, j, leave = sequences[0].size() % 100, hundred = sequences[0].size() / 100, hi = 0, score100;
    std::cout << "all_length = " << sequences[0].size() << "\n";
    ai = 0;
    bool tag = false;


    std::vector<float> sum_score;
    std::vector<int> sum_index;
    std::vector<int> index;
    std::vector<int> nscore(sequence_number);
    for (hi = 0; hi < hundred; hi++)  //处理整百
    {
        score100 = 0;//后面可以 thresh%*score100
        for (i = 0; i < sequence_number; i++) {
            nscore[i] = 0;
            pre_si[i] = si[i];
            pre_gap_num[i] += gap_num[i];
            gap_num[i] = 0;
        }//每100个分数初始化0

        for (ai = hi * 100; ai < (hi + 1) * 100; ai++)
        {
            for (i = 0; i < sequence_number; i++) if ((int)sequences[i][ai] != 7) si[i]++;
            score100 += cs[(int)sequences[0][ai]];
            for (i = 0; i < sequence_number; i++)
            {
                if (((int)sequences[i][ai] == 7) || ((int)sequences[0][ai] == 7))
                {
                    if ((int)sequences[i][ai] == 7) gap_num[i]++;
                    if (gap_tag[i]) { nscore[i] -= d; gap_tag[i] = false; }
                    else  nscore[i] -= e;
                }
                else { nscore[i] += HOXD70[(int)sequences[0][ai]][(int)sequences[i][ai]]; gap_tag[i] = true; }
            }
        }
        //std::cout << score100 << "\n";
        std::vector<int>().swap(index);
        index = index_100(nscore, score100 * thresh3 / 100);
        if (index.size() != 0)
        {
            ave_score = 0;
            for (i = 0; i < index.size(); i++) ave_score += nscore[index[i]];
            ave_score /= index.size();
            if (sum_score.size() == 0)
            {
                //std::cout << hi << "\n";

                for (i = 0; i < sequence_number; i++) { start[i] = pre_si[i]; pre_gap_num[i] = 0; }
                sum_index.assign(index.begin(), index.end());
                sum_score.push_back(ave_score);
                continue;
            }
            else if (sum_index.size() == index.size())
            {
                for (i = 0; i < sum_index.size(); i++) if (sum_index[i] != index[i]) break;
                if (i == sum_index.size()) { sum_score.push_back(ave_score);  continue; }
            }
            //else
        }
        //else

        if (sum_score.size() == 0)//当前为废100，
        {
            for (int i = 0; i < sequence_number; i++) pre_gap_num[i] = gap_num[i] = 0;//把当前100gap数清掉;
            continue;
        }
        else
        {
            score = 0; length = sum_score.size() * 100;
            for (i = 0; i < sum_score.size(); i++)score += sum_score[i];
            score /= sum_score.size();

            if ((length >= thresh1) && (sum_index.size() >= thresh2))//长度，条数，分数
            {
                os << "a score=" << score << "\n";
                for (j = 0; j < sum_index.size(); j++)
                {
                    i = sum_index[j];
                    //if ((length - pre_gap_num[i]) < thresh1) continue;
                    if (sign[i])zf = " + "; else zf = " - ";
                    os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << start[i] << " " << length - pre_gap_num[i] << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                    for (int k = ai - length - 100; k < ai - 100; k++) os << chars[(int)(sequences[i][k]) - 1];
                    os << "\n";
                }
                os << "\n";
            }
            //写完了
            for (int i = 0; i < sequence_number; i++) pre_gap_num[i] = 0;//写完一组，清空pre gap数
            std::vector<float>().swap(sum_score);
            std::vector<int>().swap(sum_index);
            if (index.size() != 0)//从上面筛选下来，连不上，但是分数合格的首100
            {
                for (i = 0; i < sequence_number; i++) start[i] = pre_si[i];
                sum_index.assign(index.begin(), index.end());
                sum_score.push_back(ave_score);
                continue;
            }
        }
    }
    //处理后续的leave
    score100 = 0;//-d * leave / 100;
    for (i = 0; i < sequence_number; i++)
    {
        nscore[i] = 0;
        pre_si[i] = si[i];
        pre_gap_num[i] += gap_num[i];
        gap_num[i] = 0;
    }//预处理
    for (; ai < sequences[0].size(); ai++)
    {
        for (i = 0; i < sequence_number; i++) if ((int)sequences[i][ai] != 7) si[i]++;
        score100 += cs[(int)sequences[0][ai]];
        for (i = 0; i < sequence_number; i++)
        {
            if (((int)sequences[i][ai] == 7) || ((int)sequences[0][ai] == 7))
            {
                if ((int)sequences[i][ai] == 7) gap_num[i]++;
                if (gap_tag[i]) { nscore[i] -= d; gap_tag[i] = false; }
                else  nscore[i] -= e;
            }
            else { nscore[i] += HOXD70[(int)sequences[0][ai]][(int)sequences[i][ai]]; gap_tag[i] = true; }
        }
    }
    std::vector<int>().swap(index);
    index = index_100(nscore, score100 * thresh3 / 100);
    if (index.size() != 0)
    {
        ave_score = 0;
        for (i = 0; i < index.size(); i++) ave_score += nscore[index[i]]; //std::cout << "???" << index[i] << "\n"; 
        ave_score /= index.size();
        new_tag = false;
        if (sum_index.size() == index.size())
        {
            for (i = 0; i < sum_index.size(); i++) if (sum_index[i] != index[i]) { new_tag = true; break; }
            if (i == sum_index.size()) //qian100 + leave
            {
                sum_score.push_back(ave_score);
                score = 0; length = (sum_score.size() - 1) * 100 + leave;
                for (i = 0; i < sum_score.size(); i++)score += sum_score[i];
                score /= sum_score.size();

                if ((length >= thresh1) && (sum_index.size() >= thresh2) && (score >= thresh3))//长度，条数，分数
                {
                    os << "a score=" << score << "\n";
                    for (j = 0; j < sum_index.size(); j++)
                    {
                        i = sum_index[j];
                        //if ((length - pre_gap_num[i] - gap_num[i]) < thresh1) continue;
                        if (sign[i])zf = " + "; else zf = " - ";
                        os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << start[i] << " " << (length - pre_gap_num[i] - gap_num[i]) << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                        for (int k = ai - length; k < ai; k++) os << chars[(int)(sequences[i][k]) - 1];
                        os << "\n";
                    }
                    os << "\n";
                }
            }
        }
        if ((sum_index.size() != index.size()) || new_tag)//qian100,  ,leave
        {
            if (sum_score.size() > 0)
            {
                score = 0; length = sum_score.size() * 100;                    //qian100,
                for (i = 0; i < sum_score.size(); i++)score += sum_score[i];
                score /= sum_score.size();
                if ((length >= thresh1) && (sum_index.size() >= thresh2))//长度，条数，分数
                {
                    os << "a score=" << score << "\n";
                    for (j = 0; j < sum_index.size(); j++)
                    {
                        i = sum_index[j];
                        //if ((length - pre_gap_num[i]) < thresh1) continue;
                        if (sign[i])zf = " + "; else zf = " - ";
                        os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << start[i] << " " << length - pre_gap_num[i] << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                        for (int k = ai - length - leave; k < ai - leave; k++) os << chars[(int)(sequences[i][k]) - 1];
                        os << "\n";
                    }
                    os << "\n";
                }
            }

            length = leave;                                                //leave
            score = ave_score;
            if ((length >= thresh1) && (index.size() >= thresh2))//长度，条数，分数
            {
                os << "a score=" << score << "\n";
                for (j = 0; j < index.size(); j++)
                {
                    i = index[j];
                    //if ((length - gap_num[i]) < thresh1) continue;
                    if (sign[i])zf = " + "; else zf = " - ";
                    os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << pre_si[i] << " " << length - gap_num[i] << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                    for (int k = ai - length; k < ai; k++) os << chars[(int)(sequences[i][k]) - 1];
                    os << "\n";
                }
                os << "\n";
            }
        }
        //else
    }
    else if (sum_score.size() > 0) //qian 100
    {
        score = 0; length = sum_score.size() * 100;
        for (i = 0; i < sum_score.size(); i++)score += sum_score[i];
        score /= sum_score.size();

        if ((length >= thresh1) && (sum_index.size() >= thresh2))//长度，条数，分数
        {
            os << "a score=" << score << "\n";
            for (j = 0; j < sum_index.size(); j++)
            {
                i = sum_index[j];
                //if ((length - pre_gap_num[i]) < thresh1) continue;
                if (sign[i])zf = " + "; else zf = " - ";
                os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << start[i] << " " << length - pre_gap_num[i] << zf << std::setw(9) << std::right << len_sequences[i] << " ";
                for (int k = ai - length - leave; k < ai - leave; k++) os << chars[(int)(sequences[i][k]) - 1];
                os << "\n";
            }
            os << "\n";
        }
    }
}
*/

//向前延伸
utils::in_block* qiansu(std::vector<std::vector<unsigned char>>& sequences, utils::in_block* pre, utils::in_block* now, 
    bool* gap_tag, int* gap_numi, int thresh)
{
    int pre_end = 0, ai, stopi = 0, gap_num=0;
    float scorei100 = 0.0;
    float scorei = 0.0;
    float p, n;
    bool tag, do_tag = false;
    utils::in_block* tmp = now;
    if (pre != NULL) pre_end = pre->end;    //确定上界
    for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
    {
        gap_numi[now->name[i]] = 0;
        if ((int)(sequences[now->name[i]][now->start]) != 7) gap_tag[now->name[i]] = true;
        else gap_tag[now->name[i]] = false;
    }

    ai = now->start - 1;
    while (ai >= pre_end)
    {
        scorei = 0.0;
        //score100 += cs[(int)sequences[0][ai]];
        for (int i = 0; i < now->name.size(); i++)
        {
            if (((int)sequences[now->name[i]][ai] == 7) || ((int)sequences[0][ai] == 7))
            {
                if (gap_tag[now->name[i]]) { scorei -= d; gap_tag[now->name[i]] = false; }
                else  scorei -= e;
            }
            else { scorei += HOXD70[(int)sequences[0][ai]][(int)sequences[now->name[i]][ai]]; gap_tag[i] = true; }
        }
        p = 100.0 / (now->end - now->start + 1);
        n = 1.0 * (now->end - now->start) / (now->end - now->start + 1);
        scorei = scorei / now->name.size() * p + now->score * n;
        //std::cout << scorei << " " << scorei + now->score << " " << (now->score_100 + cs[(int)sequences[0][ai]]) * thresh / 100.0 << " **********\n";
        scorei100 = cs[(int)sequences[0][ai]] * p + now->score_100 * n;
        if (scorei >= (scorei100 * thresh / 100.0))
        {
            for (int i = 0; i < now->name.size(); i++)
                if ((int)sequences[now->name[i]][ai] != 7) { now->si[i]--; now->length[i]++; }
            now->start--;
            now->score = scorei;
            now->score_100 = scorei100;
        }
        else break;

        tag = true;
        for (int i = 0; i < now->name.size(); i++)
            if ((int)sequences[now->name[i]][ai] == 7)
            {
                tag = false;
                gap_num++;
                gap_numi[now->name[i]]++;
                if (gap_numi[now->name[i]] >= stop_g)
                    do_tag = true;
            }
        if (tag)
        {
            gap_num = 0;
            stopi = 0;
            for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
                gap_numi[now->name[i]] = 0;
        }//重置
        else if (stopi == 0)
        {
            tmp = new utils::in_block();
            tmp->score = now->score;
            tmp->score_100 = now->score_100;
            tmp->length = now->length;
            tmp->si = now->si;
            tmp->start = now->start;
            stopi++;
        }
        else stopi++;

        if (do_tag || ((stopi >= stop_g) && (gap_num >= stopi * now->name.size() / 2)))
        {
            now->score = tmp->score;
            now->score_100 = tmp->score_100;
            now->length = tmp->length;
            now->si = tmp->si;
            now->start = tmp->start;
            break;
        }
        ai--;
    }
    if ((pre != NULL) && (ai == pre_end - 1))
    {
        if (now->name.size() == pre->name.size())
        {
            tag = true;
            for (int i = 0; i < now->name.size(); i++)
                if (now->name[i] != pre->name[i])tag = false;
            if (tag) //可以连接
            {
                for (int i = 0; i < now->name.size(); i++)
                    pre->length[i] = now->length[i] + pre->length[i];
                p = 1.0 * (pre->end - pre->start) / (pre->end - pre->start + now->end - now->start);
                n = 1.0 * (now->end - now->start) / (pre->end - pre->start + now->end - now->start);
                pre->score = now->score * n + pre->score * p;
                pre->score_100 = now->score_100 * n + pre->score_100 * p;
                pre->end = now->end;
                pre->next = now->next;
                delete now;
                now = pre;
            }
        }
    }

    return now;
}
//向后延伸
void housu(std::vector<std::vector<unsigned char>>& sequences, utils::in_block* now, bool* gap_tag, int* gap_numi, int thresh)
{
    int end = sequences[0].size(), ai, gap_num = 0, stopi = 0;
    utils::in_block* tmp = now;
    int next_end = end;
    float scorei = 0.0, scorei100 = 0.0, n, p;
    bool tag, do_tag = false;
    if (now->next != NULL)  next_end = now->next->start;    //确定上界
    for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
    {
        gap_numi[now->name[i]] = 0;
        if ((int)(sequences[now->name[i]][now->end - 1]) != 7) gap_tag[now->name[i]] = true;
        else gap_tag[now->name[i]] = false;
    }

    ai = now->end;
    while (ai < end)
    {
        if (ai == next_end)
        {
            tmp = now->next;
            if (now->name.size() == tmp->name.size())
            {
                tag = true;
                for (int i = 0; i < now->name.size(); i++)
                    if (now->name[i] != tmp->name[i])tag = false;
                if (tag) //可以连接
                {
                    for (int i = 0; i < now->name.size(); i++)
                        now->length[i] = now->length[i] + tmp->length[i];
                    p = 1.0 * (tmp->end - tmp->start) / (tmp->end - tmp->start + now->end - now->start);
                    n = 1.0 * (now->end - now->start) / (tmp->end - tmp->start + now->end - now->start);
                    now->score = now->score * n + tmp->score * p;
                    now->score_100 = now->score_100 * n + tmp->score_100 * p;

                    now->end = tmp->end;
                    now->next = tmp->next;
                    delete tmp;

                    for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
                        if ((int)(sequences[now->name[i]][now->end - 1]) != 7) gap_tag[now->name[i]] = true;
                        else gap_tag[now->name[i]] = false;
                    ai = now->end;

                    if (now->next != NULL) next_end = now->next->start;
                    else next_end = end;
                    gap_num = stopi = 0;
                    do_tag = false;
                    for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
                        gap_numi[now->name[i]] = 0;
                    continue;
                }
                else break;//不能连接
            }
            else break;//不能连接
        }
        scorei = 0.0;
        for (int i = 0; i < now->name.size(); i++)
        {
            if (((int)sequences[now->name[i]][ai] == 7) || ((int)sequences[0][ai] == 7))
            {
                if (gap_tag[now->name[i]]) { scorei -= d; gap_tag[now->name[i]] = false; }
                else  scorei -= e;
            }
            else { scorei += HOXD70[(int)sequences[0][ai]][(int)sequences[now->name[i]][ai]]; gap_tag[i] = true; }
        }
        p = 100.0 / (now->end - now->start + 1);
        n = 1.0 * (now->end - now->start) / (now->end - now->start + 1);
        scorei = scorei / now->name.size() * p + now->score * n;
        //std::cout << scorei << " " << scorei + now->score << " " << (now->score_100 + cs[(int)sequences[0][ai]]) * thresh / 100.0 << " **********\n";
        scorei100 = cs[(int)sequences[0][ai]] * p + now->score_100 * n;
        if (scorei >= (scorei100 * thresh / 100.0))
        {
            for (int i = 0; i < now->name.size(); i++)
                if ((int)sequences[now->name[i]][ai] != 7)  now->length[i]++;
            now->end++;
            now->score = scorei;
            now->score_100 = scorei100;
        }
        else break;

        tag = true;
        for (int i = 0; i < now->name.size(); i++)
            if ((int)sequences[now->name[i]][ai] == 7)
            {
                gap_num++;
                gap_numi[now->name[i]]++;
                if (gap_numi[now->name[i]] >= stop_g)
                    do_tag = true;
                tag = false;
            }
        if (tag)
        {
            gap_num = 0;
            stopi = 0;
            for (int i = 0; i < now->name.size(); i++)  //初始化gap_tag
                gap_numi[now->name[i]] = 0;
        }//重置
        else if (stopi == 0)
        {
            tmp = new utils::in_block();
            tmp->score = now->score;
            tmp->score_100 = now->score_100;
            tmp->length = now->length;
            tmp->end = now->end;
            stopi++;
        }
        else stopi++;

        if (do_tag || ((stopi >= stop_g) && (gap_num >= stopi * now->name.size() / 2)))
        {
            now->score = tmp->score;
            now->score_100 = tmp->score_100;
            now->length = tmp->length;
            now->end = tmp->end;
            break;
        }
        ai++;
    }
}


//前后延伸的maf100,全部输出，供后续筛选。
void utils::new_insert_and_write_100_maf(std::string out_file_name, std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
    std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name,
    std::vector<bool>& sign, std::vector<utils::more_block>&more_gap, int thresh1, int thresh2, int thresh3)//长度，条数，分数
{
    int thresh_extend = 98;//thresh3;
    const size_t sequence_number = insertions.size();
    std::cout << "insertions " << sequence_number << " " << insertions[0].size() << " " << insertions[1].size() << "\n";
    int* si = new int[sequence_number];
    int* pre_si = new int[sequence_number];
    int* start = new int[sequence_number];
    bool* gap_tag = new bool[sequence_number];
    bool new_tag = true;
    int* len_sequences;
    int* gap_num = new int[sequence_number];//当前block gap数量
    int* pre_gap_num = new int[sequence_number];//前面几个block gap数量
    for (int i = 0; i < sequence_number; i++) { si[i] = gap_num[i] = pre_gap_num[i] = 0; gap_tag[i] = true; }
    int  length = 0, ai = 0;
    float score = 0;
    int* score_two = NULL;
    float ave_score = 0, ave_score_100 = 0;
    int name_len = 0, pre = 0, diff;
    char chars[7] = { 'A','C','G','T','N','N','-' };
    for (int i = 0; i < name.size(); i++)
    {
        //replace(name[i].begin(), name[i].end(), ' ', '-');
        if (name_len < name[i].size())
            name_len = name[i].size();
    }
        
    const auto align_start1 = std::chrono::high_resolution_clock::now(); //插入前
if (thresh2<=2)// 逆补二次查找maf
    insertion_gap_more(os, sequences, N_insertions, name, sign, more_gap, thresh1);//thresh1

    //插入
    len_sequences = vector_insertion_gap_N(sequences, insertions, N_insertions);
    /*
    #if defined(_WIN32)//内存小
        len_sequences = vector_insertion_gap_N(sequences, insertions, N_insertions);
    #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
        len_sequences = insertion_gap_N(sequences, insertions, N_insertions);
    #endif
    */

    const auto align_start2 = std::chrono::high_resolution_clock::now(); //插入后
    std::cout << align_start2 - align_start1 << "  insert\n";

    //求分数
     //求分数
    std::ofstream tmpo(arguments::score_file, std::ios::binary | std::ios::out);
    tmpo << "I" << "\t" << "Match_num" << "\t" << "Len_sequence" << "\t" << "M/L" << "\t" << "Name" << "\n";
    tmpo << 0 << "\t" << "--" << "\t" << len_sequences[0] << "\t" << "---" << "\t" << name[0] << "\n";
    for (int i = 0; i < sequence_number; i++)
    {
        if (i != 0)
        {
            score_two = Compare_two(sequences[0], sequences[i], 0, 0, 1, 0);
            tmpo << i << "\t" << score_two[0] << "\t" << len_sequences[i] << "\t" << score_two[0] * 1.0 / len_sequences[i] << "\t" << name[i] << "\n";
        }
    }
    tmpo.close();
    const auto align_start3 = std::chrono::high_resolution_clock::now(); //求分后
    std::cout << align_start3 - align_start2 << "  score\n";
    /*int* score_two = Compare_two(sequences[0], sequences[1]);
    const auto align_start3 = std::chrono::high_resolution_clock::now(); //求分后
    std::cout << align_start3 - align_start2 << "  score\n";
    std::cout << "score: " << score_two[0] << " / " << score_two[1] << " = " << score_two[0] * 1.0 / score_two[1] << "\n";
    */

    std::string zf;
    int i = 0, j, leave = sequences[0].size() % 100, hundred = sequences[0].size() / 100, hi = 0;
    float score100;
    std::cout << "all_length = " << sequences[0].size() << "\n";
    ai = 0;
    bool tag = false;

    std::vector<float> sum_score;
    std::vector<float> sum_score_100;
    std::vector<int> sum_index;
    std::vector<int> index;
    std::vector<int> nscore(sequence_number);
    bool* out_gap;
    in_block* new_block;
    in_block* all_block = new in_block(); //所有的block 先用结构体存储,头节点。
    all_block->next = NULL;
    in_block* now = all_block;

    for (hi = 0; hi < hundred; hi++)  //处理整百
    {
        score100 = 0;//后面可以 thresh%*score100
        for (i = 0; i < sequence_number; i++) {
            nscore[i] = 0;
            pre_si[i] = si[i];
            pre_gap_num[i] += gap_num[i];
            gap_num[i] = 0;
        }//每100个分数初始化0

        for (ai = hi * 100; ai < (hi + 1) * 100; ai++)
        {
            for (i = 0; i < sequence_number; i++) if ((int)sequences[i][ai] != 7) si[i]++;
            score100 += cs[(int)sequences[0][ai]];
            for (i = 0; i < sequence_number; i++)
            {
                if (((int)sequences[i][ai] == 7) || ((int)sequences[0][ai] == 7))
                {
                    if ((int)sequences[i][ai] == 7) gap_num[i]++;
                    if (gap_tag[i]) { nscore[i] -= d; gap_tag[i] = false; }
                    else  nscore[i] -= e;
                }
                else { nscore[i] += HOXD70[(int)sequences[0][ai]][(int)sequences[i][ai]]; gap_tag[i] = true; }
            }
        }

        std::vector<int>().swap(index);
        index = index_100(nscore, score100 * thresh3 / 100);
        ave_score_100 = score100;
        if (index.size() != 0)
        {
            ave_score = 0;
            for (i = 0; i < index.size(); i++) ave_score += nscore[index[i]];
            ave_score /= index.size();
            if (sum_score.size() == 0)
            {
                for (i = 0; i < sequence_number; i++) { start[i] = pre_si[i]; pre_gap_num[i] = 0; }
                sum_index.assign(index.begin(), index.end());
                sum_score.push_back(ave_score);
                sum_score_100.push_back(score100);

                continue;
            }
            else if (sum_index.size() == index.size())
            {
                for (i = 0; i < sum_index.size(); i++) if (sum_index[i] != index[i]) break;
                if (i == sum_index.size()) { sum_score.push_back(ave_score); sum_score_100.push_back(score100); continue; }
            }
            //else
        }
        //else

        if (sum_score.size() == 0)//当前为废100，
        {
            for (int i = 0; i < sequence_number; i++) pre_gap_num[i] = gap_num[i] = 0;//把当前100gap数清掉;
            continue;
        }
        else
        {
            score = 0; score100 = 0; length = sum_score.size() * 100;
            for (i = 0; i < sum_score.size(); i++) { score += sum_score[i]; score100 += sum_score_100[i]; }
            score /= sum_score.size();
            score100 /= sum_score_100.size();

            new_block = new in_block();
            new_block->score = score;
            new_block->score_100 = score100;
            for (j = 0; j < sum_index.size(); j++)
            {
                i = sum_index[j];
                new_block->name.push_back(i);
                new_block->si.push_back(start[i]);
                new_block->length.push_back(length - pre_gap_num[i]);
                new_block->start = ai - length - 100;
                new_block->end = ai - 100;
                new_block->next = NULL;
            }
            now->next = new_block;
            now = new_block;

            //写完了
            for (int i = 0; i < sequence_number; i++) pre_gap_num[i] = 0;//写完一组，清空pre gap数
            std::vector<float>().swap(sum_score);
            std::vector<float>().swap(sum_score_100);
            std::vector<int>().swap(sum_index);
            if (index.size() != 0)//从上面筛选下来，连不上，但是分数合格的首100
            {
                for (i = 0; i < sequence_number; i++) start[i] = pre_si[i];
                sum_index.assign(index.begin(), index.end());
                sum_score.push_back(ave_score);
                sum_score_100.push_back(ave_score_100);

                continue;
            }
        }
    }
    //处理后续的leave
    score100 = 0;//-d * leave / 100;
    for (i = 0; i < sequence_number; i++)
    {
        nscore[i] = 0;
        pre_si[i] = si[i];
        pre_gap_num[i] += gap_num[i];
        gap_num[i] = 0;
    }//预处理
    for (; ai < sequences[0].size(); ai++)
    {
        for (i = 0; i < sequence_number; i++) if ((int)sequences[i][ai] != 7) si[i]++;
        score100 += cs[(int)sequences[0][ai]];
        for (i = 0; i < sequence_number; i++)
        {
            if (((int)sequences[i][ai] == 7) || ((int)sequences[0][ai] == 7))
            {
                if ((int)sequences[i][ai] == 7) gap_num[i]++;
                if (gap_tag[i]) { nscore[i] -= d; gap_tag[i] = false; }
                else  nscore[i] -= e;
            }
            else { nscore[i] += HOXD70[(int)sequences[0][ai]][(int)sequences[i][ai]]; gap_tag[i] = true; }
        }
    }
    std::vector<int>().swap(index);
    index = index_100(nscore, score100 * thresh3 / 100);
    ave_score_100 = score100;
    if (index.size() != 0)
    {
        ave_score = 0;

        for (i = 0; i < index.size(); i++) ave_score += nscore[index[i]]; //std::cout << "???" << index[i] << "\n"; 
        ave_score /= index.size();
        new_tag = false;

        if (sum_index.size() == index.size())
        {
            for (i = 0; i < sum_index.size(); i++) if (sum_index[i] != index[i]) { new_tag = true; break; }

            if (i == sum_index.size()) //qian100 + leave
            {
                sum_score.push_back(ave_score);
                sum_score_100.push_back(score100);
                score = 0; score100 = 0; length = (sum_score.size() - 1) * 100 + leave;
                for (int i = 0; i < sum_score.size(); i++) { score += sum_score[i]; score100 += sum_score_100[i]; }

                score /= sum_score.size();
                score100 /= sum_score_100.size();

                new_block = new in_block();
                new_block->score = score;
                new_block->score_100 = score100;
                for (j = 0; j < sum_index.size(); j++)
                {
                    i = sum_index[j];
                    new_block->name.push_back(i);
                    new_block->si.push_back(start[i]);
                    new_block->length.push_back(length - pre_gap_num[i] - gap_num[i]);
                    new_block->start = ai - length;
                    new_block->end = ai;
                    new_block->next = NULL;
                }
                now->next = new_block;
                now = new_block;
            }
        }

        if ((sum_index.size() != index.size()) || new_tag)//qian100,  ,leave
        {
            if (sum_score.size() > 0)
            {
                score = 0; score100 = 0; length = sum_score.size() * 100;                    //qian100,
                for (i = 0; i < sum_score.size(); i++) { score += sum_score[i]; score100 += sum_score_100[i]; }
                score /= sum_score.size();
                score100 /= sum_score_100.size();

                new_block = new in_block();
                new_block->score = score;
                new_block->score_100 = score100;
                for (j = 0; j < sum_index.size(); j++)
                {
                    i = sum_index[j];
                    new_block->name.push_back(i);
                    new_block->si.push_back(start[i]);
                    new_block->length.push_back(length - pre_gap_num[i]);
                    new_block->start = ai - length - leave;
                    new_block->end = ai - leave;
                    new_block->next = NULL;
                }
                now->next = new_block;
                now = new_block;


            }
            if (leave > 0)
            {
                length = leave;                                                //leave
                score = ave_score;

                new_block = new in_block();
                new_block->score = score;
                new_block->score_100 = ave_score_100;
                for (j = 0; j < index.size(); j++)
                {
                    i = index[j];
                    new_block->name.push_back(i);
                    new_block->si.push_back(pre_si[i]);
                    new_block->length.push_back(length - gap_num[i]);
                    new_block->start = ai - length;
                    new_block->end = ai;
                    new_block->next = NULL;
                }
                now->next = new_block;
                now = new_block;
            }
        }
        //else
    }
    else if (sum_score.size() > 0) //qian 100
    {
        score = 0; score100 = 0; length = sum_score.size() * 100;
        for (i = 0; i < sum_score.size(); i++) { score += sum_score[i]; score100 += sum_score_100[i]; }
        score /= sum_score.size();
        score100 /= sum_score_100.size();

        new_block = new in_block();
        new_block->score = score;
        new_block->score_100 = score100;
        for (j = 0; j < sum_index.size(); j++)
        {
            i = sum_index[j];
            new_block->name.push_back(i);
            new_block->si.push_back(start[i]);
            new_block->length.push_back(length - pre_gap_num[i]);
            new_block->start = ai - length - leave;
            new_block->end = ai - leave;
            new_block->next = NULL;
        }
        now->next = new_block;

    }
    
    /*
    now = all_block->next;
    while (now != NULL)
    {

        std::cout << "score: " << now->score << "\nscore100: " << now->score_100 << "\nstart: " << now->start << "\nend: " << now->end << "\n";
        std::cout << now->name.size() << "\n";
        for (i = 0; i < now->name.size(); i++)
            std::cout << "name: " << now->name[i] << "\tsi: " << now->si[i] << "\tlength: " << now->length[i] << "\n";
        std::cout << "\n";

        now = now->next;
    }*/

    if (all_block->next == NULL);
    else
    {
        //首个
        now = all_block->next;
        qiansu(sequences, NULL, now, gap_tag, gap_num, thresh_extend);
        housu(sequences, now, gap_tag, gap_num, thresh_extend);

        //非首个
        new_block = now;
        now = now->next;
        while (now != NULL)
        {
            now = qiansu(sequences, new_block, now, gap_tag, gap_num, thresh_extend);
            housu(sequences, now, gap_tag, gap_num, thresh_extend);
            new_block = now;
            now = now->next;
        }
        /*std::cout << "**********************************\n";

        now = all_block->next;
        while (now != NULL)
        {

            std::cout << "score: " << now->score << "\nscore100: " << now->score_100 << "\nstart: " << now->start << "\nend: " << now->end << "\n";
            std::cout << now->name.size() << "\n";
            for (i = 0; i < now->name.size(); i++)
                std::cout << "name: " << now->name[i] << "\tsi: " << now->si[i] << "\tlength: " << now->length[i] << "\n";
            std::cout << "\n";
            now = now->next;
        }*/

    }

    now = all_block->next;
    while (now != NULL)
    {
        //写出
        os << "a score=" << now->score << "\n";
        for (j = 0; j < now->name.size(); j++)
        {
            i = now->name[j];
            if (sign[i])zf = " + "; else zf = " - ";
            os << "s " << std::setw(name_len + 1) << std::left << name[i] << std::setw(9) << std::right << now->si[j] << " " << std::setw(9) << std::right << now->length[j] << zf << std::setw(9) << std::right << len_sequences[i] << " ";
            for (int k = now->start; k < now->end; k++) os << chars[(int)(sequences[i][k]) - 1];
            os << "\n";
        }
        os << "\n";
        //写完
        now = now->next;
    }

    //输出 全部未筛选maf 到tmp文件
    std::string tmp_file_name = out_file_name.substr(0, out_file_name.rfind(".")) + ".tmp";
    std::ofstream tmp_os(tmp_file_name, std::ios::binary | std::ios::out);

    tmp_os << sequence_number << "\n";                                 //输出序列条数
    for (i = 0; i < sequence_number; i++)                              //输出序列名
        tmp_os << name[i] << "\n";
    for (i = 0; i < sequence_number; i++)tmp_os << len_sequences[i] << " ";//输出序列长度
    tmp_os << "\n";
    for (i = 0; i < sequence_number; i++) if (sign[i]) tmp_os << 1 << " "; else  tmp_os << 0 << " ";//输出序列符号
    tmp_os << "\n\n";

    now = all_block->next;
    while (now != NULL)
    {
        //写出
        tmp_os << now->score << " " << now->name.size() << " " << (now->end - now->start) << "\n";  //score,num,all_len
        out_gap = new bool[now->end - now->start];
        memset(out_gap, 0, (now->end - now->start));
        diff = 0;
        for (j = 0; j < now->name.size(); j++)
        {
            i = now->name[j];
            tmp_os << i << " " << now->si[j] << " " << now->length[j] << " "; //si,len
            for (int k = now->start; k < now->end; k++)
            {
                if ((int)(sequences[i][k]) == 7) out_gap[k - now->start] = true;
                tmp_os << chars[(int)(sequences[i][k]) - 1];//seq
            }
            tmp_os << "\n";
        }
        for (int k = 0; k < (now->end - now->start); k++) if (out_gap[k])diff++;
        tmp_os << diff << "\n";
        delete[] out_gap;
        //写完
        now = now->next;
    }
    tmp_os.close();

    while (all_block)
    {
        now = all_block;
        all_block = all_block->next;
        delete now;
    }
    delete[] si;
    delete[] pre_si;
    delete[] start;
    delete[] gap_tag;
    delete[] len_sequences;
    delete[] gap_num;
    delete[] pre_gap_num;
    delete[] score_two;
    std::vector<float>().swap(sum_score);
    std::vector<float>().swap(sum_score_100);
    std::vector<int>().swap(sum_index);
    std::vector<int>().swap(index);
    std::vector<int>().swap(nscore);
    const auto align_start4 = std::chrono::high_resolution_clock::now(); //插入后
    std::cout << align_start4 - align_start3 << "  choose-output\n";

}

//读取tmp文件筛选
void utils::read_tmp_filter_out_maf(std::istream& is, std::ostream& os, int thresh1, int thresh2, int thresh3)//1-len,2-num,3-Continuity.
{
    int sequence_number, j; //序列数量
    unsigned char dir;
    int name_len = 0, num, ailen;
    float score;

    std::string zf;
    MAF_block deal_block;
    block iblock;
    //std::string str;
    is >> sequence_number;
    getline(is, zf);
    //std::cout << sequence_number << "\n";
    std::vector<std::string> name(sequence_number); //名字
    std::vector<int> sign(sequence_number);        //符号
    std::vector<int> len_sequences(sequence_number);//长度
    for (int i = 0; i < sequence_number; i++)//名字
        getline(is, name[i]);
    for (int i = 0; i < sequence_number; i++)//长度 
        is >> len_sequences[i];
    for (int i = 0; i < sequence_number; i++)//符号
        is >> sign[i];

    for (int i = 0; i < name.size(); i++)
        if (name_len < name[i].size())
            name_len = name[i].size();
    while (is >> deal_block.score >> num >> ailen)
    {
        //is >> deal_block.score >> num >> ailen;
        //std::cout << num<<" ??\n";
        std::vector<block>().swap(deal_block.seq);
        iblock.seqi.resize(ailen);
        for (int i = 0; i < num; i++)
        {
            is >> iblock.name >> iblock.start >> iblock.length;
            for (j = 0; j < ailen; j++)
            {
                is >> iblock.seqi[j];
                //std::cout << iblock.seqi[j];
            }
            //std::cout <<"\n" << iblock.seqi.size() << "\n";
            deal_block.seq.push_back(iblock);
        }
        is >> deal_block.tag_num;
        //std::cout << deal_block.tag_num << "--\n";
        //filter
        //std::cout << ailen << " " << num << " " << (100 - 100.0 * deal_block.tag_num / ailen) << "\n";
        if ((ailen >= thresh1) && (num >= thresh2) && ((100 - 100.0 * deal_block.tag_num / ailen) > thresh3))
        {
            //写出
            os << "a score=" << deal_block.score << "\n";
            for (int i = 0; i < num; i++)
            {
                j = deal_block.seq[i].name;
                if (sign[j])zf = " + "; else zf = " - ";
                os << "s " << std::setw(name_len + 1) << std::left << name[j] << std::setw(9) << std::right << deal_block.seq[i].start << " " << std::setw(9) << std::right << deal_block.seq[i].length << zf << std::setw(9) << std::right << len_sequences[j] << " ";
                for (int k = 0; k < ailen; k++) os << deal_block.seq[i].seqi[k];
                os << "\n";
            }
            os << "\n";
            //写完
        }
    }

}


//处理逆补串maf -------|
int sc[8][8] = { {},
    {0,1,-2,-2,-2,0,0,-2},
    {0,-2,1,-2,-2,0,0,-2},
    {0,-2,-2,1,-2,0,0,-2},
    {0,-2,-2,-2,1,0,0,-2},
    {0,0,0,0,0,0,0,-2},
    {0,0,0,0,0,0,0,-2},
    {0,-2,-2,-2,-2,-2,0,-2} };//评分矩阵
//获得最大得分子序列
bool get_max_begin_end(std::vector<unsigned char>& seq1, std::vector<unsigned char>& seq2, int& begin1, int& begin2, int& clen1, int& clen2, int& ss, int& ee)

{
    if (seq1.size() != seq1.size())
    {
        std::cout << "error:The length of the different\n";
        exit(-1);
    }
    int n = seq1.size();
    int f = -1, s, mi = -10;
    for (int i = 0; i < n; i++)//寻找的一个正数
    {
        if (sc[seq1[i]][seq2[i]] > 0)
        {
            f = i;
            break;
        }
        else if (sc[seq1[i]][seq2[i]] > mi)
        {
            mi = sc[seq1[i]][seq2[i]];
            s = i;
        }
    }
    if (f == -1)//如果没有正数，那肯定最大的就是0了
    {
        return false;
    }
    int add = sc[seq1[f]][seq2[f]], sum = sc[seq1[f]][seq2[f]], l = f, r = f, sl = f, sr = f;//把第一个看成当前最大子序列，a也从第一个开始
    for (f++; f < n; f++)
    {
        add += sc[seq1[f]][seq2[f]];//添加元素
        if (add > sum)//更新最大子序列，连同始末位置一起更新
        {
            sum = add;
            l = sl;
            r = f;
        }
        if (add < 0)//成为了累赘，舍弃变为零
        {
            add = 0;
            sl = f + 1;//把始位置也更新掉
        }
    }
    sl = 0;
    for (f = 0; f < l; f++)
    {
        if (seq1[f] != '\7')
            sl++;
    }
    sr = 0;
    for (f = seq2.size()-1; f > r; f--)
    {
        if (seq2[f] != '\7')
            sr++;
    }

    begin1 += sl;
    begin2 += sr;
    clen1 = 0; clen2 = 0;
    for (f = l; f <= r; f++)
    {
        if (seq1[f] != '\7')
            clen1++;
        if (seq2[f] != '\7')
            clen2++; 
    }
    ss = l, ee = r;
    return true;
}
//插入并输出
void utils::insertion_gap_out(std::ostream& os, int seqi, int* final_sequences, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::string>& name,
    std::vector<bool>& sign, utils::m_block* more_block, int nsum1, int nsum2, std::vector<Insertion>& N_insertion1, std::vector<Insertion>& N_insertion2, int thresh1)

{
    std::vector<std::vector<Insertion2>> all_insertions;
    std::vector<Insertion2> i_all_insertions;
    std::vector<unsigned char> tmp_vector;
    int i = 0, j = 0, pre = 0, mi, g_num, k;
    int* ii = new int[2];
    std::vector<std::vector<Insertion2>> more_insertions(2);
    int a_start = more_block->start1, b_start = more_block->start2;
    int a_len = more_block->end1 - more_block->start1;
    int b_len = more_block->end2 - more_block->start2;
    
    std::vector<unsigned char> A(a_len);
    std::vector<unsigned char> B(b_len);
    i = 0; j = more_block->start1;
    while (j < more_block->end1) A[i++] = sequences[0][j++];
    i = 0; j = more_block->end2 - 1;
    while (j > more_block->start2) B[i++] = 5 - sequences[seqi][j--];
    B[i] = 5 - sequences[seqi][j];

    for (i = 0; i < N_insertion1.size(); i++)
        N_insertion1[i].index = N_insertion1[i].index - a_start;
    for (i = 0; i < N_insertion2.size(); i++)
        N_insertion2[i].index = more_block->end2 - N_insertion2[i].index;
    std::reverse(N_insertion2.begin(), N_insertion2.end()); //逆置
    //00000000000000
    i = 0; j = 0; g_num = 0, ii[0] = 0; ii[1] = 0; k = 0;
    std::vector<Insertion2>().swap(i_all_insertions);//清空 i_all_insertions
    while (i < more_block->gap1.size() && j < N_insertion1.size())
    {
        if (std::get<0>(more_block->gap1[i]) == N_insertion1[j].index)//相等变number
        {
            g_num += std::get<1>(more_block->gap1[i]);
            if (N_insertion1[j].number > std::get<1>(more_block->gap1[i]))
            {
                i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap1[i]), (size_t)std::get<1>(more_block->gap1[i]), 0 }));  //gap 个 N
                insert_others(Insertion2({ (size_t)std::get<0>(more_block->gap1[i]) + g_num ,N_insertion1[j].number - std::get<1>(more_block->gap1[i]), 0 }), utils::Insertion2({ (size_t)std::get<0>(more_block->gap1[i]) + g_num ,0, N_insertion1[j].number - std::get<1>(more_block->gap1[i]) }), more_insertions, k, 2, ii);
            }
            else
            {
                i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap1[i]) ,N_insertion1[j].number, (size_t)std::get<1>(more_block->gap1[i]) - N_insertion1[j].number }));
            }
            i++;
            j++;
        }
        else if (std::get<0>(more_block->gap1[i]) > N_insertion1[j].index) //不等插新
        {
            insert_others(Insertion2({ N_insertion1[j].index + g_num ,N_insertion1[j].number, 0 }), utils::Insertion2({ N_insertion1[j].index + g_num ,0, N_insertion1[j].number }), more_insertions, k, 2, ii);
            j++;
        }
        else
        {
            g_num += std::get<1>(more_block->gap1[i]);
            i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap1[i]), 0, (size_t)std::get<1>(more_block->gap1[i]) }));
            i++;
        }
    }
    while (j < N_insertion1.size()) //后续多余的N，继续插入
    {
        insert_others(Insertion2({ N_insertion1[j].index + g_num ,N_insertion1[j].number, 0 }), utils::Insertion2({ N_insertion1[j].index + g_num ,0, N_insertion1[j].number }), more_insertions, k, 2, ii);
        j++;
    }
    while (i < more_block->gap1.size()) //后续多余的gap，继续插入
    {
        g_num += std::get<1>(more_block->gap1[i]);
        i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap1[i]), 0, (size_t)std::get<1>(more_block->gap1[i]) }));
        i++;
    }
    all_insertions.push_back(i_all_insertions);
    //11111111111111
    i = 0; j = 0; g_num = 0, ii[0] = 0; ii[1] = 0; k = 1;
    std::vector<Insertion2>().swap(i_all_insertions);//清空 i_all_insertions
    while (i < more_block->gap2.size() && j < N_insertion2.size())
    {
        if (std::get<0>(more_block->gap2[i]) == N_insertion2[j].index)//相等变number
        {
            g_num += std::get<1>(more_block->gap2[i]);
            if (N_insertion2[j].number > std::get<1>(more_block->gap2[i]))
            {
                i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap2[i]), (size_t)std::get<1>(more_block->gap2[i]), 0 }));  //gap 个 N
                insert_others(Insertion2({ (size_t)std::get<0>(more_block->gap2[i]) + g_num ,N_insertion2[j].number - std::get<1>(more_block->gap2[i]), 0 }), utils::Insertion2({ (size_t)std::get<0>(more_block->gap2[i]) + g_num ,0, N_insertion2[j].number - std::get<1>(more_block->gap2[i]) }), more_insertions, k, 2, ii);
            }
            else
            {
                i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap2[i]) ,N_insertion2[j].number, (size_t)std::get<1>(more_block->gap2[i]) - N_insertion2[j].number }));
            }
            i++;
            j++;
        }
        else if (std::get<0>(more_block->gap2[i]) > N_insertion2[j].index) //不等插新
        {
            insert_others(Insertion2({ N_insertion2[j].index + g_num ,N_insertion2[j].number, 0 }), utils::Insertion2({ N_insertion2[j].index + g_num ,0, N_insertion2[j].number }), more_insertions, k, 2, ii);
            j++;
        }
        else
        {
            g_num += std::get<1>(more_block->gap2[i]);
            i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap2[i]), 0, (size_t)std::get<1>(more_block->gap2[i]) }));
            i++;
        }
    }
    while (j < N_insertion2.size()) //后续多余的N，继续插入
    {
        insert_others(Insertion2({ N_insertion2[j].index + g_num ,N_insertion2[j].number, 0 }), utils::Insertion2({ N_insertion2[j].index + g_num ,0, N_insertion2[j].number }), more_insertions, k, 2, ii);
        j++;
    }
    while (i < more_block->gap2.size()) //后续多余的gap，继续插入
    {
        g_num += std::get<1>(more_block->gap2[i]);
        i_all_insertions.push_back(Insertion2({ (size_t)std::get<0>(more_block->gap2[i]), 0, (size_t)std::get<1>(more_block->gap2[i]) }));
        i++;
    }
    all_insertions.push_back(i_all_insertions);

    //消除  一列多个N导致de多个插空
    std::vector<Insertion> multi;
    while (i < more_insertions[0].size() && j < more_insertions[1].size())
    {
        if (more_insertions[0][i].index == more_insertions[1][j].index)
        {
            if (more_insertions[0][i].gap_num == 0 || more_insertions[1][j].gap_num == 0);
            else if (more_insertions[0][i].gap_num > more_insertions[1][j].gap_num)
                multi.push_back(Insertion({ more_insertions[0][i].index ,more_insertions[1][j].gap_num }));
            else multi.push_back(Insertion({ more_insertions[0][i].index ,more_insertions[0][i].gap_num }));
            i++; j++;
        }
        else if (more_insertions[0][i].index > more_insertions[1][j].index) j++;
        else i++;
    }

    //插入到原串
    int more = 0, all_size = 0, ti = 0;
    k = 0;
    //00000000000000
    i = 0;
    more = A.size();
    for (j = 0; j < all_insertions[i].size(); j++)
        more += (all_insertions[i][j].n_num + all_insertions[i][j].gap_num);
    all_size = more;
    mi = 0;
    for (j = 0; j < more_insertions[i].size(); j++)
    {
        all_size += more_insertions[i][j].n_num;
        if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
            all_size -= multi[mi++].number;
        all_size += more_insertions[i][j].gap_num;
    }
    tmp_vector.resize(more);

    ti = 0;
    k = 0;
    for (j = 0; j < all_insertions[i].size(); j++)
    {
        while (k < all_insertions[i][j].index)
            tmp_vector[ti++] = A[k++];
        for (int p = 0; p < all_insertions[i][j].n_num; p++)
            tmp_vector[ti++] = '\5';
        for (int p = 0; p < all_insertions[i][j].gap_num; p++)
            tmp_vector[ti++] = '\7';
    }while (k < A.size())tmp_vector[ti++] = A[k++];

    A.resize(all_size);
    mi = 0;
    k = 0;
    ti = 0;
    for (j = 0; j < more_insertions[i].size(); j++)
    {
        while (ti < more_insertions[i][j].index)
            A[k++] = tmp_vector[ti++];
        for (int p = 0; p < more_insertions[i][j].n_num; p++)
            A[k++] = '\5';
        if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
        {
            for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                A[k++] = '\7';
            mi++;
        }
        else
            for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                A[k++] = '\7';
    }while (ti < more) A[k++] = tmp_vector[ti++];
    //111111111111111111
    i = 1;
    ti = 0;
    k = 0;
    for (j = 0; j < all_insertions[i].size(); j++)
    {
        while (k < all_insertions[i][j].index)
            tmp_vector[ti++] = B[k++];
        for (int p = 0; p < all_insertions[i][j].n_num; p++)
            tmp_vector[ti++] = '\5';
        for (int p = 0; p < all_insertions[i][j].gap_num; p++)
            tmp_vector[ti++] = '\7';
    }while (k < B.size())tmp_vector[ti++] = B[k++];

    B.resize(all_size);
    mi = 0;
    k = 0;
    ti = 0;
    for (j = 0; j < more_insertions[i].size(); j++)
    {
        while (ti < more_insertions[i][j].index)
            B[k++] = tmp_vector[ti++];
        for (int p = 0; p < more_insertions[i][j].n_num; p++)
            B[k++] = '\5';
        if (mi < multi.size() && more_insertions[i][j].index == multi[mi].index)
        {
            for (int p = 0; p < (more_insertions[i][j].gap_num - multi[mi].number); p++)
                B[k++] = '\7';
            mi++;
        }
        else
            for (int p = 0; p < more_insertions[i][j].gap_num; p++)
                B[k++] = '\7';
    }while (ti < more) B[k++] = tmp_vector[ti++];

    std::vector<unsigned char>().swap(tmp_vector);
    std::vector<std::vector<Insertion2>>().swap(all_insertions);
    std::vector<Insertion2>().swap(i_all_insertions);
    std::vector<std::vector<Insertion2>>().swap(more_insertions);
    std::vector<Insertion>().swap(multi);
    delete[] ii;

    a_start += nsum1;
    b_start += nsum2;

    if (get_max_begin_end(A, B, a_start, b_start, a_len, b_len, i, j) && (j-i)>=thresh1)
    {
        int name_len = name[0].size() > name[seqi].size() ? name[0].size() : name[seqi].size();
        char chars[7] = { 'A','C','G','T','N','N','-' };
        os << "a score= 0 "<< "\n";
        os << "s " << std::setw(name_len + 1) << std::left << name[0] << std::setw(9) << std::right << a_start << " " << std::setw(9) << std::right << a_len << " + " << std::setw(9) << std::right << final_sequences[0] << " ";
        for (k = i; k <= j; k++) os << chars[(int)A[k] - 1];
        os <<"\n";
        if (sign[seqi])
            os << "s " << std::setw(name_len + 1) << std::left << name[seqi] << std::setw(9) << std::right << b_start << " " << std::setw(9) << std::right << b_len << " - " << std::setw(9) << std::right << final_sequences[seqi] << " ";
        else
            os << "s " << std::setw(name_len + 1) << std::left << name[seqi] << std::setw(9) << std::right << b_start << " " << std::setw(9) << std::right << b_len << " + " << std::setw(9) << std::right << final_sequences[seqi] << " ";

        for (k = i; k <= j; k++) os << chars[(int)B[k] - 1];
        
        os << "\n\n";
    }

}
//处理逆补串maf
void utils::insertion_gap_more(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences,
    const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name,
    std::vector<bool>& sign, std::vector<utils::more_block>& more_gap, int thresh1) //写回比对结果
{
    int i, j, k, jm, jn, j0, nsum1, nsum2;
    const size_t sequence_number = sequences.size();
    int* len_sequences = new int[sequence_number];
    int* final_sequences = new int[sequence_number];
    std::vector<Insertion> gapn1, gapn2;
    for (k = 0; k < sequence_number; k++)  final_sequences[k] = len_sequences[k] = sequences[k].size();
    for (i = 0; i < sequence_number; i++)
        for (j = 0; j < N_insertions[i].size(); j++)
            final_sequences[i] += N_insertions[i][j].number;

    for (i = 1; i < sequence_number; i++)
    {
        jm = jn = j0 = nsum1 = nsum2 = 0;
        gapn1.clear();
        gapn2.clear();
        while (jm < more_gap[i].size())
        {
            while (j0 < N_insertions[0].size() && jm < more_gap[i].size() && N_insertions[0][j0].index <= more_gap[i][jm]->start1)nsum1 += N_insertions[0][j0++].number;
            while (jn < N_insertions[i].size() && jm < more_gap[i].size() && N_insertions[i][jn].index <= more_gap[i][jm]->start2)nsum2 += N_insertions[i][jn++].number;
            while (j0 < N_insertions[0].size() && jm < more_gap[i].size() && N_insertions[0][j0].index < more_gap[i][jm]->end1)gapn1.push_back(N_insertions[0][j0++]);
            while (jn < N_insertions[i].size() && jm < more_gap[i].size() && N_insertions[i][jn].index < more_gap[i][jm]->end2)gapn2.push_back(N_insertions[i][jn++]);
            insertion_gap_out(os, i, final_sequences, sequences, name, sign, more_gap[i][jm], nsum1, nsum2, gapn1, gapn2, thresh1);
            jm++;
        }
    }
}


