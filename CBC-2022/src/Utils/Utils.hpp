#pragma once
//读入，写回，预处理，后处理，统一化
#include "Fasta.hpp"
#include "Pseudo.hpp"
#include "Insertion.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>


#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#include <process.h>
#include <io.h>
inline void EmptySet()
{
    EmptyWorkingSet(GetCurrentProcess());
}
void getFiles_win(std::string path, std::vector<std::string>& files);
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/types.h>
#include <dirent.h>
#include <malloc.h>
#include <unistd.h>
#include <sys/resource.h>
#include <pthread.h>
inline void EmptySet()
{
    malloc_trim(0);
}
void getFiles_linux(std::string path, std::vector<std::string>& filenames);
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

#undef max
#undef min

size_t getPeakRSS();

//测试内存峰值
inline void GetMemoryUsage()
{
    /*uint64_t mem = 0;
    int pid = getpid();
    PROCESS_MEMORY_COUNTERS pmc;
    HANDLE process = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pid);
    if (GetProcessMemoryInfo(process, &pmc, sizeof(pmc)))
        mem = pmc.PeakWorkingSetSize;
    CloseHandle(process);
    // convert mem from B to MB
    mem = mem / 1024.0 / 1024.0;*/
    int mem = getPeakRSS() / 1024.0 / 1024.0;
    std::cout << "****process mem****" << std::endl;
    std::cout << "current pid: " << getpid() << std::endl;
    std::cout << "memory usage: " << mem << "MB" << std::endl;
}

namespace utils
{ 
    struct block
    {
        int name;
        size_t start;
        size_t length;
        std::vector<unsigned char> seqi;
    };
    struct reads
    {
        std::string name;
        size_t length;
        size_t all_length;
        std::vector<unsigned char> read;
        std::vector<unsigned char> read_ni;
        std::vector<utils::Insertion> n_gap;
    };
    struct in_block
    {
        float score;
        float score_100;
        size_t start;
        size_t end;
        std::vector<size_t> name;
        std::vector<size_t> length;
        std::vector<size_t> si;
        in_block* next;
    };
    struct m_block
    {
        size_t start1;
        size_t end1;
        size_t start2;
        size_t end2;
        std::vector<std::tuple<int, int>> gap1;
        std::vector<std::tuple<int, int>> gap2;
    };
    using more_block = std::vector<m_block*>;
    struct MAF_block
    {
        float score;
        int tag_num;  //记录有tag的列数，连续性
        std::vector<block> seq;
    };
    int* Compare_two(std::vector<unsigned char>& s1, std::vector<unsigned char>& s2, int d = 0, int e = 0, int m = 1, int mis = 0);//求比对分数int d = 3, int e = 1, int m = 1, int mis = -2
    std::string remove_white_spaces(const std::string &str); //去掉空格

    unsigned char to_pseudo(char c); //预处理，char->int

    std::vector<unsigned char> to_pseudo(const std::string &str);
    std::string from_pseudo(const std::vector<unsigned char> &pseu);

    template<typename InputIterator, typename OutputIterator>
    void transform_to_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)   //预处理，char->int
    {
        std::vector<unsigned char> (*op)(const std::string &) = &to_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator, typename OutputIterator>
    void transform_from_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des) //后处理，int->char
    {
        std::string (*op)(const std::vector<unsigned char> &) = &from_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator>
    InputIterator iter_of_max(InputIterator first, InputIterator last) //找到最大的元素
    {
        auto result = first;

        for (; first != last; ++first) if (*result < *first) result = first;
        return result;
    }
    bool next_reads(utils::reads& read, std::istream& is);
    std::vector<std::vector<unsigned char>> read_to_pseudo(std::istream &is, std::vector<std::string>& name, std::vector<std::vector<utils::Insertion>> & N_gap);
    void read_to_pseudo(std::vector<std::vector<unsigned char>>& zh_sequences,
        std::istream& is, std::vector<std::size_t>& all_Length, std::vector<std::size_t>& Length, std::vector<std::string>& name,
        std::vector<std::vector<utils::Insertion>>& N_gap);
    void insertion_gap_out(std::ostream& os, int seqi, int* final_sequences, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::string>& name, std::vector<bool>& sign, utils::m_block* more_block, int nsum1, int nsum2, std::vector<Insertion>& N_insertion1, std::vector<Insertion>& N_insertion2, int thresh1);
    void insertion_gap_more(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, std::vector<utils::more_block>& more_gap, int thresh1);
    void insert_and_write(std::ostream &os, std::istream &is, const std::vector<std::vector<Insertion>> &insertions); //原文件写回比对结果
    void insert_and_write_file(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign); //写回比对结果
    int* vector_insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions);
    int* insertion_gap_N(std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions);
    void insert_and_write_fasta(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign);//写回fasta
    void insert_and_write_maf(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, int thresh);//写回maf比对结果
    // 按100分割，但不会前后回溯 
    //void insert_and_write_100_maf(std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, int thresh1=100, int thresh2=1, int thresh3=95);//长度，条数，分数
    // 按100分割，可前后回溯，同时可以筛选 
    void new_insert_and_write_100_maf(std::string out_file_name, std::ostream& os, std::vector<std::vector<unsigned char>>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign, std::vector<utils::more_block>& more_gap, int thresh1 = 100, int thresh2 = 1, int thresh3 = 95);//长度，条数，分数
    // 离线筛选 maf
    void read_tmp_filter_out_maf(std::istream& is, std::ostream& os, int thresh1, int thresh2, int thresh3);
    //用全空的列来分割，bug
    //void insert_and_write_all_maf(std::ostream& os, std::vector<utils::seqs>& sequences, std::vector<std::vector<Insertion>>& insertions, const std::vector<std::vector<Insertion>>& N_insertions, std::vector<std::string>& name, std::vector<bool>& sign);//长度，条数，分数
    template<typename InputIterator>
    static void cut_and_write(std::ostream &os, InputIterator first, InputIterator last) //一条长序列分多行写入
    {
        const size_t sequence_length = std::distance(first, last);

        for (size_t i = 0; i < sequence_length; i += Fasta::max_line_length)
        {
            if (i) os << '\n';

            size_t write_length = sequence_length - i;
            if (write_length > Fasta::max_line_length) write_length = Fasta::max_line_length;

            const auto begin = first;
            std::advance(first, write_length);
            std::copy(begin, first, std::ostream_iterator<decltype(*first)>(os));
        }
    }

}

template<typename Representation, typename Period>
std::ostream &operator<<(std::ostream &os, std::chrono::duration<Representation, Period> duration) //时间消耗
{
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms";
    return os;
}
