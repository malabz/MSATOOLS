#include "StarAligner.hpp"
#include "../Utils/Pseudo.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Arguments.hpp"
#include "../PairwiseAlignment/NeedlemanWunshReusable.hpp"


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits.h>


//先调用类初始化，再调用_get_gaps； 传入int序列数据； 返回最后结果，长度为n的vector，每个vec_set_centretor存储若干个Insertion，有index+number
void star_alignment::StarAligner::main_pairwise_align(
    std::vector<std::vector<unsigned char>>& _a_sequence,
    std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>*>& _bwta,
    std::ifstream& _is,
    std::ofstream& _os,
    std::vector<std::vector<utils::Insertion>>& _N_gap,
    std::vector<std::string>& _name,
    std::vector<std::size_t>& _all_Length,
    std::vector<std::size_t>& _Length,
    int threshold1, int threshold2, int threshold3) 
{
    StarAligner(_a_sequence, _bwta, _is, _os, _N_gap, _name, _all_Length, _Length)
        ._pairwise_align(threshold1, threshold2, threshold3);//传入参数后，调用_get_gaps第二层的返回结果
}
//类初始化
star_alignment::StarAligner::StarAligner(
    std::vector<std::vector<unsigned char>>& _a_sequence,
    std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>*>& _bwta,
    std::ifstream& _is,
    std::ofstream& _os,
    std::vector<std::vector<utils::Insertion>>& _N_gap,
    std::vector<std::string>& _name,
    std::vector<std::size_t>& _all_Length,
    std::vector<std::size_t>& _Length)
        : a_sequence(_a_sequence) //序列传入
        , bwta(_bwta)
        , is(_is)
        , os(_os)
        , N_gap(_N_gap)  
        , name(_name)
        , all_Length(_all_Length)
        , Length(_Length)
{}

bool star_alignment::StarAligner::cmp_Substrings(Substrings x, Substrings y) {
    return x.add_len > y.add_len;
}

bool star_alignment::StarAligner::cmp_Sub_align(Sub_align x, Sub_align y) {
    return x.b_start < y.b_start;
}

void star_alignment::StarAligner::NW2(std::vector<Sub_align>& align, Sub_align& align_tmp, utils::reads& now_read)
{
    //已经裁剪完的align_tmp，对他的收尾nw操作，注意正负，区间都是闭区间。
    //做完后更新 a,b的start end。
    std::vector<unsigned char>A_ans, B_ans;
    unsigned char* A = Kband::B;
    unsigned char* B = Kband::A;
    int a_len, b_len, i, j, ansa = 0, ansb = 0;
    if (align_tmp.sign)
    {
        b_len = align_tmp.b_start - align_tmp.nw_start;
        b_len = b_len > 1000 ? 1000 : b_len;
        if (b_len != 0 && align_tmp.a_start > 0)//qian b_len==0不操作.
        {
            a_len = (int)(b_len * 1.2);
            a_len = align_tmp.a_start > a_len ? a_len : align_tmp.a_start;
            for (i = 0, j = align_tmp.b_start - 1; i < b_len; i++, j--)
                B[i] = now_read.read[j];
            for (i = 0, j = align_tmp.a_start - 1; i < a_len; i++, j--)
                A[i] = a_sequence[align_tmp.id][j];
            //Align(A_ans, B_ans, A, B);
            Kband::PSA_AGP_Kband_SW(A_ans, B_ans, a_len, b_len, ansa, ansb);
            align_tmp.a_start -= ansa;
            align_tmp.b_start -= ansb;
        }

        b_len = align_tmp.nw_end - align_tmp.b_end;
        b_len = b_len > 1000 ? 1000 : b_len;
        if (b_len != 0 && align_tmp.a_end < a_sequence[align_tmp.id].size())//hou
        {
            a_len = (int)(b_len * 1.2);
            a_len = (a_sequence[align_tmp.id].size() - align_tmp.a_end) > a_len ? a_len : (a_sequence[align_tmp.id].size() - align_tmp.a_end);
            for (i = 0, j = align_tmp.b_end + 1; i < b_len; i++, j++)
                B[i] = now_read.read[j];
            for (i = 0, j = align_tmp.a_end + 1; i < a_len; i++, j++)
                A[i] = a_sequence[align_tmp.id][j];
            //Align(A_ans, B_ans, A, B);
            Kband::PSA_AGP_Kband_SW(A_ans, B_ans, a_len, b_len, ansa, ansb);
            align_tmp.a_end += ansa;
            align_tmp.b_end += ansb;
        }
    }
    else
    {
        b_len = align_tmp.b_start - align_tmp.nw_start;
        b_len = b_len > 1000 ? 1000 : b_len;
        if (b_len != 0 && align_tmp.a_end < a_sequence[align_tmp.id].size())//b-qian a-hou
        {
            a_len = (int)(b_len * 1.2);
            a_len = (a_sequence[align_tmp.id].size() - align_tmp.a_end) > a_len ? a_len : (a_sequence[align_tmp.id].size() - align_tmp.a_end);
            for (i = 0, j = now_read.length - align_tmp.b_start; i < b_len; i++, j++)
                B[i] = now_read.read_ni[j];
            for (i = 0, j = align_tmp.a_end + 1; i < a_len; i++, j++)
                A[i] = a_sequence[align_tmp.id][j];
            //Align(A_ans, B_ans, A, B);
            Kband::PSA_AGP_Kband_SW(A_ans, B_ans, a_len, b_len, ansa, ansb);
            align_tmp.a_end += ansa;
            align_tmp.b_start -= ansb;
        }

        b_len = align_tmp.nw_end - align_tmp.b_end;
        b_len = b_len > 1000 ? 1000 : b_len;
        if (b_len != 0 && align_tmp.a_start > 0)//b-hou a-qian
        {
            a_len = (int)(b_len * 1.2);
            a_len = align_tmp.a_start > a_len ? a_len : align_tmp.a_start;
            for (i = 0, j = now_read.length - align_tmp.b_end - 2; i < b_len; i++, j--)
                B[i] = now_read.read_ni[j];
            for (i = 0, j = align_tmp.a_start - 1; i < a_len; i++, j--)
                A[i] = a_sequence[align_tmp.id][j];
            //Align(A_ans, B_ans, A, B);
            Kband::PSA_AGP_Kband_SW(A_ans, B_ans, a_len, b_len, ansa, ansb);
            align_tmp.a_start -= ansa;
            align_tmp.b_end += ansb;
        }
    }
    align.push_back(align_tmp);
}

void star_alignment::StarAligner::_pairwise_align(int threshold1,int threshold2, int threshold3) 
{
    int  tmp_int;// name_len = 0;
    insert a_gap;
    insert b_gap;
    utils::reads now_read;
    std::vector<unsigned char> read_ni();
    std::vector<Substrings> C_Strings;
    std::vector<Sub_align> align;
    std::vector<triple> common_substrings;
    Sub_align align_tmp;
    while (utils::next_reads(now_read, is))
    {
        std::vector<Substrings>().swap(C_Strings);
        //std::cout << "@ " << now_read.name << " " << now_read.all_length<< " " << now_read.length << "*********************************\n";
        for (int i = 0; i < name.size(); i++)
        {
            common_substrings = _optimal_path(bwta[i]->get_common_substrings(now_read.read.cbegin(), now_read.read.cend(), threshold1, true));
            get_common_substrings_vector(C_Strings, common_substrings, i, true, now_read.length);
            
            common_substrings = _optimal_path(bwta[i]->get_common_substrings(now_read.read_ni.cbegin(), now_read.read_ni.cend(), threshold1, true));
            get_common_substrings_vector(C_Strings, common_substrings, i, false, now_read.length); 
        }
        sort(C_Strings.begin(), C_Strings.end(), cmp_Substrings); //sort C_Strings
        /*for (int i = 0; i < C_Strings.size(); i++)
        {
            std::cout << ">>> " << name[C_Strings[i].id] << " "<< C_Strings[i].sign <<" "<< C_Strings[i].add_len << "   \n";
            for (int j = 0; j < C_Strings[i].Subs.size(); j++)
            {
                std::cout << "  ??? "<<j<<" \n";
                for (int k = 0; k < C_Strings[i].Subs[j].size(); k++)
                {
                    std::cout << C_Strings[i].Subs[j][k][1] << "\t" << C_Strings[i].Subs[j][k][0] << "\t" << C_Strings[i].Subs[j][k][2] << "\n";
                }
            }
        }
        std::cout << "\n\n";*/
        std::vector<Sub_align>().swap(align);
        int k, sc, seek1, seek2;
        //std::cout << "zt1\n";
        for (int i = 0; i < C_Strings.size(); i++)
        {
            for (int j = 0; j < C_Strings[i].Subs.size(); j++)
            {
                align_tmp.id = C_Strings[i].id;
                align_tmp.sign = C_Strings[i].sign;
                align_tmp.a_start = C_Strings[i].Subs[j][0][0];
                align_tmp.a_end = C_Strings[i].Subs[j].back()[0] + C_Strings[i].Subs[j].back()[2] - 1; //闭区间
                align_tmp.b_start = C_Strings[i].Subs[j][0][1];
                align_tmp.b_end = C_Strings[i].Subs[j].back()[1] + C_Strings[i].Subs[j].back()[2] - 1;
                if (!align_tmp.sign)
                {
                    tmp_int = now_read.length - 1 - align_tmp.b_start; // tmp <- end
                    align_tmp.b_start = now_read.length - 1 - align_tmp.b_end; // start <- start
                    align_tmp.b_end = tmp_int;// end <- tmp 
                }
                align_tmp.nw_start = 0;
                align_tmp.nw_end = now_read.length - 1;
                /// /////////////////
                if (i == 0 && j == 0)
                {
                    NW2(align, align_tmp, now_read);
                }
                else
                {
                    for (k = 0; k < align.size(); k++)
                    {
                        if ((align_tmp.b_start >= align[k].b_start && align_tmp.b_end <= align[k].b_end) || (align_tmp.b_start <= align[k].b_start && align_tmp.b_end >= align[k].b_end)) // ||后面情况基本不会发生
                            break;//直接舍弃该tmp 
                        else if (align_tmp.b_start > align[k].b_end)
                        {
                            if (align_tmp.nw_start < (align[k].b_end + 1))
                                align_tmp.nw_start = align[k].b_end + 1;
                        }
                        else if (align_tmp.b_end < align[k].b_start)
                        {
                            if(align_tmp.nw_end > (align[k].b_start - 1))
                                align_tmp.nw_end = align[k].b_start - 1;
                        }
                        else if ((align_tmp.b_start < align[k].b_start) && (align_tmp.b_end >= align[k].b_start)) //tmp.end 擦 ans.start了,从尾删
                        {
                            if (align_tmp.sign)//从尾删
                            {
                                C_Strings[i].Subs[j].pop_back();
                                while (C_Strings[i].Subs[j].size() > 0)
                                {
                                    align_tmp.b_end = C_Strings[i].Subs[j].back()[1] + C_Strings[i].Subs[j].back()[2] - 1;
                                    if (align_tmp.b_end < align[k].b_start) break;
                                    C_Strings[i].Subs[j].pop_back();
                                }
                                if (C_Strings[i].Subs[j].size() == 0)//删没了
                                    break;//直接舍弃该tmp 
                                align_tmp.nw_end = align[k].b_start - 1;
                            }
                            else //从头删
                            {
                                C_Strings[i].Subs[j].erase(C_Strings[i].Subs[j].begin());
                                while (C_Strings[i].Subs[j].size() > 0)
                                {
                                    align_tmp.b_end = now_read.length - 1 - C_Strings[i].Subs[j][0][1];
                                    if (align_tmp.b_end < align[k].b_start) break;
                                    C_Strings[i].Subs[j].erase(C_Strings[i].Subs[j].begin());
                                }
                                if (C_Strings[i].Subs[j].size() == 0)//删没了
                                    break;//直接舍弃该tmp 
                                align_tmp.nw_end = align[k].b_start - 1;
                            }
                            
                        }
                        else if ((align_tmp.b_start <= align[k].b_end) && (align_tmp.b_end > align[k].b_end)) //tmp.start 擦 ans.end了，从头删
                        {
                            if (align_tmp.sign)//从头删
                            {
                                C_Strings[i].Subs[j].erase(C_Strings[i].Subs[j].begin());
                                while (C_Strings[i].Subs[j].size() > 0)
                                {
                                    align_tmp.b_start = C_Strings[i].Subs[j][0][1];
                                    if (align_tmp.b_start > align[k].b_end) break;
                                    C_Strings[i].Subs[j].erase(C_Strings[i].Subs[j].begin());
                                }
                                if (C_Strings[i].Subs[j].size() == 0)//删没了
                                    break;//直接舍弃该tmp 
                                align_tmp.nw_start = align[k].b_end + 1;
                            }
                            else //从尾删
                            {
                                C_Strings[i].Subs[j].pop_back();
                                while (C_Strings[i].Subs[j].size() > 0)
                                {
                                    align_tmp.b_start = now_read.length - C_Strings[i].Subs[j].back()[1]- C_Strings[i].Subs[j].back()[2];
                                    if (align_tmp.b_start > align[k].b_end) break;
                                    C_Strings[i].Subs[j].pop_back();
                                }
                                if (C_Strings[i].Subs[j].size() == 0)//删没了
                                    break;//直接舍弃该tmp 
                                align_tmp.nw_start = align[k].b_end + 1;
                            }
                        }
                        else
                            break; //此种情况不会发生
                    } 
                    if(k == align.size() && ((align_tmp.b_end - align_tmp.b_start)> now_read.length*0.02)) NW2(align, align_tmp, now_read);//按thresh筛选一下
                }
            }
        }
        sort(align.begin(), align.end(), cmp_Sub_align); //sort align
        //改a start end 的坐标
        for (int i = 0; i < align.size(); i++)
        {
            os << now_read.name << "\t" << name[align[i].id] << "\t";
            if (align[i].sign)os << "+\t";
            else os << "-\t";
            
            seek1=0; seek2=0;
            for (k = 0; k < N_gap[align[i].id].size(); k++)
            {
                if (align[i].a_start > N_gap[align[i].id][k].index)
                    seek1 += N_gap[align[i].id][k].number;
                if (align[i].a_end > N_gap[align[i].id][k].index)
                    seek2 += N_gap[align[i].id][k].number;
                else
                    break;
            }

            os << align[i].a_start+ seek1 << "\t" << align[i].a_end+ seek2 << "\n";
        }
    }
}

void star_alignment::StarAligner::get_common_substrings_vector(std::vector<Substrings>& C_Strings, std::vector<triple>& substrings, int id, bool sign, size_t read_len)
{
    Substrings tmp;
    tmp.id = id;
    tmp.sign = sign;
    size_t thresh1 = read_len * 0.05; //间隔
    size_t thresh2 = read_len * 0.02; //单独一个,或者小集合
    float thresh3 = 0.2;
    thresh2 = thresh2 > 50 ? thresh2 : 50;
    std::vector<std::vector<triple>> Subs;
    if (substrings.size() < 1)    //[A.index，B.index，length] , A is center
        return;
    bool tag = false; //是否已经开头
    for (int i = 1; i < substrings.size(); i++)
    {
        if(tag)
        {
            if (std::abs((int)(substrings[i][0] - substrings[i-1][0] - substrings[i][1] + substrings[i - 1][1])) < thresh1)
            {
                Subs.back().push_back(substrings[i]);
            }
            else
            {
                tag = false;
            }
        }
        else
        {
            if (std::abs((int)(substrings[i][0] - substrings[i - 1][0] - substrings[i][1] + substrings[i - 1][1])) < thresh1)
            {
                Subs.push_back(std::vector<triple>());
                Subs.back().push_back(substrings[i-1]);
                Subs.back().push_back(substrings[i]);
                tag = true;
            }
            else
            {
                if (substrings[i - 1][2] > thresh2)
                {
                    Subs.push_back(std::vector<triple>());
                    Subs.back().push_back(substrings[i - 1]);
                }
            }
        }
    }
    int x, y, z;
    for (int i = 0; i < Subs.size(); i++)
    { 
        x = 0; y = 0; z = 0;
        for (int j = 0; j < Subs[i].size(); j++)
            x += Subs[i][j][2];
        y = Subs[i].back()[1] - Subs[i][0][1];
        z = Subs[i].back()[0] - Subs[i][0][0];
        if(Subs[i].size()>3)
            tmp.Subs.push_back(Subs[i]);
        else if (y * thresh3 <= x && z * thresh3 <= x)
            tmp.Subs.push_back(Subs[i]);
    }
    if (tmp.Subs.size() == 0)
        return;
    tmp.add_len = 0;
    for (int i = 0; i < tmp.Subs.size(); i++)
        tmp.add_len += (tmp.Subs[i].back()[1] + tmp.Subs[i].back()[2] - tmp.Subs[i][0][1]);
    if (tmp.add_len>=thresh2)
        C_Strings.push_back(tmp);
}

//分两步新改,第一步，从传入全部的同源区段，[[A.index，B.index，length]...]，每个B.index选出一个A.index，相对距离最近
auto star_alignment::StarAligner::_MultiReg(const std::vector<triple>& common_substrings) //最优路径
    -> std::vector<triple> 
{
    std::cout << "***common_substrings" << common_substrings.size() << "\n";
    /*for (int i = 0; i < common_substrings.size(); ++i)
        std::cout << common_substrings[i][0] << " " << common_substrings[i][1] << " " << common_substrings[i][2] << "\n";
    */
    std::vector<triple> optimal_common_substrings;
    if (common_substrings.empty()) return optimal_common_substrings;
    int start = common_substrings[0][0];  //初始化
    int end = common_substrings[0][0] + common_substrings[0][2];//
    int i;
    float a_length,tmp, pre_tmp,b_length = common_substrings.rbegin()[0][1] + common_substrings.rbegin()[0][2]; //B的长度
    for (i = 1; i < common_substrings.size(); i++) // 找出所有后缀号列表中，最小后缀号做  start起点，最大后缀号 + length_i做end终点
    {
        if (common_substrings[i][0] < start) start = common_substrings[i][0];
        if ((common_substrings[i][0]+ common_substrings[i][2]) > end) end = common_substrings[i][0]+ common_substrings[i][2];
    }
    a_length = end - start;
    //std::cout << a_length << " "<< b_length <<"\n";
    i = 0;
    //while ((common_substrings[i][2] < threshold)) i++;
    optimal_common_substrings.push_back(common_substrings[i]);//first
    pre_tmp = common_substrings[i][0] / a_length - common_substrings[i][1] / b_length;//加入第一个元素
    pre_tmp = fabs(pre_tmp);
    for (i=i+1; i < common_substrings.size(); i++)
    {
        //if (common_substrings[i][2] < threshold)continue;
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
    
    std::cout << "***optimal_common_substrings" << optimal_common_substrings.size() << "\n";
    //for (int i = 0; i < optimal_common_substrings.size(); ++i)
        //std::cout << optimal_common_substrings[i][0] << " " << optimal_common_substrings[i][1] << " " << optimal_common_substrings[i][2] << "\n";
    return optimal_common_substrings;

}

//回溯
std::vector<int> star_alignment::StarAligner::_trace_back(const std::vector<triple>& common_substrings, int* p)
{
    std::vector<int> ansi;
    int* tmp = p;
    for (int i = 0; i < common_substrings.size(); i++)
    {
        if (*p < tmp[i])
            p = &tmp[i];
    }
    int j,i = p - tmp;
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

//分两步新改,第二步，依据动态规划，选出合适的不重叠的同源区段
auto star_alignment::StarAligner::_optimal_path(const std::vector<triple>& common_substrings) //最优路径
-> std::vector<triple>   // [A.index，B.index，length]
{
    std::vector<triple> optimal_common_substrings = common_substrings;//调用第一步结果
    //std::cout << "***_optimal_path" << optimal_common_substrings.size() << "\n";
    if (optimal_common_substrings.empty()) return optimal_common_substrings;

    int m = optimal_common_substrings.size();
    if (m<=1) return optimal_common_substrings;

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
    //std::cout << "***_optimal_path" << ans_common_substrings.size() << "\n";
    //for (int i = 0; i < ans_common_substrings.size(); ++i)
        //std::cout << ans_common_substrings[i][0] << " " << ans_common_substrings[i][1] << " " << ans_common_substrings[i][2] << "\n";

    return ans_common_substrings;
}
