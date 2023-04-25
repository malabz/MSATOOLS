/*
 ************************************************************************************************
 *                                                                                              *
 *                                   比较序列相似性                                             *
 *                                                                                              *
 * Author: ZhouTong                                                                             *
 * Date: 2023/4/25                                                                              *
 *                                                                                              *
 * 1.将序列 seq_center 拆分为长度为 k 的子串，并计算每个子串的哈希值，存储到 Bloom Filter 中；  *
 * 2.遍历序列 seq_i，将其拆分为长度为 k 的子串，并计算每个子串的哈希值，                        *
 *   然后在 Bloom Filter 中查找是否存在相同的哈希值，若存在，则说明存在相同的子串，统计数量；   *
 * 3.比较统计出的相同子串数量。                                                                 *
 * 4.这里使用的 MurmurHash 算法计算哈希值，可以快速高效地生成哈希值。而 Bloom Filter            *
 *   则可以用较小的内存空间对大规模数据进行快速查询和去重。                                     *
 * 5.总体来说，该算法属于基于哈希的近似匹配算法，适用于大规模数据集的相似度查询问题。           *
 *                                                                                              *
 * 注：返回数值意义不大，只作比较用                                                             *
 *                                                                                              *
 ************************************************************************************************
*/

#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <iostream>

#ifndef _MURMURHASH3_H_
#define _MURMURHASH3_H_

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

// Microsoft Visual Studio

#if defined(_MSC_VER) && (_MSC_VER < 1600)

typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;

// Other compilers

#else	// defined(_MSC_VER)

#include <stdint.h>

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32(const void* key, int len, uint32_t seed, void* out);

void MurmurHash3_x86_128(const void* key, int len, uint32_t seed, void* out);

void MurmurHash3_x64_128(const void* key, int len, uint32_t seed, void* out);

//-----------------------------------------------------------------------------

#endif // _MURMURHASH3_H_

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline uint32_t rotl32(uint32_t x, int8_t r)
{
    return (x << r) | (x >> (32 - r));
}

inline uint64_t rotl64(uint64_t x, int8_t r)
{
    return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

FORCE_INLINE uint32_t getblock32(const uint32_t* p, int i)
{
    return p[i];
}

FORCE_INLINE uint64_t getblock64(const uint64_t* p, int i)
{
    return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

FORCE_INLINE uint32_t fmix32(uint32_t h)
{
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    return h;
}

//----------

FORCE_INLINE uint64_t fmix64(uint64_t k)
{
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32(const void* key, int len,
    uint32_t seed, void* out)
{
    const uint8_t* data = (const uint8_t*)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t* blocks = (const uint32_t*)(data + nblocks * 4);

    for (int i = -nblocks; i; i++)
    {
        uint32_t k1 = getblock32(blocks, i);

        k1 *= c1;
        k1 = ROTL32(k1, 15);
        k1 *= c2;

        h1 ^= k1;
        h1 = ROTL32(h1, 13);
        h1 = h1 * 5 + 0xe6546b64;
    }

    //----------
    // tail

    const uint8_t* tail = (const uint8_t*)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3)
    {
    case 3: k1 ^= tail[2] << 16;
    case 2: k1 ^= tail[1] << 8;
    case 1: k1 ^= tail[0];
        k1 *= c1; k1 = ROTL32(k1, 15); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len;

    h1 = fmix32(h1);

    *(uint32_t*)out = h1;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x86_128(const void* key, const int len,
    uint32_t seed, void* out)
{
    const uint8_t* data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint32_t h1 = seed;
    uint32_t h2 = seed;
    uint32_t h3 = seed;
    uint32_t h4 = seed;

    const uint32_t c1 = 0x239b961b;
    const uint32_t c2 = 0xab0e9789;
    const uint32_t c3 = 0x38b34ae5;
    const uint32_t c4 = 0xa1e38b93;

    //----------
    // body

    const uint32_t* blocks = (const uint32_t*)(data + nblocks * 16);

    for (int i = -nblocks; i; i++)
    {
        uint32_t k1 = getblock32(blocks, i * 4 + 0);
        uint32_t k2 = getblock32(blocks, i * 4 + 1);
        uint32_t k3 = getblock32(blocks, i * 4 + 2);
        uint32_t k4 = getblock32(blocks, i * 4 + 3);

        k1 *= c1; k1 = ROTL32(k1, 15); k1 *= c2; h1 ^= k1;

        h1 = ROTL32(h1, 19); h1 += h2; h1 = h1 * 5 + 0x561ccd1b;

        k2 *= c2; k2 = ROTL32(k2, 16); k2 *= c3; h2 ^= k2;

        h2 = ROTL32(h2, 17); h2 += h3; h2 = h2 * 5 + 0x0bcaa747;

        k3 *= c3; k3 = ROTL32(k3, 17); k3 *= c4; h3 ^= k3;

        h3 = ROTL32(h3, 15); h3 += h4; h3 = h3 * 5 + 0x96cd1c35;

        k4 *= c4; k4 = ROTL32(k4, 18); k4 *= c1; h4 ^= k4;

        h4 = ROTL32(h4, 13); h4 += h1; h4 = h4 * 5 + 0x32ac3b17;
    }

    //----------
    // tail

    const uint8_t* tail = (const uint8_t*)(data + nblocks * 16);

    uint32_t k1 = 0;
    uint32_t k2 = 0;
    uint32_t k3 = 0;
    uint32_t k4 = 0;

    switch (len & 15)
    {
    case 15: k4 ^= tail[14] << 16;
    case 14: k4 ^= tail[13] << 8;
    case 13: k4 ^= tail[12] << 0;
        k4 *= c4; k4 = ROTL32(k4, 18); k4 *= c1; h4 ^= k4;

    case 12: k3 ^= tail[11] << 24;
    case 11: k3 ^= tail[10] << 16;
    case 10: k3 ^= tail[9] << 8;
    case  9: k3 ^= tail[8] << 0;
        k3 *= c3; k3 = ROTL32(k3, 17); k3 *= c4; h3 ^= k3;

    case  8: k2 ^= tail[7] << 24;
    case  7: k2 ^= tail[6] << 16;
    case  6: k2 ^= tail[5] << 8;
    case  5: k2 ^= tail[4] << 0;
        k2 *= c2; k2 = ROTL32(k2, 16); k2 *= c3; h2 ^= k2;

    case  4: k1 ^= tail[3] << 24;
    case  3: k1 ^= tail[2] << 16;
    case  2: k1 ^= tail[1] << 8;
    case  1: k1 ^= tail[0] << 0;
        k1 *= c1; k1 = ROTL32(k1, 15); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len; h3 ^= len; h4 ^= len;

    h1 += h2; h1 += h3; h1 += h4;
    h2 += h1; h3 += h1; h4 += h1;

    h1 = fmix32(h1);
    h2 = fmix32(h2);
    h3 = fmix32(h3);
    h4 = fmix32(h4);

    h1 += h2; h1 += h3; h1 += h4;
    h2 += h1; h3 += h1; h4 += h1;

    ((uint32_t*)out)[0] = h1;
    ((uint32_t*)out)[1] = h2;
    ((uint32_t*)out)[2] = h3;
    ((uint32_t*)out)[3] = h4;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x64_128(const void* key, const int len,
    const uint32_t seed, void* out)
{
    const uint8_t* data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
    const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

    //----------
    // body

    const uint64_t* blocks = (const uint64_t*)(data);

    for (int i = 0; i < nblocks; i++)
    {
        uint64_t k1 = getblock64(blocks, i * 2 + 0);
        uint64_t k2 = getblock64(blocks, i * 2 + 1);

        k1 *= c1; k1 = ROTL64(k1, 31); k1 *= c2; h1 ^= k1;

        h1 = ROTL64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;

        k2 *= c2; k2 = ROTL64(k2, 33); k2 *= c1; h2 ^= k2;

        h2 = ROTL64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
    }

    //----------
    // tail

    const uint8_t* tail = (const uint8_t*)(data + nblocks * 16);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch (len & 15)
    {
    case 15: k2 ^= ((uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((uint64_t)tail[9]) << 8;
    case  9: k2 ^= ((uint64_t)tail[8]) << 0;
        k2 *= c2; k2 = ROTL64(k2, 33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((uint64_t)tail[7]) << 56;
    case  7: k1 ^= ((uint64_t)tail[6]) << 48;
    case  6: k1 ^= ((uint64_t)tail[5]) << 40;
    case  5: k1 ^= ((uint64_t)tail[4]) << 32;
    case  4: k1 ^= ((uint64_t)tail[3]) << 24;
    case  3: k1 ^= ((uint64_t)tail[2]) << 16;
    case  2: k1 ^= ((uint64_t)tail[1]) << 8;
    case  1: k1 ^= ((uint64_t)tail[0]) << 0;
        k1 *= c1; k1 = ROTL64(k1, 31); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}

//-----------------------------------------------------------------------------





#ifndef INCLUDE_BLOOM_FILTER_HPP
#define INCLUDE_BLOOM_FILTER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <string>
#include <vector>
#include <immintrin.h>// 使用SIMD指令，需要包含此头文件

static const std::size_t bits_per_char = 0x08;    // 8 bits in 1 char(unsigned)

static const unsigned char bit_mask[bits_per_char] = {
                                                       0x01,  //00000001
                                                       0x02,  //00000010
                                                       0x04,  //00000100
                                                       0x08,  //00001000
                                                       0x10,  //00010000
                                                       0x20,  //00100000
                                                       0x40,  //01000000
                                                       0x80   //10000000
};

class bloom_parameters
{
public:

    bloom_parameters()
        : minimum_size(1),
        maximum_size(std::numeric_limits<unsigned long long int>::max()),
        minimum_number_of_hashes(1),
        maximum_number_of_hashes(std::numeric_limits<unsigned int>::max()),
        projected_element_count(10000),
        false_positive_probability(1.0 / projected_element_count),
        random_seed(0xA5A5A5A55A5A5A5AULL)
    {}

    virtual ~bloom_parameters()
    {}

    inline bool operator!()
    {
        return (minimum_size > maximum_size) ||
            (minimum_number_of_hashes > maximum_number_of_hashes) ||
            (minimum_number_of_hashes < 1) ||
            (0 == maximum_number_of_hashes) ||
            (0 == projected_element_count) ||
            (false_positive_probability < 0.0) ||
            (std::numeric_limits<double>::infinity() == std::abs(false_positive_probability)) ||
            (0 == random_seed) ||
            (0xFFFFFFFFFFFFFFFFULL == random_seed);
    }

    // Allowable min/max size of the bloom filter in bits
    unsigned long long int minimum_size;
    unsigned long long int maximum_size;

    // Allowable min/max number of hash functions
    unsigned int minimum_number_of_hashes;
    unsigned int maximum_number_of_hashes;

    // The approximate number of elements to be inserted
    // into the bloom filter, should be within one order
    // of magnitude. The default is 10000.
    unsigned long long int projected_element_count;

    // The approximate false positive probability expected
    // from the bloom filter. The default is assumed to be
    // the reciprocal of the projected_element_count.
    double false_positive_probability;

    unsigned long long int random_seed;

    struct optimal_parameters_t
    {
        optimal_parameters_t()
            : number_of_hashes(0),
            table_size(0)
        {}

        unsigned int number_of_hashes;
        unsigned long long int table_size;
    };

    optimal_parameters_t optimal_parameters;

    virtual bool compute_optimal_parameters()
    {
        /*
          Note:
          The following will attempt to find the number of hash functions
          and minimum amount of storage bits required to construct a bloom
          filter consistent with the user defined false positive probability
          and estimated element insertion count.
        */

        if (!(*this))
            return false;

        double min_m = std::numeric_limits<double>::infinity();
        double min_k = 0.0;
        double k = 1.0;

        while (k < 1000.0)
        {
            const double numerator = (-k * projected_element_count);
            const double denominator = std::log(1.0 - std::pow(false_positive_probability, 1.0 / k));

            const double curr_m = numerator / denominator;

            if (curr_m < min_m)
            {
                min_m = curr_m;
                min_k = k;
            }

            k += 1.0;
        }

        optimal_parameters_t& optp = optimal_parameters;

        optp.number_of_hashes = static_cast<unsigned int>(min_k);

        optp.table_size = static_cast<unsigned long long int>(min_m);

        optp.table_size += (((optp.table_size % bits_per_char) != 0) ? (bits_per_char - (optp.table_size % bits_per_char)) : 0);

        if (optp.number_of_hashes < minimum_number_of_hashes)
            optp.number_of_hashes = minimum_number_of_hashes;
        else if (optp.number_of_hashes > maximum_number_of_hashes)
            optp.number_of_hashes = maximum_number_of_hashes;

        if (optp.table_size < minimum_size)
            optp.table_size = minimum_size;
        else if (optp.table_size > maximum_size)
            optp.table_size = maximum_size;

        return true;
    }

};

class bloom_filter
{
protected:

    typedef unsigned int bloom_type;
    typedef unsigned char cell_type;
    typedef std::vector<unsigned char> table_type;

public:

    bloom_filter()
        : salt_count_(0),
        table_size_(0),
        projected_element_count_(0),
        inserted_element_count_(0),
        random_seed_(0),
        desired_false_positive_probability_(0.0)
    {}

    bloom_filter(const bloom_parameters& p)
        : projected_element_count_(p.projected_element_count),
        inserted_element_count_(0),
        random_seed_((p.random_seed * 0xA5A5A5A5) + 1),
        desired_false_positive_probability_(p.false_positive_probability)
    {
        salt_count_ = p.optimal_parameters.number_of_hashes;
        table_size_ = p.optimal_parameters.table_size;

        generate_unique_salt();

        bit_table_.resize(table_size_ / bits_per_char, static_cast<unsigned char>(0x00));
    }

    bloom_filter(const bloom_filter& filter)
    {
        this->operator=(filter);
    }

    inline bool operator == (const bloom_filter& f) const
    {
        if (this != &f)
        {
            return
                (salt_count_ == f.salt_count_) &&
                (table_size_ == f.table_size_) &&
                (bit_table_.size() == f.bit_table_.size()) &&
                (projected_element_count_ == f.projected_element_count_) &&
                (inserted_element_count_ == f.inserted_element_count_) &&
                (random_seed_ == f.random_seed_) &&
                (desired_false_positive_probability_ == f.desired_false_positive_probability_) &&
                (salt_ == f.salt_) &&
                (bit_table_ == f.bit_table_);
        }
        else
            return true;
    }

    inline bool operator != (const bloom_filter& f) const
    {
        return !operator==(f);
    }

    inline bloom_filter& operator = (const bloom_filter& f)
    {
        if (this != &f)
        {
            salt_count_ = f.salt_count_;
            table_size_ = f.table_size_;
            bit_table_ = f.bit_table_;
            salt_ = f.salt_;

            projected_element_count_ = f.projected_element_count_;
            inserted_element_count_ = f.inserted_element_count_;

            random_seed_ = f.random_seed_;

            desired_false_positive_probability_ = f.desired_false_positive_probability_;
        }

        return *this;
    }

    virtual ~bloom_filter()
    {}

    inline bool operator!() const
    {
        return (0 == table_size_);
    }

    inline void clear()
    {
        std::fill(bit_table_.begin(), bit_table_.end(), static_cast<unsigned char>(0x00));
        inserted_element_count_ = 0;
    }

    inline void insert(const unsigned char* key_begin, const std::size_t& length)
    {
        std::size_t bit_index = 0;
        std::size_t bit = 0;

        for (std::size_t i = 0; i < salt_.size(); ++i)
        {
            compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);

            bit_table_[bit_index / bits_per_char] |= bit_mask[bit];
        }

        ++inserted_element_count_;
    }

    template <typename T>
    inline void insert(const T& t)
    {
        // Note: T must be a C++ POD type.
        insert(reinterpret_cast<const unsigned char*>(&t), sizeof(T));
    }

    inline void insert(const std::string& key)
    {
        insert(reinterpret_cast<const unsigned char*>(key.data()), key.size());
    }

    inline void insert(const char* data, const std::size_t& length)
    {
        insert(reinterpret_cast<const unsigned char*>(data), length);
    }

    template <typename InputIterator>
    inline void insert(const InputIterator begin, const InputIterator end)
    {
        InputIterator itr = begin;

        while (end != itr)
        {
            insert(*(itr++));
        }
    }

    inline virtual bool contains(const unsigned char* key_begin, const std::size_t length) const
    {
        std::size_t bit_index = 0;
        std::size_t bit = 0;

        for (std::size_t i = 0; i < salt_.size(); ++i)
        {
            compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);

            if ((bit_table_[bit_index / bits_per_char] & bit_mask[bit]) != bit_mask[bit])
            {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    inline bool contains(const T& t) const
    {
        return contains(reinterpret_cast<const unsigned char*>(&t), static_cast<std::size_t>(sizeof(T)));
    }

    inline bool contains(const std::string& key) const
    {
        return contains(reinterpret_cast<const unsigned char*>(key.c_str()), key.size());
    }

    inline bool contains(const char* data, const std::size_t& length) const
    {
        return contains(reinterpret_cast<const unsigned char*>(data), length);
    }

    template <typename InputIterator>
    inline InputIterator contains_all(const InputIterator begin, const InputIterator end) const
    {
        InputIterator itr = begin;

        while (end != itr)
        {
            if (!contains(*itr))
            {
                return itr;
            }

            ++itr;
        }

        return end;
    }

    template <typename InputIterator>
    inline InputIterator contains_none(const InputIterator begin, const InputIterator end) const
    {
        InputIterator itr = begin;

        while (end != itr)
        {
            if (contains(*itr))
            {
                return itr;
            }

            ++itr;
        }

        return end;
    }

    inline virtual unsigned long long int size() const
    {
        return table_size_;
    }

    inline unsigned long long int element_count() const
    {
        return inserted_element_count_;
    }

    inline double effective_fpp() const
    {
        /*
          Note:
          The effective false positive probability is calculated using the
          designated table size and hash function count in conjunction with
          the current number of inserted elements - not the user defined
          predicated/expected number of inserted elements.
        */
        return std::pow(1.0 - std::exp(-1.0 * salt_.size() * inserted_element_count_ / size()), 1.0 * salt_.size());
    }

    inline bloom_filter& operator &= (const bloom_filter& f)
    {
        /* intersection */
        if (
            (salt_count_ == f.salt_count_) &&
            (table_size_ == f.table_size_) &&
            (random_seed_ == f.random_seed_)
            )
        {
            for (std::size_t i = 0; i < bit_table_.size(); ++i)
            {
                bit_table_[i] &= f.bit_table_[i];
            }
        }

        return *this;
    }

    inline bloom_filter& operator |= (const bloom_filter& f)
    {
        /* union */
        if (
            (salt_count_ == f.salt_count_) &&
            (table_size_ == f.table_size_) &&
            (random_seed_ == f.random_seed_)
            )
        {
            for (std::size_t i = 0; i < bit_table_.size(); ++i)
            {
                bit_table_[i] |= f.bit_table_[i];
            }
        }

        return *this;
    }

    inline bloom_filter& operator ^= (const bloom_filter& f)
    {
        /* difference */
        if (
            (salt_count_ == f.salt_count_) &&
            (table_size_ == f.table_size_) &&
            (random_seed_ == f.random_seed_)
            )
        {
            for (std::size_t i = 0; i < bit_table_.size(); ++i)
            {
                bit_table_[i] ^= f.bit_table_[i];
            }
        }

        return *this;
    }

    inline const cell_type* table() const
    {
        return bit_table_.data();
    }

    inline std::size_t hash_count()
    {
        return salt_.size();
    }

protected:

    inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
    {
        bit_index = hash % table_size_;
        bit = bit_index % bits_per_char;
    }

    void generate_unique_salt()
    {
        /*
          Note:
          A distinct hash function need not be implementation-wise
          distinct. In the current implementation "seeding" a common
          hash function with different values seems to be adequate.
        */
        const unsigned int predef_salt_count = 128;

        static const bloom_type predef_salt[predef_salt_count] =
        {
           0xAAAAAAAA, 0x55555555, 0x33333333, 0xCCCCCCCC,
           0x66666666, 0x99999999, 0xB5B5B5B5, 0x4B4B4B4B,
           0xAA55AA55, 0x55335533, 0x33CC33CC, 0xCC66CC66,
           0x66996699, 0x99B599B5, 0xB54BB54B, 0x4BAA4BAA,
           0xAA33AA33, 0x55CC55CC, 0x33663366, 0xCC99CC99,
           0x66B566B5, 0x994B994B, 0xB5AAB5AA, 0xAAAAAA33,
           0x555555CC, 0x33333366, 0xCCCCCC99, 0x666666B5,
           0x9999994B, 0xB5B5B5AA, 0xFFFFFFFF, 0xFFFF0000,
           0xB823D5EB, 0xC1191CDF, 0xF623AEB3, 0xDB58499F,
           0xC8D42E70, 0xB173F616, 0xA91A5967, 0xDA427D63,
           0xB1E8A2EA, 0xF6C0D155, 0x4909FEA3, 0xA68CC6A7,
           0xC395E782, 0xA26057EB, 0x0CD5DA28, 0x467C5492,
           0xF15E6982, 0x61C6FAD3, 0x9615E352, 0x6E9E355A,
           0x689B563E, 0x0C9831A8, 0x6753C18B, 0xA622689B,
           0x8CA63C47, 0x42CC2884, 0x8E89919B, 0x6EDBD7D3,
           0x15B6796C, 0x1D6FDFE4, 0x63FF9092, 0xE7401432,
           0xEFFE9412, 0xAEAEDF79, 0x9F245A31, 0x83C136FC,
           0xC3DA4A8C, 0xA5112C8C, 0x5271F491, 0x9A948DAB,
           0xCEE59A8D, 0xB5F525AB, 0x59D13217, 0x24E7C331,
           0x697C2103, 0x84B0A460, 0x86156DA9, 0xAEF2AC68,
           0x23243DA5, 0x3F649643, 0x5FA495A8, 0x67710DF8,
           0x9A6C499E, 0xDCFB0227, 0x46A43433, 0x1832B07A,
           0xC46AFF3C, 0xB9C8FFF0, 0xC9500467, 0x34431BDF,
           0xB652432B, 0xE367F12B, 0x427F4C1B, 0x224C006E,
           0x2E7E5A89, 0x96F99AA5, 0x0BEB452A, 0x2FD87C39,
           0x74B2E1FB, 0x222EFD24, 0xF357F60C, 0x440FCB1E,
           0x8BBE030F, 0x6704DC29, 0x1144D12F, 0x948B1355,
           0x6D8FD7E9, 0x1C11A014, 0xADD1592F, 0xFB3C712E,
           0xFC77642F, 0xF9C4CE8C, 0x31312FB9, 0x08B0DD79,
           0x318FA6E7, 0xC040D23D, 0xC0589AA7, 0x0CA5C075,
           0xF874B172, 0x0CF914D5, 0x784D3280, 0x4E8CFEBC,
           0xC569F575, 0xCDB2A091, 0x2CC016B4, 0x5C5F4421
        };

        if (salt_count_ <= predef_salt_count)
        {
            std::copy(predef_salt,
                predef_salt + salt_count_,
                std::back_inserter(salt_));

            for (std::size_t i = 0; i < salt_.size(); ++i)
            {
                /*
                   Note:
                   This is done to integrate the user defined random seed,
                   so as to allow for the generation of unique bloom filter
                   instances.
                */
                salt_[i] = salt_[i] * salt_[(i + 3) % salt_.size()] + static_cast<bloom_type>(random_seed_);
            }
        }
        else
        {
            std::copy(predef_salt, predef_salt + predef_salt_count, std::back_inserter(salt_));

            srand(static_cast<unsigned int>(random_seed_));

            while (salt_.size() < salt_count_)
            {
                bloom_type current_salt = static_cast<bloom_type>(rand()) * static_cast<bloom_type>(rand());

                if (0 == current_salt)
                    continue;

                if (salt_.end() == std::find(salt_.begin(), salt_.end(), current_salt))
                {
                    salt_.push_back(current_salt);
                }
            }
        }
    }

    inline bloom_type hash_ap(const unsigned char* begin, std::size_t remaining_length, bloom_type hash) const
    {
        const unsigned char* itr = begin;
        unsigned int loop = 0;

        while (remaining_length >= 8)
        {
            const unsigned int& i1 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);
            const unsigned int& i2 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);

            hash ^= (hash << 7) ^ i1 * (hash >> 3) ^
                (~((hash << 11) + (i2 ^ (hash >> 5))));

            remaining_length -= 8;
        }

        if (remaining_length)
        {
            if (remaining_length >= 4)
            {
                const unsigned int& i = *(reinterpret_cast<const unsigned int*>(itr));

                if (loop & 0x01)
                    hash ^= (hash << 7) ^ i * (hash >> 3);
                else
                    hash ^= (~((hash << 11) + (i ^ (hash >> 5))));

                ++loop;

                remaining_length -= 4;

                itr += sizeof(unsigned int);
            }

            if (remaining_length >= 2)
            {
                const unsigned short& i = *(reinterpret_cast<const unsigned short*>(itr));

                if (loop & 0x01)
                    hash ^= (hash << 7) ^ i * (hash >> 3);
                else
                    hash ^= (~((hash << 11) + (i ^ (hash >> 5))));

                ++loop;

                remaining_length -= 2;

                itr += sizeof(unsigned short);
            }

            if (remaining_length)
            {
                hash += ((*itr) ^ (hash * 0xA5A5A5A5)) + loop;
            }
        }

        return hash;
    }

    std::vector<bloom_type>    salt_;
    std::vector<unsigned char> bit_table_;
    unsigned int               salt_count_;
    unsigned long long int     table_size_;
    unsigned long long int     projected_element_count_;
    unsigned long long int     inserted_element_count_;
    unsigned long long int     random_seed_;
    double                     desired_false_positive_probability_;
};

inline bloom_filter operator & (const bloom_filter& a, const bloom_filter& b)
{
    bloom_filter result = a;
    result &= b;
    return result;
}

inline bloom_filter operator | (const bloom_filter& a, const bloom_filter& b)
{
    bloom_filter result = a;
    result |= b;
    return result;
}

inline bloom_filter operator ^ (const bloom_filter& a, const bloom_filter& b)
{
    bloom_filter result = a;
    result ^= b;
    return result;
}

class compressible_bloom_filter : public bloom_filter
{
public:

    compressible_bloom_filter(const bloom_parameters& p)
        : bloom_filter(p)
    {
        size_list.push_back(table_size_);
    }

    inline unsigned long long int size() const
    {
        return size_list.back();
    }

    inline bool compress(const double& percentage)
    {
        if (
            (percentage < 0.0) ||
            (percentage >= 100.0)
            )
        {
            return false;
        }

        unsigned long long int original_table_size = size_list.back();
        unsigned long long int new_table_size = static_cast<unsigned long long int>((size_list.back() * (1.0 - (percentage / 100.0))));

        new_table_size -= new_table_size % bits_per_char;

        if (
            (bits_per_char > new_table_size) ||
            (new_table_size >= original_table_size)
            )
        {
            return false;
        }

        desired_false_positive_probability_ = effective_fpp();

        const unsigned long long int new_tbl_raw_size = new_table_size / bits_per_char;

        table_type tmp(new_tbl_raw_size);

        std::copy(bit_table_.begin(), bit_table_.begin() + new_tbl_raw_size, tmp.begin());

        typedef table_type::iterator itr_t;

        itr_t itr = bit_table_.begin() + (new_table_size / bits_per_char);
        itr_t end = bit_table_.begin() + (original_table_size / bits_per_char);
        itr_t itr_tmp = tmp.begin();

        while (end != itr)
        {
            *(itr_tmp++) |= (*itr++);
        }

        std::swap(bit_table_, tmp);

        size_list.push_back(new_table_size);

        return true;
    }

private:

    inline void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
    {
        bit_index = hash;

        for (std::size_t i = 0; i < size_list.size(); ++i)
        {
            bit_index %= size_list[i];
        }

        bit = bit_index % bits_per_char;
    }

    std::vector<unsigned long long int> size_list;
};

#endif

using namespace std;


double CaculateSimilarity(const int k, unsigned char* seq0, unsigned char* seq1, size_t len0, size_t len1) {
    //len0 = len0 / 2;
    //len1 = len1 / 2;
    //const int k = 20;                     // 设置子串长度为 30（可根据实际情况调整）
    const uint32_t seed = 0;
    bloom_parameters parameters;         // Bloom Filter 参数
    parameters.projected_element_count = len0 / 2;   // 预计定位元素个数为序列 seq0 的一半
    parameters.false_positive_probability = 0.001;    // 误判率为 1%
    parameters.random_seed = 0xA5A5B6B6;             // 随机数种子
    parameters.compute_optimal_parameters();         // 计算出最优参数
    bloom_filter bloom(parameters);       // 创建 Bloom Filter

    uint32_t tmp_hash;
    //std::vector<uint32_t> hashes0(len0 - k + 1);  // 存储序列 seq0 中所有长度为 k 的子串的哈希值
    for (size_t i = 0; i < len0 - k + 1; ++i) {
        tmp_hash = 0;
        MurmurHash3_x86_32(seq0 + i, k, seed, &tmp_hash); // 使用 MurmurHash 进行哈希计算，并将哈希值保存至 hashes0[i] 中
        bloom.insert(tmp_hash);    // 将哈希值插入 Bloom Filter 中
    }

    int count = 0;
    for (size_t i = 0; i < len1 - k + 1; ++i) {
        tmp_hash = 0;
        MurmurHash3_x86_32(seq1 + i, k, seed, &tmp_hash); // 使用 MurmurHash 进行哈希计算，并将哈希值保存至 hashes0[i] 中
        if (bloom.contains(tmp_hash)) {  // 如果 Bloom Filter 中存在 tmp_hash，则说明存在相同子序列
            count++;
        }
    }

    //double similarity = 2.0 * count / (len0 + len1);
    //std::cout << count << " count\n";
    return count;
}



