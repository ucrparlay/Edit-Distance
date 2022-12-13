#ifndef hash_block_parallel_hpp
#define hash_block_parallel_hpp

#include <stdio.h>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <math.h>
#include <vector>
#include <string>
#include <stdint.h>
#include <cstdint>
using namespace std;
// const int PRIME_BASE = 139;

// auxiliary function for log
static size_t mylog2(size_t val)
{
    if (val == 1)
        return 0;
    unsigned int ret = 0;
    while (val > 1)
    {
        val >>= 1;
        ret++;
    }
    return ret;
}

// auxiliary function for power x^p
int mypower(int x, int p)
{
    int res = 1;
    for (int i = 0; i < p; i++)
    {
        res *= x;
    }
    return res;
}

// Function to compute the exact hash value of
// sequence, from a, to b (S[a], S[a + 1], ... S[b])
template <typename T>
int hash_value(T seq, int a, int b)
{
    int res = int(seq[a]);
    for (int pos = a + 1; pos <= b; pos++)
    {
        res *= PRIME_BASE;
        res += seq[pos];
    }
    return res;
}

// build a table
template <typename T>
void build(const T seq, vector<vector<int>> &table_seq, size_t block_size)
{
    size_t k = seq.size() / block_size;
    // pre-computed power table [p^(block_size), p^(2 * block_size), ... p^(logk * block_size)]
    vector<int> block_power_table;
    for (int i = 0; i < mylog2(k); i++)
    {
        block_power_table.push_back(mypower(PRIME_BASE, int(block_size * (i + 1))));
    }
    table_seq.resize(mylog2(k) + 1);
    // pre-computed log table [1, 2, 4, ..., logk]
    vector<int> block_log_table;
    for (int i = 0; i < mylog2(k) + 1; i++)
    {
        block_log_table.push_back(1 << i);
    }

    // for the first dim of the table
    table_seq[0].resize(k);
    for (size_t j = 0; j < k; j++)
    {
        table_seq[0][j] = hash_value(seq, j * block_size, (j + 1) * block_size - 1);
    }
    // for the remaining dims of the table
    for (int i = 1; i < mylog2(k) + 1; i++)
    {
        table_seq[i].resize(k - block_log_table[i] + 1);
        // ! can be parallelized
        for (int j = 0; j < k - block_log_table[i] + 1; j++)
        {
            table_seq[i][j] = table_seq[i - 1][j] * block_power_table[i - 1] + table_seq[i - 1][j + block_log_table[i - 1]];
        }
    }
}

// function to build the two tables. `n` is
// the length of the initial sequence size.
template <typename T>
void construct_table(T A, T B, vector<vector<int>> &table_A,
                     vector<vector<int>> &table_B,
                     vector<int> &prefix_a,
                     vector<int> &prefix_b,
                     vector<int> &suffix_a,
                     vector<int> &suffix_b,
                     vector<int> &auxiliary_single_power_table,
                     size_t n)
{
    // logn
    const size_t BLOCK_SIZE = mylog2(n);
    // build the powertable

    build(A, table_A, BLOCK_SIZE);
    build(B, table_B, BLOCK_SIZE);

    // build the prefix and suffix
    // auxiliary power table
    vector<int> block_power_table;
    for (int i = 0; i < BLOCK_SIZE; i++)
    {
        block_power_table.push_back(mypower(PRIME_BASE, i));
    }

    auxiliary_single_power_table.resize(n);
    for (int i = 0; i < n; i++)
    {
        auxiliary_single_power_table[i] = mypower(PRIME_BASE, i);
    }
    // ! can be parallelized
    int a_actual_size = BLOCK_SIZE * int(A.size() / BLOCK_SIZE);
    int b_actual_size = BLOCK_SIZE * int(B.size() / BLOCK_SIZE);
    prefix_a.resize(a_actual_size);
    prefix_b.resize(b_actual_size);
    suffix_a.resize(a_actual_size);
    suffix_b.resize(b_actual_size);
    for (int i = 0; i < a_actual_size; i++)
    {
        prefix_a[i] = 0;
        // ! can be parallelized
        for (int j = int(BLOCK_SIZE) * (i / BLOCK_SIZE); j <= i; j++)
        {
            prefix_a[i] += auxiliary_single_power_table[i - j] * int(A[j]);
        }
        suffix_a[i] = 0;
        // ! can be parallelized
        for (int k = i; k < (i / BLOCK_SIZE + 1) * int(BLOCK_SIZE); k++)
        {
            suffix_a[i] += auxiliary_single_power_table[(i / BLOCK_SIZE + 1) * int(BLOCK_SIZE) - k - 1] * int(A[k]);
        }
    }

    for (int i = 0; i < b_actual_size; i++)
    {
        prefix_b[i] = 0;
        // ! can be parallelized
        for (int j = int(BLOCK_SIZE) * (i / BLOCK_SIZE); j <= i; j++)
        {
            prefix_b[i] += auxiliary_single_power_table[i - j] * int(B[j]);
        }
        suffix_b[i] = 0;
        // ! can be parallelized
        for (int k = i; k < (i / BLOCK_SIZE + 1) * int(BLOCK_SIZE); k++)
        {
            suffix_b[i] += auxiliary_single_power_table[(i / BLOCK_SIZE + 1) * int(BLOCK_SIZE) - k - 1] * int(B[k]);
        }
    }
}

bool compare_lcp(int p, int q, int z,
                 vector<vector<int>> &table_A,
                 vector<vector<int>> &table_B,
                 vector<int> &S_A, vector<int> &P_A,
                 vector<int> &S_B, vector<int> &P_B,
                 vector<int> &aux_power_table)
{
    size_t t = S_A.size() / table_A[0].size(); // block_size
    if (p + t * mypower(2, z) >= table_A[0].size() * t || q + t * mypower(2, z) >= table_B[0].size() * t)
    {
        return false;
    }
    int next_block_A = (p / t) + 1;
    int next_block_B = (q / t) + 1;
    size_t rest_A_size = next_block_A * t - p;
    size_t rest_B_size = next_block_B * t - q;
    int hash_a_v;
    int hash_b_v;
    // A_hash_v = S_A[p] * p^{?} + Table_A[z][next_A] * p^{?} + P_A[nextA * t + rest_B]
    // need pre-computed prefix and suffix

    if (p % t == 0)
    {
        hash_a_v = table_A[z][p / t];
    }
    else if (z != 0)
    {
        hash_a_v = S_A[p] * aux_power_table[(1 << z) * t - rest_A_size] + table_A[z - 1][next_block_A] * aux_power_table[t - rest_A_size] + P_A[p + (1 << z) * t - 1];
    }
    else
    {
        hash_a_v = S_A[p] * aux_power_table[t - rest_A_size] + P_A[p + t - 1];
    }

    if (q % t == 0)
    {
        hash_b_v = table_B[z][q / t];
    }
    else if (z != 0)
    {
        hash_b_v = S_B[q] * aux_power_table[(1 << z) * t - rest_B_size] + table_B[z - 1][next_block_B] * aux_power_table[t - rest_B_size] + P_B[1 + (1 << z) * t - 1];
    }
    else
    {
        hash_b_v = S_B[q] * aux_power_table[t - rest_B_size] + P_B[q + t - 1];
    }
    if (hash_a_v == hash_b_v)
        return true;
    return false;
}

// function for query the lcp from A[p] and B[q]
template <typename T>
int query_lcp(int p, int q,
              T A, T B,
              vector<vector<int>> &table_A,
              vector<vector<int>> &table_B,
              vector<int> &S_A, vector<int> &P_A,
              vector<int> &S_B, vector<int> &P_B,
              vector<int> &aux_power_table)
{
    // find the possible block range point (omit the offset first)
    int x = 0;
    size_t t = S_A.size() / table_A[0].size(); // block_size
    while (p + t * mypower(2, x) < table_A[0].size() * t && (q + t * mypower(2, x) < table_B[0].size() * t))
    {
        if (!compare_lcp(p, q, x, table_A, table_B, S_A, P_A, S_B, P_B, aux_power_table))
        {
            break;
        }
        x++;
    }
    int pp, qq;
    //    x -= 1; // now the x is the minimum LCP,
    // and the real lcp should locate in x and x + 1 block range
    if (x > 0)
    {
        pp = p + t * mypower(2, x - 1);
        qq = q + t * mypower(2, x - 1);
        int y = x - 1;
        while (y >= 0)
        {
            if (compare_lcp(pp, qq, y, table_A, table_B, S_A, P_A, S_B, P_B, aux_power_table))
            {
                pp += t * mypower(2, y);
                qq += t * mypower(2, y);
            }
            y -= 1;
        }
    }
    else
    {
        pp = p;
        qq = q;
    }

    // if it is the beginning of the offset
    if (pp == table_A[0].size() || (qq == table_B[0].size()))
    {
        while (pp < A.size() && qq < B.size())
        {
            if (int(A[pp++]) != int(B[qq++]))
            {
                break;
            }
        }
    }
    else
    { // if not
        while (int(A[pp] == int(B[qq])))
        {
            pp++;
            qq++;
        }
    }
    return (pp - p);
}

#endif