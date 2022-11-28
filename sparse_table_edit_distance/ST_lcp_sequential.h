#ifndef lcp_h
#define lcp_h

#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#define ull unsigned long long
using namespace std;


const int PRIME_MOD = 100000;
const int PRIME_BASE = 128;

void build_hash_table(const char* s1, const char* s2,
                      vector<vector<int>>& table_s1, vector<vector<int>>& table_s2,
                      vector<int>& powerN1, vector<int>& powerN2,
                      vector<int>& logN1, vector<int>& logN2) {
    // build the first leaves layer
    ull table1_d2 = strlen(s1);
    ull table2_d2 = strlen(s2);
    ull table1_d1 = 0; int upp1 = 1;
    ull table2_d1 = 0; int upp2 = 1;
    while (upp1 <= table1_d2) {
        table1_d1 ++;
        upp1 *= 2;
    }
    
    while (upp2 < table2_d2) {
        table2_d1 ++;
        upp2 *= 2;
    }
    
    //initialize size
    table_s1.resize(table1_d1); table_s2.resize(table2_d1);
    for (int i = 0; i < table1_d1; i ++) table_s1[i].resize(table1_d2);
    for (int i = 0; i < table2_d1; i ++) table_s2[i].resize(table2_d2);
    powerN1.resize(table1_d1);
    powerN2.resize(table2_d1);
    logN1.resize(table1_d1);
    logN2.resize(table2_d1);
    for (int i = 0; i < table1_d2; i ++){
        table_s1[0][i] = int(s1[i]);
    }
    for (int i = 0; i < table2_d2; i ++) {
        table_s2[0][i] = int(s2[i]);
    }
    int len = 1;
    powerN1[0] = 1;
    for (int i = 0; i < table1_d1; i ++) {
        logN1[i] = len;
        len *= 2;
        powerN1[i + 1] = (PRIME_BASE << i) % PRIME_MOD;
    }
    len = 1;
    powerN2[0] = 1;
    for (int i = 0; i < table2_d1; i ++) {
        logN2[i] = len;
        len *= 2;
        powerN2[i + 1] = (PRIME_BASE << i) % PRIME_MOD;
    }
    

    
    // build the second to k-th layer
    for (int i = 1; i < table1_d1; i ++) {
        for (int j = 0; j < table1_d2 - logN1[i] + 1; j ++) {
            for (int pw = 0; pw < logN1[i - 1]; pw ++) {
                table_s1[i][j] = table_s1[i - 1][j] * PRIME_BASE % PRIME_MOD;
            }
            table_s1[i][j] += table_s1[i - 1][j + logN1[i - 1]];
        }
    }
    for (int i = 1; i < table2_d1; i ++) {
        for (int j = 0; j < table2_d2 - logN2[i] + 1; j ++) {
            for (int pw = 0; pw < logN2[i - 1]; pw ++) {
                table_s2[i][j] = table_s2[i - 1][j] * PRIME_BASE % PRIME_MOD;
            }
            table_s2[i][j] += table_s2[i - 1][j + logN2[i - 1]];
        }
    }
    
}


///**
//    Input: spare table, range [i,j]
//    Output: hash value in the given range
// */
//int query_hash(vector<vector<int>>& table, int l, int r, vector<int>& pw_table) {
//    if(l == r) {
//        return table[0][l];
//    }
//    // the gap should be r - l for each element
//    // in the sparse table
//    int v = 0;
//    for (int j = table.size(); j >= 0; j --) {
//        if (l + (1 << j) <= r) {
//            int gap = r - l - (1 << j);
//            v = v + table[j][l] * pw_table[gap] % PRIME_MOD;;
//            l += 1 << j;
//        }
//    }
//    return v;
//}


/**
    Input: Two sparse tables for the two sequences
    Output: LCP starting at i and j respectively
 */
int query_lcp(vector<vector<int>>& table1, vector<vector<int>>& table2,
              vector<int>& logN1, vector<int>& logN2,
              int i, int j) {
    int v1;
    int v2;
    // possible value is from 0 to the smaller one of the remaining sequences
    int pos_i = i;
    int pos_j = j;
    int l = 0;
    int r = min(table1[0].size() - i, table2[0].size() - j);
    // search from the power first
    int mid_try = 1;
    int power_id = 0;
    int power_try;
    while (l != r) {
        for (int i = 0; i < logN1.size(); i ++) {
            if (l + logN1[i] <= r) {
                power_id = i;
                power_try = logN1[i];
                mid_try = l + power_try;
            }
        }
        // test the hash value (whether equal)
        v1 = table1[power_id][pos_i];
        v2 = table2[power_id][pos_j];
        if (v1 == v2) {
            pos_i += power_try;
            pos_j += power_try;
            l = mid_try;
        } else {
            r = mid_try - 1;
        }
    }
    return l;
}


#endif /* lcp_h */
