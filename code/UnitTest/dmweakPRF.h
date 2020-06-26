//Header file weak23prf.h

#ifndef UNIT_TEST_H
#define UNIT_TEST_H


//generating input
void generate_rand_input(uint64_t input[4], std::mt19937 &generator);
void generate_rand_key(uint64_t key[4][256], std::mt19937 &generator);
void generate_rand_matrix(uint64_t randMat1[2][256], uint64_t randMat2[2][256], uint64_t randMatZ3[81][256], std::mt19937 &generator);

//phase 1 testing
void compute(uint64_t key[4][256], uint64_t input[4], uint64_t z_final[4]);
void unpackOutput(uint64_t output[4], char p2output[256]);

//phase 3-method 1 testing
void multMod3(uint64_t outM[2], uint64_t outL[2], uint64_t msbs[2][256], uint64_t lsbs[2][256], uint64_t in[4]);

//PRF alternate method testing
void InnerProdMul(uint64_t outVec[84], uint64_t randMatZ3[81][256], uint64_t in[4]);
void InnerProdMul2(uint64_t outVec[84], uint64_t randMatZ3[81][256], uint64_t in[4]);

#endif