/*
 * normal_add.c
 * Projective座標系での加算を行う(P != Q)
 *
 * 更新履歴
 * 2014/10/17 新規作成
 * 2014/10/19 バグ修正
 * 2014/11/01 各宣言をマクロに修正
 * 2014/11/02 引数の一部をconstに変更
 *
 */

#include <gmp.h>
#include "point.h"

/*
 * 異なる点の加算を行う関数
 * PROJECTIVE_POINT R:P + Qの結果を格納する点
 * PROJECTIVE_POINT P:演算される点
 * PROJECTIVE_POINT Q:演算される点
 * const mpz_t N           :mod N
 */
void normal_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, PROJECTIVE_POINT Q, const mpz_t N)
{
	mpz_t u;
	mpz_t v;
	mpz_t A;
	PROJECTIVE_POINT tP;
	PROJECTIVE_POINT tQ;

	/*一時変数*/
	mpz_t tmp;

	/*初期化*/
	mpz_init(u);
	mpz_init(v);
	mpz_init(A);
	mpz_init(tmp);
	projective_point_init(tP);
	projective_point_init(tQ);

	projective_point_set(tP, P);
	projective_point_set(tQ, Q);

	/*uの計算*/
	mpz_mul(u, tQ->Y, tP->Z);
	mpz_mul(tmp, tP->Y, tQ->Z);
	mpz_mod(tmp, tmp, N);
	mpz_sub(u, u, tmp);
	mpz_mod(u, u, N);

	/*vの計算*/
	mpz_mul(v, tQ->X, tP->Z);
	mpz_mul(tmp, tP->X, tQ->Z);
	mpz_sub(v, v, tmp);
	mpz_mod(v, v, N);

	/*Aの計算*/
	mpz_pow_ui(A, u, 2);
	mpz_mul(A, A, tP->Z);
	mpz_mul(A, A, tQ->Z);
	mpz_pow_ui(tmp, v, 3);
	mpz_mod(tmp, tmp, N);
	mpz_sub(A, A, tmp);
	mpz_pow_ui(tmp, v, 2);
	mpz_mul_ui(tmp, tmp, 2);
	mpz_mul(tmp, tmp, tP->X);
	mpz_mul(tmp, tmp, tQ->Z);
	mpz_sub(A, A, tmp);
	mpz_mod(A, A, N);

	/*X3の計算*/
	//X3=vA
	mpz_mul(R->X, v, A);
	mpz_mod(R->X, R->X, N);

	/*Y3の計算*/
	mpz_pow_ui(tmp,v,2);
	mpz_mul(tmp,tmp,tP->X);
	mpz_mul(tmp,tmp,tQ->Z);
	mpz_sub(tmp,tmp,A);
	mpz_mul(tmp,u,tmp);
	mpz_set(R->Y,tmp);
	mpz_pow_ui(tmp,v,3);
	mpz_mul(tmp,tP->Y,tmp);
	mpz_mul(tmp,tmp,tQ->Z);
	mpz_sub(R->Y,R->Y,tmp);

	mpz_mod(R->Y,R->Y,N);

	/*Z3の計算*/
	mpz_pow_ui(R->Z, v, 3);
	mpz_mul(R->Z, R->Z, tP->Z);
	mpz_mul(R->Z, R->Z, tQ->Z);
	mpz_mod(R->Z, R->Z, N);

	/*使用した変数の開放*/
	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(A);
	mpz_clear(tmp);
	projective_point_clear(tP);
	projective_point_clear(tQ);
}
