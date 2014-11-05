/*
 * ecm.c
 * 楕円曲線法によって因数分解を行う
 *
 * 更新履歴
 * 2014/10/29 新規作成
 * 2014/10/30 mpz_t f追加
 *            gcd処理追加
 * 2014/11/02 引数の一部をconstに変更
 * testtesttest
 */

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "point.h"

/* logメモ */
/* logab= logb/loga */

/*
 * 楕円曲線法にて因数分解を行う関数
 * mpz_t f              :発見した素因数
 * const mpz_t N              :素因数分解したい合成数
 * const unsigned long int A  :楕円曲線の係数a
 * const unsigned long int k  :スカラー倍
 */
void ecm(mpz_t f, const mpz_t N, const unsigned long int A, const unsigned long int k)
{
	/* 使用変数・構造体の宣言 */
	AFFINE_POINT aP;
	PROJECTIVE_POINT pP;
	int e;
	int i;

	/* Pの初期化 */
	affine_point_init(aP);
	projective_point_init(pP);

	/* Pの点の座標を指定 */
	mpz_set_ui(aP->x, 1);
	mpz_set_ui(aP->y, 1);

	/* aの決定 */
	mpz_t a;
	mpz_init(a);
	mpz_set_ui(a, A);

	/* 素数の決定 */
	unsigned long int p = 2;
	mpz_t prime;
	mpz_init(prime);
	mpz_set_ui(prime, p);

	/* Affine -> Projective 変換 */
	afftopro(pP, aP, N);

	/* 内部計算 */
	while (p <= k) {
		/* e = log p kを決める */
		e = (int)(log(k) / log(p));
		for (i = 1; i <= e; i++) {
			scalar(pP, pP, p, a, N);
			mpz_gcd(f, pP->Z, N);
			gmp_printf("gcd(%Zd, %Zd) = %Zd\n", pP->Z, N, f);
			if ( mpz_cmp_ui(f, 1) != 0) {
				goto FOUND;
			}
		}
		/* pを次の素数に */
		mpz_nextprime(prime, prime);
		p = mpz_get_ui(prime);
	}
FOUND:
	printf("Stage1: A = %ld\n", A);

	/* 使用変数・関数の開放 */
	affine_point_clear(aP);
	projective_point_clear(pP);
	mpz_clear(a);
	mpz_clear(prime);
}
