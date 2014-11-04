/*
 * double_add.c
 * Projective座標系の2倍算を行う
 *
 * R=P+Pをするメソッド
 * P:とある点
 * R:点Pを二倍算した点
 * a:楕円曲線の係数
 *
 * 更新履歴
 * 2014/10/17 新規作成
 * 2014/10/19 一部バグ修正
 *            mod処理を追加
 * 2014/10/30 コメント追加
 * 2014/11/01 宣言をマクロに変更
 * 2014/11/02 引数の一部をconstに変更
 *            メモリリーク修正
 */
#include <gmp.h>
#include "point.h"

/*
 * Projective座標系の2倍算を行う
 * PROJECTIVE_POINT R :Pの2倍算の結果を格納する点
 * PROJECTIVE_POINT P :2倍算を行う点
 * const mpz_t a            :楕円曲線の係数
 * const mpz_t N            :mod N
 */
void double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t a, const mpz_t N)
{
	mpz_t w;
	mpz_t s;
	mpz_t B;
	mpz_t h;

	PROJECTIVE_POINT tP;

	/* 一時変数 */
	mpz_t tmp;
	mpz_t tmp2;

	/* 初期化 */
	mpz_init(w);
	mpz_init(s);
	mpz_init(B);
	mpz_init(h);
	mpz_init(tmp);
	mpz_init(tmp2);
	projective_point_init(tP);

	projective_point_set(tP, P);

	/* wの計算 */
	mpz_pow_ui(w, tP->Z, 2);
	mpz_mul(w, w, a);
	mpz_pow_ui(tmp, tP->X, 2);
	mpz_mul_ui(tmp, tmp, 3);
	mpz_add(w, w, tmp);
	mpz_mod(w, w, N);

	/* sの計算 */
	mpz_mul(s, tP->Y, tP->Z);
	mpz_mod(s, s, N);

	/* Bの計算 */
	mpz_mul(B, tP->X, tP->Y);
	mpz_mul(B, B, s);
	mpz_mod(B, B, N);

	/* hの計算 */
	mpz_pow_ui(h, w, 2);
	mpz_mul_ui(tmp, B, 8);
	mpz_sub(h, h, tmp);
	mpz_mod(h, h, N);

	/* Xrの計算 */
	mpz_mul_ui(R->X, h, 2);
	mpz_mul(R->X, R->X, s);
	mpz_mod(R->X, R->X, N);

	/* Yrの計算 */
	mpz_mul_ui(R->Y, B, 4);
	mpz_sub(R->Y, R->Y, h);
	mpz_mul(R->Y, R->Y, w);
	mpz_pow_ui(tmp, tP->Y, 2);
	mpz_pow_ui(tmp2, s, 2);
	mpz_mul(tmp, tmp, tmp2);
	mpz_mul_ui(tmp, tmp, 8);
	mpz_sub(R->Y, R->Y, tmp);
	mpz_mod(R->Y, R->Y, N);

	/* Zrの計算 */
	mpz_pow_ui(R->Z, s, 3);
	mpz_mul_ui(R->Z, R->Z, 8);
	mpz_mod(R->Z, R->Z, N);

	/* 使用した変数の開放 */
	mpz_clear(w);
	mpz_clear(s);
	mpz_clear(B);
	mpz_clear(h);
	mpz_clear(tmp);
	mpz_clear(tmp2);

	projective_point_clear(tP);
}
