/*
 * main.c
 *
 * 更新履歴
 * 2014/10/30 新規作成
 * 2014/11/01 誤字修正(null, unsinedなど)
 *            define修正
 */

#include <stdio.h>
#include <gmp.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "point.h"

#define A_LOOP 10000

/* !	適当です	!
 * 素因数が見つからなかった    : 0
 * エラー終了                  : 1
 * 因数が素数で余因数が素数    : 2
 * 因数が素数で余因数が合成数  : 4
 * 因数が合成数で余因数が素数  : 8
 * 因数が合成数で余因数が合成数: 16
 */
int main (int argc, char *argv[])
{
	if (argc <= 1) {
		fprintf (stderr, "Error: Need two Argument\n");
		fprintf (stderr, "Usage: funecm [options] [k]\n");
		return 1;
	}

	/* オプション処理 */
	int opt;
	int loop = 0;
	while ((opt = getopt (argc, argv, "hl")) != -1) {
		switch (opt) {
			case 'h':
				fprintf(stdout, "Usage: funecm [options] [composite number] [k]\n");
				fprintf(stdout, "-h: help\n");
				return 0;
				break;
			case 'l':
				loop = 1;
				break;
			default:
				fprintf(stderr, "No such option\n");
				fprintf(stdout, "Usage: funecm [options] [composite number] [k]\n");
				return 1;
		}
	}

	/* ARGUMENT CONVERSION */
	mpz_t N;
	mpz_init_set_str(N, argv[optind++], 10);
	unsigned long int k;
	k = (unsigned long int)strtol(argv[optind++], NULL, 10);

	/* 修正予定 */
	if (k <= 2)
		return 0;

	switch (mpz_probab_prime_p (N, 25)) {
		case 2:
			gmp_printf("%Zd is definitely prime\n", N);
			return 0;
			break;
		case 1:
			gmp_printf("%Zd is probably prime\n", N);
			return 0;
			break;
		case 0:
			gmp_printf("%Zd is definitely composite\n", N);
			break;
		default:
			break;
	}

	AFFINE_POINT P;
	affine_point_init(P);

	mpz_t factor;
	mpz_t cofactor;

	mpz_init(factor);
	mpz_init(cofactor);

	char digits[1000];

	gmp_printf("Input number: %Zd  ", N);
	mpz_get_str(digits, 10, N);
	printf("digits: %d\n", strlen(digits));
	gmp_printf("k: %ld\n", k);

	clock_t A_start;
	clock_t total_start;
	clock_t A_end;
	clock_t total_end;

RESTART:

	total_start = clock();
	unsigned long int A;
	for (A = 1; A < A_LOOP; A++) {
		A_start = clock();

		ecm(factor, N, A, k);
		mpz_divexact(cofactor, N, factor);
		/* 因数が1又はNだった場合係数を変えてやり直す */
		if (mpz_cmp_ui(factor, 1) == 0 || mpz_cmp(factor, N) == 0) {
			A_end = clock();
			printf("stage1 time: %.3f seconds\n", (double)(A_end - A_start) / CLOCKS_PER_SEC);
			printf("factor not found\n");
			printf("--------------------------------------------------\n");
			continue;
		}
		mpz_get_str(digits, 10, factor);
		A_end = clock();
		total_end = clock();
		printf("stage1 time: %.3f seconds\n", (double)(A_end - A_start) / CLOCKS_PER_SEC);
		printf("total: %.3f seconds\n", (double)(total_end - total_start) / CLOCKS_PER_SEC);
		/* 終了ステータス */
		switch (mpz_probab_prime_p(factor, 25)) {
			case 2:
				gmp_printf("definite prime factor found: %Zd  ", factor);
				printf("digits: %d\n", strlen(digits));
				gmp_printf("cofactor: %Zd\n", cofactor);
				goto END;
			case 1:
				gmp_printf("probable prime factor found: %Zd  ", factor);
				printf("digits: %d\n", strlen(digits));
				gmp_printf("cofactor: %Zd\n", cofactor);
				goto END;
			case 0:
				gmp_printf("composite factor found: %Zd  ", factor);
				printf("digits: %d\n", strlen(digits));
				gmp_printf("cofactor: %Zd\n", cofactor);
				goto END;
			default:
				goto END;
		}
	}
END:
	if ((loop == 1) && mpz_probab_prime_p(cofactor, 25) == 0) {
		mpz_set(N, cofactor);
		printf("-------------------------RESTART-------------------------\n");
		goto RESTART;
	}
	/* メモリの解放*/
	affine_point_clear(P);
	mpz_clear(factor);
	mpz_clear(cofactor);

	return 0;
}
