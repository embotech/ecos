#include "com_verizon_jecos_NativeECOS.h"
#define CORE_PACKAGE "com/verizon/jecos/"

#include "ecos.h"

JNIEXPORT jint JNICALL Java_com_verizon_jecos_NativeECOS_EcosSolve(
		JNIEnv *env, jclass this, jint n, jint m, jint p, jint l, jint ncones,
		jintArray q, jdoubleArray Gpr, jintArray Gjc, jintArray Gir,
		jdoubleArray Apr, jintArray Ajc, jintArray Air, jdoubleArray c,
		jdoubleArray h, jdoubleArray b, jdoubleArray result)
{
	extern pwork* ECOS_setup(idxint, idxint, idxint, idxint, idxint, idxint*,
			pfloat *, idxint *, idxint *, pfloat *, idxint *, idxint *,
			pfloat *, pfloat *, pfloat *);

	extern idxint ECOS_solve(pwork *);

	extern void ECOS_cleanup(pwork*, idxint);

	printf("Native EcosSolve\n");

	idxint * qPtr = (idxint *) (*env)->GetLongArrayElements(env, q, NULL);
	pfloat *GprPtr = (*env)->GetDoubleArrayElements(env, Gpr, NULL);
	idxint *GjcPtr = (idxint *) (*env)->GetLongArrayElements(env, Gjc, NULL);
	idxint *GirPtr = (idxint *) (*env)->GetLongArrayElements(env, Gir, NULL);

	pfloat *AprPtr = (*env)->GetDoubleArrayElements(env, Apr, NULL);
	idxint *AjcPtr = (idxint *) (*env)->GetLongArrayElements(env, Ajc, NULL);
	idxint *AirPtr = (idxint *) (*env)->GetLongArrayElements(env, Air, NULL);

	pfloat *cPtr = (*env)->GetDoubleArrayElements(env, c, NULL);
	pfloat *hPtr = (*env)->GetDoubleArrayElements(env, h, NULL);
	pfloat *bPtr = (*env)->GetDoubleArrayElements(env, b, NULL);

	idxint exitflag = ECOS_FATAL;

#if PROFILING > 0
	double ttotal, tsolve, tsetup;
#endif
#if PROFILING > 1
	double torder, tkktcreate, ttranspose, tfactor, tkktsolve;
#endif

	/* set up data */
	pwork * wspace = ECOS_setup(n, m, p, l, ncones, qPtr, GprPtr, GjcPtr,
			GirPtr, AprPtr, AjcPtr, AirPtr, cPtr, hPtr, bPtr);

	if (wspace != NULL) {
		/* solve */
		exitflag = ECOS_solve(wspace);

		/* some statistics in milliseconds */
#if PROFILING > 0
		tsolve = wspace->info->tsolve * 1000;
		tsetup = wspace->info->tsetup * 1000;
		ttotal = tsetup + tsolve;
#endif
#if PROFILING > 1
		torder = wspace->info->torder * 1000;
		tkktcreate = wspace->info->tkktcreate * 1000;
		ttranspose = wspace->info->ttranspose * 1000;
		tfactor = wspace->info->tfactor * 1000;
		tkktsolve = wspace->info->tkktsolve * 1000;
#endif

#if PRINTLEVEL > 2
#if PROFILING > 1
		printf("ECOS timings (all times in milliseconds):\n\n");
		printf("1. Setup: %7.3f (%4.1f%%)\n", tsetup, tsetup / ttotal*100);
		printf("2. Solve: %7.3f (%4.1f%%)\n", tsolve, tsolve / ttotal*100);
		printf("----------------------------------\n");
		printf(" Total solve time: %7.3f ms\n\n", ttotal);

		printf("Detailed timings in SETUP:\n");
		printf("Create transposes: %7.3f (%4.1f%%)\n", ttranspose, ttranspose / tsetup*100);
		printf("Create KKT Matrix: %7.3f (%4.1f%%)\n", tkktcreate, tkktcreate / tsetup*100);
		printf(" Compute ordering: %7.3f (%4.1f%%)\n", torder, torder / tsetup*100);
		printf("            Other: %7.3f (%4.1f%%)\n", tsetup-torder-tkktcreate-ttranspose, (tsetup-torder-tkktcreate-ttranspose) / tsetup*100);
		printf("\n");

		printf("Detailed timings in SOLVE:\n");
		printf("   Factorizations: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tfactor, tfactor / tsolve*100, tfactor / ttotal*100);
		printf("       KKT solves: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tkktsolve, tkktsolve / tsolve*100, tfactor / ttotal*100);
		printf("            Other: %7.3f (%4.1f%% of tsolve / %4.1f%% of ttotal)\n", tsolve-tkktsolve-tfactor, (tsolve-tkktsolve-tfactor) / tsolve*100, (tsolve-tkktsolve-tfactor) / ttotal*100);
		printf("\n");
#endif
#endif
		/* Extract the solution */
		double *xPtr = (*env)->GetDoubleArrayElements(env, result, NULL);
		int i = 0;
		while(i < n) {
			xPtr[i] = wspace->x[i];
			i++;
		}
		(*env)->ReleaseDoubleArrayElements(env, result, xPtr, 0);

		/* clean up memory */
		ECOS_cleanup(wspace, 0);
	}

	return exitflag;
}

JNIEXPORT jstring JNICALL Java_com_verizon_jecos_NativeECOS_EcosVer(JNIEnv *env, jclass this)
{
	extern const char * ECOS_ver();
	jstring ecosVer = (*env)->NewStringUTF(env, ECOS_ver());
	return ecosVer;
}
