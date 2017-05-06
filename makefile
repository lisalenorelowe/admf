code:
	gcc -lm SpectralIntegrator.c AllocateStructs.c utility.c \
	InitializeMetric.c InitializeLapse.c InitializeShift.c ginv.c\
        FunctionDerivs.c  InitializeKij.c IntegrateStep_adm.c\
	fourierfit.c fouriereval.c fourierder.c FunctionFirstDerivs.c\
	GetRicci.c HamCon.c RKstep_adm.c Derivs_admg.c Derivs_admK.c\
	OutputStep_gij.c OutputStep_T.c OutputStep_V.c Derivs_admlapse.c\
        OutputStep_D.c Christoffel.c InitializePhiandEks.c\
        Rxx.c Rzz.c Rxz.c Rxy.c Ryy.c Ryz.c ResetLapse.c ReadRestartDump.c\
	WriteRestartDump.c MomCon.c GetDetg.c Filter.c ConCon.c OutputStep_S.c\
        Esolve.c Laplacian.c Loperator.c Atimes.c dotproduct.c Deltaop.c






