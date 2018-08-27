import SBCcode.MonteCarlo.PICOcalGlobalLikelihood as pcgl
import numpy as np

testparams = np.zeros(pcgl.n_params)
testparams[:pcgl.n_Epts] = pcgl.Epts_initialguess.ravel()

ll = pcgl.PICOcalLL(testparams)
print(ll)

residuals = pcgl.PICOcalResiduals(testparams)
