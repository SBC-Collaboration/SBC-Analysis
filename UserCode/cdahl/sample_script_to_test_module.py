import SBCcode as sbc
import SBCcode.AnalysisModules.AnalyzeDytran as ad


ev = sbc.get_event('/Users/cdahl/Desktop', 24)

out = ad.dytranAnalysis(ev)
