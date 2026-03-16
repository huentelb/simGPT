* Author: BH
* Date: 12.03.2026
* AME of cluster membership by data source

use "N:\durable\Data20\Project_JW_FamWealth\simgpt\gp_bench.dta", clear

* store cluster labels in locals
global cl1 "C1 - Late 3-gen"
global cl2 "C2 - 4-gen"								
global cl3 "C3 - 2-gen"								
global cl4 "C4 - Non-parent late"		
global cl5 "C5 - Early 3-gen"							
global cl6 "C6 - Non-parent early"	


mlogit chi i.source
est store reg
	margins, dydx(source) post


foreach n of numlist 1/6 {
	est restore reg

	margins, dydx(source) predict(outcome(`n')) post
	est store m_`n'
}

set scheme white_tableau

graph set window fontface "Arial"
	grstyle init
	grstyle set size 8pt: key_label
	grstyle set symbolsize 4pt 
	*grstyle set color #bdbdbd #252525
	*grstyle set color #bdbdbd #252525:  p#markline
	*grstyle set color #bdbdbd #252525:  p#marklab
	grstyle set symbol O


coefplot (	m_1, asequation($cl1) \ ///
			m_2, asequation($cl2) \ ///
			m_3, asequation($cl3) \ ///
			m_4, asequation($cl4) \ ///
			m_5, asequation($cl5) \ ///
			m_6, asequation($cl6) \ ///
			), drop(_cons) 									///
			xline(0, lpattern(solid) lcolor("135 31 36")) 	///
			ciopts(recast(rcap) color(black)) mcolor(black)	///
			ylab("") yscale(noline) xscale(noline) ///
			name(ame, replace)  
			
	gr export "N:\durable\Data20\Project_JW_FamWealth\simgpt\sim_results_260311\benchmark_joint\ame.png", replace	
	gr export "N:\durable\Data20\Project_JW_FamWealth\simgpt\sim_results_260311\benchmark_joint\ame.pdf", replace		
	
