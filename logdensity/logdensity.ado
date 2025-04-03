program define logdensity, eclass

	syntax varname(max=1) [if] [in], x(namelist max = 1) h(real) [DEgree(integer 1) Kernel(namelist max = 1) MInx(string) MAxx(string) g(string) dg(string)]
	marksample touse 
	
	*Drop missing data
	qui count if `touse' & missing(`varlist')
	local missing =`r(N)'
	qui count if `touse'
	local nonmissing = `r(N)' - `missing'
	
	*Default kernel is epanechnikov
	if ("`kernel'" == ""){
		local kernel "epanechnikov"
		di as txt `"Default epanechnikov kernel used"'
	}

	*Check if kernel is permitted
	local kernel_list = `" "gaussian", "epanechnikov", "triangle", "uniform", "quartic", "triweight", "cosinus" "'
	capt assert inlist( "`kernel'", `kernel_list')
	if _rc{
		di as err `"Error: Valid args to option kernel(): `kernel_list' "'
		exit
	} 

	*Check if degree is positive
	if (`degree' <= 0){
		di as err `"Error: degree() has to be a positive integer"'
		exit
	}
		
	*Check if h is positive
	if (`h'<=0){
		di as err `"Error: h() has to be positive"'
		exit
	}

	*Check if g or dg is missing
	if (("`g'" != "" & "`dg'" == "")|("`g'" == "" & "`dg'" != "")){
		di as err `"Error: You cannot have either g or dg missing"'
		exit
	}
	
	*Force non-default function to be named g2 and dg2
	if (("`g'" != "g2" & "`g'" != "")|("`dg'"!="dg2" & "`dg'" != "")){
		di as err `"Error: Please rename your functions as g2 and dg2 respectively"'
		exit
	}


	if (("`g'" == "" & "`dg'" == "")){
		di as txt`"Default g and dg functions used"'
		local g "mf_g"
		local dg "mf_dg"
	}

	else{
		di as txt`"Custom g2 and dg2 functions used"'
		local g "mf_g2"
		local dg "mf_dg2"
	}
	
	capt confirm variable `x'
	if !_rc{
	    mkmat `x', matrix(tmp)
		matrix tmp = tmp'
		local x tmp
	}
	
	*Main function
	local maxx = real("`maxx'")
	local minx = real("`minx'")
	mata:PS_est("`varlist'", "`touse'", "`kernel'","`x'",`degree',`h',`minx',`maxx',&`g'(),&`dg'())
	display("Number of non-missing observations: `nonmissing'")
	display("Order of Polynomial Approximation: `degree'")
	display("Bandwidth: `h'")

	numlist "-1/`degree'"
	matrix colnames L = `r(numlist)'
	matrix coleq L = :`varlist' :Log-density Derivative:
	ereturn scalar N = `nonmissing'
	ereturn scalar h = `h'
	ereturn matrix L = L
	matrix list e(L), noheader
end
/*==========================================================================================================================================================================*/

/*********************************************/
/* Mata Stuff                                */
/*********************************************/

mata
/*=============================*/
/*Setting up m_kernel function */
/*=============================*/
function kernel_f(u,kernel){
	
	if (kernel=="gaussian"){
		return(1:/sqrt(2*pi()):*exp(-(u:^2):/2))	
	}
	
	else if (kernel=="triweight"){
		return((35/32):*((1:-u:^2):^3):*(abs(u):<=1))
	}
	
	else if (kernel=="uniform"){
		return(1/2:*(abs(u):<=1))
	}
		
	else if (kernel=="triangle"){
		return((1:-abs(u)):*(abs(u):<=1))
	}
	
	else if (kernel=="cosinus"){
		return(((pi()/4) :* cos(pi()/2:*u)):*(abs(u):<=1))
	}
	
	else if (kernel=="quartic"){
		return((15/16:*((1:-u:^2):^2)):*(abs(u):<=1))
	}
	
	else{ /*epanechnikov*/
		return(0.75:*(1:-u:^2):*(abs(u):<=1))
	}
}

/*=============================*/
/*Setting up g,dg function     */
/*=============================*/
function mf_g(u,zl,zr,S)
{
	power = (1::S)'
	output = (((colnonmissing(power))'*(u:+zl)'):^power')':*((colnonmissing(power))'*(zr:-u)')'
	return(output)
}

function mf_dg(u,zl,zr,S)
{
	power = (1::S)'
	output = (((colnonmissing(power:-1))'*(u:+zl)'):^(power:-1)')':*((colnonmissing(power))'*(zr:-u)')':*power - (((colnonmissing(power))'*(u:+zl)'):^power')'
	return(output)
}

function mf_g2(u,zl,zr,S) return(g2(u,zl,zr,S))
function mf_dg2(u,zl,zr,S) return(dg2(u,zl,zr,S))

/*==================================================================*/
/*Setting up the term inside the integral in equation 3 of PS(2020) */
/*==================================================================*/
real rowvector expsum(real rowvector u, real rowvector store,string scalar kernel){
	h = store[1]
	S = store[2]
	beta = store[3..length(store)]'
	
	j = (1..S)
	u2 = (u'*colnonmissing(j)):^j
	h2 = (h:^j)
	huj = u2:*h2:/factorial(j)
	output = exp(huj*beta):*kernel_f(u,kernel)'
	return(output')
}

/*=============================*/
/*Carry out integration.       */
/*=============================*/
real rowvector f_integrate(real rowvector store,string scalar kernel,real scalar zl, real scalar zr){
	f = integrate(&expsum(),-1*zl,zr,60,store,kernel)
	return(f)
}

/*=============================*/
/*Implemennts the estimator    */
/*=============================*/
void PS_est(string scalar varname, string scalar touse, string scalar kernel, string scalar points, real scalar deriv,real scalar h, real scalar minx, real scalar maxx, pointer(real scalar function) scalar f, pointer(real scalar function) scalar df)
{
	vecx = st_matrix(points)
	vecx = select(vecx, colmissing(vecx):==0)
	/*=====================*/	
	/*Check if evaluation point is within support*/
	if (minx==.) minx = min(vecx):-h
	if (maxx==.) maxx = max(vecx):+h
	
	if ((rowsum(vecx:<minx)>0)|(rowsum(vecx:>maxx)>0)){
		_error("Evaluation points must be within support [minx,maxx]")
	}
	/*=====================*/
	
	cols = length(vecx)
	rows = deriv + 1
	finalresults = J(rows,cols,0)
	
	for(i=1;i<=length(vecx);i++){
		/*Setting up the data and the relevant parameters*/
		real colvector data
		st_view(data,.,varname,touse)
		n = length(data)
		S = deriv
		data = data :- minx
		x = vecx[,i] :- minx
		maxx = maxx :- minx
		zl = min((x/h,1))
		zr = min(((maxx-x)/h,1))
		u = (data:-x):/h
		u = select(u, u:<=zr :& u:>=-zl)
		gu = (*f)(u,zl,zr,S)
		dgu = (*df)(u,zl,zr,S)
	
		/*=====================*/
		/*quick check*/
		check_zl = (*f)(-zl,zl,zr,S)
		check_zr = (*f)(zr,zl,zr,S)
		if ((check_zl!=J(1,S,0))|(check_zr!=J(1,S,0))){
			_error("Please make sure both g(-zl) and g(zr) are zero vectors")
		}

		/*=====================*/
	
		store = J(S,S,0)
		for(s=1;s<=S;s++){
			R = ((u:*h):^(s:-1)):/factorial(s-1)
			R = mean(gu:*R):/h
			store[s,.] = -R
		}

		/*Backout beta*/
		dgu = mean(dgu):/(h:^2)
		beta = lusolve(store',dgu')

		/*Rosenblatt Parzen Kernel Estimator*/
		fhat_m = colsum(kernel_f(u,kernel)):/(n:*h)
		logfhat = log(fhat_m/f_integrate((h,S,beta'),kernel,zl,zr))
		ans = (logfhat, beta')'
		finalresults[.,i] = ans
	}
	st_matrix("L",(vecx', finalresults'))
}


end
