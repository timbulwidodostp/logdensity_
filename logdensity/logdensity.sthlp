{smcl}
{* *! version 0.1.0  06jun2020}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install integrate" "ssc install"}{...}
{vieweralsosee "Help integrate (if installed)" "help integrate"}{...}
{viewerjumpto "Syntax" "logdensity##syntax"}{...}
{viewerjumpto "Description" "logdensity##description"}{...}
{viewerjumpto "Options" "logdensity##options"}{...}
{viewerjumpto "Remarks" "logdensity##remarks"}{...}
{viewerjumpto "Examples" "logdensity##examples"}{...}
{viewerjumpto "Stored results" "logdensity##storedresults"}{...}
{viewerjumpto "Contact" "logdensity##contact"}{...}
{viewerjumpto "References" "logdensity##reference"}{...}
{title:Title}

{phang}
{bf:logdensity} {hline 2} Univariate density estimation


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:logdensity}
[{varname}]
{ifin}
{cmd:, x(namelist) h(#)} [{it: more options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt x(namelist)}}points at which to estimate the logarithm of the 
density and its derivative(s); can be a variable name or a row vector; missing
entries are silently dropped{p_end}
{synopt:{opt h(#)}}bandwidth to be used{p_end}
{synopt:{opt k:ernel(namelist)}}kernel function; default is 
{it:epanechnikov}{p_end}
{synopt:{opt de:gree(#)}}order of the local polynomial approximation;
        default is {it:1} (local linear){p_end}
{synopt:{opt mi:nx(#)}}left end of support of {it:varname}; default is minimum 
value of {it:varname} - {it:h}, effectively negative infinity{p_end}
{synopt:{opt ma:xx(#)}}right end of support of {it:varnae}; default is maximum 
value of {it:varname} + {it:h}, effectively infinity{p_end}
{synopt:{opt g(string)}}function used to estimate derivative(s) of the 
log-density; see below for description of default behavior{p_end}
{synopt:{opt dg(string)}}derivative of {it:g}; see below for description of 
default behavior{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:logdensity} estimates the (natural) logarithm of the density of the 
variable in {varname} using a local polynomial approximation. This command is an
implementation of the estimator developed in Pinkse, J. and Schurter, K. (2020)
"Estimates of derivatives of (log) densities and related objects"
{browse "https://arxiv.org/abs/2006.01328":available on arXiv}. An advantage of 
this approach relative to alternative kernel-based estimators for the density is
that the estimated density, defined as exp() of the log-density, is continuous
and nonnegative whenever the kernel function is nonnegative and continuous, even
when the density is estimated near the boundary of the support of the data. See
the paper for details.

{pstd}
Estimates are stored in the matrix e(L). The number of nonmissing observations
used to estimate the log-density and its derivative(s) is returned in the scalar
e(N) and the bandwidth is stored in the scalar e(h).

{pstd}
This command uses {cmd: integrate} to perform numerical integration. Please run
{stata ssc install integrate, replace} to install the latest version of
{cmd: integrate}. You may also need to run {stata integrate, installmata} to 
finish the installation (and make sure the directory in which the
command attempts to save integrate.mlib exists and is writeable).


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt x(namelist)} is a required {it:namelist} of values of {it:varname} at which
estimates of the logarithm of the density and its derivative(s) will be 
returned. The option can be specified as a variable name or a row matrix. Any 
missing entries of the variable or matrix will be silently dropped before
computing the estimates, i.e. the stored matrix e(L) will have as many rows as 
there are nonmissing observations/elements of x.

{phang}
{opt h(#)} is the bandwidth that will be used in the estimation.

{phang}
{opt k:ernel(namelist)} the kernel function that will be used to estimate the 
log-density. The default is the {it:epanechnikov} kernel. The other options are
{it:gaussian}, {it:quartic}, {it:triweight}, {it:triangle}, {it:uniform}, and 
{it:cosinus}.

{phang}
{opt de:gree(#)} specifies the order of the local polynomial approximation to 
the log-density that the estimation will use. The default is to use a local
linear approximation ({cmd:degree(1)}).

{phang}
{opt mi:nx(#)}, {opt ma:xx(#)} specify the support of the variable {it:varname}.
The default assumes the support is the entire real line. If the support is 
bounded, these options must be changed in order for the estimates to be 
consistent at the boundaries.

{phang}
{opt g(string)}, {opt dg(string)} specify the function and its derivative
that will be used to estimate the derivative(s) of the log-density. The function
g must return a vector of length {it:degree} and it must be zero outside of the
compact interval [-zl,zr], where zl = min((x-minx)/h,1) and
zr = min((maxx-x)/h,1). The function should be differentiable and dg should
return its derivative. The default value of g is a function of the function
whose {it:j}-th element is equal to g(u,zl,zr,S)=(u+zl)^j * (u-zr) for 
j=1,...,S and dg is the corresponding derivative with respect to u. The
user can replace these defaults by defining Mata functions called g2 and dg2 
that accept the same arguments (u,zl,zr,S). See the above paper on arXiv for a
discussion of alternative choices.

{marker remarks}{...}
{title:Remarks}

{pstd}
This estimation method is relatively customizable compared to alternative 
kernel-based density estimators. The current version of this Stata command  
accepts arbitrary choices of the function g, but does not directly support
alternative choices of the kernel function. Users can modify the Mata function
{it:kernel_f} and the list of available kernels in the Stata program if they
wish to explore these alternatives. All Mata and Stata code is in 
logdensity.ado.

{marker examples}{...}
{title:Examples}
{phang}{cmd:. preserve}{p_end}
{phang}{cmd:. drop _all}{p_end}
{phang}{cmd:. matrix x = (0,0.1,0.2)}{p_end}
{phang}{cmd:. set obs 1000}{p_end}
{phang}{cmd:. gen A = rchi2(2)}{p_end}
{phang}{cmd:. logdensity A, x(x) h(0.5) minx(0)}{p_end}
{phang}{cmd:. logdensity A, x(x) h(0.5) minx(0) degree(2)}{p_end}
{phang}{cmd:. mata}{p_end}
{phang}{cmd: function g2(u,zl,zr,S)}{p_end}
{phang}{{p_end}
{phang}{cmd:	power = (1::S)'}{p_end}
{phang}{cmd:	output = (((colnonmissing(power))'*(u:+zl)'):^power')':*((u:-zr):^2)}{p_end}
{phang}{cmd:	return(output)}{p_end}
{phang}{cmd: }}{p_end}
{phang}{cmd: mata mosave g2(), replace}{p_end}

{phang}{cmd:function dg2(u,zl,zr,S) }{p_end}
{phang}{{p_end}
{phang}{cmd:	power = (1::S)' }{p_end}
{phang}{cmd:	output = (((colnonmissing(power:-1))'*(u:+zl)'):^(power:-1)')':*((u:-zr):^2):*power + (((colnonmissing(power:-1))'*(u:+zl)'):^(power)')':*(2*(u:-zr))}{p_end}
{phang}{cmd:	return(output)}{p_end}
{phang}{cmd: }}{p_end}

{phang}{cmd:mata mosave dg2(), replace}{p_end}
{phang}{cmd:end}{p_end}
{phang}{cmd:. logdensity A, x(x) h(0.5) minx(0) g(g2) dg(dg2)}{p_end}
{phang}{cmd:. restore}

{marker storedresults}{...}
{title:Stored results}
{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}  The number of nonmissing observations used in estimation {p_end}
{synopt:{cmd:e(h)}}  The bandwidth used in estimation {p_end}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(L)}}   The estimated log-density and its derivatives{p_end}

{marker contact}{...}
{title:Contact}
{pstd}
Karl Schurter {browse "mailto:kschurter@psu.edu":kschurter@psu.edu}

{marker reference}{...}
{title:References}
{phang}
Pinkse, J. and Schurter, K. 2020.
{browse "https://arxiv.org/abs/2006.01328":Estimates of derivatives of (log) densities and related objects}.{p_end}




