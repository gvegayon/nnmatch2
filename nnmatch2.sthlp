{smcl}
{* *! version 0.13.04 14apr2013}{...}

{cmd:help nnmatch2}
{hline}

{title:Title}

{phang}
{bf:nnmatch2} {hline 2} Nearest Neighbor Matching Estimation for Average Treatment Effects

{title:Syntax}

{p 8 17 2}
{cmdab:nnmatch2}
{it:{help varlist: y w x}} [, {opt met:ric}({it:{help string:maha|inverse}}) {opt m}({it:{help integer:#}}) {opt l:evel}({it:{help level}}) {opt tc}({it:{help string:ate|att|atc}})]

{synoptset 13 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt f:orce}}In order to protect the users' computers, setting more than 8 
clusters it is not permitted (see the {err:WARNING} in description). With this
option the user can skip this restriction.{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:parallel} allows to implement data parallelism algorithm in order to substantially
improve speed performance of it.
{p_end}

{pstd}
In order to use {cmd:parallel} it is necesary to set the number of desired clusters
with which the user wants to work with. To do this the user should use {cmd:parallel setclusters #}
syntaxis, replacing {it:#} with the desired number of clusters. Setting more clusters
than physical cores the user's computer has it is not recommended (see the {err:WARNING}
in description).
{p_end}

{pstd}
{cmd:parallel do} is the equivalent to {cmd:do}. By using this syntax, the
loaded dataset will be splitted in the number of clusters specified by {cmd:parallel setclusters}
and the {mansection U 16Do-files:do-file} will be executed independently over
each and every one of the data clusters. After all the parallel-instances stops, the datasets will be appended
> .
{p_end}

{pstd}
{cmd:parallel :} (as a prefix) allows to, after spliting the loaded dataset,
execute a {it:stata_cmd} over the specified number of data clusters in order to
speed up computations. Like {cmd:parallel do}, after all the parallel-instances
stops, the datasets will be appended.
{p_end}

{pstd}
Every time that {cmd:parallel} runs several auxiliary files are generated which,
after finishing, are automatically deleted. In the case that the user sets {opt k:eep}
or {opt keepl:last} the auxiliar files are kept, thus the syntax {cmd:parallel clean}
becomes handy. With {cmd:parallel clean} the user can remove the last generated
auxiliar files (default option), an specific parallel instance files (using
{it:#pll_id} number), or every kept auxiliar file (with {opt all}).
{p_end}

{pstd}
Giving {it:N} clusters, within each cluster {cmd:parallel} creates the local 
macros {it:pll_id} (equal for all the clusters) and {it:pll_instance} (ranging
from 1 up to {it:N}, equalling 1 inside the first and {it:N} inside the last).
Also the global macros {it:PLL_CLUSTERS} and {it:PLL_DIR} are available within
each cluster.
{p_end}

{pstd}
As by now, {cmd:parallel} by default automatically identifies stata's
executable file path. This is necessary as it is used to run stata in batch mode
(the core of the command). Either way, after some reports, that file path is not
always correctly identified; where {cmd:parallel setstatadir} can be used to
manually set it.
{p_end}

{pstd}
{err:WARNINGS} {it:(a)} For each cluster {cmd:parallel} starts a new stata instance (thus
running as many processes as clusters), this way, if the user sets more clusters
than cores his computer has, it is possible that his computer collapses.
{it:(b)} By this time {cmd:parallel} can not stop running the clusters by itself, what
implies that, in the case of any of the clusters starts a endless loop, stoping the
clusters should be done manually by the user by killing them from the OS's tasks manager.
{p_end}

{marker details}{...}
{title:Details}

{pstd}
Inspired by the R library ``snow'' and to be used in multicore CPUs
, {cmd:parallel} implements parallel computing methods through an OS's shell 
scripting (using Stata in batch mode) to speedup computations by splitting the
dataset into a determined number of clusters in such a way to implement a 
{browse "http://en.wikipedia.org/wiki/Data_parallelism":data parallelism} algorithm.
{p_end}

{pstd}
The number of efficient computing clusters depends upon the number of physical
cores (CPUs) with which your computer is built, e.g. if you have a quad-core
computer, the correct cluster setting should be four. In the case of simultaneous
multithreading, such as that from
{browse "http://www.intel.com/content/www/us/en/architecture-and-technology/hyper-threading/hyper-threading-te
> chnology.html":Intel's hyper-threading technology (HTT)},
setting {cmd:parallel} following the number of processors threads, as it was expected,
hardly results into a perfect speedup scaling. In spite of it, after several tests
on HTT capable architectures, the results of implementing {cmd:parallel} according
to the machines physical cores versus it logical's shows small though significant differences.
{p_end}

{pstd}
{cmd:parallel} is especially handy when it comes to implementing loop-based
simulation models (or simply loops), Stata commands such as reshape, or any job
that (a) can be repeated through data-blocks, and (b) routines that processes big
datasets.
{p_end}

{pstd}
At this time {cmd:parallel} has been successfully tested in Windows and Unix
machines.Tests using Mac OS are still pending.
{p_end}

{pstd}
After several tests, it has been proven that--thanks to how {cmd:parallel} has been
written--it is possible to use the algorithm in such a way that other techniques
of parallel computing can be implemented; such as Monte Carlo Simulations, 
simultaneously running models, etc.. An extensive
example through Monte Carlo Simulations is provided {browse "http://fmwww.bc.edu/repec/bocode/p/parallel.pdf":
> here}
{p_end}

{marker caveats}{...}
{title:Caveats}

{pstd}
If the {it:stata_cmd} or {it:do-file} {help saved_results:saves results},
as {cmd: parallel} runs stata in {browse "http://www.stata.com/support/faqs/windows/batch-mode/":batch mode},
none of the results will be keept. This is also true for {help matrix:matrices},
{help scalar:scalars} and {help mata:mata objects}.
{p_end}

{pstd}
Inspite {cmd:parallel} passes-through {help program list:programs}, {help macro:macros}
and {help mata:mata objects}, by the time it is not capable of doing the same with
{help matrix:matrices} or {help scalar:scalars}.
{p_end}

{pstd}
Including {cmd:parallel} within ado-files which contains locally-defined programs it
is not recommended due to knwon issues. By now parallel can not correctly
interpret this kind of programs during the loading process causing erros.
{p_end}

{marker examples}{...}
{title:Example 1: using prefix syntax}

{pstd}In this example we'll generate a variable containing the maximum 
blood-pressure measurement ({it:bp}) by patient.{p_end}

{pstd}Setup for a quad-core computer{p_end}
        {cmd:. sysuse bplong.dta}
        {cmd:. sort patient}
        
        {cmd:. parallel setclusters 4}

{pstd}Computes the maximum of {it:bp} for each patient. We add the option {opt by(patient)}
to tell parallel not to splitt stories.{p_end}
        {cmd:. parallel, by(patient): by patient: egen max_bp = max(bp)}
        
{pstd}Which is the ``parallel way'' to do:{p_end}

        {cmd:. by patient: egen max_bp = max(bp)}
        
{pstd}Giving you the same result.{p_end}
        
{title: Example 2: using {cmd:parallel do} syntax}

{pstd}Another usage that may get big benefits from it is implementing loop-base
simulations. Imagine that we have a model that requires looping over each and
every record of a panel-data dataset.
{p_end}

{pstd}
Using {cmd:parallel}, the proper way to do this would be using the ``parallel do''
syntax
{p_end}

        {cmd:. use mybigpanel.dta, clear}

        {cmd:. parallel setclusters 4}
        {cmd:. parallel do mymodel.do}
        
        {cmd:. collapse ...}

{pstd}where {it:mymodel.do} would look something like this{p_end}
        
        {hline 35} begin of do-file {hline 12}
        {cmd:local maxiter = _N}
        {cmd:forval i = 1/`maxiter'} {cmd:{c -(}}
                        {it:...some routine...}
        {cmd:{c )-}}
        {hline 35} end of the do-file {hline 10}

{marker examples}{...}
{title:Example 3: setting the right path}

{pstd}In the case of {cmd:parallel} setting the stata.exe's path wrongly, using
{cmd:parallel setstatadir} will correct the situation. So, if 
{it:"C:\Archivos de programa\Stata12/stata.exe"} is the right path we only have
to write:

        {cmd:. parallel setstatadir "C:\Archivos de programa\Stata12/stata.exe"}

{marker saved_results}{...}
{title:Saved results}

{pstd}
{cmd:parallel} saves the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(pll_n)}}Number of parallel clusters last used{p_end}
{synopt:{cmd:r(pll_id)}}Id of the last parallel instance executed (needed to use {cmd:parallel clean}){p_end}
{synopt:{cmd:r(pll_t_setu)}}Time took to setup (before the parallelization) and to finish the job (after the p
> arallelization){p_end}
{synopt:{cmd:r(pll_t_calc)}}Time took to complete the parallel job{p_end}
{synopt:{cmd:r(pll_t_fini)}}Time took to appending and cleaning{p_end}
{synopt:{cmd:r(pll_t_reps)}}In the case of {opt keeptime}, the number of time measures performed.{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(pll_dir)}}Directory where parallel ran and stored the auxiliary files.{p_end}

{marker references}{...}        

{title:References}

{phang}Luke Tierney, A. J. Rossini, Na Li and H. Sevcikova (2012). {it:snow: Simple Network of Workstations}. 
> R package version 0.3-9. {browse "http://CRAN.R-project.org/package=snow"}{p_end}
{phang}R Core Team (2012). {it:R: A language and environment for statistical computing}. R Foundation for Stat
> istical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL {browse "http://www.R-project.org/"}.{p_end}
{phang}George Vega Y (2012). {it:Introducing PARALLEL: Stata Module for Parallel Computing}. Chilean Pension S
> upervisor, Santiago de Chile, URL {browse "http://fmwww.bc.edu/repec/bocode/p/parallel.pdf"}.{p_end}

{title:Author}

{pstd}
George Vega Yon, Superindentencia de Pensiones. {browse "mailto:gvega@spensiones.cl"}
{p_end}

{title:Contributors}

{pstd}
Damian C. Clarke (Oxford University, England), Felix Villatoro (Superintendencia de Pensiones, Chile),
Eduardo Fajnzylber (Universidad Adolfo Ibanez, Chile), Eric Melse (CAREM, Netherlands),
Research Division (Superindentendia de Pensiones, Chile)
{p_end}

{title:Also see}

{psee}
Manual: {mansection "GSM CAdvancedStatausage":{bf:[GSM] Advanced Stata usage (Mac)}},
        {mansection "GSU CAdvancedStatausage":{bf:[GSU] Advanced Stata usage (Unix)}},
        {mansection "GSW CAdvancedStatausage":{bf:[GSW] Advanced Stata usage (Windows)}}

{psee}
Online: Running Stata batch-mode in {browse "http://www.stata.com/support/faqs/mac/advanced-topics/#batch": Ma
> c},
{browse "http://www.stata.com/support/faqs/unix/batch-mode/":Unix} and 
{browse "http://www.stata.com/support/faqs/windows/batch-mode/":Windows}
{p_end}
