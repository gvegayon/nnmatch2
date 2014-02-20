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
{cmd:nnmatch2} Implements a faster version of -nnmatch- algorithm which, without 
using parallel computing, allows improving speed up to x 10. Its implementation stills
in a beta status.
{p_end}

{title:Examples}

use ldw_exper, clear

// One match
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(1) bias(bias) 
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(1) bias(bias)

{title:Author}

{pstd}
George Vega Yon, Superindentencia de Pensiones. {browse "mailto:gvega@spensiones.cl"}
Eduardo Fajnzylber Reyes, Universidad Adolfo Ib{c a'}{c n~}ez. {browse "mailto:eduardo.fajnzylber@uai.cl"}
{p_end}

