include("Application.jl");

type_prob = type_problème()

nv,nc = nbr_vars_cont()

sym_variables = positivité_variables(nv)
sym_contraintes = type_inégalité(nc)

M = saisirMatrice(nv,nc,type_prob)

S = standardiserSimplexe(M,type_prob,sym_variables,sym_contraintes)
T = copy(S)
resoudre!(S,type_prob)
resultats(S,T,type_prob,sym_variables,sym_contraintes)
