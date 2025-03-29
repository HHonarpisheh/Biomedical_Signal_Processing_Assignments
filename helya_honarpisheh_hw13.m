close all;
clear all;

load granger_example.mat
% causal graphs are common cause
Q = 3;
pval_abc = granger_test(a,b,Q,c,'delayed') % <0.05 (significant)
% Consistent: c->b , a->b
% False discovery: c->a 
pval_acb = granger_test(a,c,Q,b,'delayed')
% Missing link: a->c
pval_bac = granger_test(b,a,Q,c,'delayed')
% Missing link: c->b
pval_bca = granger_test(b,c,Q,a,'delayed')
% Missing link: a->c , a->b
pval_cab = granger_test(c,a,Q,b,'delayed')
% Missing link: None
pval_cba = granger_test(c,b,Q,a,'delayed') % <0.05 (significant)
% Consistent: c->b , a->b , a->c
% False discovery: None