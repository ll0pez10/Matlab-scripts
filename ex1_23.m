syms la u0 ur A l1 l2 N ia ib N1 i1 g

%calculo das relutancias
ra = la/(u0*ur*A);
rb = (l1 + l2)/(u0*ur*A) + g/(u0*A);
%pretty(ra);
%pretty(simplify(rb));

%FMM
fmma = N*ia;
fmmb = N*ib;
fmmc = N1*i1;

%superposicao 1

rtot = ra * rb/(ra + rb) + ra;
%pretty(simplify(rtot));

fluxa = fmma / rtot;
%pretty(simplify(flux));

flux_esqa = (fmma - ra * flux)/rb;
flux_dira = (fmma - ra * flux)/ra;

%superposicao 2

rtot = ra * ra/(ra + ra) + rb;
%pretty(simplify(rtot));

fluxb = fmmb / rtot;
%pretty(simplify(flux));

flux_esqb = (fmmb - rb * flux)/ra;
flux_dirb = (fmmb - rb * flux)/ra;

%superposicao 3

rtot = ra * rb/(ra + rb) + ra;
%pretty(simplify(rtot));

fluxc = fmmc / rtot;
%pretty(simplify(flux));

flux_esqc = (fmmc - ra * flux)/rb;
flux_dirc = (fmmc - ra * flux)/ra;

% indutancia a

mutua_inducta = N*(flux_esqb - flux_dirc)/ib*ic;
self_inducta = N*fluxa/ia;
pretty(simplify(self_inducta));

%indutancia b

mutua_inductb = N*(flux_esqa - flux_dirc)/ia*ic;
self_inductb = N*fluxb/ib;
pretty(simplify(self_inductb));

%indutancia c

mutua_inductc = N1*(flux_dira + flux_dirb)/ib*ia;
pretty(simplify(mutua_inductc));
self_inductc = N1*fluxc/i1;
pretty(simplify(self_inductc));



