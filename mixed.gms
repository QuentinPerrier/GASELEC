*-------------------------------------------------------------------------------
*                                Defining the sets
*-------------------------------------------------------------------------------
sets     h               'all hours'                     /0*8783/
         first(h)        'first hour'
         last(h)         'last hour'
         m               'month'                         /jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec/
         vec             'energy vectors'                /electricity, gas/
         tec             'all the technologies'          /offshore, onshore, PV, river, lake, biogas, PHS, battery, hydrogen, methanation, G2P, methanisation, pyrogasification, tank, P2G/
         gen             'the generation technologies'   /offshore, onshore, PV, river, lake, biogas,methanisation, pyrogasification/
         elec(tec)       'electricity technologies'      /offshore, onshore, PV, river, lake, biogas, PHS, battery, hydrogen, methanation, G2P/
         gas(tec)        'gas technologiess'             /methanisation, pyrogasification, tank, P2G/
         elec_gen(tec)   'electricity generation tecs'   /offshore, onshore, PV, river, lake, biogas/
         gas_gen(tec)    'gas generation technologies'   /methanisation, pyrogasification/
         str(tec)        'storage technologies'          /battery, phs, hydrogen, methanation, tank/
         str_elec(str)   'electric storage technologies' /battery, phs, hydrogen, methanation/
         str_gas(str)    'gas storage technologies'      /tank/
         conv(tec)       'vector change technologies'    /P2G, G2P/
         vre(tec)        'variable tecs'                 /offshore, onshore, PV/
         frr(tec)        'technologies for upward FRR'   /lake, PHS, battery, G2P, biogas/
;
first(h) = ord(h)=1;
last(h) = ord(h)=card(h);
alias(h,hh);
*-------------------------------------------------------------------------------
*                                Inputs
*-------------------------------------------------------------------------------
$ontext
2016 had 366 days, and the hours of each month is as below:
January from 0 to 743, February from 744 to 1439, March from 1440 to 2183,
April from 2184 to 2903, May from 2904 to 3647, June from 3648 to 4367,
July from 4368 to 5111, August from 5112 to 5855, September from 5856 to 6575,
October from 6576 to 7319, November from 7320 to 8039 and December from 8040 to 8783.
$offtext
parameter month(h)  /0*743 1, 744*1439 2, 1440*2183 3, 2184*2903 4
                    2904*3647 5, 3648*4367 6, 4368*5111 7, 5112*5855 8
                    5856*6575 9, 6576*7319 10, 7320*8039 11, 8040*8783 12/
$Offlisting
parameter load_factor(vre,h) 'Production profiles of VRE'
/
$ondelim
$include  inputs/vre_inputs.csv
$offdelim
/;
parameter demand_elec(h) 'electricity demand profile in each hour in MW'
/
$ondelim
$include inputs/dem_electricity.csv
$offdelim
/;
parameter demand_gas(h) 'gas demand profile in each hour in MW'
*resource: GRTgas + TIGF
/
$ondelim
$include inputs/dem_gas.csv
$offdelim
/;
Parameter lake_inflows(m) 'monthly lake inflows in GWh'
*Resource: RTE - Hourly nationwide electricity generation by sectors in 2016 for France
/
$ondelim
$include  inputs/lake_inflows.csv
$offdelim
/ ;
parameter gene_river(h) 'hourly run of river power generation in GWh'
*Resource: RTE - Hourly nationwide electricity generation by sectors in 2016 for France
/
$ondelim
$include  inputs/run_of_river.csv
$offdelim
/ ;
parameter epsilon(vre) 'additional FRR requirement for variable renewable energies because of forecast errors'
/
$ondelim
$include  inputs/reserve_requirements.csv
$offdelim
/ ;

parameter capa_ex(tec) 'existing capacities of the technologies by December 2017 in GW'
*Resource: RTE
/
$ondelim
$include  inputs/existing_capas.csv
$offdelim
/ ;
parameter capex(tec) 'annualized capex cost in M€/GW/year'
/
$ondelim
$include  inputs/annuities_mix.csv
$offdelim
/ ;
parameter capex_en(str)  'annualized capex cost in M€/GWh/year'
/
$ondelim
$include inputs/str_annuities.csv
$offdelim
/ ;
parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW'
/
$ondelim
$include  inputs/fO&M_mix.csv
$offdelim
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&M_mix.csv
$offdelim
/ ;
$Onlisting
parameter eta_in(str) 'charging efifciency of storage technologies' /PHS 0.95, battery 0.9, hydrogen 0.85, methanation 0.675, tank 1/;
parameter eta_out(str) 'discharging efficiency of storage technolgoies' /PHS 0.9, battery 0.95, hydrogen 0.5, methanation 0.43, tank 1/;
parameter eta_conv(conv) 'the vector conversion efficiency of vector change options' /G2P 0.6, P2G 0.8/;
scalar pump_capa 'pumping capacity in GW' /9.3/;
scalar max_phs 'maximum volume of energy can be stored in PHS reservoir in TWh' /0.18/;
scalar max_hydrogene 'maximum energy that can be stored in the form of hydrogene in TWh' /5/;
scalar max_biogas 'maxium energy can be generated by biogas in TWh' /15/;
scalar load_uncertainty 'uncertainty coefficient for hourly demand' /0.01/;
scalar delta 'load variation factor'     /0.1/;
scalar tank_max 'maximum gas storage volume in TWh'    /160/
scalar max_metha 'maximum yearly energy can be produced from methanisation in TWh'       /140/;
scalar max_pyro 'maximum yearly energy can be produced from pyrogasification in TWh'       /170/ ;
*-------------------------------------------------------------------------------
*                                Model
*-------------------------------------------------------------------------------
variables        GENE(tec,h)     'hourly energy generation in TWh'
                 CAPA(tec)       'overal yearly installed capacity in GW'
                 STORAGE(str,h)  'hourly electricity input of battery storage GW'
                 STORED(str,h)   'energy stored in each storage technology in GWh'
                 CONVERT(conv,h) 'vector conversion from electricity to gas and the opposite in GW'
                 CAPACITY(str)   'energy volume of storage technologies in GWh'
                 RSV(frr,h)      'required upward frequency restoration reserve in GW'
                 COST            'final investment cost in b€'
positive variables GENE(tec,h),CAPA(tec),STORAGE(str,h),STORED(str,h),CAPACITY(str),RSV(frr,h),CONVERT(conv,h) ;
equations        gene_vre                'variables renewable profiles generation'
                 gene_capa               'capacity and genration relation for technologies'
                 capa_frr                'capacity needed for the secondary reserve requirements'
                 storing_const           'the definition of stored energy in the storage options'
                 conversion              'the mechanism of conversion'
                 storage_const           'storage in the first hour is equal to the storage in the last hour'
                 lake_res                'constraint on water for lake reservoirs'
                 stored_cap              'maximum energy that is stored in storage units'
                 biogas_const            'maximum energy can be produced by biogas'
                 methanisation_max       'the maximum yearly possible energy output from methanisation'
                 pyrogasification_max    'the maximum yearly possible energy output from pyrogasification'
                 reserves                'FRR requirement'
                 adequacy_elec           'supply/demand relation for the electricity vector'
                 adequacy_gas            'supply/demand relation for the gas vector'
                 obj                     'the final objective function which is COST';
gene_vre(vre,h)..                GENE(vre,h)                     =e=     CAPA(vre)*load_factor(vre,h);
gene_capa(tec,h)..               CAPA(tec)                       =g=     GENE(tec,h);
capa_frr(frr,h)..                CAPA(frr)                       =g=     GENE(frr,h) + RSV(frr,h);
storing_const(h,h+1,str)..       STORED(str,h+1)                 =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);
storage_const(str,first,last)..  STORED(str,first)               =e=     STORED(str,last);
lake_res(m)..                    lake_inflows(m)                 =g=     sum(h$(month(h) = ord(m)),GENE('lake',h));
stored_cap(str,h)..              STORED(str,h)                   =l=     CAPACITY(str);
conversion(h,conv)..             GENE(conv,h)                    =e=     CONVERT(conv,h)*eta_conv(conv);
biogas_const..                   sum(h,GENE('biogas',h))         =l=     max_biogas*1000;
methanisation_max..              sum(h,GENE('methanisation',h))  =l=     max_metha*1000;
pyrogasification_max..           sum(h,GENE('pyrogasification',h))=l=    max_pyro*1000;
reserves(h)..                    sum(frr, RSV(frr,h))            =e=     sum(vre,epsilon(vre)*CAPA(vre))+ demand_elec(h)*load_uncertainty*(1+delta);
adequacy_elec(h)..               sum(elec,GENE(elec,h))          =g=     demand_elec(h) + sum(str_elec,STORAGE(str_elec,h)) + CONVERT('P2G',h);
adequacy_gas(h)..                sum(gas,GENE(gas,h))            =g=     demand_gas(h) + sum(str_gas,STORAGE(str_gas,h)) + CONVERT('G2P',h);
obj..                            COST                            =e=     (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,CAPACITY(str)*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec))) +sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;
*-------------------------------------------------------------------------------
*                                Initial and fixed values
*-------------------------------------------------------------------------------
CAPA.up('onshore') = 173;
CAPA.lo(tec) = capa_ex(tec);
CAPA.fx('phs') = pump_capa;
CAPA.fx('river')= capa_ex('river');
CAPA.fx('lake') = 13;
GENE.up('river',h) = gene_river(h)*capa_ex('river');
STORAGE.up('phs',h) = pump_capa;
CAPACITY.fx('phs') = max_phs*1000;
CAPACITY.fx('hydrogen') = max_hydrogene*1000;
CAPACITY.fx('tank') = tank_max*1000;
*-------------------------------------------------------------------------------
*                                Model options
*-------------------------------------------------------------------------------
model RESEG /all/;
option solvelink=2;
option RESLIM = 1000000;
option lp=cplex;
option Savepoint=1;
option solveopt = replace;
option limcol = 0;
option limrow = 0;
option SOLPRINT = OFF;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
$If exist RESEG_p1.gdx execute_loadpoint 'RESEG_p1';
Solve RESEG using lp minimizing COST;
display cost.l; display capa.l; display gene.l;
parameter sumgene_elec;
sumgene_elec = sum((elec_gen,h),GENE.l(elec_gen,h))/1000;
parameter gene_elec(elec);
gene_elec(elec) = sum(h,GENE.l(elec,h))/1000;
parameter sumgene_gas;
sumgene_gas = sum((gas_gen,h),GENE.l(gas_gen,h))/1000;
parameter gene_gas(gas);
gene_gas(gas) = sum(h,GENE.l(gas,h))/1000;
parameter dem_gas;
dem_gas = sum(h,demand_gas(h))/1000;
parameter dem_elec;
dem_elec = sum(h,demand_elec(h))/1000;
display sumgene_elec; display sumgene_gas;
display gene_elec;
display gene_gas;
display dem_elec; display dem_gas;
display capacity.l;
parameter lcoe;
lcoe = 1000*COST.l/(sumgene_elec+sumgene_gas) ;
parameter lc_elec;
lc_elec = (sumgene_elec - dem_elec)/sumgene_elec;
display lcoe;
display lc_elec;
parameter vec_change(conv);
vec_change(conv) = smax(h,CONVERT.l(conv,h));
display CAPACITY.l;
display vec_change;
*-------------------------------------------------------------------------------
*                                Outputs
*-------------------------------------------------------------------------------
$Ontext
two main output files;
The .txt file just to have a summary and general idea of the key numbers
The .csv file to have a fine output with hourly data for final data processing and analysis
$Offtext

file results_mix /results_mix.txt/ ;
*the .txt file
put results_mix;
put '                            the main results' //
//
'I)Overall investment cost is' cost.l 'b€' //
//
'II)the Renewable capacity ' //
'PV              'capa.l('PV')'  GW'//
'Offshore        'capa.l('offshore')'    GW'//
'onsore          'capa.l('onshore')'     GW' //
'run of river    'CAPA.l('river') 'GW' //
'lake            'CAPA.l('lake') 'GW' //
'biogas          'CAPA.l('biogas')' GW'//
'methanisation   'CAPA.l('methanisation')' GW'//
'pyrogasification 'CAPA.l('pyrogasification')' GW'//
//
//
'III)Needed storage volume' //
'Battery Storage         'CAPACITY.l('battery')'       GWh' //
'PHS Storage             'CAPACITY.l('phs')'       GWh'//
'hydrogen storage        'CAPACITY.l('hydrogen')' GWh'//
'methane storage         'CAPACITY.l('methanation')'   GWh'//
'tank                    'CAPACITY.l('tank')' GWh'//
//
//
'IV) Conversion units'   //
'P2G                     'vec_change('P2G')' GWh'//
'G2P                     'vec_change('G2P')' GWh'//
;

file results_mixx /results_mix.csv / ;
*the .csv file
parameter nSTORAGE(str,h);
nSTORAGE(str,h) = 0 - STORAGE.l(str,h);
put results_mixx;
results_mixx.pc=5;
put 'hour'; loop(elec, put elec.tl;); put 'elec_demand' ;put 'ElecStr' ;put 'Pump' ; put 'H2' ; put 'CH4' ; loop(gas, put gas.tl;); put 'gas_demand'; put 'tank'/ ;
loop (h,
put h.tl; loop(elec, put gene.l(elec,h);) ;put demand_elec(h); put nSTORAGE('PHS',h) ; put nSTORAGE('battery',h) ; put nSTORAGE('hydrogen',h) ; put nSTORAGE('methanation',h); loop(gas, put gene.l(gas,h);) ; put demand_gas(h); put nSTORAGE('tank',h) /
;);
*-------------------------------------------------------------------------------
*                                The End :D
*-------------------------------------------------------------------------------