clearvars

global p G T kB n0 R Gas_species Ion_species

kB = 1.38e-23; % Boltzmann constant [J/K]
R = 8.31; % Universal gas constant [J/mol/K]

p0 = 1e+5; % 2e+5; %  pressure [Pa] (partial pressure of the micture components under consideration, not including inert components)

B_N_ratio = 1e-4; % ratio of boron atoms to nitrogen atoms

condensation_flag = 0; % 1 - take carbon vapor condansation into account; 2 - don't take carbon vapor condansation into account

G_B1  =        [ 410633   263277  134165   16214       0       0 ]; 
G_B1i =        [1200701  1014080  834828  659674  581016  514406 ]; % B1- ions
G_B2 =         [ 638546   454380  307858  184461  265042  378396 ];
G_BN =         [ 366780   261719  175557  100657  127368  170113 ];
G_BN = G_BN + 1.73 * 1.60217662e-19 * 6.022140857e+23; % BN bond - 1.73 eV to agree with Darwent's reference

%G_BN_solid =   [-162652   -75347   24912    1e+6    1e+6    1e+6 ]; % At temperateratures above 3500K high values of 1e+6 are put just to exclude this species from the mixture; another option is extrapolation
G_BN_solid = ones(1,6) * 1e+6; % Prohibitively high energies

G_N1 =         [ 412171   346339  278946  210695  141670   71742 ]; 
G_N2 =         [      0        0       0       0       0       0 ];
G_N2 = G_N2 - 0.004 * 1.60217662e-19 * 6.022140857e+23; % N2 bond + 0.004 eV to agree with Darwent's reference
G_N3 =         [ 474048   527795  577999  626144  672874  718553 ];

G_BH =         [ 343045   249451  174733  111509  150303  205595 ];

G_BH2 =        [ 160806   127442  111035  104990  200384  312049 ];
G_BH3 =        [ 130510   170109  227544  295407  464745  650453 ];
G_H =          [ 165485   106760   46007  -15541  -77412  -139368];
%G_Hi =        [1457958  1350840 1230818 1102938  969484  831758 ];
G_H2 =         [      0        0       0       0       0       0 ];

G_NH =         [ 356503   336578  316493  296153  275467  254357 ];

G_NH2 =        [ 226491   268467  310569  351870  392177  431688 ];
G_N2H2 =       [ 324012   440882  555479  668356  780168  891379 ];
G_NH3 =        [  61910   179447  295689  410385  524233  637889 ];
G_N2H4 =       [ 320130   549740  772760  990928 1205835 1418579 ];

Temperatures = [ 1000     2000    3000    4000    5000     6000  ];

%{
%Gibbs energies according to J.Ragic-Peric:
for i=1:1:length(Temperatures)
    T = Temperatures(1,i)
    G_B1(i) = -196.36 * T + 597029;
    G_N1(i) = -196.24 * T + 509750;
    G_N2(i) = -258.12 * T + 59516;
    %G_BN(i) = -282.63 * T + 539708;
    G_B2(i) = -281.12 * T + 899234;
end
%G_BN = G_BN + 0.73 * 1.60217662e-19 * 6.022140857e+23;
%}

for i=1:1:length(Temperatures)
    T = Temperatures(1,i)
    
    G_B2N(i) = (267.78+0*20) * T - 1121.3e+3 + G_N1(i) + 2*G_B1(i); % Using Martin's 1989 association energy
    G_B3(i) = (274.89+0*20) * T - 849352 + 3*G_B1(i); % Using Martin's 1989 association energy
    G_BN(i)  = (145.18+0*10) * T - 430952  + G_N1(i) + G_B1(i); % Using Martin's 1989 association energy
end

n_init = ones(1,21);
A = -eye(21);
b = zeros(21,1); 


%                   G_B1     G_B1i    G_B2     G_B3     G_BN     G_B2N    G_N1     G_N2     G_N3     G_BH     G_BH2    G_BH3    G_H      G_H2     G_NH     G_NH2   G_N2H2   G_NH3    G_N2H4  G_BN_solid  G_B_ref


if (condensation_flag == 0) % don't take carbon vapor condansation into account
    Aeq0   =      [[1;0;0], [1;0;0], [2;0;0], [3;0;0], [1;1;0], [2;1;0], [0;1;0], [0;2;0], [0;3;0], [1;0;1], [1;0;2], [1;0;3], [0;0;1], [0;0;2], [0;1;1], [0;1;2], [0;2;2], [0;1;3], [0;2;4],  [0;0;0],  [0;0;0] ];
    G_BN_solid = G_BN_solid * 0;
else % take carbon vapor condansation into account
    Aeq0   =      [[1;0;0], [1;0;0], [2;0;0], [3;0;0], [1;1;0], [2;1;0], [0;1;0], [0;2;0], [0;3;0], [1;0;1], [1;0;2], [1;0;3], [0;0;1], [0;0;2], [0;1;1], [0;1;2], [0;2;2], [0;1;3], [0;2;4],  [1;1;0],  [1;0;0] ];
end

Gas_species =     [ 1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,       1,        0,        0      ];
Ion_species =     [ 0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,        0,        0      ];

fun = @(n) system_G(n);
res = zeros(length(n_init),1);

options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-50,'StepTolerance',1e-50,'OptimalityTolerance',1e-50,'MaxFunctionEvaluations',50000,'Algorithm','interior-point','HessianApproximation','bfgs','FiniteDifferenceType','central','TypicalX',ones(1,21)*1e-15);
p = p0;
for i=1:1:length(Temperatures)   
%for i=6:1:6
    Aeq = Aeq0;
    values_check = 0;
    G = [G_B1(1,i) G_B1i(1,i) G_B2(1,i) G_B3(1,i) G_BN(i) G_B2N(i) G_N1(1,i) G_N2(1,i) G_N3(1,i) G_BH(1,i) G_BH2(1,i) G_BH3(1,i) G_H(1,i) G_H2(1,i) G_NH(1,i) G_NH2(1,i) G_N2H2(1,i) G_NH3(1,i) G_N2H4(1,i) G_BN_solid(1,i) 0]; % Gibbs free energies of the components [J/mol], last zero is for boron reference state (liquid/solid)
    T = Temperatures(1,i); % Temperature [K]
    
    H_fract = 0; % fraction of Hydrogen
    B_fract = 0.5; % fraction of boron in the mixture (including solid/liquid boron)
    beq = [B_fract; 1-H_fract-B_fract; H_fract];

    n = fmincon(fun,n_init,A,b,Aeq,beq,[],[],[],options);

    n0 = p/kB/T; % total density (of all the species considered) [m^-3]
    res = [res, transpose(n) * n0 / dot(n,Gas_species)];
    
    %res = [res, transpose(n) * 2 * 1e+23];  % To compare with Dutouquet's results
end

function sG = system_G(n)
    global p G R T Gas_species Ion_species % kB
    n = max(n, 1e-16);
    n_e = dot(n, Ion_species);
    sum_n_gas = dot(n, Gas_species) + n_e;
    %p = sum_n_gas * 2 * 1e+23 * kB * T;  % To compare with Dutouquet (constant volume, or, in other words, constant number of species)     
    sG = dot(n,G) + R*T*dot(n.*Gas_species,log(n/sum_n_gas)) + R*T*n_e.*log(n_e/sum_n_gas) + R*T*log(p/1e+5)*sum_n_gas;
end

function pB = p_satur_B() % Boron saturation pressure (rather rough approximation)
    global T
    pB = 10^(-1.3516E-06 * T^2 + 1.1591E-02 * T - 2.0278E+01);
end