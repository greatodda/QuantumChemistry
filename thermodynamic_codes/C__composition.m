global p G T kB n0 R C_He_ratio Gas_species

kB = 1.38e-23; % Boltzmann constant [J/K]
R = 8.31; % Universal gas constant [J/mol/K]

p = 0.68e+5; %  pressure [Pa] (partial pressure of the micture components under consideration, not including inert components)
%p = 1e+8; %  pressure [Pa] (partial pressure of the micture components under consideration, not including inert components)
C_He_ratio = 0.01; % ratio of overal quantitity of carbon atom in all species of the mixture to helium atoms

Gas_species = [1 1 1 1 1 0];

condensation_flag = 0; % 1 - take carbon vapor condansation into account; 2 - don't take carbon vapor condansation into account


G_C5 =         [ 761595   653774  548360  445142   343949   244660  185962  147185  108687  70464  32513   -5168   -42584   -134974  -225751 ];
G_C4 =         [ 762052   658939  558096  459291   362341   267112  210758  173505  136500  99738  63216   26930    -9122    -98250  -185976 ];
G_C3 =         [ 601138   497428  397586  300796   206531   114437   60116   24266  -11307 -46611 -81655 -116447  -150944   -236329  -320256 ];
G_C2 =         [ 644008   546878  451481  357523   264759   173016  118405   82166   46055  10066 -25806  -61564   -97214   -185887  -273959 ];
G_C1 =         [ 560654   481427  402694  324474   264723   169389  123169   92427   61740  31105    512  -30012   -60497   -136499  -212211 ];

Temperatures = [   1000     1500    2000    2500     3000     3500    3800    4000    4200   4400   4600    4800     5000      5500     6000 ];

n_init = ones(1,6);

A = -eye(6);
b = zeros(6,1);

if (condensation_flag == 0) % don't take carbon vapor condansation into account
    Aeq = [1,2,3,4,5,0]; % last coefficient stands for solid state carbon
else % take carbon vapor condansation into account
    Aeq = [1,2,3,4,5,1]; % last coefficient stands for solid state carbon
end

fun = @(n) system_G(n);
res = zeros(length(n_init),1);
options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-16,'StepTolerance',1e-16,'OptimalityTolerance',1e-16); % ); % 'Algorithm','interior-point');

for i=1:1:length(Temperatures)   
    values_check = 0;
    G = [G_C1(1,i) G_C2(1,i) G_C3(1,i) G_C4(1,i) G_C5(1,i) 0]; % Gibbs free energies of the components, solid carbon is a reference state (G=0) [J/mol]
    T = Temperatures(1,i); % Temperature [K]
    beq = 1;
    n = fmincon(fun,n_init,A,b,Aeq,beq,[],[],[],options)

    n0 = p/kB/T; % total density (of all the species considered) [m^-3]
    res = [res, transpose(n) * n0*C_He_ratio/(1/dot(n,Gas_species) + C_He_ratio) / dot(n,Gas_species)];
end


function sG = system_G(n)
    global p G R T C_He_ratio Gas_species
    n = max(n,1e-16);
    sum_n_gas = dot(n,Gas_species);
    n_He = 1/C_He_ratio; % 1/C_He_ratio is normalized n_He
    sG = dot(n,G) + R*T*dot(n.*Gas_species,log(n/(sum_n_gas+n_He))) + R*T*log(p/1e+5)*sum_n_gas;
    sG = sG + R*T*n_He * log(n_He/(sum_n_gas+n_He)); 
end
