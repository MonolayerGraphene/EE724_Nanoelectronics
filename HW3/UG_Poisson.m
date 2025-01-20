
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%                                                   %%  
         %%     1D Drift Diffusion Model for pn Diodes        %%  
         %%     Equilibrium and Non Equilibrium Solver        %%
         %%                                                   %%  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
clear all;
close all;
tic;

% Defining the Fundamental and Material Constants %
nm    = 1e-7;             %cm
mfactor=1;              % arbitrary reduction in mobility

q     = 1.602E-19;        % C or [J/eV]
kb    = 1.38E-23;         % [J/K]
eps   = 1.05E-12;         % This includes the eps  = 11.7 for Si [F/cm]
T     = 300;              % [K]
ni    = 1.5E10;           % Intrinsic carrier concentration [1/cm^3]
Vt    = kb*T/q;           % [eV]
RNc   = 2.8E19;           % This is 2.8e20 in the FORTRAN file
TAUN0 = 0.1E-6;           % Electron SRH life time
TAUP0 = 0.1E-6;           % Hole SRH life time
mun0   = 1500/mfactor;            % Electron Mobility in cm2/V-s
mup0   = 1000/mfactor;             % Hole Mobility in cm2/V-s
                        
                        
dEc = Vt*log(RNc/ni);

% % MODIFY HERE- Define Doping Values 

Na = 1E17;             % [1/cm^3]
Nd = 1E19;             % [1/cm^3]


% % MODIFY BEGIN Define npn lengths % <<<<<<<<<<<<<<<<<<<<<<<NPN<<<<<<<<<<<<<<<<
profile=1; % set profile=0 for PN junction or profile =1 for NPN 
LLn=20*nm;                 %n region length
LLp=40*nm;                 %n region length
% % MODIFY END Define npn lengths % <<<<<<<<<<<<<<<<<<<<<<<NPN<<<<<<<<<<<<<<<<



% Define Max Voltage 
Vmax=0.6;

% MODIFY BEGIN filename for saving
fname='run_pn_5';
% MODIFY BEGIN filename for saving

% Calculate relevant parameters for the simulation %

Vbi = Vt*log(Na*Nd/(ni*ni));
W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd))     % [cm]
Wn  = W*sqrt(Na/(Na+Nd))                    % [cm]
Wp  = W*sqrt(Nd/(Na+Nd))                    % [cm]
Wone = sqrt(2*eps*Vbi/(q*Na))               % [cm]
E_p = q*Nd*Wn/eps                           % [V/cm]
Ldn = sqrt(eps*Vt/(q*Nd));              % extrinsic debye length
Ldp = sqrt(eps*Vt/(q*Na));                  % extrinsic debye length
Ldi = sqrt(eps*Vt/(q*ni));                  % intrinsic debye length


delta_acc = 1E-8;               % Preset the Tolerance

% Calculate relevant parameters in an input file %

% Write to a file
save input_params.txt Na Nd Vbi W Wn Wp E_p Ldn Ldp

%Material_Constants    %Define some material constants

% Setting the grid size based on the extrinsic Debye lengths %

% selection of doping profile
if profile==1    
x_max= 2*LLn+LLp; % npn% <<<<<<<<<<<<<<<<<<<<<<<NPN<<<<<<<<<<<<<<<<
else
    nfactor=3; % scaling for simulation space
    x_max= nfactor*(Wn+ Wp); % for pn % <XXXXXXXXXXXXX<<PN XXXXXXXXXXXXXXXX
end
%

dxfactor=40; %reduction of factor for Ldmin
dx=min([Ldp Ldn Ldi])/dxfactor;
n_max=round(x_max/dx);
dx=x_max/(n_max-1);


% Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni %

x=[0:dx:x_max];

if profile==1
yn=(x<LLn|x>LLn+LLp);       % n region is defined for npn% <<<<<<<<<<<<<<<<<<<<<<<NPN<<<<<<<<<<<<<<<<
else
    yn=(x<nfactor*Wn);       % n region is defined % <<<<<<<<<<<<<<<<<<<<<<<PN<<<<<<<<<<<<<<<<
end
yp=ones(1,n_max)-yn;      % p region is defined
dop=yn*Nd/ni- yp*ni/ni;     % defining doping at grid points

figure (1); plot(x, dop);

dx = dx/Ldi;    % Renormalize lengths with Ldi
% 
% Initialize the potential based on the requirement of charge
% neutrality throughout the whole structure


% use n=ND and p=NA to set initial guess for fi

fi=sign(dop).*log(abs(dop));

%figure (10); plot(x, Vt.*fi); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                      %%
%%               EQUILIBRIUM  SOLUTION PART BEGINS                      %%
%%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %(A) Define the elements of the coefficient matrix for the internal nodes and
%    initialize the forcing function


dx2=dx.*dx;

p= exp(-fi); 
n= exp(fi);
fi=UGEstatics(fi, dop, n, p , n_max, dx, delta_acc, x);
     figure (10); plot(x, Vt*fi, 'r'); hold on;
 

xx1(1) = dx*1e4;

    Ec = dEc - Vt*fi;     %Values from the second Node%
    ro = -ni*(exp(fi) - exp(-fi) - dop);
    el_field1 = -([fi 0 0] - [0 fi 0])*Vt/(dx*Ldi);
    el_field2 = -([fi 0 0] - [ 0 0 fi])*Vt/(2*dx*Ldi);
    n = exp(fi);
    p = exp(-fi);
    xx1 = x.*1e4/Ldi;
    
    el_field1 = el_field1(2: n_max+1);
    el_field2 = el_field2(2: n_max+1);    
    
Ec(1) = Ec(2);
Ec(n_max) = Ec(n_max-1);
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field1(1) = el_field1(2);
el_field2(1) = el_field2(2);
el_field1(n_max) = el_field1(n_max-1);
el_field2(n_max) = el_field2(n_max-1);
n(n_max)=n(n_max-1);
p(n_max)=p(n_max-1);

nf = n*ni;
pf = p*ni;
ro(1) = ro(2);
ro(n_max) = ro(n_max-1);

figure (11); semilogy(x, abs(n), 'r',x, abs(p), 'b' ); hold on;
    

figure(1)
plot(xx1, Vt*fi,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Potential [eV]');
title('Potential vs Position - at Equilibrium');

figure(2)
plot(xx1, el_field1,'r','LineWidth',2)
hold on;
plot(xx1, el_field2,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Electric Field [V/cm]');
title('Field Profile vs Position - at Equilibrium');
hold off;

figure(3); 
%plot(xx1, nf,'g','LineWidth',2)
plot(xx1, nf,'g','LineWidth',2)
hold on;
%plot(xx1, pf,'r','LineWidth',2)
plot(xx1, pf,'r','LineWidth',2)
hold on;
xlabel('x [um]');
ylabel('Electron & Hole Densities [1/cm^3]');
title('Electron & Hole Densities vs Position - at Equilibrium');
legend('n','p');
%axis([0 6.75 0 10.2e17])
hold off;

figure(4)
%plot(xx1, ro,'r','LineWidth',2)
plot(xx1, q*ro,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Total Charge Density [C/cm^3]');
title('Total Charge Density vs Position - at Equilibrium');
%axis([0.5 5 -3e17 8e17])

figure(5)
plot(xx1, Ec,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Conduction Band Energy (eV)');
title('Conduction Band vs Position - at Equilibrium');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 END OF EQUILIBRIUM  SOLUTION PART                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%                                                                      %%
% %%               NON-EQUILIBRIUM  SOLUTION PART BEGINS                  %%
% %%                                                                      %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%              1. Calculate Low filed mobility                         %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%  Prameters for Low field mobility calculation %%
% 
%     TL = 300;                    % Temp in Kelvin
%     N  = Na + Nd;                % Local (total) impurity concentration
% 
%     MU1N.CAUG   = 55.24;         % cm2/(V.s)
%     MU2N.CAUG   = 1429.23;       % cm2/(V.s)
%     ALPHAN.CAUG = 0.0;           % unitless
%     BETAN.CAUG  = -2.3;          % unitless
%     GAMMAN.CAUG = -3.8;          % unitless
%     DELTAN.CAUG = 0.73;          % unitless
%     NCRITN.CAUG = 1.072*10^17;   % cm-3
% 
%     MU1P.CAUG   = 49.7;          % cm2/(V.s)
%     MU2P.CAUG   = 479.37;        % cm2/(V.s)
%     ALPHAP.CAUG = 0.0;           % unitless
%     BETAP.CAUG  = -2.2;          % unitless
%     GAMMAP.CAUG = 13.7;          % unitless
%     DELTAP.CAUG = 0.70;          % unitless
%     NCRITP.CAUG = 1.606*10^17;   % cm-3
%     BETAN = 2.0;
%     BETAP = 1.0;
% 
%     
% 
%     VSATN = (2.4*10^7) / (1 + 0.8*exp(TL/600));  % Saturation Velocity of Electrons
%     VSATP = VSATN                              % Saturation Velocity of Holes
% 
% %%%%%%%%%%%%%%%%%%% END of Low Field Mobility Calculation %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2. Start the main Loop to increment the Anode voltage by Vt=KbT/q  %% 
% %%      till it reaches 0.625V.                                         %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% vindex=0;
% Vstepf=0.05;
% 
% 
% 
% 
%     Each_Step   = Vstepf*Vt;  
%     Total_Steps = round(Vmax/(Vstepf*Vt));
% 
%     Jnim1by2=[];%zeros(Total_Steps, n_max);
%     Jnip1by2=[];
%     Jpim1by2=[];
%     Jpip1by2=[];
%     Jhole=[];
%     Jelec=[];
%    Jelec1=[];
%    Jhole1=[];
%     Jtotal=[];
% for VA = 0:Each_Step:Vmax                % Start VA increment loop
%     VA 
%     
% 
%     vindex = vindex +1;
%     Vplot(vindex) = VA;
%     
%     fi(1) = fi(1) - Each_Step./Vt;            % Apply potential to Anode (1st node)  
%     
%     flag_conv2 = 0;		           % Convergence of the Poisson loop
%     k_itern= 0;
% %% Initialize the First and Last Node for Poisson's eqn
%     a(1) = 0;
%     c(1) = 0;
%     b(1) = 1;
%     f(1) = fi(1);
%     a(n_max) = 0;
%     c(n_max) = 0;
%     b(n_max) = 1;
%     f(n_max) = fi(n_max);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% 3. Start the Poisson equation solver loop to calculate the        %%
%     %%    potential for each Anode voltage increase                      %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% while(~flag_conv2)             % Start Poisson's eqn
%         
%         k_itern = k_itern + 1;
% 
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   %% 
%         %%       at each node point of the PN diode.                         %%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
% 
% 
%         %% Calculate the Electric Field at each Node
%         
%             Ef = abs([fi 0] - [0 fi]).*Vt/(dx*Ldi); %add a point to 
%                                   %calculate difference in matrix form
%             i;
%             
%             Ef(1)     = Ef(3); %assume that field is low at the edge
%             Ef(2)     = Ef(3); 
% 
%             Ef(n_max) = Ef(n_max-1);
%             Ef=Ef(1:n_max);               %reject the last point
% 
%   %         figure (16); plot(Ef); title('Ef')
%         %% Calculate the Field Dependant Mobility at each Node
%         
% 
%               pdeno  = (mup0 .* Ef / VSATP).^ BETAP;
%               mup   = mup0 * ( (1./(1 + pdeno)).^(1/BETAP)); 
%                             
%               ndeno  = (mun0 * Ef / VSATN).^ BETAN;
%               mun = mun0 * ( (1./(1 + ndeno)).^(1/BETAN)); 
%                             
%             mup(1)     = mup(2);
%         mup(n_max) = mup(n_max-1);
%         
%         mun(1)     = mun(2);
%         mun(n_max) = mun(n_max-1);
% 
%  %             figure (17); plot(xx1, mun, xx1, mup); title('mu')
% 
%    
% 
%  %%%%%%%%%%% END of FIELD Dependant Mobility Calculation %%%%%%%%%%% 
%       
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%         %% 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition %%                                    
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
% 
% 
% 
%         %(B) Define the elements of the coefficient matrix for the internal nodes and
%         %    initialize the forcing function
% 
% 
%         %for i = 2: n_max-1
%             munip1by2 = ([mun 0 0]+[0 mun 0])/2;        % i plus and i minus 
%             munim1by2 = ([0 mun 0]+[0 0 mun])/2;
%             mupip1by2 = ([mup 0 0]+[0 mup 0])/2;     
%             mupim1by2 = ([0 mup 0]+[0 0 mup])/2; 
% 
%             
%             %% Co-efficients for HOLE Continuity eqn
%                 cp = mupip1by2 .* BER([0 fi 0] - [fi 0 0 ]);
%                 ap = mupim1by2.*BER([0 fi 0] - [ 0 0 fi]);
%                 bp = -(mupim1by2.*BER([ 0 0 fi]-[0 fi 0]) + mupip1by2 .* BER([fi 0 0 ]-[0 fi 0] ));
%  
%             %% Co-efficients for ELECTRON Continuity eqn
%                 cn = munip1by2 .* BER([fi 0 0 ] -[0 fi 0]);
%                 an = munim1by2.*BER([ 0 0 fi] - [0 fi 0]);
%                 bn = -(munim1by2.*BER([0 fi 0]-[ 0 0 fi]) + munip1by2 .* BER([0 fi 0]-[fi 0 0 ]));
% 
%             %% Forcing Function for ELECTRON and HOLE Continuity eqns
%        fn = 0;%(Ldi*Ldi*dx2/Vt) .* ( p.*n - 1 ) ./ ( TAUP0*(n + 1) + TAUN0*(p+1)); % OK
%        fp = 0;%(Ldi*Ldi*dx2/Vt) .* ( p.*n - 1 )./ ( TAUP0*(n + 1) + TAUN0*(p+1));  % OK
%         %end
%             an = an(2:n_max+1);        % i plus and i minus 
%             bn = bn(2:n_max+1);        % i plus and i minus 
%             cn = cn(2:n_max+1);        % i plus and i minus 
%             %fn = fn(2:n_max+1);        % i plus and i minus 
%            
%             ap = ap(2:n_max+1);        % i plus and i minus 
%             bp = bp(2:n_max+1);        % i plus and i minus 
%             cp = cp(2:n_max+1);        % i plus and i minus 
%             %fp = fp(2:n_max+1);        % i plus and i minus 
% 
%             munip1by2 = munip1by2(2:n_max+1);        % i plus and i minus 
%             munim1by2 = munim1by2(2:n_max+1);
%             mupip1by2 = mupip1by2(2:n_max+1);     
%             mupim1by2 = mupim1by2(2:n_max+1); 
% 
% 
% 
% 
%             % at the ohmic contacts for ELECTRON and HOLE Continuity Eqns 
% 
%                     %(A) Overwrite the elements of the coefficient matrix at contacts and initialize the forcing
%         an(1) = 0;              %Co-ef for electron at Anode
%         bn(1) = 1;              %Co-ef for electron at Anode
%         cn(1) = 0;              %Co-ef for electron at Anode
%         ap(1) = 0;              %Co-ef for hole     at Anode
%         bp(1) = 1;              %Co-ef for hole     at Anode
%         cp(1) = 0;              %Co-ef for hole     at Anode
%         %fnp(1) = (Ldi*Ldi*dx2/Vt) * ( p(1)*n(1) - 1 ) / ( TAUP0*(n(1) + 1 ) + TAUN0*(p(1) + 1 ) );
%         fn(1) = n(1);
%         fp(1) = p(1);
%         
%         an(n_max) = 0;          %Co-ef for electron at Cathode
%         bn(n_max) = 1;          %Co-ef for electron at Cathode
%         cn(n_max) = 0;          %Co-ef for electron at Cathode
%         ap(n_max) = 0;          %Co-ef for hole     at Cathode
%         bp(n_max) = 1;          %Co-ef for hole     at Cathode
%         cp(n_max) = 0;          %Co-ef for hole     at Cathode
%         %fnp(n_max) = (Ldi*Ldi*dx2/Vt) * ( p(n_max)*n(n_max) - 1 ) / ( TAUP0*(n(n_max) + 1) + TAUN0*(p(n_max) + 1) );
%         fn(n_max) = n(n_max);
%         fp(n_max) = p(n_max);
% 
% % figure (102); plot((an-an1)./an); 
% % figure (103); plot((bn-bn1)./bn); 
% % figure (104); plot((cn-cn1)./cn); 
% % figure (105); plot((ap-ap1)./ap); 
% % figure (106); plot((bp-bp1)./bp); 
% % figure (107); plot((cp-cp1)./cp); 
% % 
% % pause
%             
%         %(C)  Start the iterative procedure for the solution of the linearized Continuity
%         %     equation for "ELECTRONS" using LU decomposition method:
%  
%         n=LUdecomp3(an,bn,cn,fn, n_max); % calling LU decomposition function
%                     %figure (18); semilogy(xx1, n);hold on
% 
%      
%             %%%%%%%%%%%%%%%%%%%%%%% END of ELECTRON Continuty Solver %%%%%%%%%%%  
%             
%             
%         %(D)  Start the iterative procedure for the solution of the linearized Continuity
%         %     equation for "HOLES" using LU decomposition method:
%             
%        p=LUdecomp3(ap,bp,cp,fp, n_max); % calling LU decomposition function
%                     %figure (18); semilogy(xx1, p);hold on;
% %figure (11); semilogy(x, abs(n), 'r',x, abs(p), 'b'); hold on;
% 
%         
%        %%%%%%%%%%%%%%%%%%%%%%% END of HOLE Continuty Solver %%%%%%%%%%%  
%        
%       % figure (11); semilogy(x, abs(n), 'r',x, abs(p), 'b'); hold on;
%        
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%        %% 3.3 Calculate potential fi again with new values of "n" and "p"%%
%        %%     and check for convergence                                  %%     
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        
%        % Recalculate forcing function and central coefficient b for fi 
%         a = ones(1,n_max)/dx2;
%         c = ones(1,n_max)/dx2;
%         b = -(2/dx2 + n + p);
%         f = n - p - dop - fi.*(n + p);
%             
% %% Initialize the First and Last Node for Poisson's eqn
%     a(1) = 0;
%     c(1) = 0;
%     b(1) = 1;
%     f(1) = fi(1);
%     a(n_max) = 0;
%     c(n_max) = 0;
%     b(n_max) = 1;
%     f(n_max) = fi(n_max);
%             
%             % Solve for Updated potential given the new value of Forcing 
%         % Function using LU decomposition 
%         
%         fiold=fi;  % storing previous fi to compare error
% fi=LUdecomp3(a,b,c,f, n_max); % calling LU decomposition function
% delta = fi - fiold;
%     delta_max=max(abs(delta));
%     delta_max (k_itern) = max(abs(delta));
%     
%     
% %sprintf('delta_max = %d',delta_max);      %'k_iter = %d',k_iter,'
%         
%         if(delta_max < delta_acc)
%             flag_conv2 = 1;
%         end
%  %       flag_conv2 = 1;
%  %       figure (12); plot(x, fi, 'r'); hold on;
% 
% end    
% 
% 
%   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%                        CALCULATE CURRENT                             %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%           % Electron Current       
%     
%         Jnim1by2a  = (q.*[0 mun 0]*Vt/(dx*Ldi)) * ni.*( [0 n 0] .*BER([0 fi 0]-[0 0 fi]) - [0 0 n].*BER([0 0 fi]-[0 fi 0]) );
%         Jnip1by2a  = (q.*[0 mun 0]*Vt/(dx*Ldi)) * ni.*( [n 0 0] .*BER([fi 0 0]-[0 fi 0]) - [0 n 0].*BER([0 fi 0]-[fi 0 0]) );
%         Jnim1by2  =[Jnim1by2 ; Jnim1by2a (2: n_max+1)]; 
%         Jnip1by2 =[Jnip1by2; Jnip1by2a (2: n_max+1)];
%         
%         Jelec1= (Jnim1by2a (2: n_max+1)+ Jnip1by2a (2: n_max+1))/2;
%         Jelec = [Jelec; Jelec1 ];
%         
% %       figure (101); semilogy( abs(Jelec1)); hold on;
%  
%         Jpim1by2a = (q*[0 mup 0].*Vt/(dx*Ldi)) * ni.*( [0 p 0].*BER(([ 0 0 fi]-[0 fi 0])) - [0 0 p].*BER(([0 fi 0]-[0 0 fi])) );
%         Jpip1by2a = (q*[0 mup 0]*Vt/(dx*Ldi)) * ni.*( [p 0 0].*BER(([0 fi 0]-[fi 0 0])) - [0 p 0].*BER(([fi  0 0 ]-[0 fi 0])) );
% 
%         Jpim1by2  =[Jpim1by2 ; Jpim1by2a (2: n_max+1)]; 
%         Jpip1by2 =[Jpip1by2; Jpip1by2a (2: n_max+1)];
%         Jhole1=(Jpim1by2a (2: n_max+1) + Jpip1by2a (2: n_max+1))/2;
%         Jhole  = [Jhole ; Jhole1];
% 
% 
% end
% 
% 
% 
% % Write the results of the simulation in files %
% 
% xx1(1) = dx*1e4;
% 
% nf = n*ni;
% pf = p*ni;
% 
%     Ec = dEc- Vt*fi;     %Values from the second Node%
%     ro = ni.*(-n + p + dop);  %ro = -ni*(exp(fi) - exp(-fi) - dop);
%     el_field1 = -([fi 0 0] - [0 fi 0])*Vt/(dx*Ldi);
%     el_field2 = -([fi 0 0] - [ 0 0 fi])*Vt/(2*dx*Ldi);
%     xx1 = x.*1e4/Ldi;
%     
%     el_field1 = el_field1(2: n_max+1);
%     el_field2 = el_field2(2: n_max+1);    
% 
% %Jtotal(:,1) = Jtotal(:,2);
% Jelec(:,1) = Jelec(:,2);
% Jhole(:,1) = Jhole(:,2);
% %Jtotal(:,n_max) = Jtotal(:,(n_max-1));
% Jelec(:,n_max) = Jelec(:,(n_max-1));
% Jhole(:,n_max) = Jhole(:,(n_max-1));
% Jtotal=Jhole+Jelec;
% 
% Ec(1) = Ec(2);
% Ec(n_max) = Ec(n_max-1);
% xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
% el_field1(1) = el_field1(2);
% el_field2(1) = el_field2(2);
% el_field1(n_max) = el_field1(n_max-1);
% el_field2(n_max) = el_field2(n_max-1);
% 
% ro(1) = ro(2);
% ro(n_max) = ro(n_max-1);
% 
% %% Calculate Quasi Fermi Level - Efn Efp
% %for i = 1:n_max
%     Ei   = Ec - 0.56;
%     Efn  = Ei + Vt*log(n);
%     Efp  = Ei - Vt*log(p);
% %end
%     Ev = Ec - 1.12;
% 
% figure(14)
% plot(xx1, Ec,'black','LineWidth',2.5);
% hold on;
% plot(xx1, Ev,'black','LineWidth',2.5);
% hold on;
% plot(xx1, Ei,'--black','LineWidth',2.5);
% hold on;
% plot(xx1, Efn,'r','LineWidth',2.5);
% hold on;
% plot(xx1, Efp,'b','LineWidth',2.5);
% xlabel('x [um]');
% ylabel('Energy [eV]');
% title('Quasi Fermi Levels (Efn & Efp) vs Position - at Applied Bias(0.625V)');
% legend('Ec','Ev','Ei','Efn','Efp');
% hold off;
% %axis([0 7 -1 1]);
% 
% 
% figure(6)
% plot(xx1, Ec,'b','LineWidth',2)
% xlabel('x [um]');
% ylabel('Conduction Band Energy (eV)');
% title('Conduction Band vs Position - at Applied Bias (0.625)');
% 
% figure(7)
% plot(xx1, Vt*fi,'b','LineWidth',2)
% xlabel('x [um]');
% ylabel('Potential [eV]');
% title('Potential vs Position - at Applied Bias(0.625V)');
% 
% figure(8)
% plot(xx1, el_field1,'b','LineWidth',2)
% hold on;
% plot(xx1, el_field2,'b','LineWidth',2)
% xlabel('x [um]');
% ylabel('Electric Field [V/cm]');
% title('Field Profile vs Position - at Applied Bias(0.625V)');
% hold off;
% 
% figure(9)
% %plot(xx1, nf,'g','LineWidth',2)
% semilogy(xx1, ni*abs(dop),'g','LineWidth',2)
% hold on;
% semilogy(xx1, nf,'g','LineWidth',2)
% hold on;
% %plot(xx1, pf,'b','LineWidth',2)
% semilogy(xx1, pf,'b','LineWidth',2)
% xlabel('x [um]');
% ylabel('Electron & Hole Densities [1/cm^3]');
% title('Electron & Hole Densities vs Position - at Applied Bias(0.625V)');
% legend('dop', 'n','p');
% %axis([0 6.75 0 10.2e17])
% hold off;
% 
% figure(10)
% %plot(xx1, ro,'b','LineWidth',2)
% plot(xx1, ro,'b','LineWidth',2)
% xlabel('x [um]');
% %ylabel('Total Charge Density [1/cm^3]');
% ylabel('Total Charge Density [C/cm^3]');
% title('Total Charge Density vs Position - at Applied Bias(0.625V)');
% %axis([0.5 5 -3e17 8e17])
% 
% 
% 
% figure(11)
% semilogy(Vplot, abs(Jtotal(:,round(n_max/2))),'r','LineWidth',2)
% hold on
% semilogy(Vplot, abs(Jhole(:,round(n_max/2))),'g','LineWidth',2)
% hold on
% semilogy(Vplot, abs(Jelec(:,round(n_max/2))),'b','LineWidth',2)
% xlabel('VA [V]');
% ylabel('Total Current Density [Amp/cm^2]');
% title('I vs V Plot');
% legend('Jtotal','Jhole','Jelec','2');
% hold off;
% 
% figure(12)
% semilogy(Vplot, abs(Jtotal(:,round(n_max/2))),'r','LineWidth',2)
% xlabel('VA [V]');
% ylabel('Total Current Density [Amp/cm^2]');
% title('I vs V Plot');
% %legend('Jtotal','Jhole','Jelec','2');
% 
% % figure(13)
% % plot(xx1,Jtotal((round((Total_Steps)-1)),:),'b','LineWidth',2)
% % xlabel('x [um]');
% % ylabel('Total Current Density [A/cm^2]');
% % title('Total Current Density vs Position - at Applied Bias(0.625V)');
% % %axis([0 7 0 6]);
% 
% 
% %figure(5)
% %plot(xx1, n)
% %hold all
% %plot(xx1, p)
% Iplot=abs(Jtotal(:,round(n_max/2)));
% 
% 
% save (strcat(fname,'IV.dat'), 'Vplot',  'Iplot');
% save (strcat(fname,'cond_band.dat'), 'xx1',  'Ec');
% save (strcat(fname,'tot_charge.dat'), 'xx1',  'ro');
% save (strcat(fname,'el_field.dat'), 'xx1',  'el_field1', 'el_field2');
% save (strcat(fname,'np_data.dat'), 'xx1',  'nf', 'pf');
% save (strcat(fname,'pot_1.dat'), 'xx1',  'fi');
% 

toc;