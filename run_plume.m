% ICE SHELF/TIDEWATER GLACIER PLUME MODEL
% FOR ARBITRARY ICE-OCEAN BOUNDARY GEOMETRY
% 26/08/2020 adding nitrate as a passive tracer
% Donald Slater, uploaded as part of Slater 2022, GRL
function sol = run_plume(zi,xi,Ta,Sa,Na,Q0,alpha);

% INPUTS
% zi: depths at which xi, Ta and Sa are defined
% xi: x-position of ice-ocean boundary at depths zi
% Ta: ocean temperature at depths zi
% Sa: ocean salinity at depths zi
% Na: ocean nitrate at depths zi
% Q0: subglacial discharge (NOTE UNITS m2/s)
% alpha: entrainment coefficient

% USEFUL NOTES
% 1. model assumes ocean surface at z=0 with z negative below surface
% 2. model assumes water emerges at the minimum value in zi
% 3. model assumes the ice-ocean boundary is oriented bottom left to top
% right
% 4. model cannot cope with complex in-and-out geometries; i.e.
% the gradient of zi wrt xi cannot be negative anywhere
% 5. model linearly interpolates the shape of the ice-ocean boundary
% and the ocean conditions between the supplied points

% OUTPUTS
% z - depth vector for all other outputs (m)
% b - plume width (m)
% u - plume velocity (m/s)
% T - plume temperature (degC)
% S - plume salinity (psu)
% mdot - submarine melt rate (m/d)
% total_melt - total plume melt flux (m2/d)

%% PHYSICAL PARAMETERS
par.alpha = alpha;
par.g = 9.81;
par.rho0 = 1020;
par.l1 = -5.73e-2;
par.l2 = 8.32e-2;
par.l3 = 7.53e-4;
par.cw = 3974;
par.ci = 2009;
par.Lm = 334000;
par.Cd = 0.0097;
par.GT = 1.1e-2;
par.GS = 3.1e-4;
par.Ti = -10;
par.betaS = 7.86e-4;
par.betaT = 3.87e-5;
% technical options
par.Gamma0 = 1; % defines source froude number
par.meltdragfeedback = 1; % defines whether to have melt/drag feedback on
par.EoS = 1; % 0 for linear equation of state, 1 for full non-linear

%% INPUTS AND INITIAL CONDITIONS

% a couple of basic input checks
if ~isempty(find(zi>0)), disp('WARNING: z-input should all be <=0'); end
% if ~isempty(find(gradient(zi,xi)<0)), disp('WARNING: gradient of zi wrt xi cannot be negative anywhere'); end
% order input vectors so they start at the grounding line
[zi,sort_ind] = sort(zi);
xi = xi(sort_ind);
Ta = Ta(sort_ind);
Sa = Sa(sort_ind);

% INITIAL CONDITIONS
T0 = par.l2+par.l3*min(zi); % temperature
S0 = 0; % salinity
N0 = 0; % nitrate
% reduced gravity
if par.EoS == 0,
	g0p = par.g*(par.betaS*(Sa(1)-S0)-par.betaT*(Ta(1)-T0));
elseif par.EoS == 1,
	g0p = (par.g/par.rho0)*(rho(Ta(1),Sa(1),0)-rho(T0,S0,0));
end
b0 = (par.alpha*Q0^2*par.Gamma0/g0p)^(1/3);
u0 = Q0/b0;

% MODEL SOLVES IN TERMS OF FLUXES SO NEED INITIAL FLUXES
QFLUX0 = u0*b0;
MFLUX0 = u0^2*b0;
TFLUX0 = b0*u0*T0;
SFLUX0 = b0*u0*S0;
NFLUX0 = b0*u0*N0;

% TRANSFORM DEPTH INPUT INTO AN ALONG-ICE VARIABLE
tantheta = gradient(zi,xi);
li = cumtrapz(zi,sqrt(1+1./tantheta.^2));
% is also useful to have sintheta
for ii = 1:length(tantheta),
    if tantheta(ii) == Inf, sintheta(ii) = 1;
    else sintheta(ii) = tantheta(ii) / sqrt(1+tantheta(ii)^2);
    end
end

%% SOLVE EQUATIONS
% equations defined in "equations_line.m"
clearvars l A
options = odeset('AbsTol',1e-10,'RelTol',1e-5,'Events',@eventfc_line,'refine',4);
[l,A,d1,d2,d3] = ode23(@(L,a) equations_line(L,a,par,zi,li,sintheta,Ta,Sa,Na),[li(1),li(end)],[QFLUX0,MFLUX0,TFLUX0,SFLUX0,NFLUX0],options);

%% BACK OUT SOLUTION
for i=1:length(l),

    sol.z(i) = interp1(li,zi,l(i));
	sol.b(i) = A(i,1)^2/A(i,2);
	sol.u(i) = A(i,2)/A(i,1);
	sol.T(i) = A(i,3)/A(i,1);
	sol.S(i) = A(i,4)/A(i,1);
	sol.N(i) = A(i,5)/A(i,1);

	quad1 = -par.l1*par.cw*par.GT+par.l1*par.ci*par.GS;
	quad2 = par.cw*par.GT*(sol.T(i)-par.l2-par.l3*sol.z(i))+par.GS*(par.ci*(par.l2+par.l3*sol.z(i)-par.l1*sol.S(i)-par.Ti)+par.Lm);
	quad3 = -par.GS*sol.S(i)*(par.ci*(par.l2+par.l3*sol.z(i)-par.Ti)+par.Lm);
	sol.Sb(i) = (-quad2+sqrt(quad2^2-4*quad1*quad3))/(2*quad1);
	sol.Tb(i) = par.l1*sol.Sb(i)+par.l2+par.l3*sol.z(i);
	sol.mdot(i) = 86400*par.cw*par.Cd^(1/2)*par.GT*sol.u(i)*(sol.T(i)-sol.Tb(i))/(par.Lm+par.ci*(sol.Tb(i)-par.Ti));
    
    % density and ambient values
    sol.rho(i) = rho(sol.T(i),sol.S(i),sol.z(i));
    sol.Ta(i) = interp1(zi,Ta,sol.z(i));
    sol.Sa(i) = interp1(zi,Sa,sol.z(i));
    sol.Na(i) = interp1(zi,Na,sol.z(i));
    sol.rhoa(i) = rho(sol.Ta(i),sol.Sa(i),sol.z(i));

end

% neutral buoyancy properties
id = max(find(sol.rho<sol.rhoa));
sol.zNB = sol.z(id);
sol.TNB = sol.T(id);
sol.SNB = sol.S(id);
sol.NNB = sol.N(id);
sol.QNB = sol.b(id)*sol.u(id);
sol.HNB = par.cw*sol.rho(id)*sol.QNB*(sol.TNB-(-2)); % ref to -2

end

%% line plume equations
function diff = equations_line(l,A,par,zi,li,sintheta,Ta,Sa,Na)

    % INITIALISE OUTPUT
    diff=zeros(5,1);
    
    % CALCULATE PLUME VARIABLES
    b = A(1)^2/A(2);
    u = A(2)/A(1);
    T = A(3)/A(1);
    S = A(4)/A(1);
    N = A(5)/A(1);
    
    % CALCULATE REDUCED GRAVITY
    if par.EoS == 0,
	    gp = par.g*(par.betaS*(interp1(li,Sa,l)-S)-par.betaT*(interp1(li,Ta,l)-T));
    elseif par.EoS == 1,
	    gp = (par.g/par.rho0)*(rho(interp1(li,Ta,l),interp1(li,Sa,l),0)-rho(T,S,0));
    end
    
    % NEED REAL TIME VALUE OF Z
    z = interp1(li,zi,l);
    
    % CALCULATE MELT RATE
    if par.meltdragfeedback == 1,
    
	    quad1 = -par.l1*par.cw*par.GT+par.l1*par.ci*par.GS;
	    quad2 = par.cw*par.GT*(T-par.l2-par.l3*z)+par.GS*(par.ci*(par.l2+par.l3*z-par.l1*S-par.Ti)+par.Lm);
	    quad3 = -par.GS*S*(par.ci*(par.l2+par.l3*z-par.Ti)+par.Lm);
	    Sb = (-quad2+sqrt(quad2^2-4*quad1*quad3))/(2*quad1);
	    Tb = par.l1*Sb + par.l2 + par.l3*z;
	    mdot = par.cw*par.Cd^(1/2)*par.GT*u*(T-Tb)/(par.Lm+par.ci*(Tb-par.Ti));
	    meltterm_vol = mdot;
	    meltterm_temp = mdot*Tb - par.Cd^(1/2)*par.GT*u*(T-Tb);
	    meltterm_sal = mdot*Sb - par.Cd^(1/2)*par.GS*u*(S-Sb);
    
	    dragterm = -par.Cd*u^2;
    
    elseif par.meltdragfeedback == 0,
        
	    mdot = 0; Tb=0; Sb=0; dragterm=0; meltterm_vol=0; meltterm_temp=0; meltterm_sal=0;
    
    end
    
    % ENTRAINMENT
    E = par.alpha*interp1(li,sintheta,l);
    
    % DEFINE EQUATIONS
    diff(1) = E*u + meltterm_vol;
    diff(2) = b*gp*interp1(li,sintheta,l) + dragterm;
    diff(3) = E*u*interp1(li,Ta,l) + meltterm_temp;
    diff(4) = E*u*interp1(li,Sa,l) + meltterm_sal;
    diff(5) = E*u*interp1(li,Na,l);

end

%% solver event function
function [b,isterminal,direction] = eventfc_line(t,y)

    % event function for matlab ode solver which stops integration when
    % plume width becomes very large
    % this saves time when the plume approaches max height
    
    plumewidth = y(1)^2/y(2);
    b = plumewidth - 500;  % halt integration when width = 500
    isterminal = 1;        % halt integration 
    direction = 0;         % the zero can be approached from either direction

end

%% non-linear equation of state
%**********************************************************************
% Equation of State of Sea Water At High Pressure
%**********************************************************************
function [density_seawater]=rho(temperature,salinity,depth)

    % temperature in degrees celsius, salinity in ppt, depth in metres
    t= temperature;
    S= salinity;
    rho_0=1027;
    g=9.81;
    P= rho_0*g*abs(depth)*1*10^(-5); % get pressure from depth and convert to bars
    
    %======================================================================
    % Calculating Secant Bulk Modulus
    %======================================================================
    
    kw= 19652.21+ 148.4206*t- 2.327105*t.^2+ 1.360477e-2*(t.^3)-5.155288e-5*(t.^4);
    Aw= 3.239908+ 1.43713e-3*t+ 1.16092e-4*t.^2- 5.77905e-7*t.^3;
    Bw= 8.50935e-5- 6.12293e-6*t + 5.2787e-8*(t.^2);
    k0= kw + (54.6746- 0.603459*t+ 1.09987e-2*(t.^2)- 6.1670e-5*(t.^3)).*S +(7.944e-2 + 1.6483e-2*t- 5.3009e-4*(t.^2)).*(S.^1.5);
    A= Aw+ (2.2838e-3- 1.0981e-5*t- 1.6078e-6*(t.^2)).*S+ 1.91075e-4*(S.^1.5);
    B= Bw+ (-9.9348e-7+ 2.0816e-8*t+ 9.1697e-10*t.^2).*S;
    bulk_modulus= k0+ A*P+ B*P.^2;
    
    %======================================================================
    % One atmoSphere International Equation of State [1980]
    %======================================================================
    
    A= 8.24493e-1- 4.0899e-3*t+ 7.6438e-5*t.^2- 8.2467e-7*t.^3+5.3875e-9*t.^4;
    B= -5.72466e-3 + 1.0227e-4*t- 1.6546e-6*t.^2;
    C= 4.8314e-4;
    rho_w = 999.842594 + 6.793952e-2*t - 9.095290e-3*t.^2 + 1.001685e-4*t.^3 - 1.120083e-6*t.^4 + 6.536336e-9*t.^5;
    rho_zero= rho_w+ A.*S + B.*(S.^1.5)+ C.*(S.^2);
    
    %======================================================================
    % The High Pressure International Equation of State of
    % Seawater,1980
    %======================================================================
    density_seawater = rho_zero./(1-(P./bulk_modulus));

end