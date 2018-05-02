% This code simulates Landau damping of Ion Acoustic using
% two species: 
% species 0: electron - FLUID
% species 1: ion - Kinetic

%clear all; close all;

% simulation box length
L=4*pi;
% time step
DT=0.025;
% time decentering parameter: 0.5 or 0
theta = 0.5;
% number of computational cycles
NT=600; 
% number of grid points
NG=64;
% grid spacing
dx=L/NG;
% Number of particles per cell for each species
NPC = 10000; % 2000 fluid
% number of particles: species 0 and and species 1
N0=NG*NPC; 
N1=NG*NPC;
% plasma frequency
WP=1;
% charge to mass ratio
QM0 = -1;
QM1 = 1/1836;
% no drift velocity
V0=0.0;
% thermal velocity
VT=sqrt(1/3);
% perturbation amplitude and mode
XP1=0.1;
V1=0.0;
mode=1;
% charge
Q0=WP^2/(QM0*N0/L);
Q1=-Q0*N0/N1;
% electron mass
M0=Q0/QM0;
% ion mass
M1=Q1/QM1;
% define wheter a particle kinetic or fluid
% setting a flag kfp
% kinetic if kinetic is 1

% fluid if kinetic 0
kfp0 = 0.0*ones(N0,1); %
kfp1 = 0.0*ones(N1,1); % 

% flag to turn on when a particle becomes kinetic
switchF2Kp0 = 0.0*ones(N0,1);
switchF2Kp1 = 0.0*ones(N1,1);
delta_vp0 = 0.0*ones(N0,1);
delta_vp1 = 0.0*ones(N1,1);
% threshold velocity variation to switch to kinetic particle
deltaVthreshold = VT/10;


%initial particle positions
xp0 = linspace(0,L-L/N0,N0)';
xp1 = linspace(0,L-L/N1,N1)';
% Perturbation in particle positions
% Perturbation in particle positions
%xp0=xp0+0.1*sin(2*pi*xp0/L*mode);
xp1=xp1 + 0.2*sin(2*pi*xp1/L*mode);


% This particle velocity initialization is needed for both kinetic
% and fluid since we calculate pressure from therma velocity
% initial particle velocities: e- species 0 and 1
vp0 = V0 + sqrt(abs(QM0)*5)*VT*randn(N0,1); 
vp0_old = vp0;
vp1 = -V0 + VT*randn(N1,1);
vp1_old = vp1;

%VT0 = sqrt(cov(vp0))
%VT1 = sqrt(cov(vp1))
% specific heat ratio for EOS
gamma = 5/3; %5/3

% artificial viscosity
c1 = 10.0;% linear parameter
c2 = 10.0;% quadratic parameter

% Laplacian kernel needed for switching from kinetic to fluid
NSMOOTH = 2;

% auxiliary vectors for fast interpolation
p0=1:N0;p0=[p0 p0];
p1=1:N1;p1=[p1 p1];

% ----------------------
% species 0 (one e- beam)

% calculate mass density
g1=floor(xp0/dx-.5)+1;
g=[g1;g1+1];
fraz1=1-abs(xp0/dx-g1+.5);
fraz=[fraz1;1-fraz1];
out=(g<1);g(out)=g(out)+NG; out=(g>NG);g(out)=g(out)-NG;
mat0=sparse(p0,g,fraz,N0,NG);
RHO0_C = full((M0/dx)*sum(mat0))';
% calculate pressure
fraz=[fraz1.*var(vp0);(1-fraz1).*var(vp0)];
mat0=sparse(p0,g,fraz,N0,NG);
pressure0 = full((M0/dx)*sum(mat0))';
% fix pressure to constant value to avoid noise
pressure0 = pressure0(1)*ones(NG,1);
% internal energy on the grid defined using EOS
I0_C = pressure0./(RHO0_C*(gamma-1));
% calculate particle internal energy
fraz=[fraz1;(1-fraz1)];
mat0=sparse(p0,g,fraz,N0,NG);
ep0 = M0.*mat0*I0_C;

% ----------------------
% species 0 (one e- beam)

% calculate mass density
g1=floor(xp1/dx-.5)+1;
g=[g1;g1+1];
fraz1=1-abs(xp1/dx-g1+.5);
fraz=[fraz1;1-fraz1];
out=(g<1);g(out)=g(out)+NG;
out=(g>NG);g(out)=g(out)-NG;
mat1=sparse(p1,g,fraz,N1,NG);
RHO1_C = full((M1/dx)*sum(mat1))';
% calculate pressure
fraz=[fraz1.*var(vp1);(1-fraz1).*var(vp1)];
mat1=sparse(p1,g,fraz,N1,NG);
pressure1 = full((M1/dx)*sum(mat1))';
% fix pressure to constant value to avoid noise
pressure1 = pressure1(1)*ones(NG,1);

% internal energy on the grid defined using EOS
I1_C = pressure1./(RHO1_C*(gamma-1));
% calculate particle internal energy
fraz=[fraz1;(1-fraz1)];
mat1=sparse(p1,g,fraz,N1,NG);
ep1 = M1.*mat1*I1_C;

% If particles are fluid set velocity to the fluid velocity
vp0(kfp0==0) = V0;
vp1(kfp1==0) = -V0;


% auxiliary vectors and matrices for interpolation and field solver
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);

% arrays for diagnostics
histSpectrum = [];

for it=1:NT
   
    % ----------------------
    % species 0 (first e- beam)
    g1=floor(xp0/dx-.5)+1;
    g=[g1;g1+1];
    fraz1=1-abs(xp0/dx-g1+.5);
    fraz=[fraz1;1-fraz1];
    out=(g<1);g(out)=g(out)+NG;
    out=(g>NG);g(out)=g(out)-NG;
    mat0=sparse(p0,g,fraz,N0,NG);
    % charge density species 0
    rhoq0 = full((Q0/dx)*sum(mat0))'; % charge density
    RHO0_C = rhoq0/QM0; % mass density 
    I0_C = (1./(RHO0_C*dx)).*full(mat0'*ep0);
    
    % calculate the pressure with EOS
    pressure0 = RHO0_C.*I0_C*(gamma - 1);
    for ismooth=1:NSMOOTH
            pressure0 = (0.25*[pressure0(NG); pressure0(1:NG-1)] +  0.5*pressure0 + 0.25*[pressure0(2:NG);pressure0(1)]);
            rhoq0 = (0.25*[rhoq0(NG); rhoq0(1:NG-1)] +  0.5*rhoq0 + 0.25*[rhoq0(2:NG);rhoq0(1)]);
            RHO0_C = (0.25*[RHO0_C(NG); RHO0_C(1:NG-1)] +  0.5*RHO0_C + 0.25*[RHO0_C(2:NG);RHO0_C(1)]);
            I0_C = (0.25*[I0_C(NG); I0_C(1:NG-1)] +  0.5*I0_C + 0.25*[I0_C(2:NG);I0_C(1)]);
    end 
    I0_C_old = I0_C;
    
    
    % Mass and fluid velocity on grid vertices
    g1=floor(xp0/dx)+1;
	g=[g1;g1+1];
	fraz1=1-abs(xp0/dx-g1+1);
	fraz=[fraz1;1-fraz1]	;
	out=(g<1); g(out)=1; fraz(out)=  fraz(out);
	out=(g>(NG+1));g(out)= NG+1; fraz(out)=fraz(out);
	matN=sparse(p0,g,fraz,N0,NG+1);
    % Mass on grid vertices
    M0_V = full(M0*sum(matN))';
    % Periodic Boundary condition
    M0_V(1) = M0_V(1)+ M0_V(NG+1); 
    M0_V(NG+1) = M0_V(1); 
    % calculate the fluid velocity
    % calculate fluid velocity by particle grid interpolation
    U0_V = (1./M0_V).*full(matN'*(M0*vp0));
    % Periodic Boundary condition
    U0_V(1) = U0_V(1) + U0_V(NG+1); % take care of BC
    U0_V(NG+1) = U0_V(1); 
    % Store fluid velocity for taking average later
    U0_V_old = U0_V;
    
    % calculate artificial viscosity
	cs0 =sqrt(gamma*(gamma-1)*I0_C);
    der_U0_C = U0_V(2:NG+1) - U0_V(1:NG);
    der_U0_C(der_U0_C>=0) = 0.0;
    art_viscosity0 = RHO0_C.*(0.25*c2*(gamma+1)*abs(der_U0_C) + sqrt(c2^2*(0.25*(gamma+1))^2*der_U0_C.^2 + c1^2*cs0.^2 )).*abs(der_U0_C);
    
    % ----------------------
    % ----------------------
    % species 1 (second e- beam)
    g1=floor(xp1/dx-.5)+1;
    g=[g1;g1+1];
    fraz1=1-abs(xp1/dx-g1+.5);
    fraz=[fraz1;1-fraz1];
    out=(g<1);g(out)=g(out)+NG;
    out=(g>NG);g(out)=g(out)-NG;
    mat1=sparse(p1,g,fraz,N1,NG);
    % charge density species 1
    rhoq1 = full((Q1/dx)*sum(mat1))'; % charge density
    RHO1_C = rhoq1/QM1; % mass density
    I1_C = (1./(RHO1_C*dx)).*full(mat1'*ep1);
    I1_C_old = I1_C;
    % calculate the pressure with EOS 
    pressure1 = RHO1_C.*I1_C*(gamma - 1);
    
    
    
    % Mass and fluid velocity on grid vertices
    g1=floor(xp1/dx)+1;
	g=[g1;g1+1];
	fraz1=1-abs(xp1/dx-g1+1);
	fraz=[fraz1;1-fraz1]	;
	out=(g<1); g(out)=1; fraz(out)=  fraz(out);
	out=(g>(NG+1));g(out)= NG+1; fraz(out)=fraz(out);
	matN=sparse(p1,g,fraz,N1,NG+1);
    % Mass on grid vertices
    M1_V = full(M1*sum(matN))';
    % Periodic Boundary condition
    M1_V(1) = M1_V(1)+ M1_V(NG+1); %Periodic BC
    M1_V(NG+1) = M1_V(1); %Periodic BC
    % calculate the fluid velocity
    U1_V = (1./M1_V).*full(matN'*(M1*vp1));
    % Periodic Boundary condition
    U1_V(1) = U1_V(1) + U1_V(NG+1);
    U1_V(NG+1) = U1_V(1);
    % Store fluid velocity for taking average later
    U1_V_old = U1_V;
    
    % calculate artificial viscosity
	cs1 =sqrt(gamma*(gamma-1)*I1_C);
    der_U1_C = U1_V(2:NG+1) - U1_V(1:NG);
    der_U1_C(der_U1_C>=0) = 0.0;
    art_viscosity1 = RHO1_C.*(0.25*c2*(gamma+1)*abs(der_U1_C) + sqrt(c2^2*(0.25*(gamma+1))^2*der_U1_C.^2 + c1^2*cs1.^2 )).*abs(der_U1_C);
    
    % ----------------------
    % Solve the field
    % total charge density
    rho = rhoq0 + rhoq1;
    % calculate the electrostatic potential, 0 on the right boundary
    Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
    % calculate the electric field at center cell
    Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
    % calculate the electric field on the node
    E_V = ([Phi(NG); Phi(1:NG)]-[Phi(1:NG);Phi(1)])/dx;
    % Periodic Boundary condition
    E_V(NG+1) = E_V(1);
    
    % ----------------------
    % Update Fluid Equations
    % momentum equation
    U0_V = U0_V - DT*([pressure0(1:NG);pressure0(1)] + [art_viscosity0(1:NG);art_viscosity0(1)] - [pressure0(NG); pressure0(1:NG)] - [art_viscosity0(NG); art_viscosity0(1:NG)])./M0_V + DT*QM0*E_V;
    U1_V = U1_V - DT*([pressure1(1:NG);pressure1(1)] + [art_viscosity1(1:NG);art_viscosity1(1)] - [pressure1(NG); pressure1(1:NG)] - [art_viscosity1(NG); art_viscosity1(1:NG)])./M1_V + DT*QM1*E_V;
    % Periodic boundary conditions
    U0_V(NG+1) = U0_V(1);
    U1_V(NG+1) = U1_V(1);
    % Fluid velocity at time n + 1/2
    U0_av = 0.5*(U0_V + U0_V_old);
    U1_av = 0.5*(U1_V + U1_V_old);
    U0_av_C = 0.5*(U0_av(1:NG) + U0_av(2:NG+1));
    U1_av_C = 0.5*(U1_av(1:NG) + U1_av(2:NG+1));
    % Fluid velocity difference for updating 
    U0_diff = U0_V - U0_V_old;
    U1_diff = U1_V - U1_V_old;
    U0_diff_C = 0.5*(U0_diff(1:NG) + U0_diff(2:NG+1));
    U1_diff_C = 0.5*(U1_diff(1:NG) + U1_diff(2:NG+1));
    
    % energy equation
    I0_C = I0_C - DT*((pressure0 + art_viscosity0).*(U0_av(2:NG+1) - U0_av(1:NG)))./(RHO0_C*dx);
    I1_C = I1_C - DT*((pressure1 + art_viscosity1).*(U1_av(2:NG+1) - U1_av(1:NG)))./(RHO1_C*dx);
     

    % --------------------------------------------
    % KINETIC- FLUID Switch
    % Vthermal on center from pressure - Species 0
    M0_C = full(M0*sum(mat0))';
    VT0 = sqrt(dx*pressure0./M0_C);
    VTp0 = mat0*VT0;
    % Vthermal on center from pressure - Species 1
    M1_C = full(M1*sum(mat1))';
    VT1 = sqrt(dx*pressure1./M1_C);
    VTp1 = mat1*VT1;
    
    % make the particles kinetic
    kfp0(switchF2Kp0==1) = 1; % 1 = kinetic
    kfp1(switchF2Kp1==1) = 1; % 1 = kinetic
    % add thermal nose
    VTp0 = VTp0.*randn(N0,1);
    VTp1 = VTp1.*randn(N1,1);
    vp0(switchF2Kp0==1) = vp0(switchF2Kp0==1) + VTp0(switchF2Kp0==1);
    vp1(switchF2Kp1==1) = vp1(switchF2Kp1==1) + VTp1(switchF2Kp1==1); 
    % turn off the switches
    switchF2Kp0(switchF2Kp0==1) = 0;
    switchF2Kp1(switchF2Kp1==1) = 0;
    
    % print fraction of fluid and kinetic particles
    it
    percentage_kinetic_particles = (sum(kfp0) + sum(kfp1))/(N0 + N1)
    
    % --------------------------------------------
    % advance particle velocity
    % record old particle velocities first  
    % it is important that this is after kinetic-fluid switch
    vp0_old = vp0;
    vp1_old = vp1;
    % Kinetic particles
    acc = mat0*QM0*Eg*DT;
    vp0(kfp0==1) = vp0(kfp0==1)+ acc(kfp0==1);
    acc = mat1*QM1*Eg*DT;
    vp1(kfp1==1) = vp1(kfp1==1)+ acc(kfp1==1);
    % Fluid particles
    % species 0
    acc = mat0*U0_diff_C;
    vp0(kfp0==0) = vp0(kfp0==0) + acc(kfp0==0);
    % check fluid particle fluid velocity
    delta_vp0(kfp0==0) = abs(vp0(kfp0==0) - vp0_old(kfp0==0));
    switchF2Kp0(delta_vp0 > deltaVthreshold) = 1;
    % reset delta vp0
    delta_vp0 = zeros(N0,1);
    
    % Species 1
    acc = mat1*U1_diff_C;
    vp1(kfp1==0) = vp1(kfp1==0) + acc(kfp1==0);
    % check fluid particle fluid velocity
    delta_vp1(kfp1==0) = abs(vp1(kfp1==0) - vp1_old(kfp1==0));
    switchF2Kp1(delta_vp1 > deltaVthreshold) = 1;
    % reset delta vp1
    delta_vp1 = zeros(N1,1);
    
    % --------------------------------------------
    % advance particle internal energy
    % both kinetic and fluid: use fluid advancement
    ep0 = ep0 + M0.*(mat0*(I0_C - I0_C_old));
    ep1 = ep1 + M1.*(mat1*(I1_C - I1_C_old)); 
    
    % --------------------------------------------
    % advance particle position
    xp0=xp0+((1-theta)*vp0_old + theta*vp0)*DT;
    xp1=xp1+((1-theta)*vp1_old + theta*vp1)*DT;
    % apply periodic boundary conditions
    out=(xp0<0); xp0(out)=xp0(out)+L;
    out=(xp0>=L);xp0(out)=xp0(out)-L;
    out=(xp1<0); xp1(out)=xp1(out)+L;
    out=(xp1>=L);xp1(out)=xp1(out)-L;
    
    % Diagnostics
    % movie
     plot(xp1,vp1,'.',xp0,vp0,'.')
     axis([0 L -4 4]);
     xlabel('x_p')
     ylabel('v_p')
    %plot(pressure0);
     pause(0.0001); % This is for plotting
    
    % take the Fourier transform of the electric field on the grid
    NFFT = 2^nextpow2(length(Eg)); % Next power of 2 from length of Eg
    Y = fft(Eg,NFFT)/length(Eg);
    histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];
end

time = linspace(0,NT*DT,NT);
k = 2*pi*(1/(2*dx))*linspace(0,1,NFFT/2+1);
space = linspace(0,L,NG);
semilogy(time,0.2*exp(-0.02367136103*time),'k--',time,histSpectrum(2,:),'k-*')
xlabel('time [1/\omega_{p}]')
ylabel('|E|_{k=1}')
grid on
legend('Linear Theory','Fluid e^-')
set(gca,'FontSize',18)
set(gca, 'FontName', 'Times New Roman')

