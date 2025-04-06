clc
% -------------- Model inputs ---------------------------------------------
tic
CO2_b = 2.68e19;% CO2 concentration (molecules/cm3)
BRR = 250; % Pulse packet generation frequency (s-1)
n_p =1;% Number of pulse packets
tnextp = 200e-9; % Time separation between pulse packets (s) 
t_simulation =5000e-9-665e-9; % Total simulation time (s)
t_i = 100e-12; % Initial time step, cannot be 0(s)
t_step = 50e-12; % Infinitesimally small time step (s)
Ne_i =2.5e16; % Initial electron density (cm-3)
p_n =8; % Number of packet under a single pulse packets
tp =75e-9; % Time separation between pulses in a single pulse packet(s)
rp = 185;% Plasma radius (microns)
Vp = 4/3*3.14*(rp^3)*1e-12;% Plasma volume (cm-3)
tau = 1e-9; % Pulse decay time (s)
%-------------------------Pulse packet generator---------------------------
[t,Ne,t_f,Ne_f] = Pulse_Gen(t_i,n_p,p_n,tp,tnextp,Ne_i,t_step,t_simulation,tau);
t_0 = t;
Ne_0 = Ne;
t_f_i = t_f;
Ne_i = Ne_f;
%------------------------------ ODE solver---------------------------------
s_P = size(P); %
r=ones(s_P(1),1); %rate matrix
dt=t(2)-t(1); % dt
C_all=zeros(size(Ne,2),s_P(2),n_p); % Simulation matrix for each pulse
C = [CO2_b zeros(1,s_P(2)-3) 0 0]; % initial concentration
C_in = C%
C_f =[];
N=n_p;
for i =1:N % N is number of pulses modelled as a batch reactor  
  if i==1 && N==1 
      C_out =[];
      [C_cell_conc]=ODE_solver_chunking(Flag_Ne,C_in,C_out,P,Ne_f,t_f,N,k);
  C_cell = C_cell_conc;
  elseif i==1 && N>1 
        t=t_0(i,:);
        Exp=Ne(i,:);
        dt = t(2)-t(1);
        [C_out] = ODE_solver_Ne_square_batch(P,t,k,Exp,C,Flag_Ne,dt);
        C_cell(:,:,i) = C_out;
        i
    elseif i>1 && i<N 
        t=t_0(i,:);
        Exp=Ne(i,:);
        C=C_out(end,:);
        [C_out] = ODE_solver_Ne_square_batch(P,t,k,Exp,C,Flag_Ne,dt);
        C_cell(:,:,i) = C_out;
        i
    elseif i>1 && i==N
        t=t_0(i,:);
        Exp=Ne(i,:);
        t_end = t_simulation-t(end);
        t_end = t(end):t_step:t_simulation;
        t = [t t_end];
        Exp =[Exp zeros(1,length(t_end))];
        C=C_out(end,:);
        t_f = t;
        Ne_f = Exp;
        [C_cell_conc]=ODE_solver_chunking(Flag_Ne,C_in,C_out,P,Ne_f,t_f,N,k);
        C_cell_f(:,:)= C_cell_conc;       
    end
end
if N==1
    C_all=C_cell;
else
C_all = zeros((N-1)*size(C_cell,1)+size(C_cell_f,1),s_P(2));
for i=1:size(C_cell,3)+1
    index = size(C_cell,3)+1;
    if i ==1 && i<index
    C_all(1:size(C_cell,1),:) = C_cell(:,:,i);
    elseif i>1 && i<index
    C_all(((i-1)*size(C_cell,1)+1:(i-1)*size(C_cell,1)+size(C_cell,1)),:)=C_cell(:,:,i);
    elseif i==index
        C_all(((i-1)*size(C_cell,1)+1:end),:)=C_cell_f(:,:);
    end
end
end
%--------------------Plotting----------------------------------------------   
z = C_all(1:end,2:22);
Excited = sum(z,2);
CO = C_all(:,23)%+C_all(:,30);
CO2 = C_all(:,1);
t = linspace(t_i,t_simulation,size(C_all,1));
plot(t',CO2,t',CO,t',Excited)
ylabel('Concentration (cm^-^3)')
xlabel('Time (s)')
rate = zeros(1,size(C_cell,3)+1);
%-----------------Calculation of the production rate-----------------------
if N==1
    rate= Vp*CO(end,1)*BRR/6.02e23;
else
for i=1:size(C_cell,3)+1
    if i<size(C_cell,3)+1
        rate(i) = Vp*C_cell(end,23,i)*BRR/6.02e23;
    else
        rate(i) = Vp*C_cell_f(end,23)*BRR/6.02e23;
    end
end
end
rate = 1e9*rate'; % nmols/s
toc

        
    

