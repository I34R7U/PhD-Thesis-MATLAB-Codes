function [t,Ne,t_f,Ne_f] = Pulse_Gen(t_i,n_p,p_n,tp,tnextp,Ne_i,t_step,t_simulation,tau);
t_p=t_i:t_step:tp;
Ne_0=Ne_i;
t_exp = [];
Ne_exp = [];
%-------------------------Setting up pulses--------------------------------
for i=1:p_n
    if i ==1
        t_exp =[t_exp t_p];
        Ne_exp = [Ne_exp Ne_0*(exp(-t_p/tau))];
    else
        t_exp = [t_exp  t_p+(i-1)*tp];
        Ne_exp = [Ne_exp Ne_0*(exp(-t_p/tau))];
    end
end
Tx =t_exp(end):t_step:(t_exp(end)+tnextp);
z = zeros(1,length(Tx));
t_exp_f = [t_exp Tx];
Ne_exp_f = [Ne_exp z];
%-------------------------Matching pulse packets w/t-----------------------
t = t_exp_f;
Ne = Ne_exp_f;
t_f = [];
Ne_f =t_f;
for i=1:n_p
    if i==1
    t_f = t_exp_f;
    Ne_f= Ne_exp_f;
    else
    t_f = [t_f (t+t_f(end))];
    Ne_f=[Ne_f Ne];
    end
end
%---------------------Completing Ne until the end of simulation------------
t_end = t_simulation-t_f(end);
t_end = t_f(end):t_step:t_simulation;
t_f = [t_f t_end];
Ne_f =[Ne_f zeros(1,length(t_end))];
plot(t_f,Ne_f)
%end
end