function [C_cell_conc]=ODE_solver_chunking(Flag_Ne,C_in,C_out,P,Ne_f,t_f,N,k)
dvsrs = divisors(length(Ne_f));
N_chunk =[];
for i = 1:length(dvsrs)
    if dvsrs(i) <= 25
        N_chunk = [N_chunk dvsrs(i)]
    else
    end
end
N_chunk = max(N_chunk)
if N_chunk ==1
    Ne_f = Ne_f(1,1:end-1);
    t_f= t_f(1,1:end-1); 
    dvsrs = divisors(length(Ne_f));
    for i = 1:length(dvsrs)
        if dvsrs(i) <= 25
        N_chunk = [N_chunk dvsrs(i)]
    else
    end
    end
N_chunk = max(N_chunk)
else
    N_chunk = N_chunk;
end
r_ound = round(length(Ne_f)/N_chunk);
Ne_f_chunk = zeros(N_chunk,r_ound);
t_f_chunk = Ne_f_chunk;
for i=1:N_chunk
        Ne_f_chunk(i,:) = Ne_f((i-1)*r_ound+1:(i-1)*r_ound+r_ound);
        t_f_chunk(i,:) = t_f((i-1)*r_ound+1:(i-1)*r_ound+r_ound);
end
C_cell_chunk =[];
for i=1:N_chunk
    i
    Ne = Ne_f_chunk(i,:);
    t = t_f_chunk(i,:);
    Exp=Ne;
    dt = t(2)-t(1);
    if i==1 && N ==1
    C = C_in;
    else
    C = C_out(end,:);
    end
    [C_out] = ODE_solver_Ne_square_batch(P,t,k,Exp,C,Flag_Ne,dt);
    C_cell_chunk(:,:,i) = C_out;
end
C_cell_conc = [];
for i = 1:size(C_cell_chunk,3)
C_cell_conc =[C_cell_conc; C_cell_chunk(:,:,i)];
end
end