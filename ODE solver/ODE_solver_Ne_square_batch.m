function [C_out,P,t,k,Exp,Flag_Ne,dt] = ODE_solver_Ne_square_batch(P,t,k,Exp,C,Flag_Ne,dt)
s_P = size(P);
         for i=1:length(t)
                 for l=1:s_P(1) %number of reactions (Ne_d*exp(-t(i)/40e-9))
                     r(l,1) = k(l,1)*Exp(i).^Flag_Ne(l); %Setting up rate constants/(Ne_d*exp(-t(i)/40e-9))
                     %r(l,1) = k(l,1)*Ne_d.^Flag_Ne(l)
                     for j=1:s_P(2) %numbe14r of species
                         if P(l,j) < 0 
                             r(l,1) = r(l,1).*C(end,j).^(-P(l,j)); %this is the final rate expression
                         end
                     end
                 end
             dr=sum(r.*P);
             dC = C(end,:) + dt*dr;
             if dC >=0
                 C(end+1,:) = dC;
             end
            C_out =C;
            end
end
 
