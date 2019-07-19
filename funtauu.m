function Out = funtauu(V,Vx)
    temp = (1/3.7)*(exp(-(V+Vx+22)/10.5)+28);
    mask = (V+Vx) < -80;
    temp(mask) =  (1/3.7)*exp((V(mask)+Vx(mask)+467)/66.6); 
    Out = temp;
% if (V+Vx) < -80
%     Out = (1/3.7)*exp((V+Vx+467)/66.6);
% else
%     Out = (1/3.7)*(exp(-(V+Vx+22)/10.5)+28);
% end
end