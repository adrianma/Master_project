function et = turbulence(Nsl,Nslt,et_in)

A=(et_in-1)/(1-(1/Nslt));
B=et_in-A;
et=(A/Nsl)+B;

% a=(1/(Nslt-1))*log(et_in+1);
% et=(et_in+1)-exp(a*(Nsl-1));
