function [X,Y] = Apply_linear_model(ss,U,exogenous_input)

A = ss.a;
B = ss.b(:,1);
% if exogenous input present
if(size(ss.b,2) == 2)
    E = ss.b(:,2);
end
C = ss.c;

X = cell(1,length(exogenous_input));
Y = cell(1,length(exogenous_input));

X{1} = zeros(size(A,1),1);
Y{1} = C*X{1} + D;

for ii = 2:length(exogenous_input)
    X{ii} = A*X{ii-1} + B*U(ii-1) + E*exogenous_input(ii);
    Y{ii} = C*X{ii-1};
end

X = cell2mat(X);
Y = cell2mat(Y);

end