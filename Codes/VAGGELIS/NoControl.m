function unew = NoControl(Params, step, t, xrec, urec)
%% ===================================================================
%  Control the EWH aggregation based on the internal hysterisis controllers
%
%  Evangelos Vrettos
%  October 2014
%  ===================================================================

%% Find current states
if step > 1
    x = xrec(:, step-1, :);
    u = urec(:, step-1, :);
else
    x = Params.x_init;
    u = Params.u_init;
end

unew = nan(size(u));

for i = 1:Params.n_app
    unew(:,:,i) = determine_heating_element_state_population(u(:,:,i), t, x(:,:,i), Params, i);
end

end