function [xnew, Uminew] = updateStates(step, Params, Model, xrec, unew, Umirec)
%% =====================================================================
%  Evolve the EWH states
%
%  Evangelos Vrettos
%  October 2014
%  =====================================================================

%% Find current states
if step > 1
    x = xrec(:, step-1, :);
    Umi = Umirec(:, step-1, :);
else
    x = Params.x_init;
    Umi = Params.Umi_init;
end

xnew = nan(size(x));
xtmp = nan(size(x));
Uminew = nan(size(Umi));

for i = 1:Params.n_app
    xtmp(:,:,i) = Model.Ad(:,:,i) * x(:,:,i) +  Model.Bd(:,:,i) * unew(:,:,i);
    
    [xnew(:,:,i), Uminew(:,:,i)] = buoyCor(xtmp(:,:,i), Params, Umi(:, :, i), i);
end

%% Recycle Bin
% Decide how buoyancy effect is calculated depending on the simulation time step
% U_m_i_rec(:, max([step-1 1]), i)
% if Params.t_sample <= 10
%     [x_new U_m_i] = buoyancy_correction_v1(x(:,:,i), Params, U_m_i_rec(:, max([step-1 1]), i), i);
%     x(:,:,i) = x_new;
% else
%     [x_new U_m_i] = buoyancy_correction_v2(x(:,:,i), Params, 0*U_m_i_rec(:, max([step-1 1]), i), i);
%     x(:,:,i) = x_new;
%     
%     [xnew(:,:,i), Uminew(:,:,i)] = buoyancy_correction_v3SOS(xtmp(:,:,i), Params, 0*Umi(:, :, i), i);
%     [xnew(:,:,i), Uminew(:,:,i)] = buoyancy_correction_v1(xtmp(:,:,i), Params, Umi(:, :, i), i);
% end