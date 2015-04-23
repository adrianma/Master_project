% number of experiments
N_E = 10;

% with water draws
for ii = 1:N_E
    Params{ii} = build_population(1);
    PrelimModel{ii} = precompute_system_model(Params{ii});
    WaterDrawScenarioReal{ii} = build_normalized_draw_scenario(Params{ii});
    % water draws required
    Params{ii}.bWaterDraw = 1;
    Results{ii} = simulate_population(Params{ii},PrelimModel{ii},...
        WaterDrawScenarioReal{ii},WaterDrawScenarioReal{ii},'NoControl');
end
% without water draws
for ii = 1:N_E
    Params_NW{ii} = build_population(1);
    PrelimModel_NW{ii} = precompute_system_model(Params_NW{ii});
    WaterDrawScenarioReal_NW{ii} = build_normalized_draw_scenario(...
        Params_NW{ii});
    % no water draws required
    Params_NW{ii}.bWaterDraw = 0;
    Results_NW{ii} = simulate_population(Params_NW{ii},...
        PrelimModel_NW{ii},WaterDrawScenarioReal_NW{ii},...
        WaterDrawScenarioReal_NW{ii},'NoControl');
end

% select one temperature matrix
AA_NW = Results_NW{randi(10)}.xrec;
% plot the temperatures for that particular experiment
figure;
barh(AA_NW,'DisplayName','AA_{NW}');
grid on;
ylabel('Layer number (1st is inlet, last is ambient)');
xlabel('Temperature [?C]');
title('No water draws, w_{k} = 0 \forall k');

% select one temperature matrix
AA = Results{randi(10)}.xrec;
% plot the temperatures for that particular experiment
figure;
barh(AA,'DisplayName','AA');
grid on;
ylabel('Layer number (1st is inlet, last is ambient)');
xlabel('Temperature [?C]');
title('No water draws, w_{k} \neq 0 \forall k');