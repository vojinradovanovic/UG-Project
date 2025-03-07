function sol = vitreous_model_v5(args)
arguments
    args.hasBase logical = true
    args.hasStalk logical = true
    args.hasILM logical = true
    args.hasPlot logical = true
    args.hasPrint logical = true
    args.forceNo int32 {mustBeInRange(args.forceNo, 1, 8)} = 1
    args.areaNo int32 {mustBeInRange(args.areaNo, 1, 4)} = 2
    args.relTol double = 10^(-7)
    args.absTol double = 10^(-7)
    args.angleBist double = 0
    args.eventFunc = []
    args.tSol double = 0:0.05:30
    args.commentNo double = 2
    args.scalePeriod double = 1
    args.onlyAngle logical = false
end

%% Defining variables
% disambiguation
% -------------------------
% MCC: Muller cell cone
% 
% Physical units in calculations:
% length            mm
% area              mm^2
% force             mN
% pressure          mN/mm^2 = 10^3 N/m^2 = kPa:
%                       1 Pa = 1 N/m^2 = 10^(-3) mN/mm^2
% angle             deg
% time              days:
%                       1 s = 1/3600 h = 1/86400 days


% defining constants
const.zeta = 10/86400; % d

const.vitr.stiffnessGel = 1 * 10^(-3); % kPa
const.vitr.viscosityGel = const.zeta * const.vitr.stiffnessGel; % kPa*d

const.vitr.stiffnessLiquid = 0.1 * 10^(-3); % kPa
const.vitr.viscosityLiquid = const.zeta * const.vitr.stiffnessLiquid; % kPa*d

const.vitr.l = 15.5; % mm

const.MCC.lStalk = 65 * 10^(-3); % mm
const.MCC.stiffnessInnerProcess = mean([100, 150, 200]) * 10^(-3); % kPa
const.MCC.viscosityInnerProcess = mean([50/30, 100/100, 150/200]) * 10^(-3) ...
                                        / 86400; % kPa*d
const.MCC.stiffnessSoma = mean([200, 500, 500]) * 10^(-3); % kPa
const.MCC.viscositySoma = mean([100/30, 450/100, 450/200]) * 10^(-3) ...
                                        / 86400; % kPa*d
const.MCC.stiffnessEndfoot = mean([150, 200, 250]) * 10^(-3); % kPa
const.MCC.viscosityEndfoot = mean([50/30, 100/100, 200/200]) * 10^(-3) ...
                                        / 86400; % kPa*d
const.MCC.stiffness = mean([const.MCC.stiffnessInnerProcess; ...
                            const.MCC.stiffnessSoma; ...
                            const.MCC.stiffnessEndfoot]); % kPa
const.MCC.viscosity = mean([const.MCC.viscosityInnerProcess; ...
                            const.MCC.viscositySoma; ...
                            const.MCC.viscosityEndfoot]); % kPa*d
const.MCC.angleWall = 11.5; % degree
const.MCC.rBase = 100 * 10^(-3); % mm
const.MCC.lBase = const.MCC.rBase / cosd(const.MCC.angleWall);
const.MCC.thickBase = 11 * 10^(-3); %mm
const.MCC.widthBase = const.MCC.rBase * pi / 2; % mm
const.MCC.areaBase = const.MCC.thickBase * const.MCC.widthBase; % mm^2
const.MCC.rCell = 2.2 / 2 * 10^(-3); % mm
const.MCC.areaStalk = 30 * const.MCC.rCell^2 * pi; % mm^2

const.ILM.stiffness = 227;      % kPa
const.ILM.viscosity = 0;        % kPa
const.ILM.l = const.MCC.lBase;  % mm
const.ILM.thick = 20 * 10^(-6); % mm
const.ILM.area = const.ILM.thick * const.MCC.widthBase; % mm^2



% defining variables
sol.t = sym("t");       % symbolic time variable
sol.vars = sym("vars"); % column vector of all stresses and strains
sol.eqns = [];          % column vector of all equations
sol.tSol = args.tSol;   % specify evaluation points of solution

% writing vitreous models in matrix form
commentsHomogeneous = [1 "damper" "gel" 0;
                        2 "spring" "gel" 0;
                        3 "parallel" 1 2];

commentsPartialPVD = [1 "damper" "gel" 0;
                        2 "spring" "gel" 0;
                        3 "parallel" 1 2;
                        4 "damper" "liquid" 0;
                        5 "spring" "liquid" 0;
                        6 "series" 4 5;
                        7 "spring" "liquid" 0;
                        8 "parallel" 6 7;
                        9 "damper" "gel" 0;
                        10 "spring" "gel" 0;
                        11 "parallel" 9 10;
                        12 "parallel" 8 11;
                        13 "series" 3 12;
                        14 "damper" "liquid" 0;
                        15 "spring" "liquid" 0;
                        16 "series" 14 15;
                        17 "spring" "liquid" 0;
                        18 "parallel" 16 17;
                        19 "series" 13 18];

% check if comments syntax is correct
checkComments(commentsHomogeneous);
checkComments(commentsPartialPVD);



%% Setting parameters to specify model in current code iteration
% define force signal and cross-sectional area of VMT

forceTypes = [sin(sol.t*args.scalePeriod);
                sol.t/30;
                (sol.t/30)^2;
                (1-exp(-2*sol.t))/(1-exp(-60));
                -(-exp(-4*sol.t)+2/(1+exp(2*sol.t)))*10^(-2);
                +(-exp(-4*sol.t)+2/(1+exp(2*sol.t)))*10^(-2);
                exp(-((sol.t-15)/1)^2)/10^3*2;
                0]                                       * 10^(-5); %mN

vmtRadius = [80; (80+280)/2; 280; 600] * 10^(-3) / 2; % mm
areas = vmtRadius.^2 * pi; % mm^2

const.F(sol.t) = forceTypes(args.forceNo); % mN
const.A = areas(args.areaNo); % mm^2

% define integration tolerance
const.relTol = args.relTol;
const.absTol = args.absTol;
const.rounding = round(abs(log10(const.absTol))); % disp rounding factor

% print current model setup
if args.commentNo == 1
    comments = commentsHomogeneous;
    disp("HOMOGENEOUS");
else
    comments = commentsPartialPVD;
    disp("PARTIAL PVD");
end
disp("Force: "+string(const.F))
disp("Area: "+string(const.A))

% define vitreous and MCC model
hasBase = args.hasBase;
hasStalk = args.hasStalk;
hasILM = args.hasILM;
const.simpleBase = false;    % used for testing; should be false unless testing
const.MCC.lStalk = const.MCC.lStalk - hasBase*const.MCC.thickBase;

% constants for plotting
const.xlim = [sol.tSol(1) sol.tSol(end)];

% additional data modification
const.r = 0;
const.MCC.rBase = const.MCC.rBase - const.r;
const.MCC.lBase = const.MCC.lBase - const.r;
const.ILM.l = const.MCC.lBase;

% other constants
angleBist = args.angleBist;
eventFunc = args.eventFunc;



%% Running model
% clear last warn
lastwarn('');

% set up desired model
[sol.eqns, sol.vars, sol.index] = setup(sol.t, const, comments, ...
                                            hasBase, hasStalk, hasILM);

% adding initial signal and converting sol.eqns to sym
sol.eqns = [sol.eqns;
            {sol.vars(sol.index.VMT, 1) == const.F/const.A}]; % input signal

% simplifying system
[sol.newEqns, sol.newVars, sol.R] = reduceRedundancies(sol.eqns, sol.vars);

% save values in sol for easier access
sol.cVars = sol.R.constantVariables;
sol.rVars = sol.R.replacedVariables;
sol.newVarNo = size(sol.newVars, 1);

% integrate if there are values left over after simplification
if size(sol.newEqns, 1) > 0
    % converting from symbolic to matlab function
    sol.f = daeFunction(sol.newEqns, sol.newVars);

    % guessing initial conditions
    sol.initCondGuess = zeros(sol.newVarNo, 1);
    sol.initCondDerGuess = zeros(sol.newVarNo, 1);

    % set options for ode solver
    if isempty(eventFunc)
        sol.opt = odeset('RelTol',const.relTol,'AbsTol',const.absTol);
    else
        sol.opt = odeset('RelTol',const.relTol,'AbsTol',const.absTol,...
                    'Events',eventFunc);
    end

    % adding initial conditions for desired setup
    if angleBist ~= 0
        cond = sol.vars(end, 1) == angleBist;

        % add angle initial conditions
        sol.initCondGuess = addInitCond(sol.initCondGuess, cond, sol);

    elseif hasBase
        % calculate strain initial condition
        rBase = const.MCC.rBase;
        lBase = const.MCC.lBase;
        
        % write condition; index used is index.VBA-1 to account for
        % possible ILM presence
        baseStart = sol.index.Vitreous+1;
        baseEnd = sol.index.VBA-1;
        cond = sol.vars(baseEnd, 2) == (rBase-lBase)/lBase;

        % add base strain initial conditions
        sol.initCondGuess = addInitCond(sol.initCondGuess, cond, sol);



        % select base (+ILM potentially) equations and calculate
        % stress initial condition; base starts after Vitreous, which has
        % 2*index.Vitreous-1 equations, and ends before VBA starts, which
        % has 3 equations; therefore, the equations to consider lie in:
        % 2*index.Vitreous-1+1 : (2*index.VBA-3)-1
        %
        % in terms of baseStart and baseEnd, this is equivalent to:
        % 2*baseStart-2 : 2*baseEnd-2
        baseEqns = sol.eqns(2*baseStart-2:2*baseEnd-2);

        stressCond = getBaseCond(sol.vars(baseStart:baseEnd, :), ...
                            baseEqns, cond, sol, const);

        % add base strain initial conditions
        sol.initCondGuess = addInitCond(sol.initCondGuess, stressCond, sol);
    end

    
    % checking initial conditions
    [sol.initCond, sol.initCondDer] = decic(sol.f,0,sol.initCondGuess, ...
                                        sol.initCondGuess~=0,...
                                        sol.initCondDerGuess,[], sol.opt);
    
    % print initial conditions if hasPrint
    if args.hasPrint
        fprintf("Initial conditions are: \nvar(t)\tvar(0)\tvar'(0)\n")
        disp(vpa([sol.newVars sol.initCond sol.initCondDer], const.rounding))
    end
    
    % solving the system and measuring required time
    tic
    [sol.tSol,sol.newVarSols] = ode15i(sol.f,sol.tSol,sol.initCond,...
                                    sol.initCondDer,sol.opt);
    toc

    % stop running in case of warning
    if ~isempty(lastwarn)
        error("Model convergence failure");
    end
end



%% Obtaining stress and strain data
% store names of fields in struct sol.index
sol.indexFields = fieldnames(sol.index);
sol.varSols = [];

% loop over fields of sol.index and get the respective solutions
for i = 1:size(sol.indexFields, 1)
    field = sol.indexFields{i};
    fieldVal = getfield(sol.index, field);
    
    % use different index to enlarge varSols
    j = size(sol.varSols, 2);

    if fieldVal ~= 0
        fprintf(field+" index: "+fieldVal+"\n");

        if field == "Angle"
            if args.hasPrint
                fprintf(field+": ");
            end
            sol.varSols(:, j+1) = getVarSol(sol.vars(fieldVal, 1), sol, const, args.hasPrint);

        elseif ~args.onlyAngle
            if args.hasPrint
                fprintf(field+" stress: ");
            end
            sol.varSols(:, j+1) = getVarSol(sol.vars(fieldVal, 1), sol, const, args.hasPrint);
            
            if args.hasPrint
                fprintf(field+" strain: ");
            end
            sol.varSols(:, j+2) = getVarSol(sol.vars(fieldVal, 2), sol, const, args.hasPrint);

            % if ILM, print what it's connected to
            if field=="ILM"
                if hasBase
                    field = "Base+ILM";
                elseif hasStalk
                    field = "MCC+ILM";
                end
    
                if args.hasPrint
                    fprintf(field+" stress: ");
                end
                sol.varSols(:, j+3) = getVarSol(sol.vars(fieldVal+1, 1), sol, const, args.hasPrint);
                
                if args.hasPrint
                    fprintf(field+" strain: ");
                end
                sol.varSols(:, j+4) = getVarSol(sol.vars(fieldVal+1, 2), sol, const, args.hasPrint);
            end
        end

        fprintf("\n");
    end
end

%% Additional info/tests
% % % get bistability angles from base_strain==0
% % eqc = getEquivCond(sol.vars(sol.index.MCC, 1)==1, sol)
%   [ss, dd] = getVarFormula(sol.vars(sol.index.Angle, 1), sol)
% % 
%   [angleF, angleFV] = getVarFormula(sol.vars(sol.index.Angle, 1), sol);
% % % [baseF, baseFV] = getVarFormula(sol.vars(sol.index.Base, 2), sol);
% % % [vitrF, vitrFV] = getVarFormula(sol.vars(sol.index.Vitreous, 1), sol);
% % % [MCCF, MCCFV] = getVarFormula(sol.vars(sol.index.MCC, 2), sol);
% % 
% disp("Angle condition:")
%   vpa(subs(angleF, lhs(eqc), rhs(eqc)), const.rounding)
% % % vpa(subs(baseF, baseFV, struct2cell(eqc)), const.rounding)
% % % vpa(subs(vitrF, vitrFV, struct2cell(eqc)), const.rounding)
% % % vpa(subs(MCCF, MCCFV, struct2cell(eqc)), const.rounding)
% 

% if hasBase
%     fprintf("a-b component: ")
%     % MCC
%     k1 = mean(sol.varSols(2:end, 1)./sol.varSols(2:end, 2)) ...
%             * const.MCC.areaStalk / const.MCC.lStalk;
%     % Vitreous
%     k2 = mean(sol.varSols(2:end, 3)./sol.varSols(2:end, 4)) ...
%             * const.A / const.vitr.l;
%     % Base (+ILM)
%     k3 = mean(sol.varSols(2:end, 5+4*hasILM)./sol.varSols(2:end, 6+4*hasILM)) ...
%             * const.MCC.areaBase / const.MCC.rBase * cosd(const.MCC.angleWall);
% 
%     rezz = 2*k3/((k1+k2+2*k3)*cosd(const.MCC.angleWall)) - 1
% 
%     if rezz>=0
%         disp("In mm:")
%         deltaH = sqrt(const.MCC.rBase * rezz * (rezz+2))
%         disp("In micro-m:")
%         deltaH * 10^3
%     end
% end



%% Plotting results
if args.hasPlot
    % get number of rows in plot
    plotRowNo = (size(sol.varSols, 2) - 1*(sol.index.Angle ~= 0))/2 - 1*hasILM;
    
    % create new figure and other specifspications
    fig = figure('Position', [400 410-70*plotRowNo 750 140*plotRowNo]);
    set(0,'defaultTextInterpreter','latex');
    colours = colororder("gem");
    
    % set up subplot row counter
    rowCounter = 1;
    
    % plotting time evolution
    for i = 1:size(sol.indexFields, 1)
        field = sol.indexFields{i};
    
        if field == "Angle"
            % plot MCC base angle
            subplot(plotRowNo, 2, [2*rowCounter-1, 2*rowCounter]);
            plot(sol.tSol, sol.varSols(:, end), 'LineWidth', 2, ...
                    'Color', colours(i, :));
            
            letter = char('A'+2*rowCounter-2);
            title("\textbf{"+letter+") Angle of MCC Base}");
            xlabel("\boldmath$t$ \textbf{\textbf{(days)}}");
            ylabel("\boldmath$\alpha$ $(^{\circ})$");
            xlim(const.xlim);
    
            grid on;
    
            rowCounter = rowCounter+1;
    
        elseif field ~= "VMT" && getfield(sol.index, field) ~= 0
            % plot variable stress
            subplot(plotRowNo, 2, 2*rowCounter-1);
            plot(sol.tSol, sol.varSols(:, 2*rowCounter-1+2*(field=="VBA")), 'LineWidth', 2, ...
                    'Color', colours(i, :));
            
            ax = gca;
            yExp = ax.YAxis.Exponent;
            ax.YTickMode = "manual";
            ax.YTickLabelMode = "manual";
            
            letter = char('A'+2*rowCounter-2);
            title("\textbf{"+letter+") "+field+" Stress}");
            xlabel("\boldmath$t$ \textbf{(days)}");
            if yExp == 0
                ylabel("\boldmath$\sigma$ \textbf{(kPa)}");
            else
                ylabel("\boldmath$\sigma$ $(10^{"+string(yExp)+"}$ \textbf{kPa)}");
            end
    
            xlim(const.xlim);
            
            grid on;
    
            % plot variable strain
            subplot(plotRowNo, 2, 2*rowCounter);
            plot(sol.tSol, sol.varSols(:, 2*rowCounter+2*(field=="VBA")), 'LineWidth', 2, ...
                    'Color', colours(i, :));
            
            ax = gca;
            yExp = ax.YAxis.Exponent;
            ax.YTickMode = "manual";
            ax.YTickLabelMode = "manual";
            
            letter = char(letter+1);
            title("\textbf{"+letter+") "+field+" Strain}");
            xlabel("\boldmath$t$ \textbf{(days)}");
            if yExp==0
                ylabel("\boldmath$\varepsilon$");
            else
                ylabel("\boldmath$\varepsilon$ $(10^{"+string(yExp)+"})$");
            end
    
            xlim(const.xlim);
            
            grid on;
    
            rowCounter = rowCounter + 1;
        end
    end
    
    % % plotting hysteresis graphs
    % fig = figure('Position', [400 410-70*plotRowNo 750 140*plotRowNo]);
    % rowCounter = 1;
    % 
    % for i = 1:size(sol.indexFields, 1)
    %     field = sol.indexFields{i};
    % 
    %     if field ~= "Angle" && getfield(sol.index, field) ~= 0
    %         % plot variable stress
    %         subplot(plotRowNo, 1, rowCounter);
    %         plot(sol.varSols(:, 2*rowCounter), sol.varSols(:, 2*rowCounter-1), ...
    %                 'LineWidth', 2, 'Color', colours(i, :));
    % 
    %         ax = gca;
    %         xExp = ax.XAxis.Exponent;
    %         yExp = ax.YAxis.Exponent;
    %         % ax.XTickMode = "manual";
    %         % ax.XTickLabelMode = "manual";
    %         % ax.YTickMode = "manual";
    %         % ax.YTickLabelMode = "manual";
    %         ax.XAxisLocation = "origin";
    %         ax.YAxisLocation = "origin";
    % 
    %         letter = char('A'+rowCounter-1);
    %         title("\textbf{"+letter+") "+field+" Stress-Strain graph}");
    % 
    %         if xExp==0
    %             xlabel("\boldmath$\varepsilon$");
    %         else
    %             xlabel("\boldmath$\varepsilon\times10^{"+string(xExp)+"}$");
    %         end
    % 
    %         if yExp+3 == 0
    %             ylabel("\boldmath$\sigma$ \textbf{(Pa)}");
    %         else
    %             ylabel("\boldmath$\sigma\times10^{"+string(yExp+3)+"}$ \textbf{(Pa)}");
    %         end
    % 
    %         grid on;
    % 
    %         rowCounter = rowCounter+1;
    %     end
    % end
end



end


%% Functions
function eqn = spring(var, stiffness)
    % SPRING set up a spring:
    %   var is an array of the form [stress, strain]
    %   stiffness is the stiffness constant

    syms stress strain
    stress = var(1);
    strain = var(2);

    eqn = stress == stiffness * strain;
end

function eqn = damper(var, viscosity)
    % DAMPER set up a damper:
    %   var is an array of the form [stress, strain]
    %   viscosity is the viscosity constant 

    syms stress strain
    stress = var(1);
    strain = var(2);

    eqn = stress == viscosity * diff(strain);
end

function eqns = parallel(vars)
    % PARALLEL set up a parallel bond:
    %   vars contains three elements of the form [stress, strain],
    %       where first two elements are connected in parallel
    %       and the total stress and strain are stored in the third

    syms stresses strains
    stresses = vars(:, 1);
    strains = vars(:, 2);

    eqns = [stresses(3) == stresses(1) + stresses(2);
            strains(3) == strains(1);
            strains(3) == strains(2)];
end

function eqns = series(vars)
    % SERIES set up a series bond:
    %   vars contains three elements of the form [stress, strain],
    %       where first two elements are connected in series
    %       and the total stress and strain are stored in the third

    syms stresses strains
    stresses = vars(:, 1);
    strains = vars(:, 2);

    eqns = [strains(3) == strains(1) + strains(2);
            stresses(3) == stresses(1);
            stresses(3) == stresses(2)];
end

function eqns = vmtForceDisp(vars)
    % VMTFORCEDISP set up a "VMT" bond via input forces and displacements,
    % where forces and displacements are additive and the total
    % displacement is 0:
    %   vars contains three elements of the form [forces, displacements]

    syms forces displacements

    forces = vars(:, 1);
    displacements = vars(:, 2);

    eqns = [displacements(3) == displacements(1) + displacements(2);
            displacements(3) == 0;
            forces(3) == forces(1) + forces(2)];
end

function eqns = vbaForceDisp(angleBase, vars, const)
    % VBAFORCEDISP set up a "VBA" bond with the Muller cell cone base at an
    % angle via input forces and displacements:
    %   vars contains three elements of the form [forces, displacements];
    %   vars(1, :) is the vitreous and vars(2, :) is the base

    syms forces displacements
    forces = vars(:, 1);
    displacements = vars(:, 2);

    % calculate length of base spring
    lBase = const.MCC.lBase;

    eqns = [forces(3) == forces(1) + 2*forces(2)*sind(angleBase);
            displacements(3) == displacements(1);
            displacements(2) == sqrt(displacements(1)^2 + const.MCC.rBase^2) - lBase];
end

function eqns = addEqns(vars, const, comments, index, hasILM)
    % ADDEQNS creates a system of equations from given parameters:
    %   vars - model variables
    %   const - model constants
    %   comments - instructions to set up the model
    %   index - various helpful indices
    %   hasILM - check if ILM is supposed to be included

    varNo = size(vars, 1);
    eqns  = [];

    for i = 1:varNo
        switch comments(i, 2)
            case "damper"
                if comments(i, 3) == "gel"
                    eqns = [eqns; damper(vars(i, :), const.vitr.viscosityGel)];
                elseif comments(i, 3) == "liquid"
                    eqns = [eqns; damper(vars(i, :), const.vitr.viscosityLiquid)];
                elseif comments(i, 3) == "MCC"
                    eqns = [eqns; damper(vars(i, :), const.MCC.viscosity)];
                elseif comments(i, 3) == "MCC IP"
                    eqns = [eqns; damper(vars(i, :), const.MCC.viscosityInnerProcess)];
                elseif comments(i, 3) == "MCC EF"
                    eqns = [eqns; damper(vars(i, :), const.MCC.viscosityEndfoot)];
                elseif comments(i, 3) == "MCC S"
                    eqns = [eqns; damper(vars(i, :), const.MCC.viscositySoma)];
                elseif comments(i, 3) == "ILM"
                    eqns = [eqns; spring(vars(i, :), const.ILM.viscosity)];
                end

            case "spring"
                if comments(i, 3) == "gel"
                    eqns = [eqns; spring(vars(i, :), const.vitr.stiffnessGel)];
                elseif comments(i, 3) == "liquid"
                    eqns = [eqns; spring(vars(i, :), const.vitr.stiffnessLiquid)];
                elseif comments(i, 3) == "MCC"
                    eqns = [eqns; spring(vars(i, :), const.MCC.stiffness)];
                elseif comments(i, 3) == "MCC IP"
                    eqns = [eqns; spring(vars(i, :), const.MCC.stiffnessInnerProcess)];
                elseif comments(i, 3) == "MCC EF"
                    eqns = [eqns; spring(vars(i, :), const.MCC.stiffnessEndfoot)];
                elseif comments(i, 3) == "MCC S"
                    eqns = [eqns; spring(vars(i, :), const.MCC.stiffnessSoma)];
                elseif comments(i, 3) == "ILM"
                    eqns = [eqns; spring(vars(i, :), const.ILM.stiffness)];
                end

            case "parallel"
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];
                eqns = [eqns; parallel(vars(indices, :))];

            case "series"
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];
                
                % If Base then add them via force-disp, otherwise normal
                % stress-strain series bond
                if i == index.Base
                    lengths = [const.MCC.lBase/2;
                                const.MCC.lBase/2;
                                const.MCC.lBase];
                    areas = ones(3, 1) * const.MCC.areaBase;

                    eqns = [eqns; series(vars(indices, :) .* [areas, lengths])];

                else
                    eqns = [eqns; series(vars(indices, :))];
                end

            case "parallel ILM"
                % simple parallel bond with force and disp
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];

                % create input force and displacement sym vectors
                syms forces displacements [3 1]
    
                % create lengths and areas column vector; no signs as all
                % elements are treated as positive
                lengths = [];
                areas = [];
                for j=1:3
                    switch indices(j)
                        case index.Base
                            % calculate length and cross-sectional area of
                            % MCC base
                            lengths(j, 1) = const.MCC.lBase;
                            areas(j, 1) = const.MCC.areaBase;
    
                        case index.ILM
                            % calculate length and cross-sectional area of
                            % ILM
                            lengths(j, 1) = const.ILM.l;
                            areas(j, 1) = const.ILM.area;
    
                        case index.ILM+1
                            lengths(j, 1) = const.MCC.lBase;
                            areas(j, 1) = const.MCC.areaBase+const.ILM.area;

                        otherwise
                            error("Error in comment line "+i+": "+ ...
                                    "Check indices being connected.");
                    end
                end

                % calculate input forces and displacements; input forces
                % are multiplied by signs to adjust for their direction
                forces = vars(indices, 1) .* areas;
                displacements = vars(indices, 2) .* lengths;
    
                % add new equations
                eqns = [eqns; parallel([forces displacements])];

            case "series ILM"
                % simple series bond with force and disp
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];

                % create input force and displacement sym vectors
                syms forces displacements [3 1]
    
                % create lengths and areas column vector; no signs as all
                % elements are treated as positive
                lengths = [];
                areas = [];
                for j=1:3
                    switch indices(j)
                        case index.MCC
                            % calculate length and cross-sectional area of
                            % MCC stalk
                            lengths(j, 1) = const.MCC.lStalk;
                            areas(j, 1) = const.MCC.areaStalk;
    
                        case index.ILM
                            % calculate length and cross-sectional area of
                            % ILM
                            lengths(j, 1) = const.ILM.thick;
                            areas(j, 1) = const.A;
    
                        case index.ILM+1
                            lengths(j, 1) = const.MCC.lStalk+const.ILM.thick;
                            areas(j, 1) = const.A;

                        otherwise
                            error("Error in comment line "+i+": "+ ...
                                    "Check indices being connected.");
                    end
                end

                % calculate input forces and displacements; input forces
                % are multiplied by signs to adjust for their direction
                forces = vars(indices, 1) .* areas;
                displacements = vars(indices, 2) .* lengths;
    
                % add new equations
                eqns = [eqns; series([forces displacements])];

            case "VBA"
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];

                % create input force and displacement sym vectors
                syms forces displacements [3 1]
    
                % create lengths, areas and signs column vectors
                lengths = [];
                areas = [];
                signs = [];
                for j=1:3
                    switch indices(j)
                        case index.Base
                            % calculate length and cross-sectional area of 
                            % MCC base
                            lengths(j, 1) = const.MCC.lBase;
                            areas(j, 1) = const.MCC.areaBase;
    
                            % Base deformation follows the input force
                            % positively, their sign is the same
                            signs(j, 1) = 1;
    
                        case index.Vitreous
                            % calculate length and cross-sectional area of 
                            % vitreous
                            lengths(j, 1) = const.vitr.l;
                            areas(j, 1) = const.A;
    
                            % the vitreous shrinks in response to a positive
                            % input force; the sign is negative
                            signs(j, 1) = -1;
    
                        case index.VBA
                            % the element combining the vitreous and the 
                            % MCC base is imagined as a spring with 
                            % length const.vitr.l and 
                            % cross-sectional area const.A
                            lengths(j, 1) = const.vitr.l;
                            areas(j, 1) = const.A;
    
                            % VBA is taken to be of the same orientation as
                            % the vitreous, so the sign is negative
                            signs(j, 1) = -1;

                        case hasILM*(index.ILM+1)
                            % ILM and base have already been connected;
                            % this is that element
                            lengths(j, 1) = const.MCC.lBase;
                            areas(j, 1) = const.MCC.areaBase+const.ILM.area;
    
                            % Base+ILM deformation follows the input force
                            % positively, their sign is the same
                            signs(j, 1) = 1;

                        otherwise
                            error("Error in comment line "+i+": "+ ...
                                    "Check indices being connected.");
                    end
                end

                % calculate input forces and displacements; input forces
                % are multiplied by signs to adjust for their direction
                forces = vars(indices, 1) .* areas .* signs;
                displacements = vars(indices, 2) .* lengths;
    
                % memorise angleBase variable
                angleBase = vars(index.Angle, 1);
    
                % calculate new equations
                eqns = [eqns;
                        vbaForceDisp(angleBase, [forces displacements], const)];

            case "VMT"
                indices = [str2double([comments(i, 3), comments(i, 4)]), i];

                % create input force and displacement sym vectors
                syms forces displacements [3 1]
    
                % create lengths, areas and signs column vectors
                lengths = [];
                areas = [];
                signs = [];
                for j=1:3
                    switch indices(j)
                        case index.MCC
                            % calculate length and cross-sectional area of 
                            % MCC stalk
                            lengths(j, 1) = const.MCC.lStalk;
                            areas(j, 1) = const.MCC.areaStalk;
    
                            % MCC expands for a positive input
                            signs(j, 1) = 1;
    
                        case {index.Vitreous, index.VBA}
                            % if MCC is connected to the vitreous or 
                            % to the VBA
                            lengths(j, 1) = const.vitr.l;
                            areas(j, 1) = const.A;
    
                            % the vitreous/VBA contracts for 
                            % a positive input
                            signs(j, 1) = -1;
    
                        case index.VMT
                            % index combining everything
                            lengths(j, 1) = lengths(j-1, 1) + lengths(j-2, 1);
                            areas(j, 1) = const.A;
    
                            % the resultant input force is positive
                            signs(j, 1) = 1;

                        case hasILM*(index.ILM+1)
                            % if MCC and ILM have already been connected in
                            % series
                            lengths(j, 1) = const.MCC.lStalk+const.ILM.thick;
                            areas(j, 1) = const.A;
    
                            % MCC+ILM expands for a positive input
                            signs(j, 1) = 1;

                        otherwise
                            error("Error in comment line "+i+": "+ ...
                                    "Check indices being connected.");
                    end
                end

                % calculate input forces and displacements; input forces
                % are multiplied by signs to adjust for their direction
                forces = vars(indices, 1) .* areas .* signs;
                displacements = vars(indices, 2) .* lengths;
    
                % add new equations
                eqns = [eqns; vmtForceDisp([forces displacements])];

            case "angleBase"
                % define extra variables to simplify expression;
                % a positive deltaL is viewed as an increase in MCC stalk 
                % length and a decrease in Vitreous length, but
                % index.Vitreous is used because the base is connected to
                % the Vitreous and then to the MCC stalk
                deltaL = -vars(index.Vitreous, 2)*const.vitr.l;
                rBase = const.MCC.rBase;
    
                % add equations for vars(i, 1) = angleBase(t)
                %   and vars(i, 2) = dummy(t)
                eqns = [eqns;
                        vars(i, 1) == simplify(atand(deltaL/rBase));
                        vars(i, 2) == 0];

            otherwise
                error("Error in comment ("+i+",2): incorrect word");
        end
    end
end

function [eqns, vars, index] = setup(t, const, comments, ...
                                        hasBase, hasStalk, hasILM)
    % SETUP set up the model of vitreous and MCC:
    %   t is the time variable
    %   const is all the constant characteristics
    %   comments contain all the element and bond information
    %   hasBase is a boolean to include for a bistable base   
    %   hasStalk is a logical variable about including MCC in the model

    arguments
        t sym
        const   struct
        comments (:, 4) string
        hasBase logical = false
        hasStalk logical = true
        hasILM logical = false
    end

    % check hasILM complies with hasStalk and hasBase
    if hasILM > (hasStalk || hasBase)
        error("Error: ILM cannot be turned on without both MCC Stalk and MCC Base");
    end

    % get number of variables in comments
    varNo = size(comments, 1);

    % preset index
    index.MCC = 0;
    index.Vitreous = varNo; % currently last comment is the entire vitreous
    index.Base = 0;
    index.ILM = 0;
    index.VBA = 0;  % index of vitreous-base attachment (for bistability)
    index.VMT = 0;
    index.Angle = 0;

    % check hasBase
    if hasBase
        % add MCC base element to comments; only 1 base spring is added to
        % avoid repetition; this is adjusted for in vbaForceDisp
        if const.simpleBase
            baseComments = ["spring" "MCC" 0;       % 1
                            "VBA" varNo varNo+1];   % 2
        else
            baseComments = ["damper" "MCC IP" 0;            % 1
                            "spring" "MCC IP" 0;            % 2
                            "parallel" varNo+1 varNo+2;     % 3
                            "damper" "MCC S" 0;             % 4
                            "spring" "MCC S" 0;             % 5
                            "parallel" varNo+4 varNo+5;     % 6
                            "series" varNo+3 varNo+6;       % 7
    
                            "VBA" varNo varNo+7];           % 8
        end

        baseNo = size(baseComments, 1);

        % memorize base and VBA indices
        index.Base = varNo+baseNo-1;
        index.VBA = varNo+baseNo;

        % add ILM if hasILM
        if hasILM
            % adjust baseNo to take out old VBA
            baseNo = baseNo-1;

            % add ILM to baseComments without last line about VBA
            baseComments = [baseComments(1:end-1, :);

                            "spring" "ILM" 0;                           % 1
                            "damper" "ILM" 0;                           % 2
                            "parallel" varNo+baseNo+1 varNo+baseNo+2;   % 3
                            "parallel ILM" varNo+baseNo varNo+baseNo+3; % 4
    
                            "VBA" varNo varNo+baseNo+4];                % 5

            % recalculate baseNo
            baseNo = size(baseComments, 1);

            index.ILM = varNo+baseNo-2;
            index.VBA = varNo+baseNo;
        end

        % recalculate number of variables in comments
        baseComments = [(varNo+1:varNo+baseNo)' baseComments];
        comments = [comments; baseComments];
        varNo = size(comments, 1);
    end
    
    % check hasStalk
    if hasStalk
        % determine type of variable used
        if hasBase
            type = "MCC EF";
        else
            type = "MCC";
        end

        % add MCC to comments
        MCCComments = ["damper" type 0;               % 1
                       "spring" type 0;               % 2
                       "parallel" varNo+1 varNo+2;    % 3

                       "VMT" varNo varNo+3];          % 4

        MCCNo = size(MCCComments, 1);

        % memorize MCC and vitreous indices
        index.MCC = varNo+MCCNo-1;

        % add ILM if hasILM and ILM has not been added yet
        if hasILM && index.ILM == 0
            % adjust MCCNo to take out old VMT
            MCCNo = MCCNo-1;

            % add ILM to MCCComments without last line about VMT
            MCCComments = [MCCComments(1:end-1, :);

                            "spring" "ILM" 0;                       % 1
                            "damper" "ILM" 0;                       % 2
                            "parallel" varNo+MCCNo+1 varNo+MCCNo+2; % 3
                            "series ILM" varNo+MCCNo varNo+MCCNo+3; % 4
    
                            "VMT" varNo varNo+MCCNo+4];             % 5

            % recalculate MCCNo
            MCCNo = size(MCCComments, 1);

            index.ILM = varNo+MCCNo-2;
        end

        % recalculate number of variables in comments
        MCCComments = [(varNo+1:varNo+MCCNo)' MCCComments];
        comments = [comments; MCCComments];
        varNo = size(comments, 1); 
    end

    % the VMT index is always the last variable; angle is not added yet, so
    % it will be equal to varNo
    index.VMT = varNo;

    % creating new variables
    syms stress(t) strain(t) [varNo 1]
    vars = [stress(t) strain(t)];   % matrix with all variables
    eqns = [];                      % vector of equations

    % if hasBase or hasStalk is true then add angle to vars
    if hasBase || hasStalk
        comments = [comments; [varNo+1 "angleBase" 0 0]];

        index.Angle = varNo+1;

        syms angleBase(t) dummy(t)
        vars = [vars; [angleBase(t) dummy(t)]];
    end

    % check that comments are in order
    checkComments(comments);

    % add equations based on comments
    eqns = [eqns; addEqns(vars, const, comments, index, hasILM)];
end

function checkComments(comments)
    % CHECKCOMMENTS checks if comments are properly set up

    varNo = size(comments, 1);
    if comments(varNo, 2) == "angleBase"
        elementIndices = string(1:varNo-2); % exclude angle and last connection
    else
        elementIndices = string(1:varNo-1); % no angle, exclude last connection
    end
    connected = string([]); % stores indices of connected elements
    
    for i = 1:varNo
        if comments(i, 1) ~= string(i)
            disp("Error in comment ("+i+",1): incorrect index");

        elseif ismember(comments(i, 2), ["spring" "damper"])
            if ~ismember(comments(i, 3), ...
                    ["gel" "liquid" "MCC" "MCC IP" "MCC EF" "MCC S" "ILM"])
                disp("Error in comment ("+i+",3): incorrect word");
            end

        elseif ismember(comments(i, 2), ["parallel" "parallel ILM" ...
                                            "series" "series ILM" ...
                                            "VMT" "VBA"])
            if ismember(comments(i, 3), connected)
                disp("Error in comment ("+i+",3): index already connected");
            elseif ~ismember(comments(i, 3), elementIndices)
                disp("Error in comment ("+i+",3): index should not be connected");
            else
                connected(length(connected)+1) = comments(i, 3);
            end

            if ismember(comments(i, 4), connected)
                disp("Error in comment ("+i+",4): index already connected");
            elseif ~ismember(comments(i, 4), elementIndices)
                disp("Error in comment ("+i+",4): index should not be connected");
            else
                connected(length(connected)+1) = comments(i, 4);
            end

        elseif comments(i, 2) ~= "angleBase"
            disp("Error in comment ("+i+",2): incorrect word");
        end
    end
end

function output = hasFuncs(expr, funcs)
    % HASFUNC performs function of HAS() for each function in funcs:
    %   has(expr, f) for all f in funcs

    % save original size of funcs
    funcSize = size(funcs);

    % reshape funcs
    funcs = funcs(:);
    len = size(funcs, 1);

    % define output vector
    output = logical([]);

    for i = 1:len
        output = [output, has(expr, funcs(i))];
    end

    % reshape output to original size of funcs
    output = reshape(output, funcSize);
end

function [varFormula, varFVars] = getVarFormula(var, sol)
    % GETVARFORMULA returns a formula of the desired variable

    if any(has(sol.cVars(:, 1), var))
        % get formula for var from 2nd column of sol.cVars row where var is
        varFormula = sol.cVars(has(sol.cVars(:, 1), var), 2);
        varFVars = {sol.t};

    elseif any(has(sol.rVars(:, 1), var))
        % get formula for var from 2nd column of sol.rVars row where var is
        varFormula = sol.rVars(has(sol.rVars(:, 1), var), 2);

        % get free variables mentioned in varFormula, including t
        varFVars = [{sol.t};
                    sym2cell(sol.newVars(hasFuncs(varFormula, sol.newVars)))];

    else
        varFormula = var;
        varFVars = {var};
    end
end

function varSol = getVarSol(var, sol, const, hasPrint)
    % GETVARSOL returns a solution array for a desired variable after solving
    % the system of equations:
    %   var is the function for which the numerical solution is constructed
    %   newVars is the set of independent variables in the system
    %   cVars is the set of variables which already have an expression in
    %       terms of t
    %   rVars is the set of variables which have a formula in terms of t
    %       and newVars
    %
    %   hasPrint is an optional variable for printing varFormula

    arguments
        var sym
        sol struct
        const struct
        hasPrint logical = true
    end

    [varFormula, varFVars] = getVarFormula(var, sol);

    if varFormula == var
        % in case of a free variable get solution
        if hasPrint
            fprintf(string(var)+' is a free variable.\n');
        end

        % as var is a free variable (therefore in sol.newVars), get
        %   numerical solution directly
        varSol = sol.newVarSols(:, has(sol.newVars, var));

    else
        % print formula
        if hasPrint
            fprintf(string(vpa(varFormula, const.rounding))+"\n");
        end

        % substitute into formula and get solution; varSol will stay the
        % size of tSol because stress and strain variables in varFVars are
        % treated as syms and not symfuncs
        varSol = double(subs(varFormula, varFVars, ...
                [{sol.tSol}; ...
                num2cell(sol.newVarSols(:, hasFuncs(varFormula, sol.newVars)), 1).']));
    end
end

function equivCond = getEquivCond(cond, sol)
    % GETEQUIVCOND substitutes an event condition for one which is
    %   expressed in terms of newVars after simplification:
    %
    %   cond - current condition

    % get name of condition variable
    condVar = sol.vars(hasFuncs(cond, sol.vars));

    % get condition value, which is on RHS of cond
    condRhs = rhs(cond);

    % get formula and formulaVars of condVar
    [formula, formulaVars] = getVarFormula(condVar, sol);

    if formula == condVar
        equivCond = cond;

    else    
        % substitute formulaVars, which are contain symfun data types as
        % functions of sol.t, with v, which are dummy sym variables
        syms v [size(formulaVars)]
        subFormula = subs(formula, formulaVars, v);
    
        % get equivalent equation
        equation = subFormula == condRhs;

        equivCond = [];
    
        % get equivalent condition by solving equation for v
        for i=1:size(formulaVars, 1)
            % get solution for v(i) and use first solution (arbitrary
            % choice, lowers no. of rhs params to 1)
            params = subs(solve(equation, v(i)), v, formulaVars);
            if ~isempty(params)
                equivCond = [equivCond;
                                formulaVars(i) == params(1)];
            end
        end
    end
end

function sol = clearSol(sol, lockedFields)
    % CLEARSOL replaces all sol fields with [] except lockedFields

    fieldNames = fieldnames(sol);

    for i = 1:size(fieldNames, 1)
        fName = fieldNames{i};

        if ~ismember(fName, lockedFields)
            sol.(fName) = [];
        end
    end
end

function baseCond = getBaseCond(vars, eqns, cond, mainSol, mainConst)
    % GETBASECOND calculates the initial condition for base stress from the
    %   base strain initial condition:
    %       vars - MCC base variables
    %       eqns - MCC base equations
    %       cond - MCC strain initial condition
    %       mainConst, mainSol - const and sol structs from main code

    % create local sol with same form as mainSol, but empty all fields
    % except mainSol.t, which is carried over
    sol = mainSol;
    lockedFields = ["t"; "tSol"; "opt"];
    sol = clearSol(sol, lockedFields);

    % adding initial signal and turn eqns to sym
    eqns = [eqns; cond]; % input signal

    % create local sol for base and store vars and eqns
    sol.vars = vars;
    sol.eqns = eqns;
    
    % simplifying system
    [sol.newEqns, sol.newVars, sol.R] = reduceRedundancies(eqns, vars);

    sol.cVars = sol.R.constantVariables;
    sol.rVars = sol.R.replacedVariables;
    sol.newVarNo = size(sol.newVars, 1);

    if sol.newVarNo>0
        sol.f = daeFunction(sol.newEqns, sol.newVars);

        sol.initCondGuess = zeros(sol.newVarNo, 1);
        sol.initCondDerGuess = zeros(sol.newVarNo, 1);
        
        % try running decic; otherwise inform of insufficient tolerance
        try
            [sol.initCond, sol.initCondDer] = decic(sol.f,0,sol.initCondGuess, ...
                                            [],sol.initCondDerGuess,[], ...
                                            sol.opt);
        catch e
            disp("Error: Adjust local sol.opt in getBaseCond function.");
        end
        
        % fprintf("Initial conditions are: \nvar(t)\tvar(0)\tvar'(0)\n")
        % disp(vpa([sol.newVars sol.initCond sol.initCondDer], ...
        %             mainConst.rounding))
        
        % solving the system
        [sol.tSol,sol.newVarSols] = ode15i(sol.f,sol.tSol,sol.initCond,...
                                        sol.initCondDer,sol.opt);

        % get solution for base stress, located at sol.vars(end, 1)
        baseSol = getVarSol(sol.vars(end, 1), sol, mainConst, false);

        % output stabilized baseCond as the condition
        baseCond = sol.vars(end, 1) == baseSol(end);
    end
end

function initCondGuess = addInitCond(initCondGuess, cond, sol)
    % ADDINITCOND adds respective initial conditions to my guesses, after
    %   which they will be fixed when using decic

    eqc = getEquivCond(cond, sol);

    [varF, varFVars] = getVarFormula(lhs(cond), sol);

    % for every var in varFVars add its rhs to the guesses
    for i = 1:size(varFVars, 1)
        var = varFVars{i};
        initCondGuess(hasFuncs(var, sol.newVars)) = rhs(eqc);
    end
end

function [values, terminal, direction] = getEvent(t, y, dy, evt, sol, const)
    % GETEVENT defines an ode15i event specified by evt of t, vars and
    %   dvars
    %
    % output:
    %   values - a variable that describes the event; an event in
    %       matlab is when a variable is equal to 0
    %   terminal - whether the ode solver should stop integrating
    %   direction - from which direction the event variable should
    %       approach 0

    % group all possible functions and derivatives in evt into a cell;
    % the order ({diffs; funcs; t}) is important, as it impacts how values
    % are substituted
    vars = [sym2cell(diff(sol.newVars, sol.t));
            sym2cell(sol.newVars);
            {sol.t}];

    % group all numerical values corresponding to vars into a cell
    params = [num2cell(dy);
                num2cell(y);
                {t}];

    % calculate output
    values = vpa(subs(evt, vars, params), const.rounding);
    terminal = 1;
    direction = 0;
end

function eventFunc = getEventFunc(evt, sol, const)
    % GETEVENTFUNC returns a function handle for an event specified by evt

    % transform evt into equivalent in terms of newVars
    equivEvent = getEquivCond(evt, sol);

    % rewrite equivEvent into form "lhs==0"
    equivEvent = lhs(equivEvent) - rhs(equivEvent);

    % create event function handle
    eventFunc = @(t, y, dy) getEvent(t, y, dy, equivEvent, sol, const);
end



%% UNUSED CODE
% UNUSED FUNCTIONS:
%
% function [eqns, vars] = homogeneous(t, const)
%     % HOMOGENEOUS set up homogeneous vitreous model (no MCC):
%     %   t is the time variable
%     %   const is constant characteristics
% 
%     % creating new variables
%     elementNo = 2;              % two elements: 1 spring and 1 damper
%     bondNo = elementNo - 1;     % number of bonds is always 1 less than elementNo
%     varNo = elementNo + bondNo; % number of stress and strain variables
%     syms stress(t) strain(t) [varNo 1]
% 
%     % create matrix with all stresses and strains
%     vars = [stress(t) strain(t)];
% 
%     % springs and dampers are set up in left-right, top-bottom direction;
%     % bonds are set up immediately after all the components are set up
%     eqns = [];
% 
%     % setting up the system
%     eqns = [damper(vars(1, :), const.vitr.viscosityGel);
%             spring(vars(2, :), const.vitr.stiffnessGel)];
%     eqns = [eqns;
%             parallel(vars([1 2 3], :))];
% end
%
% function comments = shiftIndices(comments)
%     % SHIFTINDICES shifts all indices in bonds by 1 upward
% 
%     varNo = size(comments, 1);
% 
%     for i = 1:varNo
%         % check first column index is correct
%         if str2double(comments(i, 1)) ~= i
%             comments(i, 1) = i;
%         end
% 
%         % shift bond indices by 1 upward
%         if ismember(comments(i, 2), ["parallel" "series" "VMT"])
%             comments(i, 3) = str2double(comments(i, 3)) + 1;
%             comments(i, 4) = str2double(comments(i, 4)) + 1;
%         end
%     end
% end

% OTHER UNUSED CODE PARTS
% 
% % OLD LEGEND WITH PA
%         if yExp+3 == 0
%             ylabel("\boldmath$\sigma$ \textbf{(Pa)}");
%         else
%             ylabel("\boldmath$\sigma\times10^{"+string(yExp+3)+"}$ \textbf{(Pa)}");
%         end
%
%
%
% % set options for ode solver
% event = sol.vars(sol.index.Angle, 1)==-10;
% eventFunc = getEventFunc(event, sol, const);
% sol.opt = odeset('RelTol',const.relTol,'AbsTol',const.absTol,...
%                     'Events',eventFunc);
%

