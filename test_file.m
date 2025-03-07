%% set up constants
hasBase = true;
hasStalk = true;
hasILM = true;
hasPlot = false;

forces = 1:8; % 1-sin, 2-lin, 3-quad, 4-exp, 5--pert, 6-+pert, 7-gaus, 8-nil
areas = 1:3; % 4th is extra
comments = 1:2;

relTol = 5*10^(-7);
absTol = 5*10^(-5);
tSol = 0:0.05:30;



%% run simple model
commentNo = 2;
hasBase = false;

for forceNo = [2 3 4]
    for areaNo = areas
        if hasBase
            disp("BISTABLE MODEL")
        else
            disp("SIMPLE MODEL")
        end
        sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
                            forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
                            commentNo=commentNo, hasPlot=hasPlot, tSol=tSol);
    end
end



% %% run bistable model
% commentNo = 2;
% hasBase = true;
% angleNegBist = 0;
% 
% for forceNo = [2 3 4]
%     for areaNo = areas
%         if hasBase
%             sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                                 forceNo=5, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                                 commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, ...
%                                 onlyAngle=true, hasPrint=false);
% 
%             angleNegBist = sol.varSols(end, end)
%         else
%             angleNegBist = 0
%         end
% 
%         sol = vitreous_model_v5(hasBase=hasBase, hasStalk=hasStalk, hasILM=hasILM, ...
%                     forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, ...
%                     angleBist = angleNegBist);
% 
%         % disp("forceNo: "+forceNo)
%         % disp("hasBase: "+hasBase)
%         % disp("Starting angle: "+angleNegBist)
%         % 
%         % disp("Final angle: "+sol.varSols(end, end))
%         % disp()
%     end
% end
% 
% 
% 
% %% get all final values
% simpleVals = [];
% bistVals = [];
% angleNegBist = 0;
% commentNo = 2;
% 
% for hasBase = [false true]
%     for forceNo = [2 3 4]
%         for areaNo = areas
%             if hasBase
%                 sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                                     forceNo=5, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                                     commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, ...
%                                     onlyAngle=true, hasPrint=false);
% 
%                 angleNegBist = sol.varSols(end, end)
% 
%                 % add forceNo, areaNo and angleNegBist to values
%                 sizeBistVals = size(bistVals);
%                 bistVals(sizeBistVals(1)+1, 1:3) = [forceNo, areaNo, angleNegBist]
%             else
%                 angleNegBist = 0
% 
%                 sizeSimpleVals = size(simpleVals);
%                 simpleVals(sizeSimpleVals(1)+1, 1:2) = [forceNo, areaNo]
%             end
% 
%             sol = vitreous_model_v5(hasBase=hasBase, hasStalk=hasStalk, hasILM=hasILM, ...
%                         forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                         commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, ...
%                         angleBist = angleNegBist);
% 
%             if hasBase
%                 % get stresses: MCC, Vitr, Base, ILM, VBA
%                 bistVals(sizeBistVals(1)+1, 4:8) = sol.varSols(end, [1 3 5 7 11]);
% 
%                 % get strains: MCC, Vitr, Base, ILM, VBA
%                 bistVals(sizeBistVals(1)+1, 9:13) = sol.varSols(end, [2 4 6 8 12]);
% 
%                 % get angle
%                 bistVals(sizeBistVals(1)+1, 14) = sol.varSols(end, end);
%             else
%                 % get stresses: MCC, Vitr, ILM
%                 simpleVals(sizeSimpleVals(1)+1, 3:5) = sol.varSols(end, [1 3 5]);
% 
%                 % get strains: MCC, Vitr, ILM
%                 simpleVals(sizeSimpleVals(1)+1, 6:8) = sol.varSols(end, [2 4 6]);
% 
%                 % get angle
%                 simpleVals(sizeSimpleVals(1)+1, 9) = sol.varSols(end, end);
%             end
%         end
%     end
% end
% 
% 
% 
% %% get angle of bistability and test +pert
% 
% forceNo = 5;
% areaNo = 2;
% commentNo = 2;
% 
% sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                     forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     commentNo=commentNo, hasPlot=hasPlot, tSol=tSol);
% 
% angleNegBist = sol.varSols(end, end);
% 
% sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                     forceNo=6, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, ...
%                     angleBist=angleNegBist);
% 
% 
% 
% %% get delta (phase lag between stress and strain)
% forceNo = 1;
% areaNo = 2;
% commentNo = 2;
% step = 0.001;
% 
% for scalePeriod = [1 10 100 1000]
%     for commentNo = comments
%         disp("scalePeriod: "+scalePeriod)
%         tSol = 0:step:30;
% 
%         sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                         forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                         commentNo=commentNo, hasPlot=hasPlot, tSol=tSol, hasPrint=false, ...
%                         scalePeriod=scalePeriod);
% 
%         rowsDelta = sol.tSol < 2*pi/scalePeriod*floor(scalePeriod*30/(2*pi));
%         length(tSol(rowsDelta))/length(tSol)*100
% 
%         mcc_stress = sol.varSols(rowsDelta, [1 3 5]) ./ ...
%                         max(sol.varSols(rowsDelta, [1 3 5]));
%         mcc_strain = sol.varSols(rowsDelta, [2 4 6]) ./ ...
%                         max(sol.varSols(rowsDelta, [2 4 6]));
% 
%         cosDelta = mean(2*mcc_stress.*mcc_strain);
%         acosd(cosDelta)
%     end
% end
% 
% 
% 
% %% change vitreous model
% % get number of rows in plot
% plotRowNo = 3+hasILM+2*hasBase;
% 
% % create new figure and other specifspications
% fig = figure('Position', [400 420-60*plotRowNo 750 120*plotRowNo]);
% set(0,'defaultTextInterpreter','latex');
% colours = colororder("gem");
% lineStyles = linestyleorder("mixedstyles");
% 
% forceNo = 1;
% areaNo = 2;
% 
% for commentNo = comments
%     % get solution
%     sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                     forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     commentNo=commentNo, hasPlot=hasPlot);
% 
%     % set up subplot row counter and xlims
%     rowCounter = 1;
%     xlims = [sol.tSol(1) sol.tSol(end)];
% 
%     % plotting time evolution
%     for i = 1:size(sol.indexFields, 1)
%         field = sol.indexFields{i};
% 
%         if field == "Angle"
%             % plot MCC base angle
%             subplot(plotRowNo, 2, [2*rowCounter-1, 2*rowCounter]);
%             plot(sol.tSol, sol.varSols(:, end), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 letter = char('A'+2*rowCounter-2);
%                 title("\textbf{"+letter+") Angle of MCC Base}");
%                 xlabel("\boldmath$t$ \textbf{\textbf{(days)}}");
%                 ylabel("\boldmath$\alpha$ $(^{\circ})$");
%                 xlim(xlims);
%             end
% 
%             if commentNo==1
%                 grid on;
%             end
%             hold on;
% 
%             angleBist = sol.varSols(end, end);
%             vpa(angleBist)
% 
%             rowCounter = rowCounter+1;
% 
%         elseif field ~= "VMT" && getfield(sol.index, field) ~= 0
%             % plot variable stress
%             subplot(plotRowNo, 2, 2*rowCounter-1);
%             plot(sol.tSol, sol.varSols(:, 2*rowCounter-1+2*(field=="Vba")), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 ax = gca;
%                 yExp = ax.YAxis.Exponent;
%                 ax.YTickMode = "manual";
%                 ax.YTickLabelMode = "manual";
% 
%                 letter = char('A'+2*rowCounter-2);
%                 title("\textbf{"+letter+") "+field+" Stress}");
%                 xlabel("\boldmath$t$ \textbf{(days)}");
%                 if yExp == 0
%                     ylabel("\boldmath$\sigma$ \textbf{(kPa)}");
%                 else
%                     ylabel("\boldmath$\sigma$ $(10^{"+string(yExp)+"}$ \textbf{kPa)}");
%                 end
% 
%                 xlim(xlims);
%             end
% 
%             if commentNo==1
%                 grid on;
%             end
%             hold on;
% 
%             % plot variable strain
%             subplot(plotRowNo, 2, 2*rowCounter);
%             plot(sol.tSol, sol.varSols(:, 2*rowCounter+2*(field=="Vba")), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 ax = gca;
%                 yExp = ax.YAxis.Exponent;
%                 ax.YTickMode = "manual";
%                 ax.YTickLabelMode = "manual";
% 
%                 letter = char(letter+1);
%                 title("\textbf{"+letter+") "+field+" Strain}");
%                 xlabel("\boldmath$t$ \textbf{(days)}");
%                 if yExp==0
%                     ylabel("\boldmath$\varepsilon$");
%                 else
%                     ylabel("\boldmath$\varepsilon$ $(10^{"+string(yExp)+"})$");
%                 end
% 
%                 xlim(xlims);
%             end
% 
%             if commentNo==1
%                 grid on;
%             end
%             hold on;
% 
%             rowCounter = rowCounter + 1;
%         end
%     end
% end
% 
% % if forceNo == 5
% %     angleBist = sol.varSols(end, end);
% %     vpa(angleBist)
% % 
% %     vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
% %                     forceNo=6, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
% %                     commentNo=commentNo, hasPlot=true, angleBist=angleBist);
% % end
% 
% 
% 
% %% force area change
% % get number of columns in plot
% plotColNo = 3; % depends on no. of forces tested: 2, 3 and 4
% endMCCStress = [];
% 
% % create new figure and other specifspications
% fig = figure('Position', [400 300 750 400]);
% set(0,'defaultTextInterpreter','latex');
% colours = colororder("gem");
% lineStyles = linestyleorder("mixedstyles");
% 
% for forceNo = forces(2:4)
%     % initialize subplot
% 
%     for areaNo = areas
%         commentNo = 2;
%         % get solution
%         sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                         forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                         commentNo=commentNo, hasPlot=hasPlot);
% 
%         endMCCStress = [endMCCStress; sol.varSols(end, 1)];
% 
%         subplot(2, 3, forceNo-1);
%         plot(sol.tSol, sol.varSols(:, 2), ...
%                'Color', colours(forceNo-1, :), 'LineStyle', lineStyles(areaNo));
% 
%         grid on;
%         hold on;
% 
%         subplot(2, 3, forceNo-1+3);
%         plot(sol.tSol, sol.varSols(:, 4), ...
%                'Color', colours(forceNo-1, :), 'LineStyle', lineStyles(areaNo));
% 
%         grid on;
%         hold on;
%     end
% end
% 
% 
% 
% %% force 5
% plotRowNo = 3+hasILM+2*hasBase;
% 
% % create new figure and other specifspications
% fig = figure('Position', [400 410-70*plotRowNo 750 140*plotRowNo]);
% set(0,'defaultTextInterpreter','latex');
% colours = colororder("gem");
% lineStyles = linestyleorder("mixedstyles");
% 
% for areaNo = areas
%     forceNo = 5;
%     commentNo = 2;
%     % get solution
%     sol = vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                     forceNo=forceNo, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     commentNo=commentNo, hasPlot=hasPlot);
% 
%     disp("Area "+areaNo+" angle: "+sol.varSols(end, end))
% 
%     % set up subplot row counter and xlims
%     rowCounter = 1;
%     xlims = [sol.tSol(1) sol.tSol(end)];
% 
%     % plotting time evolution
%     for i = 1:size(sol.indexFields, 1)
%         field = sol.indexFields{i};
% 
%         if field == "Angle"
%             % plot MCC base angle
%             subplot(plotRowNo, 2, [2*rowCounter-1, 2*rowCounter]);
%             plot(sol.tSol, sol.varSols(:, end), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 letter = char('A'+2*rowCounter-2);
%                 title("\textbf{"+letter+") Angle of MCC Base}");
%                 xlabel("\boldmath$t$ \textbf{\textbf{(days)}}");
%                 ylabel("\boldmath$\alpha$ $(^{\circ})$");
%                 xlim(xlims);
%             end
% 
%             grid on;
%             hold on;
% 
%             rowCounter = rowCounter+1;
% 
%         elseif field ~= "VMT" && getfield(sol.index, field) ~= 0
%             % plot variable stress
%             subplot(plotRowNo, 2, 2*rowCounter-1);
%             plot(sol.tSol, sol.varSols(:, 2*rowCounter-1+2*(field=="Vba")), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 ax = gca;
%                 yExp = ax.YAxis.Exponent;
%                 ax.YTickMode = "manual";
%                 ax.YTickLabelMode = "manual";
% 
%                 letter = char('A'+2*rowCounter-2);
%                 title("\textbf{"+letter+") "+field+" Stress}");
%                 xlabel("\boldmath$t$ \textbf{(days)}");
%                 if yExp == 0
%                     ylabel("\boldmath$\sigma$ \textbf{(kPa)}");
%                 else
%                     ylabel("\boldmath$\sigma$ $(10^{"+string(yExp)+"}$ \textbf{kPa)}");
%                 end
% 
%                 xlim(xlims);
%             end
% 
%             grid on;
%             hold on;
% 
%             % plot variable strain
%             subplot(plotRowNo, 2, 2*rowCounter);
%             plot(sol.tSol, sol.varSols(:, 2*rowCounter+2*(field=="Vba")), 'LineWidth', 2, ...
%                     'Color', colours(i, :), 'LineStyle', lineStyles(commentNo));
% 
%             if commentNo == 1
%                 ax = gca;
%                 yExp = ax.YAxis.Exponent;
%                 ax.YTickMode = "manual";
%                 ax.YTickLabelMode = "manual";
% 
%                 letter = char(letter+1);
%                 title("\textbf{"+letter+") "+field+" Strain}");
%                 xlabel("\boldmath$t$ \textbf{(days)}");
%                 if yExp==0
%                     ylabel("\boldmath$\varepsilon$");
%                 else
%                     ylabel("\boldmath$\varepsilon$ $(10^{"+string(yExp)+"})$");
%                 end
% 
%                 xlim(xlims);
%             end
% 
%             grid on;
%             hold on;
% 
%             rowCounter = rowCounter + 1;
%         end
%     end
% end
% 
% 
% 
% %% zoom in
% vitreous_model_v5(hasBase = hasBase, hasStalk = hasStalk, hasILM = hasILM, ...
%                     forceNo=6, areaNo=areaNo, relTol=relTol,absTol=absTol, ...
%                     angleBist=angleBist, tSol=14.7:0.00001:14.8);
% 
% 
% 
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