function report = run_paper_results(mode, outputDirectory)
% Selected paper evidence only. Default 'verify' is read-only and config-free.
% Run Python tools/verify_paper_results.py first for exact certificate checks.
% 'render' makes deterministic tables and four scientific figures, no solver.
if nargin<1, mode='verify'; end
root=fileparts(mfilename('fullpath'));
if nargin<2, outputDirectory=fullfile(root,'results','figures','paper'); end
assert(any(strcmp(mode,{'verify','render'})),'Mode must be verify or render');
v=jsondecode(fileread(fullfile(root,'results/research/R135-paper-verified.json')));
assert(v.passed && strcmp(sha256(fullfile(root,'results/research/R135-paper-fixture.json')),v.fixtureSha256));
assert(v.bonuses.S(1)>1.5*v.bonuses.U(2) && v.bonuses.U(1)>0);
assert(abs(v.ratioLower-v.bonuses.S(1)/v.bonuses.U(2))<1e-13);
assert(v.capacities.S0(2)<v.capacities.U0(1) && v.capacities.S1(2)<v.capacities.U1(1));
report=struct('passed',true,'mode',mode,'matlabVersion',version,'newOptimizerCalls',0, ...
    'fixtureSha256',v.fixtureSha256,'verification','floating presentation checks; exact proof uses Python verifier');
if strcmp(mode,'verify')
    fprintf('Paper presentation checks passed; exact bonus ratio lower %.10f. No files written.\n',v.ratioLower);return
end
if ~isfolder(outputDirectory),mkdir(outputDirectory);end
blue=[.04 .30 .53]; orange=[.76 .32 .08]; gray=[.40 .43 .48];
% Figure 1: endpoints are rigorous bounds, not confidence intervals or optima.
f=figure('Visible','off','Color','w','Position',[100 100 1120 460]);cl=onCleanup(@()close(f));
t=tiledlayout(f,1,2,'TileSpacing','loose','Padding','loose');
ax=nexttile(t);hold(ax,'on');
names={'U exact','S exact','U, \rho = 0.10','S, \rho = 0.10'};
vals=[v.capacities.U0';v.capacities.S0';v.capacities.U1';v.capacities.S1'];
for i=1:4
    color=blue;if mod(i,2)==0,color=orange;end
    plot(ax,vals(i,:),[i i],'-','Color',color,'LineWidth',5);
    plot(ax,vals(i,1),i,'o','Color',color,'MarkerFaceColor',color,'MarkerSize',7);
    plot(ax,vals(i,2),i,'o','Color',color,'MarkerFaceColor','w','MarkerSize',7);
end
set(ax,'YTick',1:4,'YTickLabel',names,'YDir','reverse','FontSize',11);ylim(ax,[.5 4.5]);xlim(ax,[0 .58]);grid(ax,'on');
xlabel(ax,'Target amplitude \beta');title(ax,'Total capacity');
ax=nexttile(t);hold(ax,'on');
for i=1:2
    val=v.bonuses.U;color=blue;if i==2,val=v.bonuses.S;color=orange;end
    plot(ax,val,[i i],'-','Color',color,'LineWidth',6);
    plot(ax,val(1),i,'o','Color',color,'MarkerFaceColor',color,'MarkerSize',8);
    plot(ax,val(2),i,'o','Color',color,'MarkerFaceColor','w','MarkerSize',8);
end
set(ax,'YTick',[1 2],'YTickLabel',{'U','S'},'YDir','reverse','FontSize',11);ylim(ax,[.5 2.8]);xlim(ax,[0 .28]);grid(ax,'on');
xlabel(ax,'Absolute bonus \Delta');title(ax,'Tolerance bonus');
text(ax,.012,2.55,'S bonus > 1.50 times U bonus','FontWeight','bold','FontSize',11);
title(t,'Capacity and absolute tolerance bonus','FontSize',13);
savefigures(f,outputDirectory,'01-capacity-and-absolute-bonus');clear cl
% Figure 2: analytical bounds; plotting grid is not a numerical experiment.
f=figure('Visible','off','Color','w','Position',[100 100 930 550]);cl=onCleanup(@()close(f));ax=axes(f);hold(ax,'on');
r=linspace(.02,.20,501);p=v.values;
lower=p.envelopeBaselineWitness./(1-p.blendKappa*r);
ii=r>.10;lower(ii)=p.envelopeTolerantWitness*.9./(1-r(ii));
lower=max(0,lower-p.envelopeBaselineUpper./(1-r));
upper=p.supportUpper./(p.pairing-p.penaltyUpper*r)-p.unrestrictedBaselineLower./(1-r);
ylim(ax,[0 .4]);iv=v.closedToleranceInterval;
patch(ax,[iv(1) iv(2) iv(2) iv(1)],[0 0 .4 .4],[.90 .94 .89],'EdgeColor','none','HandleVisibility','off');
plot(ax,r,lower,'Color',orange,'LineWidth',2.3,'DisplayName','S bonus: certified lower bound');
plot(ax,r,upper,'Color',blue,'LineWidth',2.3,'DisplayName','U bonus: certified upper bound');
xline(ax,.10,':','DisplayName','Primary \rho = 0.10');
xline(ax,iv(1),'--','HandleVisibility','off');xline(ax,iv(2),'--','HandleVisibility','off');
xlim(ax,[.02 .20]);grid(ax,'on');set(ax,'FontSize',12);xlabel(ax,'Global effect-relative allowance \rho');ylabel(ax,'Absolute bonus \Delta');
legend(ax,'Location','northwest');title(ax,'Certified ordering on [0.056736, 0.137809]','FontSize',13);
savefigures(f,outputDirectory,'02-tolerance-interval');clear cl
% Figure 3: universal actual response bounds and a feasible witness are distinct.
f=figure('Visible','off','Color','w','Position',[100 100 1120 490]);cl=onCleanup(@()close(f));t=tiledlayout(f,1,2,'TileSpacing','loose','Padding','loose');
ax=nexttile(t);data=zeros(2,3);
for j=1:2
    a=v.responseBounds(j);data(j,:)=[a.baselineXIncrementUpper,a.actualXIncrementLower,v.physical(j+1).actualChange(1)];
end
b=bar(ax,data);b(1).FaceColor=gray;b(2).FaceColor=orange;b(3).FaceColor=blue;
set(ax,'XTickLabel',{'A','FL11'},'FontSize',11);ylabel(ax,'Actual X increment (white Y = 100)');ylim(ax,[0 7.7]);grid(ax,'on');
legend(ax,{'Exact S upper','Universal S lower','Feasible witness'},'Location','northwest');title(ax,'Actual X response bounds');
ax=nexttile(t);cost=[v.responseBounds(1).YDeficitLower,-v.physical(2).actualChange(2)];
bar(ax,1:2,cost,.55,'FaceColor',orange);set(ax,'XTick',1:2,'XTickLabel',{'Universal bound','Witness'},'FontSize',11);
ylabel(ax,'A-illuminant Y deficit from 50');ylim(ax,[0 .65]);grid(ax,'on');
text(ax,1,.28,sprintf('> %.8f',floor(cost(1)*1e8)/1e8),'HorizontalAlignment','center','FontSize',11);
text(ax,2,.51,sprintf('%.8f',cost(2)),'HorizontalAlignment','center','FontSize',11);title(ax,'Necessary A-Y deficit');
title(t,'Class S at \beta_* \approx 0.4027278731','FontSize',13);
savefigures(f,outputDirectory,'03-actual-response-and-cost');clear cl
% Figure 4: compare the two prescribed directions in the SAME alpha units.
f=figure('Visible','off','Color','w','Position',[100 100 850 480]);cl=onCleanup(@()close(f));ax=axes(f);
bar(ax,[v.scope.alphaXLower,v.scope.alphaZUpper],.55,'FaceColor',orange);grid(ax,'on');
set(ax,'XTick',1:2,'XTickLabel',{'X/Y: lower bound','Z/Y: upper bound'},'FontSize',12);ylim(ax,[0 .0205]);
ylabel(ax,'Absolute bonus in \alpha = K\beta');title(ax,'Direction comparison: S, \rho = 0.10','FontSize',13);
text(ax,1,v.scope.alphaXLower+.001,sprintf('>= %.10f',floor(v.scope.alphaXLower*1e10)/1e10),'HorizontalAlignment','center');
text(ax,2,v.scope.alphaZUpper+.001,sprintf('<= %.12f',ceil(v.scope.alphaZUpper*1e12)/1e12),'HorizontalAlignment','center');
text(ax,1.5,.019,'X bonus > 45.278 times Z bonus (Z may be zero)','HorizontalAlignment','center','FontWeight','bold');
savefigures(f,outputDirectory,'04-common-alpha-direction');clear cl
% Tables include absolute targets and errors, as well as the global ratio.
bounds=boundstrings(v.printedCapacityBounds);
capacity=table(string(names'),bounds(:,1),bounds(:,2),'VariableNames',{'Condition','LowerBeta','UpperBeta'});
writetable(capacity,fullfile(outputDirectory,'capacity-bounds.csv'));
bounds=boundstrings(v.printedBonusBounds);
bonus=table(["U";"S"],bounds(:,1),bounds(:,2),'VariableNames',{'Class','DeltaLower','DeltaUpper'});
writetable(bonus,fullfile(outputDirectory,'absolute-bonus-bounds.csv'));
physical=zeros(3,11);
for j=1:3
    a=v.physical(j);physical(j,:)=[a.target(:)',a.actualChange(:)',a.error(:)',a.absoluteErrorUpper,a.globalRelativeErrorUpper];
end
response=array2table(physical,'VariableNames',{'TargetX','TargetY','TargetZ','ActualChangeX','ActualChangeY','ActualChangeZ','ErrorX','ErrorY','ErrorZ','ErrorNormApprox','GlobalRelativeRatioApprox'});
response=addvars(response,["D65";"A";"FL11"],'Before',1,'NewVariableNames','Illuminant');
writetable(response,fullfile(outputDirectory,'witness-absolute-response.csv'));
writej(fullfile(outputDirectory,'render-receipt.json'),report);
fprintf('Four deterministic figures and three evidence tables saved to %s\n',outputDirectory);
end
function bounds=boundstrings(rows)
% jsondecode preserves nested arrays of strings as nested cell arrays.
bounds=strings(numel(rows),2);
for i=1:numel(rows),bounds(i,:)=reshape(string(rows{i}),1,2);end
end
function savefigures(f,d,name)
exportgraphics(f,fullfile(d,[name '.png']),'Resolution',160);
exportgraphics(f,fullfile(d,[name '.pdf']),'ContentType','vector');
end
function h=sha256(p)
f=fopen(p,'rb');assert(f>=0);cl=onCleanup(@()fclose(f));bytes=fread(f,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function writej(p,x)
f=fopen(p,'w');assert(f>=0);cl=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(x,PrettyPrint=true));
end
