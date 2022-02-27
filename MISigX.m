
%% MI Significance Plots
% load('Mat Files//Using B//AllMIsB_XR32.mat')
% load('Mat Files//Using B//AllMIsA_XR32.mat')
% load('AllMIsB_X32J.mat')
% load('AllMIsA_X32J.mat')
load('AllMIsB_X32.mat')
load('AllMIsA_X32.mat')
%%
ConnA = squeeze(AllMIsA); 
ConnB = squeeze(AllMIsB);   
%
BaselineConn = squeeze(ConnB);
InvolvedConn = squeeze(ConnA(:,:,:,:,:,1,:));
NonInvolvedConn = squeeze(ConnA(:,:,:,:,:,2,:));
BaselineP = permute(BaselineConn, [1 2 6 3 4 5]);
BaselineP = median(BaselineP(:,:,:,:),4);
BaselineP = FillMat(BaselineP);
BaselineP = (8/7)*squeeze(mean(BaselineP,1));
InvolvedP = permute(InvolvedConn, [1 2 6 3 4 5]);
InvolvedP = median(InvolvedP(:,:,:,:),4);
InvolvedP = FillMat(InvolvedP);
InvolvedP = (8/7)*squeeze(mean(InvolvedP,1));
NonInvolvedP = permute(NonInvolvedConn, [1 2 6 3 4 5]);
NonInvolvedP = median(NonInvolvedP(:,:,:,:),4);
NonInvolvedP = FillMat(NonInvolvedP);
NonInvolvedP = (8/7)*squeeze(mean(NonInvolvedP,1));
Muscles = {'FDI', 'APB', 'ECRL', 'FDS', 'FCU', 'EDC', 'UT', 'LT'};
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.025], [0.1 0.1], [0.07 0.01]);
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
YLabel = 'Degree';
YLim = [0 1];
AllS = [];
Combos = [1 2; 2 3; 1 3];
AllDF = {};
AllCS = [];
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
for m=1:4
    subplot(1,4,m)
    Comparison2 = [];
    MuscleB = BaselineP(m,:);
    % Involved
    MuscleI = InvolvedP(m,:);
    % Non-involved
    MuscleNI = NonInvolvedP(m,:);
    Comparison2(1:length(MuscleB),1) = (MuscleB);
    Comparison2(1:length(MuscleI),2) = (MuscleI);
    Comparison2(1:length(MuscleNI),3) = (MuscleNI);
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel, YLim);
    Allps{m} = p;
    AllDF{m} = df;
    AllCS(m) = Chi_sq;
    AllS(m) = ComputeSensitivity(Comparison2);
    title(Muscles{m});
    set(gca, 'Fontsize',30);
    if m>1
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
    yticks(0:0.1:1);
end
%
figure('units','normalized','outerposition',[0 0 1 1]);
for m=5:8
    subplot(1,4,m-4)
    Comparison2 = [];
    MuscleB = BaselineP(m,:);
    % Involved
    MuscleI = InvolvedP(m,:);
    % Non-involved
    MuscleNI = NonInvolvedP(m,:);
    Comparison2(1:length(MuscleB),1) = (MuscleB);
    Comparison2(1:length(MuscleI),2) = (MuscleI);
    Comparison2(1:length(MuscleNI),3) = (MuscleNI);
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel, YLim);   
    AllDF{m} = df;
    Allps{m} = p;
    AllCS(m) = Chi_sq;
    AllS(m) = ComputeSensitivity(Comparison2);
    title(Muscles{m});
    set(gca, 'Fontsize',30);
    if m>5
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
    yticks(0:0.1:1);
end
T2 = array2table(AllS, 'VariableNames',Muscles );
% writetable(rows2vars(T2),'AllSensitivity.xlsx','Sheet',1);

%% MI Significance Plots For Conditions
%
ConnA = squeeze(AllMIsA); 
ConnB = squeeze(AllMIsB);   
%
BaselineConn = squeeze(ConnB);
InvolvedConn = squeeze(ConnA(:,:,:,:,:,1,:));
NonInvolvedConn = squeeze(ConnA(:,:,:,:,:,2,:));
BaselineP = permute(BaselineConn, [1 2 6 5 3 4]);
BaselineP = median(BaselineP(:,:,:,:,:),5);
InvolvedP = permute(InvolvedConn, [1 2 6 5 3 4]);
InvolvedP = median(InvolvedP(:,:,:,:,:),5);
NonInvolvedP = permute(NonInvolvedConn, [1 2 6 5 3 4]);
NonInvolvedP = median(NonInvolvedP(:,:,:,:,:),5);
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.05], [0.1 0.1], [0.07 0.01]);
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
YLabel = 'Mean Degree';
YLim = [0 1];
AllS = [];
Combos = [1 2; 2 3; 1 3];
AllConds = {'K';'NTVK';'VK';'V'; 'T'; 'TK'; 'TV'; 'TVK'};
figure('units','normalized','outerposition',[0 0 1 1]);
for m=1:4
    subplot(1,4,m)
    Comparison2 = [];
    MuscleB = BaselineP(:,:,:,m);
    MuscleB = FillMatN(MuscleB,8);
    MuscleB = (8/7)*squeeze(mean(MuscleB,1));
    % Involved
    MuscleI = InvolvedP(:,:,:,m);
    MuscleI = FillMatN(MuscleI,8);
    MuscleI = (8/7)*squeeze(mean(MuscleI,1));
    % Non-involved
    MuscleNI = NonInvolvedP(:,:,:,m);
    MuscleNI = FillMatN(MuscleNI,8);
    MuscleNI = (8/7)*squeeze(mean(MuscleNI,1));
    Comparison2(1:length(MuscleB),1) = mean(MuscleB,1);
    Comparison2(1:length(MuscleI),2) = mean(MuscleI,1);
    Comparison2(1:length(MuscleNI),3) = mean(MuscleNI,1);
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel, YLim);  
    AllDF{m} = df;
    AllCS(m) = Chi_sq;
    AllS(m) = ComputeSensitivity(Comparison2);
    title(AllConds{m});
    set(gca, 'Fontsize',30);
    if m>1
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
end
%
figure('units','normalized','outerposition',[0 0 1 1]);
for m=5:8
    subplot(1,4,m-4)
    Comparison2 = [];
    MuscleB = BaselineP(:,:,:,m);
    MuscleB = FillMatN(MuscleB,8);
    MuscleB = (8/7)*squeeze(mean(MuscleB,1));
    % Involved
    MuscleI = InvolvedP(:,:,:,m);
    MuscleI = FillMatN(MuscleI,8);
    MuscleI = (8/7)*squeeze(mean(MuscleI,1));
    % Non-involved
    MuscleNI = NonInvolvedP(:,:,:,m);
    MuscleNI = FillMatN(MuscleNI,8);
    MuscleNI = (8/7)*squeeze(mean(MuscleNI,1));
    Comparison2(1:length(MuscleB),1) = mean(MuscleB,1);
    Comparison2(1:length(MuscleI),2) = mean(MuscleI,1);
    Comparison2(1:length(MuscleNI),3) = mean(MuscleNI,1);
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel, YLim); 
    AllDF{m} = df;
    AllCS(m) = Chi_sq;
    AllS(m) = ComputeSensitivity(Comparison2);
    title(AllConds{m});
    set(gca, 'Fontsize',30);
    if m>5
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
end
T2 = array2table(AllS, 'VariableNames',AllConds );
% writetable(rows2vars(T2),'AllSensitivity.xlsx','Sheet',2);
%% Try to cluster conditions

%% Plot Conditions as Performance Score
ConnA = squeeze(AllMIsA); 
ConnB = squeeze(AllMIsB);   
%
BaselineConn = squeeze(ConnB);
InvolvedConn = squeeze(ConnA(:,:,:,:,:,1,:));
NonInvolvedConn = squeeze(ConnA(:,:,:,:,:,2,:));
BaselineP = permute(BaselineConn, [1 2 6 5 3 4]);
BaselineP = median(BaselineP(:,:,:,:,:),5);
InvolvedP = permute(InvolvedConn, [1 2 6 5 3 4]);
InvolvedP = median(InvolvedP(:,:,:,:,:),5);
NonInvolvedP = permute(NonInvolvedConn, [1 2 6 5 3 4]);
NonInvolvedP = median(NonInvolvedP(:,:,:,:,:),5);
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.01], [0.1 0.1], [0.07 0.01]);
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
YLabel = 'Mean Degree';
YLim = [0 1];
Combos = [1 2; 2 3; 1 3];
AllConds = {'K';'NTVK';'VK';'V'; 'T'; 'TK'; 'TV'; 'TVK'};
Comp1 = [];
Comp2 = [];
for m=1:8
    Comparison2 = [];
    MuscleB = BaselineP(:,:,:,m);
    MuscleB = FillMatN(MuscleB,8);
    MuscleB = (8/7)*squeeze(mean(MuscleB,1));
    % Involved
    MuscleI = InvolvedP(:,:,:,m);
    MuscleI = FillMatN(MuscleI,8);
    MuscleI = (8/7)*squeeze(mean(MuscleI,1));
    % Non-involved
    MuscleNI = NonInvolvedP(:,:,:,m);
    MuscleNI = FillMatN(MuscleNI,8);
    MuscleNI = (8/7)*squeeze(mean(MuscleNI,1));
    Comparison2(1:length(MuscleB),1) = mean(MuscleB,1);
    Comparison2(1:length(MuscleI),2) = mean(MuscleI,1);
    Comparison2(1:length(MuscleNI),3) = mean(MuscleNI,1);
%     Score1 = abs(Comparison2(:,3)-Comparison2(:,1))./Comparison2(:,3);
%     Score2 = abs(Comparison2(:,3)-Comparison2(:,2))./Comparison2(:,3);
      Score1 = Comparison2(:,1)./Comparison2(:,3);
      Score2 = Comparison2(:,2)./Comparison2(:,3);
      Comp1(:,m) = Score1;
      Comp2(:,m) = Score2;
end
%
% Comp1 = Comp1/max(max(Comp1));
% Comp2 = Comp2/max(max(Comp2));
AllConds = {'K';'NTVK';'VK';'V'; 'T'; 'TK'; 'TV'; 'TVK'};
AllCondsNew = {'TVK';'TV';'VK';'V';'NTVK';'T';'K'; 'TK';};
IdxsYes = [];
for ii=1:8
    IdxsYes(ii) = find(strcmp(AllCondsNew{ii},AllConds));
    
end
XLabels = AllCondsNew;
YLabel = 'Score';
YLim = [0 2];
figure('units','normalized','outerposition',[0 0 1 1]);
YCoOrd = 1.25;
PerformanceShow(Comp1(:,IdxsYes), XLabels, YLabel, YLim, YCoOrd);   
title('B:NI');
figure('units','normalized','outerposition',[0 0 1 1]);
YCoOrd = 1.5;
PerformanceShow(Comp2(:,IdxsYes), XLabels, YLabel, YLim, YCoOrd);     
title('I:NI');
%
%%
T1 = array2table((round(Comp1,4)), 'VariableNames',AllCondsNew );
writetable((T1),'AllScores1.xlsx','Sheet',1);
T2 = array2table((round(Comp2,4)), 'VariableNames',AllCondsNew );
writetable((T2),'AllScores1.xlsx','Sheet',2);
%%
IndConns = {'K', 'T', 'V'};
AllConds = {'K';'NTVK';'VK';'V'; 'T'; 'TK'; 'TV'; 'TVK'};
AllP = [];
NoValues = [];
YesValues = [];
AllYes1 = [];
AllNo1 = [];
for ii=1:length(IndConns)
    IdxsNow = find(contains(AllConds, IndConns{ii}));
    IdxsNowNot = find(~contains(AllConds, IndConns{ii}));
    
    IdxsNo = [];
    IdxsNo(1) = 2;
    for jj=1:3
        IdxsNo(jj+1) = IdxsNowNot(jj);
    end
    IdxsYes = [];
    k=1;
    for jj=1:5
        IdxsNow1 = IdxsNow(jj);
        if IdxsNow1~=2
         IdxsYes(k) = IdxsNow1;
         k=k+1;
        end
    end
    AllYes = Comp2(:,IdxsYes);
    AllYes = reshape(AllYes, [numel(AllYes) 1]);
    AllYes1(:,ii) = AllYes;
    AllNo = Comp2(:,IdxsNo);
    AllNo = reshape(AllNo, [numel(AllNo) 1]);
    AllNo1(:,ii) = AllNo;
    AllP(ii) = signrank(AllYes, AllNo);
    YesValues(:,ii) = median(AllYes,1);
    NoValues(:,ii) = median(AllNo,1);
end
%% 
figure('units','normalized','outerposition',[0 0 1 1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.025], [0.075 0.075], [0.07 0.02]);
XLabels = {'Without'; 'With'};
YLabel = 'I:NI';
Combos = [1 2];
YLim = [0.3 1.4];
for jj=1:3
% Tactile, K, V
subplot(1,3,jj)
Comparison = []; Comparison(:,1) = AllNo1(:,jj); Comparison(:,2) = AllYes1(:,jj);
SignificancePlotAllShow1CZ(Comparison, Combos, XLabels, YLabel, YLim); 
if jj>1
    ylabel('');
    yticklabels([]);
else
ylabel(YLabel, 'Fontweight', 'bold');
end
set(gca, 'TickLabelInterpreter', 'tex');
xticklabels(XLabels);
yticks((0:0.2:2));
ylim([0.3 1.6]);
title(IndConns{jj});
set(gca, 'box','off');
set(gca, 'Fontsize',24);
end
%%
T0 = array2table((AllP), 'VariableNames',IndConns );
writetable((T0),'CondYesNo.xlsx','Sheet',2);
T1 = array2table(([YesValues; NoValues]), 'VariableNames',IndConns );
writetable((T1),'CondYesNo.xlsx','Sheet',1);
%% Conditions and Conditions Not
% Define Condition Idxs
AllConds = {'K';'NTVK';'VK';'V'; 'T'; 'TK'; 'TV'; 'TVK'};
AllCondsNew = {'TVK';'NTVK';'VK';'T'; 'TV';'K'; 'V';'TK';};
IdxsNew = [];
for ii=1:8
    IdxsNew(ii) = find(strcmp(AllCondsNew{ii},AllConds));
end
% Plot unimodal comparisons
figure('units','normalized','outerposition',[0 0 1 1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.025], [0.075 0.075], [0.07 0.02]);
for jj=1:4
    subplot(1,4,jj)
    Comparison = [];
    Idx1 = ((jj-1)*2)+1;
    Idx2 = Idx1+1;
    Comparison(:,1) = Comp2(:,IdxsNew(Idx1));
    Comparison(:,2) = Comp2(:,IdxsNew(Idx2));
    XLabels = {AllCondsNew{Idx1}; AllCondsNew{Idx2}};
    YLabel = 'I:NI';
    SignificancePlotAllShow1CZ(Comparison, [1 2], XLabels, YLabel, [0 1]);
    if jj>1
        ylabel('');
        yticklabels([]);
    end
    ylim([0.25 1.5]);
    yticks(0.25:0.25:1.5);
    set(gca, 'Fontsize',30);
    set(gca,'box','off');
end
%%
figure('units','normalized','outerposition',[0 0 1 1]);
Comparison = [];
jj=1;
Idx1 = ((jj-1)*2)+1;
Idx2 = Idx1+1;
Comparison(:,1) = Comp2(:,IdxsNew(Idx1));
Comparison(:,2) = Comp2(:,IdxsNew(Idx2));
XLabels = {AllCondsNew{Idx1}; AllCondsNew{Idx2}};
YLabel = 'I:NI';
SignificancePlotAllShow1CZ(Comparison, [1 2], XLabels, YLabel, [0 1]);
if jj>1
    ylabel('');
    yticklabels([]);
end
ylim([0.25 2]);
yticks(0.25:0.25:1.5);
set(gca, 'Fontsize',30);
set(gca,'box','off');
%%
BaselineConn = squeeze(ConnB);
BaselineConn1 = [];
k=1;
Muscles = {'FDI', 'APB', 'ECRL', 'FDS', 'FCU', 'EDC', 'UT', 'LT'};
Pairings = {};
for i=1:7
    for j=i+1:8
        PairingNow = [Muscles{i} ' - '  Muscles{j}];
        Pairings{k} = PairingNow;
        BaselineConn1(k,:,:,:,:) = BaselineConn(i,j,:,:,:,:);
        k=k+1;
    end
end
 %
InvolvedConn = squeeze(ConnA(:,:,:,:,:,1,:));
InvolvedConn1 = [];
k=1;
for i=1:7
    for j=i+1:8
        InvolvedConn1(k,:,:,:,:) = InvolvedConn(i,j,:,:,:,:);
        k=k+1;
    end
end

NonInvolvedConn = squeeze(ConnA(:,:,:,:,:,2,:));
NonInvolvedConn1 = [];
k=1;
for i=1:7
    for j=i+1:8
        NonInvolvedConn1(k,:,:,:,:) = NonInvolvedConn(i,j,:,:,:,:);
        k=k+1;
    end
end

% 1. Find The Mean of Degrees
BaselineConnA = permute(BaselineConn, [1 2 6 3 4 5]);
BaselineConnA = BaselineConnA(:,:,:,:);
BaselineConnA = median(BaselineConnA,4);
BaselineConnA = FillMatN(BaselineConnA,8);
BaselineConnA = (8/7)*squeeze(mean(BaselineConnA,1));

InvolvedConnA = permute(InvolvedConn, [1 2 6 3 4 5]);
InvolvedConnA = InvolvedConnA(:,:,:,:);
InvolvedConnA = median(InvolvedConnA,4);
InvolvedConnA = FillMatN(InvolvedConnA,8);
InvolvedConnA = (8/7)*squeeze(mean(InvolvedConnA,1));

NonInvolvedConnA = permute(NonInvolvedConn, [1 2 6 3 4 5]);
NonInvolvedConnA = NonInvolvedConnA(:,:,:,:);
NonInvolvedConnA = median(NonInvolvedConnA,4);
NonInvolvedConnA = FillMatN(NonInvolvedConnA,8);
NonInvolvedConnA = (8/7)*squeeze(mean(NonInvolvedConnA,1));

Comparison1 = [];
Comparison1(1:length(BaselineConnA),1) = mean(BaselineConnA,1);
Comparison1(1:length(InvolvedConnA),2) = mean(InvolvedConnA,1);
Comparison1(1:length(NonInvolvedConnA),3) = mean(NonInvolvedConnA,1);
%%
count =1;
for a=1:32
    if (Comparison1(a,2)>Comparison1(a,1))&&(Comparison1(a,3)>Comparison1(a,2))
        count=count+1;
    end
end
count =2;
for a=1:32
    if (Comparison2(a,2)>Comparison2(a,1))&&(Comparison2(a,3)>Comparison2(a,2))
        count=count+1;
    end
end
%%
YLabel = 'Mean Degree';
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
figure('units','normalized','outerposition',[0 0 1 1]);
subplot = @(m,n,p) subtightplot (m, n, p, [0 0.075], [0.05 0.05], [0.05 0.025]);
Combos = [1 2; 1 3; 2 3;];
subplot(1,4,1)
[df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison1, Combos, XLabels, YLabel,[0 1]); 
AllS = [];
AllS(1) = ComputeSensitivity(Comparison1);

set(gca, 'Fontsize',25);  
YLim = ylim;
ylim([YLim(1) YLim(2)+0.1]);


% 3. Find the Mean Shortest Path
% figure('units','normalized','outerposition',[0 0 1 1]);
SPMI;
BaselineSPA = permute(AllSPsB, [6 1 2 3 4 5]);
BaselineSPA = BaselineSPA(:,:,:,:);
BaselineSPA(isnan(BaselineSPA)) = median(BaselineSPA(:));
BaselineSPA = median(BaselineSPA,4);

InvolvedSPA = permute(squeeze(AllSPsA(:,:,:,:,:,1,:)), [6 1 2 3 4 5]);
InvolvedSPA = InvolvedSPA(:,:,:,:);
InvolvedSPA(isnan(InvolvedSPA)) = median(InvolvedSPA(:));
InvolvedSPA = median(InvolvedSPA,4);

NonInvolvedSPA = permute(squeeze(AllSPsA(:,:,:,:,:,2,:)), [6 1 2 3 4 5]);
NonInvolvedSPA = NonInvolvedSPA(:,:,:,:);
NonInvolvedSPA(isnan(NonInvolvedSPA)) = median(NonInvolvedSPA(:));
NonInvolvedSPA = median(NonInvolvedSPA,4);

BaselineSP = [];
InvolvedSP = [];
NonInvolvedSP = [];
k=1;
for i=1:7
    for j=i+1:8
        BaselineSP(:,k) = BaselineSPA(:,i,j);
        InvolvedSP(:,k) = InvolvedSPA(:,i,j);
        NonInvolvedSP(:,k) = NonInvolvedSPA(:,i,j);
        k=k+1;
    end
end
       
Comparison3 = [];
Comparison3(1:length(BaselineSP),1) = mean(BaselineSP,2);
Comparison3(1:length(InvolvedSP),2) = mean(InvolvedSP,2);
Comparison3(1:length(NonInvolvedSP),3) = mean(NonInvolvedSP,2);
YLabel = ('Mean Shortest Path');
subplot(1,4,3)
[df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison3, Combos, XLabels, YLabel,[0 20]);
AllS(3) = ComputeSensitivity(Comparison3);
set(gca, 'Fontsize',25);  
%
% 4. Find Mean Clustering Coefficient
CC1MI;
% figure('units','normalized','outerposition',[0 0 1 1]);
%
% Average MI across all edges and conditions
BaselineCC = permute(AllCCsB, [5 1 2 3 4]);
BaselineCC = BaselineCC(:,:,:);
BaselineCC(isnan(BaselineCC)) = median(BaselineCC(:));
BaselineCC = median(BaselineCC,3);
% Find MEAN across edges
BaselineCC = mean(BaselineCC,2);
InvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,1,:)), [5 1 2 3 4]);
InvolvedCC = InvolvedCC(:,:,:);
InvolvedCC(isnan(InvolvedCC)) = median(InvolvedCC(:));
InvolvedCC = median(InvolvedCC,3);
InvolvedCC = mean(InvolvedCC,2);
NonInvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,2,:)), [5 1 2 3 4]);
NonInvolvedCC = NonInvolvedCC(:,:,:);
NonInvolvedCC(isnan(NonInvolvedCC)) = median(NonInvolvedCC(:));
NonInvolvedCC = median(NonInvolvedCC,3);
NonInvolvedCC = mean(NonInvolvedCC,2);
Comparison4 = [];
Comparison4(1:length(BaselineCC),1) = BaselineCC;
Comparison4(1:length(InvolvedCC),2) = InvolvedCC;
Comparison4(1:length(NonInvolvedCC),3) = NonInvolvedCC;
YLabel = ('Mean Clustering Coefficient');
 subplot(1,4,2)
[df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison4, Combos, XLabels, YLabel,[0 1]);
AllS(2) = ComputeSensitivity(Comparison4);
set(gca, 'Fontsize',25);  

% 6. Find Global Efficiency of Network
GE1MI;
%
% Average MI across all edges and conditions
BaselineCC = permute(squeeze(AllCCsB), [4 1 2 3]);
BaselineCC = BaselineCC(:,:);
BaselineCC(isnan(BaselineCC)) = median(BaselineCC(:));
BaselineCC = median(BaselineCC,2);
% Find MEAN across edges
InvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,1,:)), [4 1 2 3]);
InvolvedCC = InvolvedCC(:,:);
InvolvedCC(isnan(InvolvedCC)) = median(InvolvedCC(:));
InvolvedCC = median(InvolvedCC,2);
NonInvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,2,:)), [4 1 2 3]);
NonInvolvedCC = NonInvolvedCC(:,:);
NonInvolvedCC(isnan(NonInvolvedCC)) = median(NonInvolvedCC(:));
NonInvolvedCC = median(NonInvolvedCC,2);
Comparison6 = [];
Comparison6(1:length(BaselineCC),1) = BaselineCC;
Comparison6(1:length(InvolvedCC),2) = InvolvedCC;
Comparison6(1:length(NonInvolvedCC),3) = NonInvolvedCC;
YLabel = ('Global Efficiency');
subplot(1,4,4)
[df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison6, Combos, XLabels, YLabel,[0 1]);
AllS(4) = ComputeSensitivity(Comparison6);
set(gca, 'Fontsize',25);  
Metrics = {'Mean Deg'; 'Mean CC'; 'Mean SP'; 'GE'};
T2 = array2table(AllS, 'VariableNames',Metrics );
writetable(rows2vars(T2),'AllSensitivity.xlsx','Sheet',3);
%% 5. Find Path Length
% 2. Find the Mean of Edges 
% figure('units','normalized','outerposition',[0 0 1 1]);
BaselineConn = permute(BaselineConn1, [5 1 2 3 4]);
BaselineConn = BaselineConn(:,:,:);
BaselineConn(isnan(BaselineConn)) = 0;
BaselineConn = median(BaselineConn,3);
% Find MEAN across edges
BaselineConn = mean(BaselineConn,2);
InvolvedConn = permute(InvolvedConn1, [5 1 2 3 4]);
InvolvedConn = InvolvedConn(:,:,:);
InvolvedConn(isnan(InvolvedConn)) = 0;
InvolvedConn = median(InvolvedConn,3);
InvolvedConn = mean(InvolvedConn,2);
NonInvolvedConn = permute(NonInvolvedConn1, [5 1 2 3 4]);
NonInvolvedConn = NonInvolvedConn(:,:,:);
NonInvolvedConn(isnan(NonInvolvedConn)) = 0;
NonInvolvedConn = median(NonInvolvedConn,3);
NonInvolvedConn = mean(NonInvolvedConn,2);
Comparison2 = [];
Comparison2(1:length(BaselineConn),1) = BaselineConn;
Comparison2(1:length(InvolvedConn),2) = InvolvedConn;
Comparison2(1:length(NonInvolvedConn),3) = NonInvolvedConn;
YLabel = 'Mean Connectivity';
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
% subplot(1,5,1)
% [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',25);  
YLim = ylim;
ylim([YLim(1) YLim(2)+0.1]);
%%
figure('units','normalized','outerposition',[0 0 1 1]);
PathLengthRoryXMI;
%
% Average MI across all edges and conditions
BaselinePL = permute(AllPLsB, [4 3 2 1]);
BaselinePL = BaselinePL(:,:,:);
BaselinePL(isnan(BaselinePL)) = median(BaselinePL(:));
BaselinePL = median(BaselinePL,3);
% Find MEAN across edges
BaselinePL = mean(BaselinePL,2);
InvolvedPL = permute(squeeze(AllPLsA(:,:,:,1,:)), [4 3 2 1]);
InvolvedPL = InvolvedPL(:,:,:);
InvolvedPL(isnan(InvolvedPL)) = median(InvolvedPL(:));
InvolvedPL = median(InvolvedPL,3);
InvolvedPL = mean(InvolvedPL,2);
NonInvolvedPL = permute(squeeze(AllPLsA(:,:,:,2,:)), [4 3 2 1]);
NonInvolvedPL = NonInvolvedPL(:,:,:);
NonInvolvedPL(isnan(NonInvolvedPL)) = median(NonInvolvedPL(:));
NonInvolvedPL = median(NonInvolvedPL,3);
NonInvolvedPL = mean(NonInvolvedPL,2);
Comparison5 = [];
Comparison5(1:length(BaselinePL),1) = BaselinePL;
Comparison5(1:length(InvolvedPL),2) = InvolvedPL;
Comparison5(1:length(NonInvolvedPL),3) = NonInvolvedPL;
YLabel = ('Average Path Length');
% subplot(1,3,3)
Comparison5(Comparison5>1) = median(Comparison5(:));
SignificancePlotAllShow1(Comparison5, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);
%


%%
WCC1MI;
figure('units','normalized','outerposition',[0 0 1 1]);
%
% Average MI across all edges and conditions
BaselineCC = permute(AllCCsB, [5 1 2 3 4]);
BaselineCC = BaselineCC(:,:,:);
BaselineCC(isnan(BaselineCC)) = median(BaselineCC(:));
BaselineCC = median(BaselineCC,3);
% Find MEAN across edges
BaselineCC = mean(BaselineCC,2);
InvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,1,:)), [5 1 2 3 4]);
InvolvedCC = InvolvedCC(:,:,:);
InvolvedCC(isnan(InvolvedCC)) = median(InvolvedCC(:));
InvolvedCC = median(InvolvedCC,3);
InvolvedCC = mean(InvolvedCC,2);
NonInvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,2,:)), [5 1 2 3 4]);
NonInvolvedCC = NonInvolvedCC(:,:,:);
NonInvolvedCC(isnan(NonInvolvedCC)) = median(NonInvolvedCC(:));
NonInvolvedCC = median(NonInvolvedCC,3);
NonInvolvedCC = mean(NonInvolvedCC,2);
Comparison2(1:length(BaselineCC),1) = BaselineCC;
Comparison2(1:length(InvolvedCC),2) = InvolvedCC;
Comparison2(1:length(NonInvolvedCC),3) = NonInvolvedCC;
YLabel = ('MI WCC (BCT)');
% subplot(1,3,2)
SignificancePlotAllShow1(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);
% subplot(1,3,3): Path Length


%%
figure('units','normalized','outerposition',[0 0 1 1]);
BCRoryX;
%
% Average MI across all edges and conditions
BaselinePL = permute(AllPLsB, [4 3 2 1]);
BaselinePL = BaselinePL(:,:,:);
BaselinePL(isnan(BaselinePL)) = median(BaselinePL(:));
BaselinePL = median(BaselinePL,3);
% Find MEAN across edges
BaselinePL = mean(BaselinePL,2);
InvolvedPL = permute(squeeze(AllPLsA(:,:,:,1,:)), [4 3 2 1]);
InvolvedPL = InvolvedPL(:,:,:);
InvolvedPL(isnan(InvolvedPL)) = median(InvolvedPL(:));
InvolvedPL = median(InvolvedPL,3);
InvolvedPL = mean(InvolvedPL,2);
NonInvolvedPL = permute(squeeze(AllPLsA(:,:,:,2,:)), [4 3 2 1]);
NonInvolvedPL = NonInvolvedPL(:,:,:);
NonInvolvedPL(isnan(NonInvolvedPL)) = median(NonInvolvedPL(:));
NonInvolvedPL = median(NonInvolvedPL,3);
NonInvolvedPL = mean(NonInvolvedPL,2);
Comparison2(1:length(BaselinePL),1) = BaselinePL;
Comparison2(1:length(InvolvedPL),2) = InvolvedPL;
Comparison2(1:length(NonInvolvedPL),3) = NonInvolvedPL;
YLabel = ('Average Path Length');
subplot(1,3,3)
SignificancePlotAllShow1(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);

% subplot(1,3,3): Path Length
% figure('units','normalized','outerposition',[0 0 1 1]);
PathLengthRoryXMI;
%

%%
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.08], [0.11 0.03], [0.07 0.01]);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1) 
BaselineConn = permute(BaselineConn1, [5 1 2 3 4]);
BaselineConn = BaselineConn(:,:,:);
BaselineConn(isnan(BaselineConn)) = 0;
BaselineConn = median(BaselineConn,3);
% Find MEAN across edges
BaselineConn = mean(BaselineConn,2);
InvolvedConn = permute(InvolvedConn1, [5 1 2 3 4]);
InvolvedConn = InvolvedConn(:,:,:);
InvolvedConn(isnan(InvolvedConn)) = 0;
InvolvedConn = median(InvolvedConn,3);
InvolvedConn = mean(InvolvedConn,2);
NonInvolvedConn = permute(NonInvolvedConn1, [5 1 2 3 4]);
NonInvolvedConn = NonInvolvedConn(:,:,:);
NonInvolvedConn(isnan(NonInvolvedConn)) = 0;
NonInvolvedConn = median(NonInvolvedConn,3);
NonInvolvedConn = mean(NonInvolvedConn,2);
Comparison2(1:length(BaselineConn),1) = BaselineConn;
Comparison2(1:length(InvolvedConn),2) = InvolvedConn;
Comparison2(1:length(NonInvolvedConn),3) = NonInvolvedConn;
YLabel = 'MI: Degree';
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
SignificancePlotAllShow1(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);
YLim = ylim;
ylim([YLim(1) YLim(2)+0.1]);

WCC1MI;
%
% Average MI across all edges and conditions
BaselineCC = permute(AllCCsB, [5 1 2 3 4]);
BaselineCC = BaselineCC(:,:,:);
BaselineCC(isnan(BaselineCC)) = median(BaselineCC(:));
BaselineCC = median(BaselineCC,3);
% Find MEAN across edges
BaselineCC = mean(BaselineCC,2);
InvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,1,:)), [5 1 2 3 4]);
InvolvedCC = InvolvedCC(:,:,:);
InvolvedCC(isnan(InvolvedCC)) = median(InvolvedCC(:));
InvolvedCC = median(InvolvedCC,3);
InvolvedCC = mean(InvolvedCC,2);
NonInvolvedCC = permute(squeeze(AllCCsA(:,:,:,:,2,:)), [5 1 2 3 4]);
NonInvolvedCC = NonInvolvedCC(:,:,:);
NonInvolvedCC(isnan(NonInvolvedCC)) = median(NonInvolvedCC(:));
NonInvolvedCC = median(NonInvolvedCC,3);
NonInvolvedCC = mean(NonInvolvedCC,2);
Comparison2(1:length(BaselineCC),1) = BaselineCC;
Comparison2(1:length(InvolvedCC),2) = InvolvedCC;
Comparison2(1:length(NonInvolvedCC),3) = NonInvolvedCC;
YLabel = ('MI: Weighted CC');
subplot(1,3,2)
SignificancePlotAllShow1(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);


% Average MI across all edges and conditions
BaselinePL = permute(AllPLsB, [4 3 2 1]);
BaselinePL = BaselinePL(:,:,:);
BaselinePL(isnan(BaselinePL)) = median(BaselinePL(:));
BaselinePL = median(BaselinePL,3);
% Find MEAN across edges
BaselinePL = mean(BaselinePL,2);
InvolvedPL = permute(squeeze(AllPLsA(:,:,:,1,:)), [4 3 2 1]);
InvolvedPL = InvolvedPL(:,:,:);
InvolvedPL(isnan(InvolvedPL)) = median(InvolvedPL(:));
InvolvedPL = median(InvolvedPL,3);
InvolvedPL = mean(InvolvedPL,2);
NonInvolvedPL = permute(squeeze(AllPLsA(:,:,:,2,:)), [4 3 2 1]);
NonInvolvedPL = NonInvolvedPL(:,:,:);
NonInvolvedPL(isnan(NonInvolvedPL)) = median(NonInvolvedPL(:));
NonInvolvedPL = median(NonInvolvedPL,3);
NonInvolvedPL = mean(NonInvolvedPL,2);
Comparison2(1:length(BaselinePL),1) = BaselinePL;
Comparison2(1:length(InvolvedPL),2) = InvolvedPL;
Comparison2(1:length(NonInvolvedPL),3) = NonInvolvedPL;
YLabel = ('MI: Average Path Length');
subplot(1,3,3)
Comparison2(Comparison2>1) = median(Comparison2(:));
SignificancePlotAllShow1(Comparison2, Combos, XLabels, YLabel,[]); 
set(gca, 'Fontsize',30);
%%
BaselineConn = permute(BaselineConn1, [5 4 1 2 3]);
BaselineConn = BaselineConn(:,:,:,:);
BaselineConn(isnan(BaselineConn)) = median(BaselineConn(:));
% Find MEDIAN across frequency components, weights and trials
BaselineConn = median(BaselineConn,4);
% Find MEAN across edges as in the paper
BaselineConn = mean(BaselineConn,3);
InvolvedConn = permute(InvolvedConn1, [5 4 1 2 3]);
InvolvedConn = InvolvedConn(:,:,:,:);
InvolvedConn(isnan(InvolvedConn)) = median(InvolvedConn(:));
% Find MEDIAN across frequency dimensions and weights
InvolvedConn = median(InvolvedConn,4);
% Find MEAN across edges as in the paper
InvolvedConn = mean(InvolvedConn,3);
NonInvolvedConn = permute(NonInvolvedConn1, [5 4 1 2 3]);
NonInvolvedConn = NonInvolvedConn(:,:,:,:);
NonInvolvedConn(isnan(NonInvolvedConn)) = median(NonInvolvedConn(:));
% Find MEDIAN across frequency dimensions and weights
NonInvolvedConn = median(NonInvolvedConn,4);
% Find MEAN across edges as in the paper
NonInvolvedConn = mean(NonInvolvedConn,3);
%%
subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.01], [0.1 0.1], [0.07 0.01]);
XLabels = {'\bf B'; '\bf I'; '\bf NI'; };
YLabel = 'MI Mean';
Combos = [1 2; 2 3; 1 3];
figure('units','normalized','outerposition',[0 0 1 1]);
Comparison2 = [];
AllConds = {'K';'NTVK';'V';'VK'; 'T'; 'TK'; 'TV'; 'TVK'};
YLim = [0 0.9];  
for c=1:4
    subplot(1,4,c)
    Comparison2 = [];
    Condition = BaselineConn(:,c);
    BaselineC = Condition;
    % Involved
    Condition = InvolvedConn(:,c);
    InvolvedC = Condition;
    % Non-involved
    Condition = NonInvolvedConn(:,c);
    NonInvolvedC = Condition;
    Comparison2(1:length(BaselineC),1) = BaselineC;
    Comparison2(1:length(InvolvedC),2) = InvolvedC;
    Comparison2(1:length(NonInvolvedC),3) = NonInvolvedC;
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel,YLim); 
    set(gca, 'Fontsize',30);
    title(AllConds{c});
    if c>1
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
end
%
figure('units','normalized','outerposition',[0 0 1 1]);
for c=5:8
    subplot(1,4,c-4)
    Comparison2 = [];
    Condition = BaselineConn(:,c);
    BaselineC = Condition;
    % Involved
    Condition = InvolvedConn(:,c);
    InvolvedC = Condition;
    % Non-involved
    Condition = NonInvolvedConn(:,c);
    NonInvolvedC = Condition;
    Comparison2(1:length(BaselineC),1) = BaselineC;
    Comparison2(1:length(InvolvedC),2) = InvolvedC;
    Comparison2(1:length(NonInvolvedC),3) = NonInvolvedC;
    [df, Chi_sq, p] = SignificancePlotAllShow1C(Comparison2, Combos, XLabels, YLabel,YLim); 
    set(gca, 'Fontsize',30);
    title(AllConds{c});
    if c>5
%         yticks([]);
        ylabel([]);
        set(gca, 'yticklabels', []);
    end
end
