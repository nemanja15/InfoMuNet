% load('Xmas\Mat Files\With Decimated\Using\AllBaselinesDecimatedXR32.mat')
Subjects = [1248,1275,1332,1460,1506,1516,1542,1553,1577, ...
1603,1610,1208,1291,1412,1458,1466,1487,1494 ...
1499,1509,1512,1514,1517,1521,1536,1558,1562,1579,1580,1581,1591,1599];
%
T = readtable('Strategies for dexterity_demographics.xlsx');
SubsCell = table2cell(T);
for s=1:32
        InvolvedSides{s,1} = FindInvolvedSide(SubsCell, Subjects(s));
end

AllMIsB = [];
%
n_fft = 1024;
T=1;
EmptyTrials = {};
%%
for sub=1:length(Subjects)
    if length(InvolvedSides{sub,1})==length('Left')
        idx1 = 1;
        idx2 = 9;
    else
        idx1 = 9;
        idx2 = 1;
    end
        
        for cond =1:8
            for LH=1:2
                for trial = 1:7
                    % Involved
                  Session = squeeze(AllBaselinesD(:,idx1:idx1+7,trial,LH,cond,sub));
                  Session1 = Session(:,1);
                  Clipped = ClipVectorFromEnd(Session1);
                  Session = Session(1:length(Clipped),:);
                   if length(Clipped) < n_fft
                      RandomTrial = ChooseRandomTrial(trial);  
                      Session = squeeze(AllBaselinesD(:,idx1:idx1+7,RandomTrial,LH,cond,sub));
                      Session1 = Session(:,1);
                      Clipped = ClipVectorFromEnd(Session1);
                      Session = Session(1:length(Clipped),:);
                   end
                  for i=1:7
                          for j=i+1:8
                        MINow = RoryMI((Session(:,i)), (Session(:,j)), 'eqspace');
                        AllMIsB(i,j,trial, LH,cond,sub) = MINow;
                          end
                  end
                end
            end
        end
end
%
save('AllMIsB_X32.mat', 'AllMIsB', 'EmptyTrials');
%%
% load('Xmas\Mat Files\With Decimated\Using\AllAlternatingsDecimatedXR32.mat')
AllMIsA = [];
L=1;
AltLengths = [];
EmptyTrials = {};
T = 1;
%
for sub=1:length(Subjects)
    if length(InvolvedSides{sub,1})==length('Left')
        idx1 = 1;
        idx2 = 9;
    else
        idx1 = 9;
        idx2 = 1;
    end
        
    for inv=1:2
        for cond =1:8
            for LH=1:2
                for trial = 1:7
                    % Involved
                  Session = squeeze(AllAlternatingsD(:,idx1:idx1+7,trial,LH,cond,inv,sub));
                  Session1 = Session(:,1);
                  Clipped = ClipVectorFromEnd(Session1);
                  if length(Clipped) < n_fft
                      RandomTrial = ChooseRandomTrial(trial);  
                      Session = squeeze(AllAlternatingsD(:,idx1:idx1+7,RandomTrial,LH,cond,inv,sub));
                      Session1 = Session(:,1);
                      Clipped = ClipVectorFromEnd(Session1);
                      if length(Clipped) < n_fft
                      RandomTrial = ChooseRandomTrial(trial);  
                      Session = squeeze(AllAlternatingsD(:,idx1:idx1+7,RandomTrial,LH,cond,inv,sub));
                      Session1 = Session(:,1);
                      Clipped = ClipVectorFromEnd(Session1);
                      end
                      Session = Session(1:length(Clipped),:);
                  end
                    Session1 = Session(:,1);
                    Clipped = ClipVectorFromEnd(Session1);
                    if length(Clipped) < n_fft
                        EmptyTrials{T} = [sub,inv,cond,LH,trial];
                        T = T+1;
                        continue;
                    end
                  for i=1:7
                      for j=i+1:8
                        MINow = RoryMI((Session(:,i)), (Session(:,j)), 'eqspace');
                        AllMIsA(i,j, trial,LH,cond,inv,sub) = MINow;
                      end
                  end
                end
            end
        end
    end
end


%
AllMIsA = single(AllMIsA);
save('AllMIsA_X32.mat', 'AllMIsA')
%
