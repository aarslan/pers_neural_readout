function  [BSthreshold] = permutLabelClassif(Ytest,YTestPredict,BSiter,threshold) 
%convert YTestPredict en -1/ 1
if ~isint(YTestPredict(1))
YTestPredict = sign(YTestPredict);
end
%estimate a statistical threshold above chance level for binary classification
        randAcc=[];
        for perm=1:BSiter
            YtestRand=Ytest(randperm(length(Ytest)));
            perfTestRand=sum(YTestPredict==YtestRand)/length(Ytest);
            randAcc=[randAcc, perfTestRand];
        end
        BSthreshold=quantile(randAcc,threshold);
end