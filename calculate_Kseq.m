% Calculate K
% Mi & Ti:the vectors of events
% T0=0

function K = calculate_Kseq (Mi, Ti, Iri, tstart, Ml, theta, Nbo)
intLambda = 0;
No = length(Ti);
% Divided by 2 parts and first is tj~tj-1, second is tNo~tstart
% Then sum them up
for j=2:length(Ti)
    intLambda = intLambda+int3LambdaETAS (Mi, Ti, Iri, Ti(j-1), Ti(j), Ml, [theta(1);1;theta(3:end)]);
end
intLambda = intLambda+int3LambdaETAS (Mi, Ti, Iri, Ti(end), tstart, Ml, [theta(1);1;theta(3:end)]);

K = (No-Nbo)/sum(intLambda);

end
     