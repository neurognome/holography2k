function CoC = function_2DCoC(refAsk,refGet,modelterms)

FitX =  polyfitn(refAsk,refGet(:,1),modelterms);
FitY =  polyfitn(refAsk,refGet(:,2),modelterms);
% FitZ =  polyfitn(refAsk,refGet(:,3),modelterms);

CoC.FitX = FitX;
CoC.FitY = FitY;
% CoC.FitZ = FitZ;