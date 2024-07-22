function [SLMXYZP] = function_SItoSLM2D(SIXYZ,CoC)
%Function takes coordinates in SI space XY and optotuneZ (n x 3 matrix)
%and CoC as created in alignSLMtoCam
%Outputs in SLM coordinates XYZ and power normalization (n x 4 Matrix)

estSLM = function_Eval2DCoC(CoC.SItoSLM,SIXYZ);

SLMXYZP = [estSLM ones(size(estSLM,1), 1)];