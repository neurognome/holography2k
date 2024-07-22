function Get = function_Eval2DCoC(CoC,Ask)

FitX = CoC.FitX;
FitY = CoC.FitY;


GetX = polyvaln(FitX,Ask);
GetY = polyvaln(FitY,Ask);

Get = [GetX GetY]; 