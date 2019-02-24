ClearAll["Global`*"]

fileName = NotebookDirectory[] <> "eigenvalues/";
exportfileFunction[zeemanField_] := 
 fileName <> "eigenvaluesZem" <> ToString[N[zeemanField*100]] <> "png"
finalPlot = True;
exportAnimatedGif = True;
exportListOfAllPlots = False;
exportSinglePlot = False;

(* Plot Parameters *)
is = 1000;
zemPlotVal = 0.5;
zemRangeMin = 0;
zemRangeMax = 1;
numberOfPlots = 20;
labStyle = {FontFamily -> "Latex", FontSize -> 35};
kxRange = 1;
kzRange = 1;
enRange = 2;
pRangeScale = 0.1;
vPoint = {1.3, -.75, 0};
tix = {{-1, 0, 1}, Automatic, Automatic};
axLabs = {"\!\(\*SubscriptBox[\(k\), \(x\)]\)", "\!\(\*SubscriptBox[\(k\), \(z\)]\)", "E"};
boxRatios = {.125, 1, 1};
pStyle = {
   {Yellow, Opacity[.65]},
   {Pink, Opacity[.65]},
   {Blue, Opacity[1]},
   {Green, Opacity[1]}};

(* adjusts resolution of final plots based on Plotting Options *)
If[finalPlot,
  {meshPoints = 35,
   plotPoints = 100,
   perfGoal = "Quality"},
  {meshPoints = 10,
   plotPoints = 25,
   perfGoal = "Speed"},
  {meshPoints = 10,
   plotPoints = 25,
   perfGoal = "Speed"}];

(* Hamiltonian parameters *)
mo = -2;
m1 = -4;
m2 = -1/5;
\[Mu] = 1/2;
\[Eta] = 1;
\[Beta] = -1/5;
\[Gamma] = 1;
bz = 3/5;
a[kx_, ky_, kz_] := a[kx, ky, kz] = mo - m1*kz^2 - m2*(kx^2 + ky^2)
b[kx_, ky_, kz_] := b[kx, ky, kz] = \[Eta]*kx
c[kx_, ky_, kz_] := c[kx, ky, kz] = -\[Eta]*ky
d[kx_, ky_, kz_] := d[kx, ky, kz] = (\[Beta] + \[Gamma])*kz*(ky^2 + kx^2)
ee[kx_, ky_, kz_] := ee[kx, ky, kz] = -2*(\[Beta] - \[Gamma])*kx*ky*kz

(* Hamiltonian *)
NH4Lat[kx_, ky_, kz_, Bz_] := NH4Lat[kx, ky, kz, Bz] =
  a[kx, ky, kz]*tzs0 + b[kx, ky, kz]*txsz + c[kx, ky, kz]*tys0 + d[kx, ky, kz]*txsx + ee[kx, ky, kz]*txsy + Bz*t0sz


(* Calculate eigenvalues *)
eigSys[kx_, ky_, kz_, Bz_] := eigSys[kx, ky, kz, Bz] = Eigensystem[NH4Lat[kx, ky, kz, Bz]];
eigs[kx_, ky_, kz_, Bz_] = Indexed[eigSys[kx, ky, kz, Bz], 1];

(* Create the Plot *)
plotWeyl[zeemanField_] := Plot3D[{
   eigs[kx, 0, kz, zeemanField][[1]],
   eigs[kx, 0, kz, zeemanField][[2]],
   eigs[kx, 0, kz, zeemanField][[3]],
   eigs[kx, 0, kz, zeemanField][[4]]},
  {kx, -kxRange, kxRange},
  {kz, -kzRange, kzRange},
  BoxRatios -> boxRatios,
  PlotTheme -> "Detailed",
  PlotLegends -> None,
  AxesLabel -> axLabs,
  LabelStyle -> labStyle,
  Ticks -> tix,
  ViewPoint -> vPoint,
  PlotRange -> {
    {-(kxRange + pRangeScale*kxRange), (kxRange + pRangeScale*kxRange)},
    {-(kzRange + pRangeScale*kzRange), (kxRange + pRangeScale*kzRange)},
    {-enRange, enRange}},
  PlotStyle -> pStyle,
  ImageSize -> is,
  Mesh -> {meshPoints, meshPoints},
  PerformanceGoal -> perfGoal,
  PlotPoints -> plotPoints]


(* Export plot *)
exportPlots[zeemanField_] := Export[exportfileFunction[zeemanField], plotWeyl[zeemanField]];

(* Export Animated Gif *)
If[exportAnimatedGif, {
   tableOfPlots = ParallelTable[plotWeyl[iZem], {iZem, zemRangeMin, zemRangeMax, (zemRangeMax - zemRangeMin)/numberOfPlots}],
   Export[fileName <> "eigenValues.gif", tableOfPlots]}];

If[exportListOfAllPlots,
  Table[exportPlots[iZem], {iZem, zemRangeMin, zemRangeMax, (zemRangeMax - zemRangeMin)/numberOfPlots}]];

If[exportSinglePlot,
 exportPlots[zemPlotVal]]

























































