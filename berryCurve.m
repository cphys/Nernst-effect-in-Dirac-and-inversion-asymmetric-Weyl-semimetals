ClearAll["Global`*"]


(* Plot Parameters *)
kzValMin = 55/100;
kzValMax = 82/100;
kxValMin = -50/100;
kxValMax = 50/100;
fontSize = 40;
strengths = Table[2^i, {i, -4, 4, 1/2}] // N;
bounds = Partition[
   Flatten@{0, MovingAverage[strengths, 2], \[Infinity]}, 2, 1];
imageSize = 600;
perfGoal = "Speed";
streamColor = Blue;
fLabel = {"\!\(\*SubscriptBox[\(k\), \(y\)]\)", "\!\(\*SubscriptBox[\(k\), \(z\)]\)"};


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

(* Pauli matricies *)
s0 = PauliMatrix[0];
s1 = PauliMatrix[1];
s2 = PauliMatrix[2];
s3 = PauliMatrix[3];
tzs0 = KroneckerProduct[s3, s0];
txsz = KroneckerProduct[s1, s3];
tys0 = KroneckerProduct[s2, s0];
txsx = KroneckerProduct[s1, s1];
txsy = KroneckerProduct[s1, s2];
t0sz = KroneckerProduct[s0, s3];


(* Hamiltonian *)
NH4Lat[kx_, ky_, kz_, Bz_] := NH4Lat[kx, ky, kz, Bz] =
  a[kx, ky, kz]*tzs0 + b[kx, ky, kz]*txsz + c[kx, ky, kz]*tys0 + d[kx, ky, kz]*txsx + ee[kx, ky, kz]*txsy + Bz*t0sz


(* Calculate eigenvalues *)
eigSys[kx_, ky_, kz_, Bz_] = Eigensystem[NH4Lat[kx, ky, kz, Bz]];
eigs[kx_, ky_, kz_, Bz_] = Indexed[eigSys[kx, ky, kz, Bz], 1];
eigVec[kx_, ky_, kz_, Bz_] = Indexed[eigSys[kx, ky, kz, Bz], 2];


(* Calculate Derivitives for use in Berry Curvature *)
dHkx[kx_, ky_, kz_, Bz_] = D[NH4Lat[kx, ky, kz, Bz], kx];
dHky[kx_, ky_, kz_, Bz_] = D[NH4Lat[kx, ky, kz, Bz], ky];
dHkz[kx_, ky_, kz_, Bz_] = D[NH4Lat[kx, ky, kz, Bz], kz];


\[CapitalOmega]yz[kx_, ky_, kz_, Bz_] := \[CapitalOmega]yz[kx, ky, kz, Bz] = (
   Eigsys = eigSys[kx, ky, kz, Bz];
   Dhx = dHkx[kx, ky, kz, Bz];
   Dhy = dHky[kx, ky, kz, Bz];
   Dhz = dHkz[kz, ky, kz, Bz];
   Energy = Re[Indexed[Eigsys, 1]];
   eigenVec = Indexed[Eigsys, 2];
   Mx = Table[
     Conjugate[Indexed[eigenVec, i]].Dhx.Indexed[eigenVec, j], {i, 1, 4}, {j, 1, 4}];
   My = Table[
     Conjugate[Indexed[eigenVec, i]].Dhy.Indexed[eigenVec, j], {i, 1, 4}, {j, 1, 4}];
   Mz = Table[
     Conjugate[Indexed[eigenVec, i]].Dhz.Indexed[eigenVec, j], {i, 1, 4}, {j, 1, 4}];
   i\[CapitalOmega]y = Re[Table[
      I*Sum[
        If[Or[
          i != j,
          Indexed[Energy, i] != Indexed[Energy, j]],
         (Indexed[Mz, {i, j}]*Indexed[Mx, {j, i}] - Indexed[Mx, {i, j}]* Indexed[Mz, {j, i}])/((Indexed[Energy, i] - Indexed[Energy, j])^2),
         0],
        {j, 1, 4}], {i, 1, 4}]];
   i\[CapitalOmega]z = Re[Table[
      I*Sum[
        If[Or[
          i != j,
          Indexed[Energy, i] != Indexed[Energy, j]],
         (Indexed[Mx, {i, j}]*Indexed[My, {j, i}] - Indexed[My, {i, j}]* Indexed[Mx, {j, i}])/((Indexed[Energy, i] - Indexed[Energy, j])^2),
         0],
        {j, 1, 4}], {i, 1, 4}]];
   {i\[CapitalOmega]y, i\[CapitalOmega]z})


berryPlot = Show[
   MapThread[
    StreamPlot[{
       \[CapitalOmega]yz[0, ky, kz, bz][[1, 1]] + \[CapitalOmega]yz[0, ky, kz, bz][[1, 3]],
       \[CapitalOmega]yz[0, ky, kz, bz][[2, 1]] + \[CapitalOmega]yz[0, ky, kz, bz][[2, 3]]},
      {ky, -1, 1}, {kz, -1, 1},
      PlotTheme -> "Scientific",
      FrameLabel -> fLabel,
      ImageSize -> imageSize,
      LabelStyle -> {FontFamily -> "Latex", FontSize -> fontSize},
      StreamPoints -> {Automatic, 0.5/#1},
      PerformanceGoal -> perfGoal,
      FrameStyle -> Black,
      StreamStyle -> streamColor,
      Background -> Automatic,
      RegionFunction -> 
       Function[{x, y, vx, vy, n}, First@#2 <= n <= Last@#2]
      ] &, {strengths, bounds}],
   PlotRange -> {{-1, 1}, {-1, 1}}];

Export[NotebookDirectory[] <> "berryCurve.png", berryPlot];
