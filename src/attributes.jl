
include(".//Intgl.jl")
include(".//molecule.jl")
function Basis_attributes_finder(atom)
	basis_set_STO3G = Dict("H" =>  [[[3.425250914, 0.6239137298, 0.168855404], [0.1543289673, 0.5353281423, 0.4446345422], [0, 0, 0]]], 
		"He" =>  [[[6.362421394, 1.158922999, 0.3136497915], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)]],
		"Li" => [[[16.11957475, 2.936200663, 0.794650487], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		[[0.6362897469, 0.1478600533, 0.0480886784], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[0.6362897469, 0.1478600533, 0.0480886784], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[0.6362897469, 0.1478600533, 0.0480886784], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[0.6362897469, 0.1478600533, 0.0480886784], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]],
		"Be" =>  [[[30.16787069, 5.495115306, 1.487192653], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		[[1.31483311, 0.3055389383, 0.0993707456], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)],
		[[1.31483311, 0.3055389383, 0.0993707456], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)],
		[[1.31483311, 0.3055389383, 0.0993707456], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)],
		[[1.31483311, 0.3055389383, 0.0993707456], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]],
		"B" =>  [[[48.79111318, 8.887362172, 2.40526704], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		[[2.236956142, 0.5198204999, 0.16906176], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)],
		[[2.236956142, 0.5198204999, 0.16906176], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)],
		[[2.236956142, 0.5198204999, 0.16906176], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)],
		[[2.236956142, 0.5198204999, 0.16906176], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]], 
		"C" =>  [[[71.61683735, 13.04509632, 3.53051216], [0.1543289673, 0.5353281423, 0.4446345422], [0, 0, 0]],
		[[2.941249355, 0.6834830964, 0.2222899159], [-0.09996722919, 0.3995128261, 0.7001154689], [0, 0, 0]],
		[[2.941249355, 0.6834830964, 0.2222899159], [0.155916275, 0.6076837186, 0.3919573931], [1, 0, 0]], 
		[[2.941249355, 0.6834830964, 0.2222899159], [0.155916275, 0.6076837186, 0.3919573931], [0, 1, 0]], 
		[[2.941249355, 0.6834830964, 0.2222899159], [0.155916275, 0.6076837186, 0.3919573931], [0, 0, 1]]],
		"N" => [[[99.10616896, 18.05231239, 4.885660238], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		[[3.780455879, 0.8784966449, 0.2857143744], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[3.780455879, 0.8784966449, 0.2857143744], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)],
		[[3.780455879, 0.8784966449, 0.2857143744], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)],
		 [[3.780455879, 0.8784966449, 0.2857143744], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]], 
		 "O" => [[[130.7093214, 23.80886605, 6.443608313], [0.1543289673, 0.5353281423, 0.4446345422], [0, 0, 0]], 
		[[5.033151319, 1.169596125, 0.38038896], [-0.09996722919, 0.3995128261, 0.7001154689], [0, 0, 0]],
		[[5.033151319, 1.169596125, 0.38038896], [0.155916275, 0.6076837186, 0.3919573931], [1, 0, 0]],
		[[5.033151319, 1.169596125, 0.38038896], [0.155916275, 0.6076837186, 0.3919573931], [0, 1, 0]], 
		[[5.033151319, 1.169596125, 0.38038896], [0.155916275, 0.6076837186, 0.3919573931], [0, 0, 1]]],
		"F" => [[[166.679134, 30.36081233, 8.216820672], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		[[6.464803249, 1.502281245, 0.4885884864], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[6.464803249, 1.502281245, 0.4885884864], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[6.464803249, 1.502281245, 0.4885884864], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[6.464803249, 1.502281245, 0.4885884864], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]], 
		"Ne" => [[[207.015607, 37.70815124, 10.20529731], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[8.24631512, 1.916266291, 0.6232292721], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[8.24631512, 1.916266291, 0.6232292721], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[8.24631512, 1.916266291, 0.6232292721], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[8.24631512, 1.916266291, 0.6232292721], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)]], 
		"Na" => [[[250.77243, 45.67851117, 12.36238776], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[12.04019274, 2.797881859, 0.909958017], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[12.04019274, 2.797881859, 0.909958017], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[12.04019274, 2.797881859, 0.909958017], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[12.04019274, 2.797881859, 0.909958017], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[1.478740622, 0.4125648801, 0.1614750979], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"Mg" => [[[299.2374137, 54.50646845, 14.75157752], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[15.12182352, 3.513986579, 1.142857498], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[15.12182352, 3.513986579, 1.142857498], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[15.12182352, 3.513986579, 1.142857498], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[15.12182352, 3.513986579, 1.142857498], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"Al" =>  [[[351.4214767, 64.01186067, 17.32410761], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[18.89939621, 4.391813233, 1.42835397], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[18.89939621, 4.391813233, 1.42835397], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[18.89939621, 4.391813233, 1.42835397], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[18.89939621, 4.391813233, 1.42835397], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[1.395448293, 0.3893265318, 0.1523797659], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"Si" =>  [[[407.7975514, 74.28083305, 20.10329229], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)],
		 [[23.19365606, 5.389706871, 1.752899952], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		 [[23.19365606, 5.389706871, 1.752899952], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		 [[23.19365606, 5.389706871, 1.752899952], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		 [[23.19365606, 5.389706871, 1.752899952], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		 [[1.478740622, 0.4125648801, 0.1614750979], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		 [[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		 [[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		 [[1.478740622, 0.4125648801, 0.1614750979], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]],
		"P" => [[[468.3656378, 85.31338559, 23.08913156], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[28.03263958, 6.514182577, 2.118614352], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[28.03263958, 6.514182577, 2.118614352], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[28.03263958, 6.514182577, 2.118614352], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[28.03263958, 6.514182577, 2.118614352], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[1.743103231, 0.4863213771, 0.1903428909], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[1.743103231, 0.4863213771, 0.1903428909], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[1.743103231, 0.4863213771, 0.1903428909], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[1.743103231, 0.4863213771, 0.1903428909], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"S" => [[[533.1257359, 97.1095183, 26.28162542], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[33.32975173, 7.745117521, 2.518952599], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[33.32975173, 7.745117521, 2.518952599], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[33.32975173, 7.745117521, 2.518952599], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[33.32975173, 7.745117521, 2.518952599], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[2.029194274, 0.5661400518, 0.2215833792], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[2.029194274, 0.5661400518, 0.2215833792], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[2.029194274, 0.5661400518, 0.2215833792], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[2.029194274, 0.5661400518, 0.2215833792], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"Cl" =>  [[[601.3456136, 109.5358542, 29.64467686], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[38.96041889, 9.053563477, 2.944499834], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[38.96041889, 9.053563477, 2.944499834], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[38.96041889, 9.053563477, 2.944499834], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[38.96041889, 9.053563477, 2.944499834], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[2.129386495, 0.5940934274, 0.232524141], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[2.129386495, 0.5940934274, 0.232524141], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[2.129386495, 0.5940934274, 0.232524141], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[2.129386495, 0.5940934274, 0.232524141], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]], 
		"Ar" => [[[674.4465184, 122.8512753, 33.24834945], [0.1543289673, 0.5353281423, 0.4446345422], (0, 0, 0)], 
		[[45.16424392, 10.495199, 3.413364448], [-0.09996722919, 0.3995128261, 0.7001154689], (0, 0, 0)], 
		[[45.16424392, 10.495199, 3.413364448], [0.155916275, 0.6076837186, 0.3919573931], (1, 0, 0)], 
		[[45.16424392, 10.495199, 3.413364448], [0.155916275, 0.6076837186, 0.3919573931], (0, 1, 0)], 
		[[45.16424392, 10.495199, 3.413364448], [0.155916275, 0.6076837186, 0.3919573931], (0, 0, 1)], 
		[[2.621366518, 0.731354605, 0.2862472356], [-0.219620369, 0.2255954336, 0.900398426], (0, 0, 0)], 
		[[2.621366518, 0.731354605, 0.2862472356], [0.01058760429, 0.5951670053, 0.462001012], (1, 0, 0)], 
		[[2.621366518, 0.731354605, 0.2862472356], [0.01058760429, 0.5951670053, 0.462001012], (0, 1, 0)], 
		[[2.621366518, 0.731354605, 0.2862472356], [0.01058760429, 0.5951670053, 0.462001012], (0, 0, 1)]])
	attributes = []                         
	for i in 1:(lastindex(basis_set_STO3G[atom]))
		push!(attributes,basis_set_STO3G[atom][i])
    end

	return attributes
end
function Basis_attributes_finder_DZP(atom)
	basis_set_DZP =Dict(
	"H"  =>  [[[19.2406,2.89920,0.65340],[0.032828,0.231208,0.817238],[0,0,0]],
			[[0.17760],[1.000000],[0,0,0]],
			[[1.000000],[1.000000],[1,0,0]],
			[[1.000000],[1.000000],[0,1,0]],
			[[1.000000],[1.000000],[0,0,1]]],
	"O"  =>  [[[7816.54,1175.82,273.188,81.1696,27.1836,3.41360],[0.002031,0.015436,0.073771,0.247606,0.611832,0.241205],[0,0,0]],
			[[9.5322],[1.000000],[0,0,0]],
			[[0.9398],[1.000000],[0,0,0]],
			[[0.2846],[1.000000],[0,0,0]],
			[[35.1832,7.9040,2.3051,0.7171],[0.019580,0.124189,0.394727,0.627375],[1,0,0]],
			[[35.1832,7.9040,2.3051,0.7171],[0.019580,0.124189,0.394727,0.627375],[0,1,0]],
			[[35.1832,7.9040,2.3051,0.7171],[0.019580,0.124189,0.394727,0.627375],[0,0,1]],
			[[0.2137],[1.000000],[1,0,0]],
			[[0.2137],[1.000000],[0,1,0]],
			[[0.2137],[1.000000],[0,0,1]],
			[[0.850000],[1.000000],[0,1,1]],
			[[0.850000],[1.000000],[1,0,1]],
			[[0.850000],[1.000000],[1,1,0]],
			[[0.850000],[1.000000],[2,0,0]],
			[[0.850000],[1.000000],[0,2,0]],
			[[0.850000],[1.000000],[0,0,2]]])

	attributes = []
	for i in 1:lastindex(basis_set_DZP[atom])
		push!(attributes,basis_set_DZP[atom][i])
    end
	return attributes
end




function Basis_attributes_finder_cc_pVDZ(atom)
	basis_set_VDZ = Dict(
	"H"  => [[[1.301000E+01,1.962000E+00,4.446000E-01,1.220000E-01],[1.968500E-02,1.379770E-01,4.781480E-01,5.012400E-01],[0,0,0]],
		  [[1.220000E-01],[1.0],[0,0,0]],
		  [[7.270000E-01],[1.0],[1,0,0]],
		  [[7.270000E-01],[1.0],[0,1,0]],
		  [[7.270000E-01],[1.0],[0,0,1]]],
	"O"  => [[[1.172000E+04,1.759000E+03,4.008000E+02,1.137000E+02,3.703000E+01,1.327000E+01,5.025000E+00,1.013000E+00,3.023000E-01],[7.100000E-04,5.470000E-03,2.783700E-02,1.048000E-01, 2.830620E-01,4.487190E-01,2.709520E-01,1.545800E-02,-2.585000E-03],[0,0,0]],
		  [[1.172000E+04,1.759000E+03,4.008000E+02,1.137000E+02,3.703000E+01,1.327000E+01,5.025000E+00,1.013000E+00,3.023000E-01],[-1.600000E-04,-1.263000E-03,-6.267000E-03,-2.571600E-02,-7.092400E-02,-1.654110E-01,-1.169550E-01,5.573680E-01,5.727590E-01],[0,0,0]],
		  [[3.023000E-01],[1.0],[0,0,0]],
		  [[1.770000E+01, 3.854000E+00,1.046000E+00,2.753000E-01],[4.301800E-02,2.289130E-01,5.087280E-01,4.605310E-01],[1,0,0]],
		  [[1.770000E+01, 3.854000E+00,1.046000E+00,2.753000E-01],[4.301800E-02,2.289130E-01,5.087280E-01,4.605310E-01],[0,1,0]],
		  [[1.770000E+01, 3.854000E+00,1.046000E+00,2.753000E-01],[4.301800E-02,2.289130E-01,5.087280E-01,4.605310E-01],[0,0,1]],
		  [[2.753000E-01],[1.0],[1,0,0]],
		  [[2.753000E-01],[1.0],[0,1,0]],
		  [[2.753000E-01],[1.0],[0,0,1]],
		  [[1.185000E+00],[1.0],[0,1,1]],
		  [[1.185000E+00],[1.0],[1,0,1]],
		  [[1.185000E+00],[1.0],[1,1,0]],
		  [[1.185000E+00],[1.0],[2,0,0]],
		  [[1.185000E+00],[1.0],[0,2,0]]])
	attributes = []                         
	for i in 1:(lastindex(basis_set_VDZ[atom]))
		push!(attributes,basis_set_VDZ[atom][i])
    end
	return attributes
end
using LinearAlgebra
function Basis_attributes_finder_321guc_H(atom)
	basis_set_VDZ = Dict(
	"H"  => [[[5.447178, 0.82454724],[0.1562849786947539, 0.904690876669632],[0,0,0]],
		  [[0.18319158],[1.0],[0,0,0]]])
 	attributes = []                         
	for i in 1:(lastindex(basis_set_VDZ[atom]))
		push!(attributes,basis_set_VDZ[atom][i])
    end
	return attributes
end

function sort_attri(orbital_objects)
    exps_ = []
    coefs_ = []
    origins = []
    shells = []
    norms = []
    for obj in orbital_objects
        E = Float64.(obj.exps)
        C = Float64.(obj.coefs)
        N = Float64.(obj.norms)
        push!(exps_, E)
        push!(coefs_, C)
        push!(origins, (obj.origin))
        push!(shells, (obj.shell))
        push!(norms, N)
    end
    return exps_, coefs_, origins, shells, norms
end


function orbital_config(atoms, geom, basis_set)
    attributes = []
	#println(lastindex(atoms))
	for i in 1:lastindex(atoms)
		if basis_set == "sto3g"
			temp_attri = Basis_attributes_finder(atoms[i])
		elseif basis_set == "dzp"
			temp_attri = Basis_attributes_finder_DZP(atoms[i])
		elseif basis_set == "cc-pvdz"
			temp_attri = Basis_attributes_finder_cc_pVDZ(atoms[i])
		elseif basis_set=="3-21g-uc"
			temp_attri=Basis_attributes_finder_321guc_H(atoms[i])
		end

		for j in 1:lastindex(temp_attri)
			push!(temp_attri[j], geom[i])
		end

		push!(attributes, temp_attri)
	end
	#println(attributes)
	#println(size(attributes[2][1]))
	orbital_objects = [] 
	for i in 1:lastindex(attributes)
		for j in 1:lastindex(attributes[i])
			norms=BasisFunction(attributes[i][j][4], attributes[i][j][3], attributes[i][j][1], attributes[i][j][2],size(attributes[i][j][1])[1],[0.0])
			N,norms=normalize_basis!(norms)
			push!(orbital_objects, BasisFunction(attributes[i][j][4], attributes[i][j][3], attributes[i][j][1], attributes[i][j][2],size(attributes[i][j][1])[1],norms))
		end
	end
	#display(orbital_objects)
	return sort_attri(orbital_objects)
end




