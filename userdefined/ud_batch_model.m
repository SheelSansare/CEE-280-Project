function ud_batch_model
%
%  This MASTAN2 utility allows the user to create all
%  or a portion of a MASTAN2 model by numerical input
%  (instead of graphical input).  The user should start by
%  making a copy of the ud_batch_model.m file.  The copied file
%  can then be edited using Matlab's edit command and
%  eventually executed by typing the copied filename
%  at the Matlab command prompt.  For example,
%  >>!copy ud_batch_model.m tower.m
%  >>edit tower.m
%  >>tower.m

%
%  Create space for all arrays (Do not modify).
%
node_info = []; elem_info = [];
sect_info = []; sect_name = [];
mat_info = []; mat_name = [];
nload_info = []; uniload_info = [];
fixity_info = []; settle_info = [];
%
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%  Start of user defined information (Should be modified).
%
% model_title = [text string describing section]
model_title = 'Example frame prepared by RDZ';
%
% Node information
%  Requires three arrays to be constructed, including
%  node_info, support_info, nload_info
%
node_info = [...
			0	0	0;	...	%Node 1's XYZ Coordinates
			0	120	0;	...	%Node 2's XYZ Coordinates
			240	120	0;	...	%Node 3's XYZ Coordinates
			240	0	0;	...	%Node 4's XYZ Coordinates
			];
%
% support_info = [Support condition of nodes' six d.o.f. (free=NaN, fixed=0, val for prescribed displacement)]
%
support_info =	[...
				0 0 0 NaN NaN NaN;...			Node 1's support condition
				NaN NaN NaN NaN NaN NaN;...		Node 2's support condition
				NaN NaN NaN NaN NaN NaN;...		Node 3's support condition
				0 0 0 0 0 0;...			Node 4's support condition
				];
%
% nload_info = [concentrated force or moment on nodes' six d.o.f.]
%
nload_info =	[...
				0 0 0 0 0 0;...			Node 1's forces/moments
				0 -100 0 0 0 0;...		Node 2's forces/moments
				0 -100 0 0 0 0;...		Node 3's forces/moments
				0 0 0 0 0 0;...			Node 4's forces/moments
				];
%
% Element information
%  Requires six arrays to be constructed, including
%  elem_info, uniload_info, thermal_info, sect_info, sect_name, mat_info, mat_name
%
%
% elem_info = [Node i, Node j, Section #, Material #, Beta Angle (rads), ...
%				Flexure condition Node i (rigid=0,pin=1,spring=2), Flexure condition Node j (rigid=0,pin=1,spring=2), ...
%				Warping condition Node i (fixed=0,free=1,cont=2), Warping condition Node j (fixed=0,free=1,cont=2), ...
%				Major-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
%				Minor-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
%				Major-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
%				Minor-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
%				Major-axis spring moment capacity Mpz node i (val = value to inf (unlimited)),...
%				Minor-axis spring moment capacity Mpz node i (val = value to inf (unlimited)),...
%				Major-axis spring moment capacity Mpz node j (val = value to inf (unlimited)),...
%				Minor-axis spring moment capacity Mpz node j (val = value to inf (unlimited))]
elem_info = [...
			1 2 1 2 0 0 0 0 0 inf inf inf inf inf inf inf inf;...   %Element 1's information
			2 3 2 1 0 0 0 0 0 inf inf inf inf inf inf inf inf;...   %Element 2's information
			4 3 1 2 0 0 0 0 0 inf inf inf inf inf inf inf inf;...   %Element 3's information
			];
		
%
% uniload_info = [uniformly distributed loads along elements' local axes (x,y,z)]
%
uniload_info =	[...
				0 0 0 ;...			Element 1's distributed load
				0 -.25 0 ;...		Element 2's distributed load
				0 0 0;...			Element 3's distributed load
				];
%
% thermal_info = [thermal effects along elements' local axes, including
%					thermal coef, dT(centroid), Tgradient(y'), Tgradient(z')]
%
thermal_info =	[...
				0 0 0 0;...			Element 1's temp. effects
				0 0 0 0;...			Element 2's temp. effects
				0 0 0 0;...			Element 3's temp. effects
				];
%
% Section information
%  Requires two arrays to be constructed, including
%  sect_info and sect_name
%
% sect_info = [Area, Moment of inertia Izz, Moment of inertia Iyy, Torsion constant J,
%				Warping coefficient Cw, Plastic section modulus Zzz, Plastic section modulus Zyy, ...
%				Shear area Ayy, Shear area Azz YldSurf(1) YldSurf(2) YldSurf(3)]
sect_info = [...
			4.44 48 3.41 0.137 51.8 13.6 2.67 Inf Inf 1 1 1; ...	%Section 1's properties
			22.8 658 519 921 0 111 100 14 12 1 1 1;...				%Section 2's properties
			];
%
% sect_name = [brief text string describing section]
%
sect_name = {...
			'W8X15';...				Section 1's name
			'HSS14X12X1/2';...		Section 2's name
			};
%
% Material information
%  Requires two arrays to be constructed, including
%  mat_info and mat_name
%
% mat_info = [Modulus of elasticity E, Poisson Ratio v, Yield strength Fy, Weight density Wt. Dens.]
%
mat_info = [...
			29000 0.3 50 0;...	%Material 1's properties
			10500 0.28 35 0;...				%Material 2's properties
			];
%
% mat_name = [Brief text string describing material]
%
mat_name = {...
			'Steel';...				Material 1's name
			'Metal';...				Material 2's name
			};
%
% End of user defined information
%
%%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%
%
%  Clean-up/pad arrays and save function (Do not modify).
%
webdir = el_webdir(size(elem_info,1),elem_info(:,1:2),elem_info(:,5),node_info,node_info,zeros(size(node_info,1),6));
%
elem_info = [	elem_info(:,1:5) ones(size(elem_info,1),2) webdir ones(size(elem_info,1),1) ...
				elem_info(:,6:7) ones(size(elem_info,1),2) ...
 				elem_info(:,8:9) ones(size(elem_info,1),2) ...
 				elem_info(:,10:17)];
nload_info = [nload_info NaN(size(nload_info,1),6)];
fixity_info = support_info;
fixity_info(find(~isnan(support_info))) = 0;
fixity_info = [fixity_info NaN(size(fixity_info,1),6)];
settle_info = support_info;
settle_info(find(~isnan(settle_info)& settle_info==0)) = NaN;
settle_info = [settle_info NaN(size(settle_info,1),6)];
uniload_info = [uniload_info zeros(size(uniload_info,1),3) NaN(size(uniload_info,1),6) thermal_info];

analysis_info=[]; periods_info=[]; first_el_settings=[]; sec_el_settings=[];
first_inel_settings=[]; sec_inel_settings=[]; ebuck_settings=[]; ibuck_settings=[];
period_settings=[]; deflect_settings=[]; axial_settings=[]; sheary_settings=[];
shearz_settings=[]; momx_settings=[]; momy_settings=[]; momz_settings=[];
bimom_settings = [];
%
filt = '*.mat';
[newmatfile, newpath, filterindex] = uiputfile(filt, 'Save Model As');
if filterindex
	file_str = [newpath newmatfile];
	save(file_str,'model_title','node_info','elem_info', 'sect_info', 'sect_name', ...
		'mat_info','mat_name','fixity_info','nload_info','uniload_info','settle_info', ...
		'analysis_info','periods_info','first_el_settings','sec_el_settings','first_inel_settings', ...
		'sec_inel_settings','ebuck_settings','ibuck_settings','period_settings', ...
		'deflect_settings','axial_settings','sheary_settings', ...
		'shearz_settings','momx_settings','momy_settings','momz_settings','bimom_settings','-v6');
	disp('Batch model construction completed successfully');
else
	disp('Batch model construction cancelled');
end