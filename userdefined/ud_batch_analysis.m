function anal_results = ud_batch_analysis
%%
[FileName,PathName,FilterIndex] = uigetfile('*.mat', 'Select an Existing *.mat File');
if ~[PathName FileName], return; end
load([PathName FileName]);
disp([PathName FileName]);
%%
if ~exist('node_info'), node_info = []; end
if ~exist('elem_info'), elem_info = []; end
if ~exist('sect_info'), sect_info = []; end
if ~exist('mat_info'), mat_info = []; end
if ~exist('fixity_info'), fixity_info = []; end
if ~exist('nload_info'), nload_info = []; end
if ~exist('settle_info'), settle_info = []; end
if ~exist('uniload_info'), uniload_info = []; end
if ~exist('analysis_info'), analysis_info = []; end
if ~exist('periods_info'), periods_info = []; end
if ~exist('stiff'), stiff = []; end
if ~exist('ele_pldef'), ele_pldef = []; end
if ~exist('defl'), defl = []; end
if ~exist('react'), react = []; end
if ~exist('ele_for'), ele_for = []; end
if ~exist('ele_yld'), ele_yld = []; end
if ~exist('applfm'), applfm = []; end
h_stat_mes = [];
%%
%%
%Please change these variables to control the analysis solution
%
struct_type = 2; %SpaceFrame(1)|PlanarFrame(2)|SpaceTruss(3)|PlanarTruss(4)
E_t = 3; %E(1)|E_t(2)|E_tm(3)
anatype = 4;  %1st-orderElastic(1)|2nd-orderElastic(2)|1st-orderInelastic(3)|2nd-orderInelastic(4)|Elastic Critical Load(5)
sol_scheme = 2; %Simple Step(1)|Predictor-Corrector(2)|
incr_size = 0.01; %increment size
incr_nums = 10000; %number of increments or number of modes
stop_ratio = 1.00; %maximum applied load ratio
restart = 1; % start new (1) | restart from prev. analysis (2)
apratios =[];
user_def = 0; % MSA2 analysis engine (0) | User Defined analysis (1)
kftoggle = 0;
%%
w_max = 46.2958;
w_dist = w_max*[0:0.1:1];
axial_force = 800;
elem_info(:,4)=2;
for i_count = 1:size(w_dist,2)
	restart = 1;
	stop_ratio = 1.00;
	nload_info(:,:) = 0;
% peform analysis
	nload_info(1:25,2) = -w_dist(i_count);
    [analysis_info,periods_info,stiff,ele_pldef,defl,react,ele_for,ele_yld,...
	    applfm,status_mes] = ud_batch_anaprep(h_stat_mes,struct_type,E_t,...
	    anatype,sol_scheme,incr_size,incr_nums,stop_ratio,restart,apratios,...
	    user_def,node_info,elem_info,sect_info,sect_name,mat_info,mat_name,fixity_info,...
	    nload_info,settle_info,uniload_info,analysis_info,periods_info,...
	    stiff,ele_pldef,defl,react,ele_for,ele_yld,applfm,kftoggle);
	anal_results(i_count,1) = abs(applfm(13,2,end)/w_dist(1,end));
	anal_results(i_count,2) = 0;	
	if analysis_info(1,4) == 0
		restart = 2;
		apratios = analysis_info(5:end)';
		nload_info(1:25,2) = 0;
		nload_info(25,1) = -axial_force;
		stop_ratio = 10.00;
		[analysis_info,periods_info,stiff,ele_pldef,defl,react,ele_for,ele_yld,...
			applfm,status_mes] = ud_batch_anaprep(h_stat_mes,struct_type,E_t,...
			anatype,sol_scheme,incr_size,incr_nums,stop_ratio,restart,apratios,...
			user_def,node_info,elem_info,sect_info,sect_name,mat_info,mat_name,fixity_info,...
			nload_info,settle_info,uniload_info,analysis_info,periods_info,...
			stiff,ele_pldef,defl,react,ele_for,ele_yld,applfm,kftoggle);
		anal_results(i_count,2) = abs(ele_for(13,1,end));
	end
	anal_results(i_count,3) = size(analysis_info(5:end),2)-4;	
	anal_results
end
%
disp('Batch analysis completed successfully');