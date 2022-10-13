function [analysis_info,periods_info,stiff,ele_pldef,defl,react,ele_for,ele_yld,...
    applfm,status_mes] = ud_batch_anaprep(h_stat_mes,struct_type,E_t,...
    anatype,sol_scheme,incr_size,incr_nums,stop_ratio,restart,apratios,...
    user_def,node_info,elem_info,sect_info,sect_name,mat_info,mat_name,fixity_info,...
    nload_info,settle_info,uniload_info,analysis_info,periods_info,...
    stiff,ele_pldef,defl,react,ele_for,ele_yld,applfm,kftoggle)
%Function of MASTAN2
%
%  Authors :   R.D. Ziemian
%                 Bucknell University
%              W. McGuire
%                 Cornell University
%
%  Copyright (c)   
%   Permission is not granted to modify and re-distribute this 
%   code in any manner.  All standard disclaimers apply. 

h_fig = figure('visible','off');
uicontrol(h_fig,'style','pushbutton','tag','stop_analysis',...
    'UserData',0,'visible','off')
    
node_total = size(node_info,1);
elem_total = size(elem_info,1);

for i_ele = 1:elem_total
    elem_sect_name{i_ele} = sect_name{elem_info(i_ele,3)};
end

for i_ele = 1:elem_total
    elem_mat_name{i_ele} = mat_name{elem_info(i_ele,4)};
end


if struct_type < 3
	truss = 0;
else
	truss = 1;
end

flag = 0;

nele = elem_total;
nnodes = node_total;
coord = node_info(1:nnodes,1:3);
ends = elem_info(1:nele,[1 2 12 13 16 17 20 21 22 23 24 25 26 27]);
beta_ang = elem_info(1:nele,5);

for i = 1:nnodes
	for j = 1:6
		if isnan(fixity_info(i,j)) & isnan(settle_info(i,j));
			fixity(i,j) = NaN;
		elseif ~isnan(fixity_info(i,j)) & isnan(settle_info(i,j))
			fixity(i,j) = fixity_info(i,j);
		else
			fixity(i,j) = settle_info(i,j);
		end
		if truss & j>3 & isnan(fixity(i,j))
			fixity(i,j) = 0;
		end
	end
end

A = sect_info(elem_info(1:nele,3),1);
Zzz = sect_info(elem_info(1:nele,3),6);
Zyy = sect_info(elem_info(1:nele,3),7);
if ~truss
	for ele = 1:nele
		if elem_info(1:nele,12)== 1 & elem_info(1:nele,13)== 1
			Izz(ele) = 0;
			Iyy(ele) = 0;
			Ayy(ele,1) = inf;
			Azz(ele,1) = inf;
		else
			Izz(ele) = sect_info(elem_info(ele,3),2);
			Iyy(ele) = sect_info(elem_info(ele,3),3);
			Ayy(ele) = sect_info(elem_info(ele,3),8);
			Azz(ele) = sect_info(elem_info(ele,3),9);
		end
	end
	J   = sect_info(elem_info(1:nele,3),4);
	Cw  = sect_info(elem_info(1:nele,3),5);
else
	Izz = zeros(nele,1);
	Iyy = zeros(nele,1);
	J = zeros(nele,1);
	Cw = zeros(nele,1);
	Ayy(1:nele,1) = inf;
	Azz(1:nele,1) = inf;
end
YldSurf = sect_info(elem_info(1:nele,3),10:12);

E = mat_info(elem_info(1:nele,4),1);
v = mat_info(elem_info(1:nele,4),2);
Fy = mat_info(elem_info(1:nele,4),3);
Wt = mat_info(elem_info(1:nele,4),4);

webdir = elem_info(1:nele,8:10);

concen = nload_info(1:nnodes,1:6);

w = uniload_info(1:nele,1:3);
thermal = uniload_info(1:nele,13:16);
for ele = 1:nele
	  del(ele,:) = coord(ends(ele,2),:)-coord(ends(ele,1),:);
	  L(ele) = norm(del(ele,:),2);
end

%  Set-up for 2D analysis
if struct_type == 1 | struct_type == 3
	two_dim = 0;
else
	two_dim = 1;
end

if two_dim ==1
	fixity(1:nnodes,3:5) = 0;
	if anatype ~= 8 && anatype ~= 9, ends(1:nele,5:6) = 1; end

end

%
if anatype == 1
	applfm =[];
	if user_def == 0
		[stiff,defl,react,ele_for,aflag] = msa_1el(nnodes,coord,concen,fixity,nele,ends,...
					A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype);
    else
		stiff =[];
		if two_dim ~= 1
			[defl,react,ele_for,aflag] = ud_3d1el(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype);
			if aflag == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
			end
		else
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			[defl,react,ele_for,aflag] = ud_2d1el(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,anatype);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3)  defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end

				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 6])
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				else
					ele_for = zeros(nele,16);
				end
			end
        end
	end
	ele_yld = zeros(nele,2);
	ele_pldef = zeros(nele,14);
	limit_state = 0;
	if aflag == 1
	   apratios = [1.00];
	   applfm(:,:,1) = concen;
	   status_mes = sprintf('Applied Load Ratio = %5.3f  ---->  Success:  Analysis Complete',apratios(size(apratios,1)));
	elseif aflag ==0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	end

elseif anatype == 2
	if size(analysis_info,2) > 4
		if two_dim ~= analysis_info(1) | truss ~= analysis_info(2) | anatype ~= analysis_info(3)
			apratios = [];
			limit_state = 0;
		end
		if size(apratios,1) > 0
			limit_state = analysis_info(4);
		else
			limit_state = 0;
		end
	else
		limit_state = 0;
	end
%
	if restart == 1 | isempty(apratios)
		total_pre = 0;
		applfm =[];
	else
		total_pre = size(apratios,1);
	end
%
	if user_def == 0
		if sol_scheme > 22 & two_dim ~= 1 & ~truss
			status_mes = sprintf('Error:  Iterative solution not functional for space frame analysis');
			return
		end
		[stiff,defl,react,ele_for,aflag,apratios,limit_state] = msa_2el(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,stiff,defl,react,ele_for,apratios,limit_state,h_stat_mes,kftoggle);
	else
		stiff =[];
		if two_dim ~= 1
			[defl,react,ele_for,aflag,apratios,limit_state] = ud_3d2el(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,apratios,limit_state,h_stat_mes);
%
			if aflag == 1 | aflag == -1 | aflag == 3 | aflag == 99
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 6 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 6 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif all(size(ele_for) == [nele 12 size(apratios,1)])
					ele_for = [ele_for zeros(nele,4,size(apratios,1))];
				elseif all(size(ele_for) == [nele 14 size(apratios,1)])
					ele_for = [ele_for zeros(nele,2,size(apratios,1))];
				elseif all(size(ele_for) == [nele 16 size(apratios,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			if total_pre > 0
				defl =     defl(:,[1 2 6],:);
				react =   react(:,[1 2 6],:);
				ele_for = ele_for(:,[1 2 6 7 8 12],:);
			end
%
			[defl,react,ele_for,aflag,apratios,limit_state] = ud_2d2el(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,apratios,limit_state,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1 | aflag == -1 | aflag == 3 | aflag == 99
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 3 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				else
					defl_tmp = zeros(nnodes,6,size(apratios,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 3 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				else
					react_tmp = zeros(nnodes,6,size(apratios,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif any(size(ele_for) ~= [nele 6 size(apratios,1)])
					ele_for = zeros(nele,16,size(apratios,1));
				else
					ele_for_tmp = zeros(nele,16,size(apratios,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[];
			 end
%
			end
%
		end
	end
	ele_yld = zeros(nele,2,size(ele_for,3));
	ele_pldef = zeros(nele,14,size(ele_for,3));
%
	total = size(apratios,1);
	if total == 0
		load_val = 0;
	else
		load_val = apratios(total);
	end
%
	if total > 0
		for iii=total_pre+1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm = [];
	end
%
	if aflag == 1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Success:  Analysis Complete',total,load_val);
	elseif aflag == 3
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped by User',total,load_val);
	elseif aflag == 0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif aflag == -1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Halted: Limit Reached',total,load_val);
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	elseif aflag == 99
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Halted: No Convergence',total,load_val);
	end

elseif anatype == 3
	if size(analysis_info,2) > 4
		if two_dim ~= analysis_info(1) | truss ~= analysis_info(2) | anatype ~= analysis_info(3)
			apratios = [];
		end
		if size(apratios,1) > 0
			limit_state = analysis_info(4);
		else
			limit_state = 0;
		end
	else
		limit_state = 0;
	end
%
	if restart == 1 | isempty(apratios)
		total_pre = 0;
		applfm =[];
	else
		total_pre = size(apratios,1);
	end
%
	if user_def == 0
		[stiff,defl,react,ele_for,ele_yld,ele_pldef,aflag,apratios,limit_state] = msa_1in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,stiff,defl,react,ele_for,ele_yld,ele_pldef,apratios,limit_state,h_stat_mes,kftoggle);
	else
		stiff =[];
		if two_dim ~= 1
			[defl,react,ele_for,ele_yld,aflag,apratios,limit_state] = ud_3d1in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,ele_yld,apratios,limit_state,h_stat_mes);
%
			if aflag == 1 | aflag == -1 | aflag == 2 | aflag == 3
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 6 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 6 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif all(size(ele_for) == [nele 12 size(apratios,1)])
					ele_for = [ele_for zeros(nele,4,size(apratios,1))];
				elseif all(size(ele_for) == [nele 14 size(apratios,1)])
					ele_for = [ele_for zeros(nele,2,size(apratios,1))];
				elseif all(size(ele_for) == [nele 16 size(apratios,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(apratios,1));
				end
				if ndims(ele_yld) ~= 3
					ele_yld = zeros(nele,2,size(apratios,1));
				elseif any(size(ele_yld) ~= [nele 2 size(apratios,1)])
					ele_yld = zeros(nele,2,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
				if ndims(ele_yld) ~= 2
					ele_yld = zeros(nele,2);
				elseif any(size(ele_yld) ~= [nele 2])
					ele_yld = zeros(nele,2);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			if total_pre > 0
				defl =     defl(:,[1 2 6],:);
				react =   react(:,[1 2 6],:);
				ele_for = ele_for(:,[1 2 6 7 8 12],:);
			end
%
			[defl,react,ele_for,ele_yld,aflag,apratios,limit_state] = ud_2d1in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,ele_yld,apratios,limit_state,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1 | aflag == -1 | aflag == 2 | aflag == 3
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 3 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				else
					defl_tmp = zeros(nnodes,6,size(apratios,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 3 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				else
					react_tmp = zeros(nnodes,6,size(apratios,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif any(size(ele_for) ~= [nele 6 size(apratios,1)])
					ele_for = zeros(nele,16,size(apratios,1));
				else
					ele_for_tmp = zeros(nele,16,size(apratios,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
				if ndims(ele_yld) ~= 3
					ele_yld = zeros(nele,2,size(apratios,1));
				elseif any(size(ele_yld) ~= [nele 2 size(apratios,1)])
					ele_yld = zeros(nele,2,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
				if ndims(ele_yld) ~= 2
					ele_yld = zeros(nele,2);
				elseif any(size(ele_yld) ~= [nele 2])
					ele_yld = zeros(nele,2);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
		end
	end
%
	total = size(apratios,1);
	if total == 0
		load_val = 0;
	else
		load_val = apratios(total);
	end
%
	if total > 0
		for iii=total_pre+1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm =[];
	end
%
	if aflag == 1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Success:  Analysis Complete',total,load_val);
	elseif aflag == 2
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f -> Analysis Stopped, Extreme Deflections',total,load_val);
	elseif aflag == 3
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped by User',total,load_val);
	elseif aflag == 4
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped, Unload/Reload',total,load_val);
	elseif aflag == 0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif aflag == -1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Halted: Limit Reached',total,load_val);
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	end

elseif anatype == 4
	if size(analysis_info,2) > 4
		if two_dim ~= analysis_info(1) | truss ~= analysis_info(2) | anatype ~= analysis_info(3)
			apratios = [];
		end
		if size(apratios,1) > 0
			limit_state = analysis_info(4);
		else
			limit_state = 0;
		end
	else
		limit_state = 0;
	end
%
	if restart == 1 | isempty(apratios)
		total_pre = 0;
		applfm =[];
	else
		total_pre = size(apratios,1);
	end
%
	if user_def == 0
		[stiff,defl,react,ele_for,ele_yld,ele_pldef,aflag,apratios,limit_state] = msa_2in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,stiff,defl,react,ele_for,ele_yld,ele_pldef,apratios,limit_state,h_stat_mes,kftoggle);
	else
		stiff =[];
		if two_dim ~= 1
			[defl,react,ele_for,ele_yld,aflag,apratios,limit_state] = ud_3d2in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,ele_yld,apratios,limit_state,h_stat_mes);
%
			if aflag == 1 | aflag == -1 | aflag == 2 | aflag == 3 | aflag == 4 | aflag == 99
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 6 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 6 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif all(size(ele_for) == [nele 12 size(apratios,1)])
					ele_for = [ele_for zeros(nele,4,size(apratios,1))];
				elseif all(size(ele_for) == [nele 14 size(apratios,1)])
					ele_for = [ele_for zeros(nele,2,size(apratios,1))];
				elseif all(size(ele_for) == [nele 16 size(apratios,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(apratios,1));
				end
				if ndims(ele_yld) ~= 3
					ele_yld = zeros(nele,2,size(apratios,1));
				elseif any(size(ele_yld) ~= [nele 2 size(apratios,1)])
					ele_yld = zeros(nele,2,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
				if ndims(ele_yld) ~= 2
					ele_yld = zeros(nele,2);
				elseif any(size(ele_yld) ~= [nele 2])
					ele_yld = zeros(nele,2);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			if total_pre > 0
				defl =     defl(:,[1 2 6],:);
				react =   react(:,[1 2 6],:);
				ele_for = ele_for(:,[1 2 6 7 8 12],:);
			end
%
			[defl,react,ele_for,ele_yld,aflag,apratios,limit_state] = ud_2d2in(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,E_t,anatype,sol_scheme,incr_nums,incr_size,...
				stop_ratio,restart,defl,react,ele_for,ele_yld,apratios,limit_state,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1 | aflag == -1 | aflag == 2 | aflag == 3 | aflag == 4 | aflag == 99
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 3 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				else
					defl_tmp = zeros(nnodes,6,size(apratios,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 3 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				else
					react_tmp = zeros(nnodes,6,size(apratios,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif any(size(ele_for) ~= [nele 6 size(apratios,1)])
					ele_for = zeros(nele,16,size(apratios,1));
				else
					ele_for_tmp = zeros(nele,16,size(apratios,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
				if ndims(ele_yld) ~= 3
					ele_yld = zeros(nele,2,size(apratios,1));
				elseif any(size(ele_yld) ~= [nele 2 size(apratios,1)])
					ele_yld = zeros(nele,2,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
				if ndims(ele_yld) ~= 2
					ele_yld = zeros(nele,2);
				elseif any(size(ele_yld) ~= [nele 2])
					ele_yld = zeros(nele,2);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		end
	end
%
	total = size(apratios,1);
	if total == 0
		load_val = 0;
	else
		load_val = apratios(total);
	end
%
	if total > 0
		for iii=total_pre+1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm = [];
	end
%
	if aflag == 1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Success:  Analysis Complete',total,load_val);
	elseif aflag == 2
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f -> Analysis Stopped, Extreme Deflections',total,load_val);
	elseif aflag == 3
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped by User',total,load_val);
	elseif aflag == 4
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped, Unload/Reload',total,load_val);
	elseif aflag == 0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif aflag == -1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Halted: Limit Reached',total,load_val);
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	end

elseif anatype == 5
	applfm =[];
	if user_def == 0
		[defl,react,ele_for,aflag,apratios] = msa_ecl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,incr_nums,h_stat_mes);
	else
		if two_dim ~= 1
		[defl,react,ele_for,aflag,apratios] = ud_3decl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,incr_nums,h_stat_mes);
%
			if aflag == 1
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 6 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 6 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif all(size(ele_for) == [nele 12 size(apratios,1)])
					ele_for = [ele_for zeros(nele,4,size(apratios,1))];
				elseif all(size(ele_for) == [nele 14 size(apratios,1)])
					ele_for = [ele_for zeros(nele,2,size(apratios,1))];
				elseif all(size(ele_for) == [nele 16 size(apratios,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			[defl,react,ele_for,aflag,apratios] = ud_2decl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,incr_nums,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 3 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				else
					defl_tmp = zeros(nnodes,6,size(apratios,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 3 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				else
					react_tmp = zeros(nnodes,6,size(apratios,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif any(size(ele_for) ~= [nele 6 size(apratios,1)])
					ele_for = zeros(nele,16,size(apratios,1));
				else
					ele_for_tmp = zeros(nele,16,size(apratios,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[];
			 end
%
			end
%
		end
	end
%
	ele_yld = zeros(nele,2,size(ele_for,3));
	ele_pldef = zeros(nele,14,size(ele_for,3));
	limit_state = 0;
	stiff =[];
%
	total = size(apratios,1);
	if total > 0
		for iii=1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm = [];
	end
%
	if aflag == 1
	   if total > 0
	     status_mes = sprintf('# of Modes Calculated = %i  ---->  Success:  Analysis Complete',size(apratios,1));
	   else
	     status_mes = sprintf('Warning:  Analysis complete with no positive buckling load ratios detected');
	   end
	elseif aflag ==0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	end

elseif anatype == 6
	applfm =[];
	if user_def == 0
		[defl,react,ele_for,aflag,apratios] = msa_icl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,incr_nums,h_stat_mes);
	else
		if two_dim ~= 1
		[defl,react,ele_for,aflag,apratios] = ud_3dicl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,incr_nums,h_stat_mes);
%
			if aflag == 1
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 6 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 6 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif all(size(ele_for) == [nele 12 size(apratios,1)])
					ele_for = [ele_for zeros(nele,4,size(apratios,1))];
				elseif all(size(ele_for) == [nele 14 size(apratios,1)])
					ele_for = [ele_for zeros(nele,2,size(apratios,1))];
				elseif all(size(ele_for) == [nele 16 size(apratios,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(apratios,1));
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			[defl,react,ele_for,aflag,apratios] = ud_2dicl(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,incr_nums,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1
%
			 if size(apratios,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(apratios,1));
				elseif any(size(defl) ~= [nnodes 3 size(apratios,1)])
					defl = zeros(nnodes,6,size(apratios,1));
				else
					defl_tmp = zeros(nnodes,6,size(apratios,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(apratios,1));
				elseif any(size(react) ~= [nnodes 3 size(apratios,1)])
					react = zeros(nnodes,6,size(apratios,1));
				else
					react_tmp = zeros(nnodes,6,size(apratios,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(apratios,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(apratios,1));
				elseif any(size(ele_for) ~= [nele 6 size(apratios,1)])
					ele_for = zeros(nele,16,size(apratios,1));
				else
					ele_for_tmp = zeros(nele,16,size(apratios,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
			 elseif size(apratios,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
			 else
				apratios=[]; defl=[]; react=[]; ele_for=[];
			 end
%
			end
%
		end
	end
%
	ele_yld = zeros(nele,2,size(ele_for,3));
	ele_pldef = zeros(nele,14,size(ele_for,3));
	limit_state = 0;
	stiff =[];
%
	total = size(apratios,1);
	if total > 0
		for iii=1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm = [];
	end
%
	if aflag == 1
	   if total > 0
	     status_mes = sprintf('# of Modes Calculated = %i  ---->  Success:  Analysis Complete',size(apratios,1));
	   else
	     status_mes = sprintf('Warning:  Analysis complete with no positive buckling load ratios detected');
	   end
	elseif aflag ==0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	elseif aflag == 3
	   apratios = [];
	   status_mes = sprintf('Analysis Stopped by User ---> No Results Provided');
	elseif aflag == 2
	   apratios = [];
	   status_mes = sprintf('Error:  Elastic Buckling Controls');
	elseif aflag == -1
	   apratios = [];
	   status_mes = sprintf('Error:  Solution did not converge during force distribution analysis');
	elseif aflag == -2
	   apratios = [];
	   status_mes = sprintf('Error:  No convergence during iterative eigenvalue analysis');
	elseif aflag == -99
	   apratios = [];
	   status_mes = sprintf('Error:  Inelastic Critical Load Analysis does not work with thermal loads');
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	end

elseif anatype == 7
	applfm =[];
	if user_def == 0
		if sol_scheme == 1
			anatype_prev = 1;
		else
			if analysis_info(3) ~= 7
				anatype_prev = analysis_info(3);
			else
				anatype_prev = periods_info(2);
			end
		end
		[defl,react,ele_for,aflag,periods] = msa_natp(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,...
				incr_nums,stop_ratio,incr_size,sol_scheme,restart,analysis_info,stiff,h_stat_mes,mass_type,addMDir);
	else
		anatype_prev = 1;
		if two_dim ~= 1
		[defl,react,ele_for,aflag,periods] = ud_3dnatp(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,...
				incr_nums,h_stat_mes);
%
			if aflag == 1
%
			 if size(periods,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(periods,1));
				elseif any(size(defl) ~= [nnodes 6 size(periods,1)])
					defl = zeros(nnodes,6,size(periods,1));
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(periods,1));
				elseif any(size(react) ~= [nnodes 6 size(periods,1)])
					react = zeros(nnodes,6,size(periods,1));
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(periods,1));
				elseif all(size(ele_for) == [nele 12 size(periods,1)])
					ele_for = [ele_for zeros(nele,4,size(periods,1))];
				elseif all(size(ele_for) == [nele 14 size(periods,1)])
					ele_for = [ele_for zeros(nele,2,size(periods,1))];
				elseif all(size(ele_for) == [nele 16 size(periods,1)])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16,size(periods,1));
				end
			 elseif size(periods,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 6])
					defl = zeros(nnodes,6);
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 6])
					react = zeros(nnodes,6);
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif all(size(ele_for) == [nele 12])
					ele_for = [ele_for zeros(nele,4)];
				elseif all(size(ele_for) == [nele 14])
					ele_for = [ele_for zeros(nele,2)];
				elseif all(size(ele_for) == [nele 16])
					ele_for = ele_for;
				else
					ele_for = zeros(nele,16);
				end
			 else
				periods=[]; defl=[]; react=[]; ele_for=[]; ele_yld=[];
			 end
%
			end
%
		else
%
			coord = coord(:,[1 2]);
			concen = concen(:,[1 2 6]);
			fixity = fixity(:,[1 2 6]);
			w = w(:,[1 2]);
			ends = ends(:,[1:4 7 9 11 13]);
%
			[defl,react,ele_for,aflag,periods] = ud_2dnatp(nnodes,coord,concen,fixity,nele,ends,...
				A,Izz,Zzz,Ayy,E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,incr_nums,h_stat_mes);
%
			concen = [ concen(:,[1 2]) zeros(nnodes,3) concen(:,3)];
%
			if aflag == 1
%
			 if size(periods,1) > 1
				if ndims(defl) ~= 3
					defl = zeros(nnodes,6,size(periods,1));
				elseif any(size(defl) ~= [nnodes 3 size(periods,1)])
					defl = zeros(nnodes,6,size(periods,1));
				else
					defl_tmp = zeros(nnodes,6,size(periods,1));
					defl_tmp(:,[1 2 6],:) = defl;
					defl = defl_tmp;
%					defl = [defl(:,[1 2],:) zeros(nnodes,3,size(periods,1)) defl(:,3,:)];
				end
				if ndims(react) ~= 3
					react = zeros(nnodes,6,size(periods,1));
				elseif any(size(react) ~= [nnodes 3 size(periods,1)])
					react = zeros(nnodes,6,size(periods,1));
				else
					react_tmp = zeros(nnodes,6,size(periods,1));
					react_tmp(:,[1 2 6],:) = react;
					react = react_tmp;
%					react = [react(:,[1 2],:) zeros(nnodes,3,size(periods,1)) react(:,3,:)];
				end
				if ndims(ele_for) ~= 3
					ele_for = zeros(nele,16,size(periods,1));
				elseif any(size(ele_for) ~= [nele 6 size(periods,1)])
					ele_for = zeros(nele,16,size(periods,1));
				else
					ele_for_tmp = zeros(nele,16,size(periods,1));
					ele_for_tmp(:,[1 2 6 7 8 12],:) = ele_for;
					ele_for = ele_for_tmp;
				end
			 elseif size(periods,1) == 1
				if ndims(defl) ~= 2
					defl = zeros(nnodes,6);
				elseif any(size(defl) ~= [nnodes 3])
					defl = zeros(nnodes,6);
				else
					defl = [defl(:,[1 2]) zeros(nnodes,3) defl(:,3)];
				end
				if ndims(react) ~= 2
					react = zeros(nnodes,6);
				elseif any(size(react) ~= [nnodes 3])
					react = zeros(nnodes,6);
				else
					react = [react(:,[1 2]) zeros(nnodes,3) react(:,3)];
				end
				if ndims(ele_for) ~= 2
					ele_for = zeros(nele,16);
				elseif any(size(ele_for) ~= [nele 6])
					ele_for = zeros(nele,16);
				else
					ele_for = [ele_for(:,[1 2]) zeros(nele,3) ele_for(:,3) ...
						        ele_for(:,[4 5]) zeros(nele,3) ele_for(:,6) zeros(nele,4)];
				end
			 else
				periods=[]; defl=[]; react=[]; ele_for=[];ele_yld=[];
			 end
%
			end
%
		end
	end
%
	if sol_scheme == 1
		limit_state = 0;
	else
		if incr_size <= size(stiff,3)
			limit_state = 0;
		else
			limit_state= -1;
		end
	end
%
	if ~isempty(apratios)
		total = size(apratios,1);
	else
		total = 0;
	end
	if total > 0
		for iii=1:total
			if iii > 1
				applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
			else
				applfm(:,:,1) = apratios(iii)*concen;
			end
		end
	else
		applfm = [];
	end
%
	total = size(periods,2)-sum(isnan(periods));
	if aflag == 1
	   if total > 0
	     status_mes = sprintf('# of Modes Calculated = %i  ---->  Success:  Analysis Complete',total);
	   else
		defl=[]; react=[]; ele_for=[]; periods=[];
	     status_mes = sprintf('Warning:  Analysis complete with no vibrational modes calculated.');
	   end
	elseif aflag ==0
		periods = [];
		status_mes = sprintf('Error:  Unstable Structure');
	elseif isinf(aflag)
	   status_mes = sprintf('Error:  No Analysis Code Available');
    end
    
elseif anatype == 8 | anatype == 9

    limit_state = 0;
    applfm = [];
%  Confirm system is stable
    ends_prev = ends;
    ends(1:nele,5:6) = 1;
    [stiff,defl,react,ele_for,aflag] = msa_1el(nnodes,coord,concen,fixity,nele,ends,...
					A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype);
    ends = ends_prev;
    if aflag == 1 && sol_scheme == 2
%           Confirm all sections employed are in FE++2015 database (aisc.txt)
            fid = fopen('aisc.txt');
            data_base = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
                 'Headerlines',1);
            fclose(fid);
            for i_sect = 1:size(sect_name,2)
                if ~isempty(find(elem_info(:,3)==i_sect, 1))
                    i_s = find(strcmp(data_base{1},sect_name{i_sect}), 1);
                    if isempty(i_s) && aflag == 1, aflag = -4; end
                end
            end
    end
%  Confirm model is planar if 2D analysis selected
    if aflag == 1
        if two_dim == 1
            if all(coord(1,3) == coord(:,3))
                aflag = 1;
            else
                aflag = 2;
            end
        end
    end
%  Perform analysis
    if aflag == 1
        status_mes = sprintf('FE++2015 analyis is being performed...');
        set(h_stat_mes,'String',status_mes); drawnow;
%        if sol_scheme == 1, mat_info(:,3) = 1e6*mat_info(:,3); end  %Artificially increase Fy for second-order elastic analysis
        %
        [defl,react,ele_for,ele_yld,aflag,apratios,limit_state] = msa_fepp(nnodes,coord,concen,fixity,nele,ends,elem_info,mat_info,...
            A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,sect_name,elem_sect_name,mat_name,elem_mat_name,...
            truss,E_t,two_dim,sol_scheme,incr_nums,incr_size,...
            stop_ratio,restart,defl,react,ele_for,ele_yld,apratios,limit_state,h_stat_mes);
        %
        ele_yld = zeros(nele,2,size(ele_for,3));
        ele_pldef = zeros(nele,14,size(ele_for,3));
        limit_state = 0;
        stiff =[];
        %
        total = size(apratios,1);
        if total == 0
            load_val = 0;
        else
            load_val = apratios(total);
        end
        %
        if total > 0
            for iii=1:total
                if iii > 1
                    applfm(:,:,iii) = applfm(:,:,iii-1) + (apratios(iii)-apratios(iii-1))*concen;
                else
                    applfm(:,:,1) = apratios(iii)*concen;
                end
            end
        else
            applfm = [];
        end
    end
%
	if aflag == 1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Success:  Analysis Complete',total,load_val);
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == 2
	   status_mes = sprintf('Analyis not performed:  Model is not a planar frame');
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == -2
	   status_mes = sprintf('Analyis aborted by User.');
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == 3
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped by User',total,load_val);
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == -3
	   status_mes = sprintf('Analyis failed to run.');
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == 4
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Stopped, Unload/Reload',total,load_val);
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == -4
	   status_mes = sprintf('Analyis not performed:  Section(s) not included in database');
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == 0
	   apratios = [];
	   status_mes = sprintf('Error:  Unstable Structure');
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif aflag == -1
	   status_mes = sprintf('Incr # %i, Applied Load Ratio = %5.3f  -->  Analysis Halted: Limit Reached',total,load_val);
	   set(h_stat_mes,'String',status_mes); drawnow;
	elseif isinf(aflag)
	   apratios = [];
	   status_mes = sprintf('Error:  No Analysis Code Available');
	   set(h_stat_mes,'String',status_mes); drawnow;
	end
    
end

if aflag == 0, 
	analysis_info =[two_dim truss anatype 0];
	return;
end

if truss
	for  dof = 4:6
		for node = 1:nnodes
			react(node,dof) = 0;
		end
	end

	for  dof = 4:6
		for ele = 1:nele
			ele_for(ele,dof) = 0;
			ele_for(ele,dof+6) = 0;
		end
	end
end

if anatype == 7
    periods_info = [ sol_scheme anatype_prev incr_size incr_nums periods];
	if sol_scheme == 1
		analysis_info(1:2) = [two_dim truss];
	end
		analysis_info(3) = anatype;
else
	analysis_info = [two_dim truss anatype limit_state apratios'];
end
