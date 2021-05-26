%%% New ���Գ������㷽��
%%%���汾����New_opf_main_1_3_prde_activepowercompara,
%% ---- ���ɵ��ɾ������Ϣ -------------------------------------------------
%�����������
mpc=ext2int(case1354xy);
[baseMVA,bus,gen,branch,gencost]=deal(mpc.baseMVA,mpc.bus,mpc.gen,mpc.branch,mpc.gencost);
%���ɵ��ɾ���
Ybus=makeYbus(baseMVA,bus,branch);
Gbus=real(Ybus);%ʵ��
Bbus=imag(Ybus);%�鲿
%������Ϣ
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;%ĸ�߲����к�
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;%����������к�
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;%������ɱ������к�
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;%ĸ�߲����к�
[ref, pv, pq] = bustypes(bus,gen);%�ڵ�����
nbus=size(bus,1);%ĸ������
branch=branch(branch(:,BR_STATUS)==1,:);
nbranch=size(branch,1);
gen=gen(gen(:,GEN_STATUS)==1,:);
ng=size(gen,1);%���������
No_ref=[1:ref-1,ref+1:nbus];%�ǲο��ڵ�
%����ֵ
sumpload=[3959];%	3733	3280	3016	3129	3393	3959	4411	4637	4977	5033	4751	4524	4298	4977	5090	4807	5033	5116	5184	5090	5014	4901	4581];
Pload=bus(:,PD).*sumpload/sumpload(1);
Qload=bus(:,QD).*sumpload/sumpload(1);
%ʱ�������
nt=size(sumpload,2);
%��������
Pmax=gen(:,PMAX);
Pmin=gen(:,PMIN);
Qmax=gen(:,QMAX);
Qmin=gen(:,QMIN);
%ȡ��Ψһ����·���
unibranch=unique(branch(:,1:2),'rows');
nunibranch=size(unibranch,1);%ȥ���ظ������·����

%%---- ������� ------------------------------------------------------------
V_2=sdpvar(nbus,nt,'full');
theta=sdpvar(nbus,nt,'full');
theta_xi=sdpvar(nunibranch,nt,'full');
theta_vi=sdpvar(nunibranch,nt,'full');
theta_ti=sdpvar(nunibranch,nt,'full');
P=sdpvar(ng,nt,'full');
Q=sdpvar(ng,nt,'full');

%%---- �������� ------------------------------------------------------------
Vgen_bus=sparse(gen(:,1),1:ng,1,nbus,ng);
VFunibranch_bus=sparse(unibranch(:,1),1:nunibranch,1,nbus,nunibranch);
VTunibranch_bus=sparse(unibranch(:,2),1:nunibranch,1,nbus,nunibranch);
VFTunibranch_bus_sum=VFunibranch_bus+VTunibranch_bus;
VFTunibus_branch_sum=VFunibranch_bus'+VTunibranch_bus';
VFTunibus_branch_diff=VFunibranch_bus'-VTunibranch_bus';
VFbus_branch=sparse(1:nbranch,branch(:,1),1,nbranch,nbus);
VTbus_branch=sparse(1:nbranch,branch(:,2),1,nbranch,nbus);
VFTbus_branch_diff=VFbus_branch-VTbus_branch;
VFTbus_branch_sum=VFbus_branch+VTbus_branch;

%%---- ��������+���� -------------------------------------------------------
%�γ�Bbus_����
Bbus_=-Bbus;
Bbus_(logical(eye(nbus,nbus)))=sum(Bbus.*(1-eye(nbus,nbus)),2);
%�γ�Gbus_����
Gbus_=Gbus;
Gbus_(logical(eye(nbus,nbus)))=-sum(Gbus.*(1-eye(nbus,nbus)),2);
% %�γɷ�������ʼ��޲���(Ϲ���)
% alpha1=Pmax./(Qmin/4-Qmax/4);
% alpha2=Pmax./(Qmax/4-Qmin/4);
%ref,pv,�ڵ��ѹԼ��
V_2_ref=bus(ref,VM).^2;
% V_2_pv=bus(pv,VM).^2;
V_2_MAX=bus(:,VMAX).^2;
V_2_MIN=bus(:,VMIN).^2;
%�ο��ڵ����Լ��
theta_ref=bus(ref,VA)*pi/180;
%�γɷ�����Ĺ���+����Լ��
alpha_g=[-ones(ng,1),ones(ng,1),zeros(ng,1),zeros(ng,1)];%�������
beta_g=[zeros(ng,1),zeros(ng,1),-ones(ng,1),ones(ng,1)];
gamma_g=[Pmin,-Pmax,Qmin,-Qmax];
% �γɶ�Ӧ��Gbus��Bbus�ĵ�������
II=sub2ind([nbus,nbus],unibranch(:,1),unibranch(:,2));
Gbranch=Gbus(II);
Bbranch=Bbus(II);
%%---- �γɼ�С�ɳ������� -----------------------------------------------
%��ѹ����
Vmin=max(bus(:,VMIN));
Vmax=max(bus(:,VMAX));
a=-0.5;
Va2=Vmin^2:0.005:Vmax^2;
Vb2=Va2;
Vab=-sqrt(Va2'.*Vb2);
b=1.05*max([-min(a*Va2'+a*Vb2-Vab,[],'all'),0]);
%�ǶȲ���
c=0.125;

%% ---- ģ�͹��� ------------------------------------------------------------
con=[];
for t=1:nt
    %�ڵ��й�ƽ��
    con=con+[(-Vgen_bus*P(:,t)+Pload(:,t)+baseMVA*(diag(Gbus).*V_2(:,t)+Bbus_*theta(:,t)+...
        (-1)*VFTunibranch_bus_sum*(Gbranch.*theta_xi(:,t)))==0)];
    %�ڵ��޹�ƽ��
    con=con+[(-Vgen_bus*Q(:,t)+Qload(:,t)+baseMVA*(-diag(Bbus).*V_2(:,t)-Gbus_*theta(:,t)+...
        VFTunibranch_bus_sum*(Bbranch.*theta_xi(:,t)))==0)];
    %�ɳ�Լ��;
    for n=1:nunibranch
        con=con+[cone([2*theta_vi(n,t);V_2(unibranch(n,1),t)-V_2(unibranch(n,2),t)],V_2(unibranch(n,1),t)+V_2(unibranch(n,2),t))];
    end
    con=con+[theta_vi(:,t)<=a*VFTunibus_branch_sum*V_2(:,t)+b];
    for n=1:nunibranch
        con=con+[cone([theta(unibranch(n,1),t)-theta(unibranch(n,2),t);theta_ti(n,t);1],theta_ti(n,t)+1)];
    end
    con=con+[theta_ti(:,t)<=c];
    con=con+[theta_xi(:,t)==theta_vi(:,t)+theta_ti(:,t)];
    %�����������ΧԼ��
    for g=1:size(alpha_g,2)
        con=con+[(alpha_g(:,g).*P(:,t)+beta_g(:,g).*Q(:,t)+gamma_g(:,g)<=0)];
    end
    %�ڵ��ѹԼ��
    con=con+[(V_2_MIN-V_2(:,t)<=0)];
    con=con+[(V_2(:,t)-V_2_MAX<=0)];
    %�ο��ڵ����
    con=con+[(theta(ref,t)-theta_ref==0)];
end

%%---- Ŀ�꺯�����ܳɱ���С�� -----------------------------------------------
obj=sum(P(:,:),2)'*gencost(1:ng,COST+1)+sum(Q(:,:),2)'*gencost(ng+1:2*ng,COST+1);

%% ---- ��� ---------------------------------------------------------------
%���ò���
ses=sdpsettings('solver','gurobi');
%����
tic
dd=solvesdp(con,obj,ses);
if dd.problem~=0
    error('ģ�Ͳ�����!');
end
toc
%% ---- ������� ------------------------------------------------------------
% V=double(V_2).^0.5;
% theta=double(theta);
% P=double(P);
% Q=double(Q);
% obj=double(obj);
% for n=1:nunibranch
%     ess_voltage(n)=double(V_2(unibranch(n,1),t)+V_2(unibranch(n,2),t)-norm([2*theta_vi(n,t);V_2(unibranch(n,1),t)-V_2(unibranch(n,2),t)]));
%     ess_theta(n)=double(theta_ti(n,t)+1-norm([theta(unibranch(n,1),t)-theta(unibranch(n,2),t);theta_ti(n,t);1]));
% end

% lamda_P=dual(con(['lamda_P']));
% lamda_Q=dual(con(['lamda_Q']));
% condition=-Gbranch.*(lamda_P(unibranch(:,1))+lamda_P(unibranch(:,2)))+Bbranch.*(lamda_Q(unibranch(:,1))+lamda_Q(unibranch(:,2)));