function Pars = Par_Read_Struct(parfil)

% function Pars = Par_Read_Struct(parfil)

ParMATRIX = textread(parfil,'%s','delimiter','\n','whitespace','');

Pars.seqfil= '';
Pars.studyUID='';
Pars.nD=3;
Pars.pss=0;
Pars.pss0=0;
Pars.ppe2= 0;
Pars.gap=0;

Pars.lro= 0;
Pars.lpe= 0;
Pars.lpe2= 0;

Pars.np=0;
Pars.nv=0;
Pars.nv2= 0;
Pars.ns=0;
Pars.NS=0;

Pars.Slice_Thickness= 0;
Pars.thk= 0;
Pars.pro= 0;

Pars.theta=0;
Pars.phi=0;
Pars.psi=0;
Pars.orient= '';

Pars.te=0;
Pars.tr=0;
Pars.sw= 0;

Pars.p1pat= '';
Pars.FlipA= 0;
Pars.p1= 0;
Pars.p2= 0;
Pars.tpwr1= 0;
Pars.tpwr2= 0;

Pars.FC_= 0;
Pars.RO_CVRD= 0;
Pars.PE1_CVRD= 0;
Pars.PE2_CVRD= 0;
Pars.sat_type1= 0;


Pars.Double_precision= '';
Pars.RCVRS_= 0;
Pars.time_complete= '';
Pars.rfcoil= '';
Pars.sfrq= 0;
Pars.B0= 0;

Pars.acqcycles= 0;
Pars.nt=0;
Pars.weight=0;
Pars.age=0;
Pars.gender= '';
Pars.name= '';
%ParMATRIX=textread('H:\47_DAta\s_20070115_02\09.fid\procpar');

% Parameters' locations list (by line number in the procpar file)
% nv2   961
% nv    946
%RO_CVRD 34
%PE2_CVRD 22
%PE1_CVRD 31
%pss    1171
%FC_    25
%te     1607
%TR     1718
%lpe2   849


% pss    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='pss '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            pss=str2num(Line2(3:Length));
            break
        end
    end
end

Pars.pss=pss;
% gap
gap=0;
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='gap '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            gap=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.gap=gap;




% pss0    ... for 2D gems case    
pss0=0;
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>4
        if Line1(1:5)=='pss0 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            pss0=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.pss0=pss0;

% nD      Dimension of scan (2D or 3D)    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='nD '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            nD=str2num(Line2(3:Length));
            break
        end
    end
end


Pars.nD=nD;

% studyid_
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>8
        if Line1(1:9)=='studyid_ '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            studyUID=Line2(6:Length-1);
            break
        end
    end
end
studyUID=strcat(studyUID(1:8),studyUID(10:11));
Pars.studyUID=studyUID;

% np
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='np '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            np=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.np=np;



% te
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='te '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            te=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.te=te;



% tr    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='tr '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            tr=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.tr=tr;
    
% nv    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='nv '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            nv=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.nv=nv;

% ns    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='ns '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            ns=str2num(Line2(3:Length));
            break
        end
    end
end

Pars.ns=ns;


NS=1;
% NS    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='NS '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            NS=str2num(Line2(3:Length));
            break
        end
    end
end
if ns > NS
    NS=ns;
end
Pars.NS=NS;

% nt    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='nt '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            nt=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.nt=nt;    
    
% weight
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:7)=='weight '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            weight=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.weight=weight;



% age    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='age '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            age=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.age=age;    

% psi    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='psi '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            psi=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.psi=psi;

% phi    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='phi '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            phi=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.phi=phi;

% theta    
theta=0;
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>5
        if Line1(1:6)=='theta '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            theta=str2num(Line2(3:Length));
            break
        end
    end
end

Pars.theta=theta;

% % theta    
% for i=1:max(size(ParMATRIX))
%     Line1=cell2mat(ParMATRIX(i));
%     if max(size(Line1))>9
%         if Line1(1:6)=='theta '
%             Line2=cell2mat(ParMATRIX(i+1));
%             Length=max(size(Line2));
%             theta=Line2(4:Length-1); % different indexing to avoide the (((("))))s
%             break
%         end
%     end
% end
% 


% seqfil    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:7)=='seqfil '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            seqfil=Line2(4:Length-1);
            break
        end
    end
end
Pars.seqfil=seqfil;

% dp (double precission ??? )    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:3)=='dp '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            Double_precision=Line2(4);
            break
        end
    end
end
Pars.Double_precision=Double_precision;


% p1pat
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>5
        if Line1(1:6)=='p1pat '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            p1pat=Line2(4:Length-1);
            break
        end
    end
end
Pars.p1pat=p1pat;

% thk    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='thk '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            thk=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.thk=thk;

% gap    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='gap '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            gap=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.gap=gap;
    

    
% rcvrs    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>5
        if Line1(1:6)=='rcvrs '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            rcvrs=Line2(4:Length-1); % different indexing to avoide the (((("))))s
            break
        end
    end
end



RCVRS_=0;

for i=1:Length-4
    if rcvrs(i)=='y'
        RCVRS_=RCVRS_+1;
    end
end
Pars.RCVRS_=RCVRS_;


% lpe 
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='lpe '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lpe=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.lpe=lpe;


% lro    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='lro '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lro=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.lro=lro;


% pro    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>3
        if Line1(1:4)=='pro '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            pro=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.pro=pro;


% fliplist    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>8
        if Line1(1:9)=='fliplist '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            FlipA=str2num(Line2(3:5));
            break
        end
    end
end
Pars.FlipA=FlipA;



% orient    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:7)=='orient '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            orient=Line2(4:Length-1);
            break
        end
    end
end
Pars.orient=orient;

% time_complete    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>13
        if Line1(1:14)=='time_complete '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            time_complete=Line2(4:Length-1);
            break
        end
    end
end
Pars.time_complete=time_complete;

% gender    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:7)=='gender '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            gender=Line2(4:Length-1);
            break
        end
    end
end

Pars.gender=gender;

% name    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>4
        if Line1(1:5)=='name '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            name=Line2(4:Length-1);
            break
        end
    end
end

Pars.name=name;

% rfcoil    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>6
        if Line1(1:7)=='rfcoil '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            rfcoil=Line2(4:Length-1);
            break
        end
    end
end
Pars.rfcoil=rfcoil;

% sw    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='sw '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            sw=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.sw=sw;


% sfrq   
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>4
        if Line1(1:5)=='sfrq '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            sfrq=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.sfrq=sfrq;



% B0    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>2
        if Line1(1:3)=='B0 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            B0=str2num(Line2(3:Length));
            break
        end
    end
end
Pars.B0=B0;

if nD==3
	% nv2    
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>3
            if Line1(1:4)=='nv2 '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                nv2=str2num(Line2(3:Length));
                break
            end
        end
    end
    Pars.nv2=nv2;
	
	% FC_    
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>3
            if Line1(1:4)=='FC_ '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                FC_=str2num(Line2(3:Length));
                break
            end
        end
	end
	if exist('FC_')
        FC_=FC_;
	else
        FC_=0;
    end
    Pars.FC_=FC_;
    
    
	% Pars.ppe2    
    ppe2=pss;
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>4
            if Line1(1:5)=='ppe2 '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                ppe2=str2num(Line2(3:Length));
                break
            end
        end
    end
	
    Pars.ppe2=ppe2;

	% RO_CVRD
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>7
            if Line1(1:8)=='RO_CVRD '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                RO_CVRD=str2num(Line2(3:Length));
                break
            end
        end
	end
	Pars.RO_CVRD=RO_CVRD;
        
	% PE1_CVRD    
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>8
            if Line1(1:9)=='PE1_CVRD '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                PE1_CVRD=str2num(Line2(3:Length));
                break
            end
        end
	end
	Pars.PE1_CVRD=PE1_CVRD;	
        
	% PE2_CVRD    
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>8
            if Line1(1:9)=='PE2_CVRD '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                PE2_CVRD=str2num(Line2(3:Length));
                break
            end
        end
	end
	Pars.PE2_CVRD=PE2_CVRD;		
	
	
	% lpe2    
	for i=1:max(size(ParMATRIX))
        Line1=cell2mat(ParMATRIX(i));
        if max(size(Line1))>4
            if Line1(1:5)=='lpe2 '
                Line2=cell2mat(ParMATRIX(i+1));
                Length=max(size(Line2));
                lpe2=str2num(Line2(3:Length));
                break
            end
        end
    end
    Pars.lpe2=lpe2;
	Slice_Thickness=lpe2*10/nv2;  % mm
    
end
% if seqfil == 'gems'
%     nv2=1; lpe2=0; Pars.ppe2=10000; RO_CVRD=2; FC_=0; Slice_Thickness=thk; PE1_CVRD=2; PE2_CVRD=2;
%     save(strcat(Dir_,'params'),'np', 'nv', 'ns', 'te', 'tr', 'FlipA', 'lpe', 'thk', 'lro', 'seqfil', 'sw', 'RCVRS_', 'B0')  % 'pss0',
% else
%     save(strcat(Dir_,'params'),'np', 'nv', 'nv2', 'te', 'tr', 'FlipA', 'lpe', 'lpe2', 'lro', 'seqfil', 'ppe2', 'RO_CVRD', 'PE1_CVRD', 'PE2_CVRD', 'Slice_Thickness', 'sw', 'RCVRS_', 'Slice_Thickness', 'B0' , 'FC_')
% end






% SOME PARS ?
% =====================================================================
% 
% np/2        te      tr      nv      nv2     pss    FC_    
%     
% lpe2        nv2     RO_CVRD     PE1_CVRD    PE2_CVRD    
% 
% RCVRS_      lpe     lro         FlipA       Pars.ppe2    
% 
% seqfil      sw      B0          lpe2
% 
% Slice_Thickness





% UNUSED IN DCMWRITE ????????????
% 
% 
% ********************************
%      pss    FC_    
%     
%              RO_CVRD     PE1_CVRD    PE2_CVRD    
% 
% RCVRS_      lpe     lro         FlipA       Pars.ppe2    
% ********************************
% % % AZZY1=['pss   ' num2str(pss) 'np/nv/nv2/NS   ' num2str(np/NS) '/' num2str(nv) '/' num2str(nv2) '/' num2str(NS)  ]
% % % AZZY2=[  'te/tr/thk   ' num2str(te) '/' num2str(tr) '/' num2str(thk)]
% % % AZZY3=[ 'RCVRS_/lpe2/FlipA   ' num2str(RCVRS_) '/' num2str(lpe2) '/' num2str(FlipA) ]
% % % AZZY4= [ 'ppe2   ' num2str(Pars.ppe2) '   Scan Time   ' num2str(nv*nv2*tr) 'pulse shape' num2str(p1pat) ]
% % % AZZY5= [ 'RO_CVRD/sw/FC_   ' num2str(RO_CVRD) '/' num2str(sw) '/'  num2str(FC_) ]


% sprintf('%s\n%s%s\n%s%s\n%s%s\n%s%s\n', Dir_, 'np/NS = ', num2str(np/NS), 'nv = ' , num2str(nv), 'nv2 = ' , num2str(nv2), 'NS = ' , num2str(NS))
% sprintf('%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n','FlipA = ', num2str(FlipA), 'pulse shape = ', p1pat,'sw = ' , num2str(sw), 'te = ', num2str(te), 'tr = ' ,num2str(tr))
% sprintf('%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n', 'thk = ', num2str(thk), 'lpe2 = ',num2str(lpe2),'pss = ', num2str(pss), 'Pars.ppe2 = ', num2str(Pars.ppe2),'Scan Time = ', num2str(nv*nv2*tr))
% sprintf('%s%s\n%s%s\n%s%s\n', 'RO_CVRD = ',num2str(RO_CVRD),  'FC_ = ', num2str(FC_), 'RCVRS_ = ' ,num2str(RCVRS_))

% sprintf('%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n\n%s%s\n%s%s\n%s%s\n','DIR = ', Dir_, 'np/NS = ', num2str(np/NS), 'nv = ' , num2str(nv), 'nv2 = ' , num2str(nv2), 'NS = ' , num2str(NS),'FlipA = ', num2str(FlipA), 'pulse shape = ', p1pat,'sw = ' , num2str(sw), 'te = ', num2str(te), 'tr = ' ,num2str(tr), 'thk = ', num2str(thk), 'lpe2 = ',num2str(lpe2),'pss = ', num2str(pss), 'Pars.ppe2 = ', num2str(Pars.ppe2),'Scan Time = ', datestr(nv*nv2*tr/24/60/60,13), 'RO_CVRD = ',num2str(RO_CVRD),  'FC_ = ', num2str(FC_), 'RCVRS_ = ' ,num2str(RCVRS_))

% % % % sprintf('%s\n%s\n%s\n%s\n%s\n%s\n ', Dir_, AZZY1, AZZY4, AZZY3, AZZY2, AZZY5)

% ParMATRIX = textread(parfil,'%s','delimiter','\n','whitespace','');


%ParMATRIX=textread('H:\47_DAta\s_20070115_02\09.fid\procpar');

% Parameters' locations list (by line number in the procpar file)
% nv2   961
% nv    946
%RO_CVRD 34
%PE2_CVRD 22
%PE1_CVRD 31
%pss    1171
%FC_    25
%te     1607
%TR     1718
%lpe2   849


% np
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='np '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            np=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.np=np;
% ns
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='ns '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            ns=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.ns=ns;
% NS
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='NS '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            NS=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.NS=NS;
% te    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='te '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            te=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.te=te;
% tr    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='tr '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            tr=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.tr=tr;
    
% nv    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='nv '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            nv=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.nv=nv;

% p1    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='p1 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            p1=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.p1=p1;

% p2    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='p2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            p2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.p2=p2;

nv2=1;    
% nv2    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='nv2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            nv2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.nv2=nv2;
    
    
% pss    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='pss '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            pss=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.pss=pss;
    
% thk    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='thk '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            thk=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.thk=thk;


FC_=0;    
% FC_    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='FC_ '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            FC_=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.FC_=FC_;

lpe2=1;    
% lpe2    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:5)=='lpe2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lpe2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.lpe2=lpe2;

    
% % nv2    
% for i=1:max(size(ParMATRIX))
%     Line1=cell2mat(ParMATRIX(i));
%     if max(size(Line1))>9
%         if Line1(1:4)=='nv2 '
%             Line2=cell2mat(ParMATRIX(i+1));
%             Length=max(size(Line2));
%             nv2=str2double(Line2(3:Length));
%             break
%         end
%     end
% end

RO_CVRD=1;    
% RO_CVRD
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:8)=='RO_CVRD '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            RO_CVRD=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.RO_CVRD=RO_CVRD;


PE1_CVRD=1;    
% PE1_CVRD    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:9)=='PE1_CVRD '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            PE1_CVRD=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.PE1_CVRD=PE1_CVRD;


PE2_CVRD=1;
% PE2_CVRD    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:9)=='PE2_CVRD '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            PE2_CVRD=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.PE2_CVRD=PE2_CVRD;

rcvrs='y';
% rcvrs
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:6)=='rcvrs '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            rcvrs=Line2(4:Length-1); % different indexing to avoide the (((("))))s
            break
        end
    end
end

tpwr1=0;
% tpwr1
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:6)=='tpwr1 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            tpwr1=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.tpwr1=tpwr1;

tpwr2=0;
% tpwr2
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:6)=='tpwr2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            tpwr2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.tpwr2=tpwr2;

% RCVRS_=0;
% 
% for i=1:Length-4
%     if rcvrs(i)=='y'
%         RCVRS_=RCVRS_+1;
%     end
% end




% lpe 
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='lpe '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lpe=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.lpe=lpe;

% lro    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='lro '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lro=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.lro=lro;

% fliplist    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:9)=='fliplist '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            if max(size(Line2))>4
                FlipA=str2double(Line2(3:5));
                FlipA2=str2double(Line2(6:Length));
                FlipA= [FlipA FlipA2];
            else
                FlipA=str2double(Line2(3:Length));
            end
            break
        end
    end
end
Pars.FlipA=FlipA;


sat_type1=0;
%sat_type1
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:10)=='sat_type1 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            sat_type1=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.sat_type1=sat_type1;

% Pars.ppe2    
ppe2=pss;
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:5)=='ppe2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            Pars.ppe2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.ppe2=ppe2;

% seqfil    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:7)=='seqfil '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            seqfil=Line2(4:Length-1);
            break
        end
    end
end
Pars.seqfil=seqfil;


% sw    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='sw '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            sw=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.sw=sw;


% B0    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:3)=='B0 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            B0=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.B0=B0;

% lpe2    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:5)=='lpe2 '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            lpe2=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.lpe2=lpe2;




% psi
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='psi '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            psi=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.psi=psi;
% phi    
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:4)=='phi '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            phi=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.phi=phi;

if ns > 1
    lpe2=ns*thk+(ns-1)*gap;
    Pars.lpe2=lpe2;
end

% acqcycles    
acqcycles=1;
for i=1:max(size(ParMATRIX))
    Line1=cell2mat(ParMATRIX(i));
    if max(size(Line1))>9
        if Line1(1:10)=='acqcycles '
            Line2=cell2mat(ParMATRIX(i+1));
            Length=max(size(Line2));
            acqcycles=str2double(Line2(3:Length));
            break
        end
    end
end
Pars.acqcycles=acqcycles;
Slice_Thickness=lpe2*10/(nv2*NS);  % mm
Pars.Slice_Thickness=Slice_Thickness;


return
