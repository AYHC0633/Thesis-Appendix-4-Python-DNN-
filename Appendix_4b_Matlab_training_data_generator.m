% clc 
% close all
% clear
%% calculator select 
    Calculator_1 = 1; % AF array factor method
    Calculator_2 = 0; % Matlab method
    CSTPhase_lodar = 0 ;% loading phase result
    Chebyshev_Amplitude_taper = 0; %Amplitude taper for CASSIOPeiA array
   
    load2ndenable = 0;
RRtheta = 90;%80:0.5:100% from 75 to 115
RRPhi = 0:0.1:360;
% polar plot axis
    DTR = pi/180; %DegToRad
    RTD = 180/pi; %DegToRad
  
%% phase Calculator 1 

if Calculator_1 == 1
    fprintf('[%s] Calculator 1 running... \n', datestr(now,'HH:MM:SS'));

    theta = linspace(0,360,361); % 2 thread for 2 second 200 runtime =2 min
    phi = linspace(0,360,361);
    LL = -40; %minium axis limit, in power dB
    UL = 0;  %upper axis limit, in power dB

%Array factor equation 
    C= physconst('LightSpeed');
    F=2.45e9;
    WL=C/F;
    k =(2*pi)/WL; % wavelength number
    PhaShif = 0*DTR;% Phase Shift
    d = WL/4; % Spacing
    r = 100; %distance from circular array
    U = k*d*cos(theta) + PhaShif; % ou
    
    AmpCon = exp(-j*k*r)/r; % please refer to Balans: Antenna theory Ch6 p80 (6-113)
    In = 1; % phase excitation 
        
        % Unit element config.
        ElementSpacing = WL/4;
       
        Nu = 1; % number of element 
        
        % UnitCell config.
        UNu = 8; %Number of Unitcell
        dx = 0.1; % Unitcell spacing
        LayPhaShifx = 0;
        Clk = -1; % choose -1 or 1

        %Layer config 
        LNu = 8; % No. of Layer
        dz = 0.1;
        StrPhaShifz = 0;	

        Spara = readCSTSparameter();

   for Rtheta = 1:length(RRtheta)
        for RPhi =1:length(RRPhi)
        fprintf('[%s] progress %0.2f %% ... \n', datestr(now,'HH:MM:SS'),(Rtheta*length(RRPhi)+RPhi)/(length(RRtheta)*length(RRPhi))*100);
        % Angle steering 
        DStTh = RRtheta(Rtheta);
        DStPh = RRPhi(RPhi);
        StTh = DStTh; % desired beam angle Theta in radian (up and down), in azuith plane
        StPh = DStPh; % desired beam angle Phi in radian (left or right)
        Iamp = ones(1,8)*3.5;%[1:4  4:-1:1];

        Rotate = pi/4;
    % Array factor Equation/ nomalise electric field 
        ESum =0; % initial value
        [Th,Ph] = meshgrid(theta,phi);
    
        if mod(UNu,2) == 0
            rtpos = (UNu/2+0.5);
        else 
            rtpos = round(UNu/2);
        end
    
        Lay = 1:LNu;
        Nel = 1:UNu ;
        SpaceShift = 0.14;
        MidSpacing = [zeros(1,UNu/2) ones(1,UNu/2)]*0.28;
    
        LayPhaShifx = -(((Nel-rtpos)*dx+MidSpacing-SpaceShift)*k*sind(StTh)).*(cosd(StPh-(Clk*(Lay-1)*(180/LNu))))';% steering config in Layer
        StrPhaShifz = -ones(1,length(Nel)).*((Lay-1)*k*dz*cosd(StTh))'; % steering config in structure
        PhaShifTab = (LayPhaShifx + StrPhaShifz)';



    % %Power level conversion    
    %     PL= 10*log10(((Eview.*conj(Eview))/(Nu*UNu*LNu))*1.6);

    %phase shift table for above structure calculation
        
    
         Exp_AFPhaseShift = exp(j*PhaShifTab(:));
         MutalCoupling =Reflection_Coeffic(StTh,StPh,Spara);
         RealisticPhaseshift= angle(Exp_AFPhaseShift+MutalCoupling);

         PPhaShifTab = mod(RealisticPhaseshift,2*pi);
         PPPhaShifTab = PPhaShifTab(:)';

         thetaNormalise = RRtheta(Rtheta)-RRtheta(1);
         databasestore(RPhi+(length(RRPhi)*(Rtheta-1)),:) =[RRPhi(RPhi) thetaNormalise PPPhaShifTab];
         % fprintf('theta %d , phi %d \n', Rtheta,RPhi);
        end
   end

   % %----------Phase Check
phasecheck = 1;
if phasecheck == 1
    
        phiSter =38;
        thetaSter = 90;
        StPhi = find(RRPhi==phiSter);
        Sttheta = find(RRtheta==thetaSter);
        InputPhase = databasestore(StPhi+(length(RRPhi)*(Sttheta-1)),3:end);
        CSTPhase_lodar = InputPhase;
        InputPhase = reshape(InputPhase,8,8)';
        
        % InputPhase = PhaShifTab;%phase shift
        for ith = 1:length(theta)     
            for iph = 1:length(phi)
         
                 ESum = exp(j*(...
                     (((Nel-rtpos)*dx+MidSpacing-SpaceShift)*k*sind(theta(ith)).*(cosd(phi(iph)-(Clk*(Lay-1)*(180/LNu))))') + ...
                     (ones(1,length(Nel)).*((Lay-1)*k*dz*cosd(theta(ith)))') + InputPhase ...
                     ));

                 Eview(iph,ith) = sum(sum(ESum));
                 EviewN(iph,ith) = Eview(iph,ith)/(64);
             end 
        end

            %figure
           hold on
            x = phi; %x = -179:1:180
            StPhi = find(RRPhi==phiSter);
            PlotTheta = find(theta==90);
            plot(x,abs(EviewN(:,PlotTheta)))%(x,abs(Eview(:,90)))
            xlabel('\phi (deg.)')
            ylabel('Amplitude')
            title('Electric field in AF equation')
            axis([0 360 -inf inf])
            grid on
             
end

       
       Title = ["Phi" "Theta" (repelem(["Element"],64)'+num2str((1:64)'))']
       filename ='HelicalPhaseTestV3.csv';
       writematrix(Title,filename,'WriteMode','append');
       writematrix(databasestore,filename,'WriteMode','append');
       
       if load2ndenable == 1
%          %PhaseShift Table printf
            fprintf('[%s] Phase loaded in txt file... \n', datestr(now,'HH:MM:SS'));
            fileID = fopen('MTC_phaseCalibration.txt','w');
            fprintf(fileID,'Cal_1 CASSIOPeiA %i layer x %i Unitcell Array phase result  \n \n',LNu,UNu);
            fprintf(fileID,'Layer | Unitcell | Element0 | Element120 | Element240 \n');
            fprintf(fileID,'----------------------------------------------------- \n');
            for L = 1:LNu
                 for U = 1:UNu
                 fprintf(fileID,'   %i       %i',L,U);
                    for E = 1:Nu
                    fprintf(fileID,'        %0.6f',PPhaShifTab((L-1)*UNu*Nu+(U-1)*Nu+E));
                    end
                 fprintf(fileID,'\n');
                 end
            fprintf(fileID,'\n');
            end

            fclose(fileID);
            %type 'MTC_phaseCalibration.txt'
      
       
            fprintf('[%s] Phase loaded in 2nd txt file... \n', datestr(now,'HH:MM:SS'));
            fileID = fopen('Cal_1_CST_phase_explore.txt','w');
            fprintf(fileID,'Cal_1 CASSIOPeiA %i layer x %i Unitcell Array phase result  \n \n',LNu,UNu);
            fprintf(fileID,'----------------------------------------------------- \n');
            for port = 1:LNu*UNu*Nu
                fprintf(fileID,'.ExcitationPortMode "%2.0f", "1", "1", "%4.2f", "default", "True" \n',port,PPhaShifTab(port));
            end
       end
            fprintf('[%s] Calculator 1 export result completed... \n', datestr(now,'HH:MM:SS'));
         

end

function AF_Phase = AF_Calculator(phi,theta)
    Lay = 1:LNu;
    Nel = 1:UNu ;
    SpaceShift = 0.14;
    MidSpacing = [zeros(1,UNu/2) ones(1,UNu/2)]*0.28;

    LayPhaShifx =-(((Nel-rtpos)*dx+MidSpacing-SpaceShift)*k*sin(StTh)).*(cos(StPh-((Lay-1)*(180/LNu)*DTR)))';% steering config in Layer
    StrPhaShifz = -ones(1,length(Nel)).*((Lay-1)*k*dz*cos(StTh))'; % steering config in structure
    PhaShifTab = ((LayPhaShifx + StrPhaShifz)*RTD)';

    AF_Phase = PhaShifTab ;
end


function  sumofSmn = Reflection_Coeffic(theta,phi,Spar) % angle in radian

    DTR = pi/180; %DegToRad
    RTD = 180/pi; %DegToRad
    C= physconst('LightSpeed');
    F=2.45e9;
    WL=C/F;
    k =(2*pi)/WL; % wavelength number
    PhaShif = 0*DTR;% Phase Shift
    dx = 0.1;%WL/4; % Spacing
    dz = 0.1; 
    rtpos =4.5000;

    Lay = 1:8;
    Nel = 1:8 ;
    SpaceShift = 0.14;
    MidSpacing = [zeros(1,8/2) ones(1,8/2)]*0.28;
    disarraySpac = ((Nel-rtpos)*dx+MidSpacing-SpaceShift);

    phase = (-((k*disarraySpac*sind(theta)).*cosd(phi-((Lay-1)*(180/8)))')-ones(1,8).*((Lay-1)*k*dz*cosd(theta))')';
    sumofSmn = Spar*(exp(j*phase(:)'))';%sum(SParameter.*exp(-j*phase(:)),2);
    % Reflect = exp(j*k*m*d*sin(angle*DTR))*sumofSmn;
end


function SParameter = readCSTSparameter()

    Tat = readtable('SparameterData8x8inMag.txt');
    TaT = Tat(:,1:2);
    AA = table2array(TaT);
    Position=find(round(AA(:,1),4)==2.45);
    SParamete = AA(Position,2);
    SParameter = reshape(SParamete,64,64);% [S11 S12 S13.....;S21 S22 S23 .....;]
%     SParameter = 10.^(SParameter/10);
end
