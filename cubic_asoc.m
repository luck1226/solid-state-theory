% this is to calculate the spectrum of 3D cubic crystal
% which includes 4-orbital s,px,py,pz from cation and anion
% atoms

% total Hamiltonian is an 4x4 matrix
% for cubic crystal CsPbI3 or ABO3

clear all; clc;


L = 200;

N=8;


tSS = (-1)*0.25*1.9;     % Jin paper 0.25 and the minus factor is 
                         % is extremely important
                         % without the minus sign, the whole model
                         % won't work !!!

tSP = 0.4;      % Jin paper 0.4 % figure(b) uses 0.04 figure(c) uses 0.4

tPPs = 0.9;     % Jin paper 0.9

tPPp = 0.15;     % Jin paper 0.15


%on-site for I s-orbital
es = -1.5;

%on-site for Bi s-orbital % figure(a) use -2 while figure(b)(c) use 2
%eBs = 2;

%on-site for I p-orbital
ep = 4;

%on-site for Bi p-orbital
%eBp = 6;

% At I site, the split of Pz from Px and Py
% due to the crystal field

%eCFEA = -0.6;
%eCFEA = -0.2;
eCFEA = 0;

%eCFEB = 0.5;
%eCFEB = 0.1;
eCFEB = 0;


%eta = 0.5;
eta = 1;

lambda = 1;


Ham = zeros (N,N);

ee = zeros(L,N,N);

ee_dat = zeros(L,N);

vv = zeros(L,N,N);


% assigning on-site energies for A and B local orbitals
    
    
%Ham(1,1) = eAs;
    
    
%Ham(5,5) = eBs;
    
    
%Ham(2,2) = ep;
    
    
%Ham(3,3) = Ham(2,2);
    
    
%Ham(4,4) = Ham(2,2)+eCFEA;
    
    
%Ham(6,6) = eBp;
    
    
%Ham(7,7) = Ham(6,6);
    
    
%Ham(8,8) = Ham(6,6)+eCFEB;

jj = sqrt(-1);

% atom A spin up




for k =1:L

    value = pi*(k/L);

    if k<L/2
    
    
        kx = value*2;
    
    
        ky = kx;
    
    
        kz = kx;
        
    else
       
        kz = (pi-value)*2;
        
        ky = kz;
        
        kx = pi;
        
    end
    
    coX = 2*cos(kx);

    coY = 2*cos(ky);
    
    coZ = 2*cos(kz);
    
    snX = 2*sin(kx);
    
    snY = 2*sin(ky);
    
    snZ = 2*sin(kz);
    
   
    % s-orbital@A to s-orbital@B
    Ham(1,1) = es + tSS*(coX+coY+eta*coZ);
    
    %Ham(5,1) = conj(Ham(1,5));
    
    
    % s-orbital@A to px,py,pz-orbital@B
    Ham(1,2) = tSP*snX;
    
    Ham(2,1) = conj(Ham(1,2));
    
    Ham(1,3) = tSP*snY;
    
    Ham(3,1) = conj(Ham(1,3));
    
    Ham(1,4) = tSP*eta*snZ;
    
    Ham(4,1) = conj(Ham(1,4));
    
    
    % s-orbital@B to p-orbital@A
    %Ham(5,2) = tSP*snX;
    
    %Ham(2,5) = conj(Ham(5,2));
    
    %Ham(5,3) = tSP*snY;
    
    %Ham(3,5) = conj(Ham(5,3));
    
    %Ham(5,4) = tSP*eta*snZ;
    
    %Ham(4,5) = conj(Ham(5,4));
    
    
    
    
    % px-orbital@A to p-orbital@B
    
    Ham(2,2) = ep + tPPs*coX+tPPp*(coY+eta*coZ);
    
    %Ham(6,2) = conj(Ham(2,6));
    
    %Ham(2,7) = t3*snY;
    
    %Ham(7,2) = conj(Ham(2,7));
    
    %Ham(2,8) = t3*snZ;
    
    %Ham(8,2) = conj(Ham(2,8));
    
    
    
    % py-orbital@A to p-orbital@B
    
    %Ham(3,6) = t3*snX;
    
    %Ham(6,3) = conj(Ham(3,6));
    
    Ham(3,3) = ep + tPPs*coY+tPPp*(coX+eta*coZ);
    
    %Ham(7,3) = conj(Ham(3,7));
    
    %Ham(3,8) = t3*snZ;
    
    %Ham(8,3) = conj(Ham(3,8));
    
    
    
    %pz-oribtal@A to p-orbital@B
    
    %Ham(4,6) = t3*snX;
    
    %Ham(6,4) = conj(Ham(4,6));
    
    %Ham(4,7) = t3*snY;
    
    %Ham(7,4) = conj(Ham(4,7));
    
    Ham(4,4) = ep + tPPs*eta*coZ+tPPp*(coY+coX);
    
    %Ham(8,4) = conj(Ham(4,8));
    
    
    Ham(5:8,5:8) = Ham(1:4,1:4);
    
    
    Ham(2,3) = (-1)*jj*lambda/2;

    Ham(3,2) = conj(Ham(2,3));

    Ham(2,8) = lambda/2;

    Ham(8,2) = conj(Ham(2,8));

    Ham(3,8) = (-1)*jj*lambda/2;

    Ham(8,3) = conj(Ham(3,8));
    
    Ham(6,7) = jj*lambda/2;

    Ham(7,6) = conj(Ham(6,7));

    Ham(6,4) = (-1)*lambda/2;

    Ham(4,6) = conj(Ham(6,4));

    Ham(7,4) = (-1)*jj*lambda/2;

    Ham(4,7) = conj(Ham(7,4));
      
    
    [vv(k,:,:),ee(k,:,:)] = eig(Ham,'nobalance');
    
    
    ee_dat(k,:) = eig(Ham,'nobalance');

end


%figure(2);

colorMap = zeros(L,N);


for i=1:L
    
   
    for j=1:N
        
    
        % component falling on s-orbital @ I
        temp1 = vv(i,1,j);
        
        % component falling on s-oribtal @ Bi
        temp2 = vv(i,5,j);
        
        colorMap(i,j) = temp1*conj(temp1)+temp2*conj(temp2);
        
        
        
    end
    
    
end


figure(1);


for j=1:N

    
    %plot((1:L),ee(1:L,j,j),'b-'); hold on;
    scatter((1:L),ee(1:L,j,j),6,colorMap(:,j)); hold on;
    
    %for k=1:L
    
    

    
        %scatter(j,ee(j,k,k)); hold on;
    
        %scatter(k,ee(1:L,j,j)); hold on;

        
    %end
        
end



%Ea1 = (eAs+eBp)/2+sqrt(12*tSP^2+(eAs-eBp)^2/4);


%Ea2 = (eAs+eBp)/2-sqrt(12*tSP^2+(eAs-eBp)^2/4);


%Eb1 = (eBs+eAp)/2+sqrt(12*tSP^2+(eBs-eAp)^2/4);


%Eb2 = (eBs+eAp)/2-sqrt(12*tSP^2+(eBs-eAp)^2/4);


%Es = (eBs+eAs)/2+sqrt(36*tSS^2+(eBs-eAs)^2/4);


%Ep = (eBp+eAp)/2-sqrt((2*tPPs+4*tPPp)^2+(eBp-eAp)^2/4);

%scatter(1,Ea1,'*'); 


%scatter(1,Ea2,'*');


%plot(1,Eb1,'*','MarkerSize',12);


%plot(1,Eb2,'square','MarkerSize',12);


%plot(1,eAp,'diamond','MarkerSize',12);


%plot(L/2,Es,'o','MarkerSize',12);


%plot(L/2,Ep,'^','MarkerSize',12);

%scatter(1,eBp,'*');


