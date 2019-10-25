clc; clear all

sfile = importdata('dstran.txt');
[m n] = size(sfile);

%% The material properties for the UMAT

props = [210000.0 0.3 240. 1206.0];


%% What is the strain increment type to be used
stran = zeros(1,6); stress = stran';
index = sfile(1,1);
statev = zeros(1,2);

for i = 1:m

    if(sfile(i,1) == index+1)
        disp('----------Start-----------')
        
        sfile(i-1,1)
        
        dstran = sfile(i-1,2:7);
        stran = stran + dstran;
        dstran = -1*dstran';
        oldstress = stress;
        
 
        
        [ddsdde, stress, statev] = UMAT(stress, dstran, stran, props, statev);
        
        
        a = ddsdde*dstran;
        b = stress - oldstress;
        
        pstore(index) = statev(1);
        fystore(index) = statev(2) + props(3);
    
        stress_store(index,:) = stress;
        strain_store(index,:) = stran;
        
        index = index+1;
        disp('----------End-----------')
    end
end

%% Plot dat shit

plot(strain_store(:,2), stress_store(:,2))
plot(pstore, fystore)
