n = 250;
N = 1000;
snr = .05;
Trials = 50;
EC = [10:10:250];

Data = zeros(Trials,length(EC));
Data2 = zeros(Trials,length(EC));
Data3 = zeros(Trials,length(EC));

for(k=1:1:length(EC))
    
    L = [1:1:EC(k)];
    LC = setdiff(1:N,L);
    
    for(t = 1:1:Trials)
        
        F = randn(N,n);
        [F,~] = qr(F,0);
        F = sqrt(N/n)*F(:,1:n)';
        G = (n/N)*F;
        
        f = randn(n,1);
        f = f./norm(f);
        
        FC = G' * f;
        noise = randn(size(LC'));
        noise = snr * norm(FC(LC))/norm(noise) * noise;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        f_R = F*FC;
        
        M = (F(:,L)' * G(:,L))';
        C = (eye(length(L)) - M) \ eye(length(L));
        g = f_R + F(:,L) * (C * (G(:,L)' * f_R));
        
        Data(t,k) = norm(f-g);
        Data2(t,k) = norm(f-f_R);
        
    end
    
    k

end

X = repmat(EC,Trials,1);
X = reshape(X,[length(EC)*Trials,1]);
Y = reshape(Data,[length(EC)*Trials,1]);
Z = reshape(Data2,[length(EC)*Trials,1]);
plot(X,Y,'x');
hold on;
plot(EC,median(Data));
plot(X,Z,'o');
title('Erasure Set Size vs Reconstruction Error');
xlabel('Erasure Set Size');
ylabel('Reconstruction Error');
hold off;
